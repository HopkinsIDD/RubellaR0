
#' load_and_clean_R0_data.R
#' Purpose: Load and format serological survey data needed for R0 stan model
#' Author: Shaun Truelove, satruelove@gmail.com
#' Date updated: 30 Aug 2018



# SETUP -------------------------------------------------------------------


# Load packages


#if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
#if(!require('googlesheets')) install.packages('googlesheets'); 
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)





# LOAD AND SETUP DATA -----------------------------------------------------





# raw_files = NULL #paste0(data_dir,"/", new_data)


#' Title
#'
#' @param raw_files 
#'
#' @return
#' @export
#'
#' @examples
combine_raw_sero <- function(raw_files = NULL){
  
  sero_data <- lapply(raw_files, readr::read_csv) %>%
    data.table::rbindlist(fill = TRUE) %>% # rbind all the data together
    as_tibble() %>%
    mutate(row_NotNA = rowSums(apply(!is.na(.), 2, as.numeric))) %>%
    filter(row_NotNA>0) %>% # remove extra rows added during reading from csvs
    distinct() %>% # remove duplicate rows (must be exact duplicates so may not get rid of all duplicate records)
    dplyr::select(-row_NotNA)
  
  # create an ID for each row to remove duplicates 
  sero_data <- remove_sero_dups(sero_data)
  
  return(sero_data)
}


remove_sero_dups <- function(sero_data){
  
  revert_age_vars <- FALSE
  if ("min_age_years" %in% colnames(sero_data)){
    revert_age_vars <- TRUE
    #rename
    sero_data <- sero_data %>% rename(min.age.years = min_age_years, 
                                      max.age.years = max_age_years)
  }
  
  # create an ID for each row to remove duplicates 
  sero_data <- sero_data %>%
    rowwise() %>%
    mutate(row_id = paste(purrr::discard(c(reference_id, paste0("age", min.age.years,"-", max.age.years), gender, 
                                           country, city.province.district, region, urban.rural, ref_survey_number, 
                                           N, N.seropositive, seropositive.immune,
                                           start.year), is.na), collapse = "_")) %>%
    ungroup() %>%
    mutate(row_id = gsub(" ", "", row_id)) %>%
    filter(is.na(manual_skip) | !manual_skip)
  
  # check for duplicates
  # sum(duplicated(sero_data$row_id))
  dup_ids <- sero_data %>%
    mutate(dups = duplicated(row_id)) %>%
    filter(dups) %>%
    pull(row_id)
  # sero_data %>% filter(row_id %in% dup_ids) %>% arrange(row_id) %>% View()
  
  sero_data <- sero_data %>% 
    mutate(dup_remove = (row_id %in% dup_ids),
           row_nas = rowSums(apply(is.na(.), 2, as.numeric)))
  
  sero_data <- sero_data %>% 
    group_by(row_id) %>%
    arrange(row_id, row_nas) %>%
    mutate(dup_inds = 1:length(row_id)) %>%
    ungroup() %>%
    filter(dup_inds==1) %>%
    arrange(row_id) %>%
    dplyr::select(-any_of(c("dup_inds", "row_nas", "dup_remove")))
  
  if (revert_age_vars){
    revert_age_vars <- TRUE
    sero_data <- sero_data %>% rename(min_age_years = min.age.years, 
                                      max_age_years = max.age.years)
  }
  
  return(sero_data)
}







load_and_clean_sero_raw <- function(taiwan=TRUE, 
                                    countries=NULL, 
                                    new_data=NULL, 
                                    data_dir=NULL){
  
  # Source functions and data
  source('source/ISO_code_source.R') #loaded in the next line
  require(dplyr)
  
  sero_data <- combine_raw_sero(raw_files = new_data)
  
  # Put NAs in for all blanks
  sero_data <- sero_data %>% mutate_at(vars(colnames(.)), .funs = list(~ifelse(.=="", NA, as.character(.))))
  
  # FIX "TRUE" and "FALSE"
  sero_data <- sero_data %>% mutate_at(vars(colnames(.)), .funs = list(~ifelse(.=="TRUE", TRUE, .)))
  sero_data <- sero_data %>% mutate_at(vars(colnames(.)), .funs = list(~ifelse(.=="FALSE", FALSE, .)))
  
  # Remove any entries not wanted
  sero_data <- sero_data %>% filter(!further_review_needed | is.na(further_review_needed))
  sero_data <- sero_data %>% filter(!manual_skip | is.na(manual_skip))
  sero_data <- sero_data %>% filter(!(is.na(N_full_survey) & is.na(N) & is.na(seropositive.immune)))
  
  
  # ~ Fix Author and Year Published -----------------------------------------
  
  # get first author
  auth <- stringi::stri_match_first_regex(sero_data$reference, "(.*?)\\ ")[,2]
  auth <- gsub(",","", auth)
  auth <- gsub("\\.","", auth)
  auth2 <- stringi::stri_match_first_regex(sero_data$reference_id, "(.*?)\\_")[,2]
  auth[is.na(auth)] <- auth2[is.na(auth)]
  if (length(auth)==0) auth <- auth2
  sero_data$author <- tolower(auth)
  
  # Get published year
  yr <- stringi::stri_match_first_regex(sero_data$reference_id, "(\\w+)_(\\w+)", cg_missing="")[,3]
  # sero_data$year_published[is.na(sero_data$year_published)] <- yr[is.na(sero_data$year_published)]
  sero_data$year_published <- yr
  
  
  # ~ Clean Countries and Codes ---------------------------------------------
  
  # Change Taiwan to Taiwan
  if (taiwan){
    sero_data <- sero_data %>% mutate(country = ifelse(!is.na(city.province.district) & tolower(city.province.district)=="taiwan", "Taiwan", country))
  }
  
  # Add ISO3 and UN codes
  
  match_iso_country <- tibble(country = sero_data$country %>% unique())
  match_iso_country$ISO <- sapply(match_iso_country$country, get.iso)  
  match_iso_country$UNcode <- get.UNcode.from.ISO3(match_iso_country$ISO)
  
  # Fix Country Name
  match_iso_country$country2 <- globaltoolboxlite::get_country_name_ISO3(match_iso_country$ISO)
  
  # merge back in
  sero_data <- sero_data %>%
    dplyr::select(-any_of(c("ISO", "UNcode", "country2"))) %>%
    left_join(match_iso_country, by="country")
  
  # if restricting to specific countries
  if (!is.null(countries)){
    sero_data <- sero_data %>% filter(ISO %in% countries)
  }
  
  
  # ~ Fix age range variables -----------------------------------------------
  
  age.range <- sero_data$age.group.years
  age.min <- age.max <- age.mid <- rep(NA, length(age.range))
  
  range.yr.ind <- grepl('-', age.range) & !is.na(age.range)
  single.yr.ind <- !grepl('-', age.range) & !grepl('<', age.range) & !grepl('>', age.range) & age.range>=1 & !is.na(age.range)
  
  # single year age ranges
  age.min[single.yr.ind] <- as.numeric(age.range[single.yr.ind])
  age.max[single.yr.ind] <- as.numeric(age.range[single.yr.ind])+0.99
  
  # ranges with '-'
  tmp <- data.frame(strsplit(age.range[range.yr.ind], '-'), stringsAsFactors = F)
  age.min[range.yr.ind] <- as.numeric(tmp[1,])
  age.max[range.yr.ind] <- as.numeric(tmp[2,])+0.99 
  
  # ranges with '>'
  tmp <- data.frame(strsplit(age.range[grepl('>', age.range)], '>'), stringsAsFactors = F)
  age.min[grepl('>', age.range)] <- as.numeric(as.character(unlist(tmp[2,])))+1
  age.max[grepl('>', age.range) & age.min<50] <- 45
  age.max[grepl('>', age.range) & age.min>=50] <- 65
  
  # ranges with '<'
  tmp <- data.frame(strsplit(age.range[grepl('<', age.range)], '<'), stringsAsFactors = F)
  age.min[grepl('<', age.range)] <- 0
  age.max[grepl('<', age.range)] <- as.numeric(as.character(unlist(tmp[2,])))-.01
  
  # Add age min and max back into the dataset
  sero_data <- sero_data %>% rename(min_age_years = min.age.years, 
                                    max_age_years = max.age.years)
  
  sero_data$min_age_years[is.na(sero_data$min_age_years)] <- age.min[is.na(sero_data$min_age_years)]
  sero_data$max_age_years[is.na(sero_data$max_age_years)] <- age.max[is.na(sero_data$max_age_years)]
  sero_data$age.group.years.new <- paste0(sero_data$min_age_years, '-', sero_data$max_age_years,'yr')
  sero_data$min_age_years <- as.numeric(sero_data$min_age_years)
  sero_data$max_age_years <- as.numeric(sero_data$max_age_years)
  rm(age.min, age.max, age.range, range.yr.ind, single.yr.ind)
  
  # Add Mid-age variable
  sero_data$mid.age.years <- ((sero_data$max_age_years+0.01) + sero_data$min_age_years) / 2
  
  
  # ~ Clean Gender ----------------------------------------------------------
  
  # Get gender from Notes
  MandF <- ((grepl('\\bmale', tolower(sero_data$notes)) & grepl('female', tolower(sero_data$notes))) | 
              (grepl('\\bmen\\b', tolower(sero_data$notes)) & grepl('women', tolower(sero_data$notes))) |
              grepl('\\bboth\\b', tolower(sero_data$notes)))
  F_only <- ((grepl('female', tolower(sero_data$notes)) | grepl('women', tolower(sero_data$notes))) & !MandF)
  M_only <- ((grepl('\\bmale', tolower(sero_data$notes)) | grepl('\\bmen\\b', tolower(sero_data$notes))) & !MandF & !F_only)
  
  sero_data$gender[MandF & is.na(sero_data$gender)] <- "both"
  sero_data$gender[F_only & is.na(sero_data$gender)] <- "female"
  sero_data$gender[M_only & is.na(sero_data$gender)] <- "male"
  
  
  # ~ Reference and Survey IDs ----------------------------------------------
  
  # Make reference_id for those without
  sero_data$reference_id[is.na(sero_data$reference_id)] <- paste0(sero_data$author,"_",sero_data$year_published)[is.na(sero_data$reference_id)]
  
  # survey_id_2
  sero_data <- sero_data %>% mutate(survey_id_2 = paste0(reference_id, "--", touchstone))
  
  # Set up unique integer ids for country/survey (for each survey or reference but individual countries, as some refs have multiple)
  sero_data <- sero_data %>% mutate(country.year = paste0(country2, ' - ', start.year)) %>% 
    mutate(country.year = factor(country.year, levels=sort(unique(country.year))))
  sero_data <- sero_data %>% mutate(country.survey.year = paste0(country2, ' - ', 
                                                                 ifelse(!is.na(city.province.district), paste0(city.province.district,' - '),''),
                                                                 ifelse(!is.na(region), paste0(region,' - '),''),
                                                                 ifelse(!is.na(urban.rural), paste0(urban.rural,' - '),''),
                                                                 start.year,'-',start.year, 
                                                                 ifelse(!is.na(ref_survey_number),paste0(' - survey',ref_survey_number),'')))
  sero_data <- sero_data %>% mutate(country.survey.year = factor(country.survey.year, levels=unique(country.survey.year)))
  sero_data <- sero_data %>% mutate(survey.ind = as.integer(country.survey.year))
  
  
  # ~ Fix the Numbers -------------------------------------------------------
  
  # Convert numeric
  sero_data <- sero_data %>% mutate(N.seropositive = as.numeric(N.seropositive),
                                    N = as.numeric(N),
                                    N_full_survey = as.numeric(N_full_survey),
                                    seropositive.immune = as.numeric(seropositive.immune))
  
  source('source/get_demog_data.R') # needed for scaled mid ages
  
  # Surveys without Age numbers, just overall --> use age structure and assume sample follows that
  if (any(is.na(sero_data$N) & !is.na(sero_data$N_full_survey))){
    
    survey_fix_N <- as.character(unique(sero_data$country.survey.year[is.na(sero_data$N) & !is.na(sero_data$N_full_survey)]))
    
    rm(pop, popproj, popF, popM, popFprojMed, popMprojMed)
    
    for (s in 1:length(survey_fix_N)){
      
      UNcode_ <- sero_data$UNcode[sero_data$country.survey.year==survey_fix_N[s]][1]
      year_   <- sero_data$start.year[sero_data$country.survey.year==survey_fix_N[s]][1]
      N_full_survey_ <- as.numeric(sero_data$N_full_survey[sero_data$country.survey.year==survey_fix_N[s]][1])
      
      age_groups_ <- sero_data[sero_data$country.survey.year==survey_fix_N[s], c("min_age_years", "max_age_years")] %>% as.data.frame()
      age_groups_$max_age_years <- floor(age_groups_$max_age_years)
      
      N_ <- numeric(nrow(age_groups_))
      
      pop_ <- getPopAge(uncode = UNcode_, age.classes = c(1:100))
      pop_ <- pop_$pop.age.byageclasses.1950.2100 %>% as.data.frame() %>%
        mutate(age = 1:100) %>%
        dplyr::select(age, `year_`)
      
      for (a in 1:nrow(age_groups_)){
        pop_age_ <- pop_ %>% filter(age>=age_groups_[a,1] & age<=age_groups_[a,2])
        N_[a] <- sum(pop_age_[,2])
      }
      
      N_ <- round(N_/sum(N_) * N_full_survey_)
      
      sero_data$N[sero_data$country.survey.year==survey_fix_N[s]] <- N_
      
    }
  }
  
  
  # First fill in the missing %sero positive when Ns are available (this will fix any miscalculated %s)
  sero_data$seropositive.immune[!is.na(sero_data$N.seropositive) & !is.na(sero_data$N)] <- (sero_data$N.seropositive / sero_data$N)[!is.na(sero_data$N.seropositive) & !is.na(sero_data$N)]
  
  # For those missing N.seropositive, calculate
  sero_data$N.seropositive[is.na(sero_data$N.seropositive)] <- (sero_data$N * sero_data$seropositive.immune)[is.na(sero_data$N.seropositive)]
  
  #Remove rows with N=0
  sero_data <- sero_data %>% filter(N>0)
  
  
  # Add Birth rate per study year/location/age
  missing_years_surveyed <- is.na(sero_data$start.year) | tolower(sero_data$start.year)=="not mentioned"
  sero_data$start.year[missing_years_surveyed] <- as.integer(sero_data$year_published[missing_years_surveyed]) - 2
  sero_data$end.year[missing_years_surveyed] <- as.integer(sero_data$year_published[missing_years_surveyed]) - 2
  sero_data$start.year <- as.integer(sero_data$start.year)
  sero_data$end.year <- as.integer(sero_data$end.year)
  
  sero_data$mid.survey.year <- round(((sero_data$start.year + sero_data$end.year) / 2), 0)
  sero_data$country.midyears <-  as.factor(paste0(sero_data$ISO, ' - ', sero_data$mid.survey.year))
  unique_country_years <- as.character(unique(sero_data$country.midyears))
  
  
  # Check again for duplicates with the updated variables
  
  sero_data <- remove_sero_dups(sero_data)
  
  
  # Add Birth rate per study year/location/age
  
  #source('source/get_demog_data.R') # needed for scaled mid ages
  sero_data$birth.rate <- sero_data$scaled.mid.age <- NA
  for (cy in 1:length(unique_country_years)){
    dat_tmp <- sero_data %>% filter(as.character(country.midyears)==unique_country_years[cy])
    
    if (is.na(dat_tmp$UNcode[1]) | sum(is.na(dat_tmp$min_age_years))==nrow(dat_tmp) | 
        sum(is.na(dat_tmp$max_age_years))==nrow(dat_tmp) | is.na(dat_tmp$mid.survey.year[1])) next
    
    tmp <- GetSurveyDemo(min.age.year=dat_tmp$min_age_years, 
                         max.age.year=dat_tmp$max_age_years, 
                         mid.year=dat_tmp$mid.survey.year[1], UN.code=dat_tmp$UNcode[1])
    
    sero_data$birth.rate[as.character(sero_data$country.midyears)==unique_country_years[cy]] <- tmp$mu.survey
    sero_data$scaled.mid.age[as.character(sero_data$country.midyears)==unique_country_years[cy]] <- tmp$mid.age.year
  }
  
  # # Remove conflicting packages loaded
  # if("wpp2017" %in% (.packages())){ detach("package:wpp2017", unload=TRUE) }
  # if("reshape2" %in% (.packages())){ detach("package:reshape2", unload=TRUE) }
  # if("plyr" %in% (.packages())){ detach("package:plyr", unload=TRUE) }
  # if("fitdistrplus" %in% (.packages())){ detach("package:fitdistrplus", unload=TRUE) }
  # if("MASS" %in% (.packages())){ detach("package:MASS", unload=TRUE) }
  # library(dplyr)
  rm(pop, popF, popM, popFprojMed,popFTproj,popFT,popMT,popMprojMed,popMTproj,tfr,tfrprojMed,
     who,percentASFR,popproj,population, who_regions)
  
  # for age groups which span 1 year or less, mid-ages don't scale, so we do it manually
  sero_data <- sero_data %>% mutate(scaled.mid.age = ifelse(((!is.na(scaled.mid.age) & 
                                                                ((max_age_years - min_age_years) <= 1) |  scaled.mid.age==0)), 
                                                            (max_age_years + min_age_years) / 2, scaled.mid.age))
  
  # Sort the data based on country, year, city, region, survey, age
  sero_data <- arrange(sero_data, country2, start.year, city.province.district, 
                       ref_survey_number, region, urban.rural, scaled.mid.age)
  
  return(sero_data)
}



get_birthrate <- function(year=2002, UN.code=586){
  year[which(year<1950)] <- 1950
  out <- getDemography(UN.code, seq(1, 101, 1))
  cbr <- out$cbr.1950.2100[year-1950+1]/1000
  return(cbr)
}



clean_sero_for_analysis <- function(sero_data=NULL, 
                                    save.data=TRUE, 
                                    version=NULL, 
                                    clean_sero_data = NULL,
                                    pull_raw_data=TRUE, ...){
  
  if (is.null(sero_data)){
    if (pull_raw_data) {
      sero_data <- load_and_clean_sero_raw(taiwan, 
                                           countries, 
                                           new_data, 
                                           data_dir)
      dir.create(file.path("data", version), recursive = TRUE)
      
      # read clean data if present
      if (!is.null(clean_sero_data)){
        sero_data_cl <- readr::read_csv(clean_sero_data)
        sero_data <- sero_data %>% bind_rows(sero_data_cl)
      }
      readr::write_csv(sero_data, file.path("data",version,"sero_data_raw.csv"))
    } else {
      if (is.null(version)){
        sero_data <- read_csv("data/sero_data_raw.csv")
      } else {
        sero_data <- readr::read_csv(file.path("data",version,"sero_data_raw.csv")) %>% as.data.frame()
      }
    }
  }
  
  
  
  
  # Remove Manually determined skips
  # sero_data <- sero_data %>% filter(manual_skip!=1) # We now do this in `load_and_clean_sero_data`
  
  # Remove rows with missing countries or ages
  sero_data <- sero_data %>% filter(!is.na(country))
  
  # Remove those with missing age at this point
  sero_data <- sero_data %>% filter(!is.na(min_age_years) | !is.na(max_age_years))
  
  # Remove ages where max is <6 months
  sero_data <- sero_data %>% filter(!is.na(max_age_years) & max_age_years>.5)
  
  # Sort the data based on country, year, city, region, survey, age
  sero_data <- arrange(sero_data, country2, start.year, city.province.district, 
                       ref_survey_number, region, urban.rural, scaled.mid.age)
  
  
  # # Set up unique integer ids for references
  # sero_data$reference <- factor(sero_data$reference, levels = sort(unique(sero_data$reference)))
  # sero_data$refs.ind <- as.integer(sero_data$reference)
  
  # Set up unique integer ids for countries 
  sero_data$ISO <- factor(sero_data$ISO, levels=sort(unique(sero_data$ISO)))
  sero_data$country2 <- factor(sero_data$country2, levels = sort(unique(sero_data$country2)))
  sero_data$country.ind <- as.integer(sero_data$country2)
  sero_data$sub.region <- get.subregion(sero_data$ISO)
  sero_data$sub.region <- factor(sero_data$sub.region, levels = sort(unique(sero_data$sub.region)))
  sero_data$region.ind <- as.integer(sero_data$sub.region)
  
  # Check for duplicate study entries
  sero_data <- sero_data %>% filter(!(reference_id=="desinor_2004" & entered_by=="IG")) # remove known duplicate
  sero_data <- remove_sero_dups(sero_data)
  
  # Set up unique integer ids for country/survey (for each survey or reference but individual countries, as some refs have multiple)
  sero_data <- sero_data %>% mutate(country.year = paste0(country2, ' - ', start.year)) 
  sero_data <- sero_data %>% mutate(country.year = factor(country.year, levels=sort(unique(country.year))))
  sero_data <- sero_data %>% mutate(country.survey.year = paste0(country2, ' - ', 
                                                                 ifelse(!is.na(city.province.district), paste0(city.province.district,' - '),''),
                                                                 ifelse(!is.na(region), paste0(region,' - '),''),
                                                                 ifelse(!is.na(urban.rural), paste0(urban.rural,' - '),''),
                                                                 start.year,'-',start.year, 
                                                                 ifelse(!is.na(ref_survey_number),paste0(' - survey',ref_survey_number),'')))
  sero_data <- sero_data %>% mutate(country.survey.year = factor(country.survey.year, levels=unique(country.survey.year)))
  sero_data <- sero_data %>% mutate(survey.ind = as.integer(country.survey.year))
  
  
  # ~ Age Group Inclusion Criteria  -----------------------------------------
  
  # Note studies which do not have at least 2 age groups above 15 years and 2 age groups below
  sero_data <- sero_data %>% 
    group_by(survey.ind) %>% 
    mutate(meet_age_crit         = ifelse((sum(max_age_years<15)>=2 & sum(min_age_years>=15)>=2) |
                                            (sum(max_age_years<15)>=2 & sum(min_age_years<=15 & max_age_years>=15)==1 & sum(min_age_years>=15)>=1), TRUE, FALSE), 
           meet_age_crit_partial = ifelse(sum(max_age_years<15)>=1 & sum(min_age_years<=15 & max_age_years>=15)==1 & sum(min_age_years>=15)>=1, TRUE, FALSE))
  
  
  # Remove columns before joining
  cols_ <- c("VIMC_priority", "WHO_region", "gavi_elig", "who_priority_68", "vacc_in_schedule", "year_intro_full","year_intro_part")
  cols_to_remove <- cols_[cols_ %in% colnames(sero_data)]
  sero_data <- sero_data %>% dplyr::select(-cols_to_remove)
  
  
  # ~ Add VIMC Priority Countries -------------------------------------------
  
  # Note if a priority country for VIMC
  priority_countries <- read_csv(file="data/vimc_countries_of_interest.csv") %>%
    as.data.frame() %>% mutate(req_country=TRUE)
  
  sero_data <- left_join(sero_data, priority_countries %>% dplyr::select(-country), by=c("ISO"="iso3")) %>% 
    #select(-VIMC_priority, -year_intro_full, -year_intro_part, -who_priority_68, -WHO_region, -vacc_in_schedule) %>%
    rename(VIMC_priority=req_country) %>% mutate(VIMC_priority=ifelse(is.na(VIMC_priority), FALSE, VIMC_priority))
  
  
  # ~ Add Year Vaccine Intro & If Serosurvey After --------------------------
  
  year_vacc_intro <- read_csv(file="data/year_vaccine_introduction_rubella.csv") %>% as.data.frame() %>% select(-Country)
  colnames(year_vacc_intro) <- c("ISO", "WHO_region", "gavi_elig", "who_priority_68", "vacc_in_schedule", "year_intro_full","year_intro_part")
  year_vacc_intro$year_intro_full[year_vacc_intro$year_intro_full=="n/a"] <- NA
  year_vacc_intro$year_intro_part[year_vacc_intro$year_intro_part=="n/a"] <- NA
  year_vacc_intro$vacc_in_schedule <- tolower(year_vacc_intro$vacc_in_schedule) 
  
  sero_data <- left_join(sero_data, 
                         year_vacc_intro, 
                         by=c("ISO"="ISO")) %>%
    mutate(prior_to_vacc_full = (end.year<=year_intro_full), prior_to_vacc_part = (end.year<=year_intro_part)) 
  
  
  # Sort the data based on country, survey, year, age
  sero_data <- arrange(sero_data, country2, start.year, city.province.district, ref_survey_number, region, urban.rural, scaled.mid.age)
  
  
  
  # ~ Combine Male and Female -----------------------------------------------
  # Combine Male and Female when detailed separately for same age and survey
  
  sero_data$survey_age_ind <- paste0(sero_data$survey.ind,': ', sero_data$age.group.years.new)
  sero_data_clean <- sero_data # Keep an original without these combined
  
  survey.age.unique <- unique(sero_data$survey_age_ind)
  for (i in 1:length(survey.age.unique)){
    dat_keep <- NULL
    dat_tmp <- sero_data %>% filter(survey_age_ind == survey.age.unique[i])
    if (nrow(dat_tmp)==1) next
    
    # Combine Genders
    if ("both" %in% dat_tmp$gender){
      dat_keep <- dat_tmp %>% filter(gender=='both')
    } else if ("male" %in% dat_tmp$gender & "female" %in% dat_tmp$gender & nrow(dat_tmp)==2){
      dat_keep <- dat_tmp[1,]
      dat_keep <- dat_keep %>% mutate(gender='both', N=sum(dat_tmp$N), N.seropositive=sum(dat_tmp$N.seropositive))
      dat_keep$seropositive.immune <- dat_keep$N.seropositive / dat_keep$N
    }
    
    if (!is.null(dat_keep)){
      sero_data_clean <- sero_data_clean %>% filter(survey_age_ind != survey.age.unique[i])
      sero_data_clean <- rbind(sero_data_clean, dat_keep)
    } else {
      print(paste0(i, " -- ", survey.age.unique[i], " -- ", substr(dat_tmp$reference[1], 1,10)))
    }
  }
  
  # ~ Save ------------------------------------------------------------------
  
  sero_data <- sero_data %>% arrange(country, survey_id_2, gender, start.year, urban.rural, min_age_group)
  sero_data_clean <- sero_data_clean %>% arrange(country, survey_id_2, gender, start.year, urban.rural, min_age_group)
  
  if (save.data){
    sero_data_tmp  <- tibble(sero_data)
    # sero_data_tmp  <- data.frame(apply(sero_data_tmp,2,as.character))
    sero_data_tmp2 <- tibble(sero_data_clean)
    # sero_data_tmp2 <- data.frame(apply(sero_data_tmp2,2,as.character))
    
    
    if (is.null(version)){
      write.csv(sero_data_tmp, file="data/sero_data_orig.csv", row.names=FALSE)
      write.csv(sero_data_tmp2, file="data/sero_data_clean.csv", row.names=FALSE)
    } else {
      write.csv(sero_data_tmp, file=file.path("data", version, "sero_data_orig.csv"), row.names=FALSE)
      write.csv(sero_data_tmp2, file=file.path("data", version, "sero_data_clean.csv"), row.names=FALSE)
    }
  }
  
  return(list(sero_data_orig=sero_data_tmp, sero_data_clean=sero_data_tmp2))
}




# Final data cleaning for Stan ----
get_sero_data <- function(meets_crit = TRUE, 
                          adjust_low = TRUE, 
                          run_from_saved = TRUE, 
                          pull_raw_data = NULL,
                          analysis_version, 
                          min_age_group=0.6, 
                          ...) {
  
  require(dplyr)
  
  if (run_from_saved) {
    sero_data <- readr::read_csv(file.path("data", analysis_version, "sero_data_clean.csv")) %>% as.data.frame() 
    
  } else {
    sero_data_list <- clean_sero_for_analysis(sero_data=NULL, 
                                              save.data=TRUE, 
                                              version=analysis_version, 
                                              pull_raw_data = pull_raw_data,
                                              ...)
    sero_data <- sero_data_list$sero_data_clean %>% as.data.frame() 
    rm(sero_data_list)
  } 
  
  sero_data <- sero_data %>% filter(!(country.survey.year %in% 
                                        c(#"Singapore - Toa Payoh - 1986-1986", 
                                          "Japan - Sapporo - 1967-1967"#,
                                          #"US - Atlanta - 1967-1967",
                                          #"Saudi Arabia - Mult. - 1989-1989 - survey1"
                                        )))
  
  # Clean data to format for age-specific estimation
  sero_data <- sero_data %>% filter(!is.na(N) & N>0 &   # make sure no surveys with missing data & remove age groups with 0 data
                                      max_age_years>=min_age_group & # Remove any before X months of age (default is 7 months, or 0.6: 6.99/12 = 0.58)
                                      !is.na(scaled.mid.age)) %>%
    dplyr::rename(Y=N.seropositive) # rename variable
  
  # subtract 6 months from age (maternal immunity)
  sero_data <- sero_data %>% mutate(scaled.mid.age = ifelse(max_age_years>=1, scaled.mid.age-.5, (max_age_years - .5)/2))
  
  if (adjust_low) {
    # Age groups with max above Remove some of positives below age 6 mo
    #  - to attempt to fix the problem of maternal antibodies, we will try to just remove those children
    age.span <- sero_data$max_age_years - sero_data$min_age_years
    n <- with(sero_data, ifelse(min_age_years<0.5, round(N*((.5-min_age_years)/age.span)), 0))
    sero_data$N <- sero_data$N - n
    sero_data$Y <-sero_data$Y - n
    sero_data$Y[sero_data$Y<0] <- 0
  }
  
  # Restrict to those meeting the specified age criteria
  if (tolower(meets_crit)=="full"){
    sero_data <- sero_data %>% filter(meet_age_crit==TRUE)
  } else if (tolower(meets_crit)=="partial"){
    sero_data <- sero_data %>% filter(meet_age_crit_partial==TRUE | meet_age_crit==TRUE)
  } 
  
  # Restrict to studies occurring before vaccination was introduced
  sero_data <- sero_data %>% mutate(prior_to_vacc_full = ifelse(is.na(prior_to_vacc_full), TRUE, prior_to_vacc_full),
                                    prior_to_vacc_part = ifelse(is.na(prior_to_vacc_part), TRUE, prior_to_vacc_part)) %>%
    filter(prior_to_vacc_full==TRUE & prior_to_vacc_part==TRUE)
  
  # sort before we assign factor indicators
  sero_data <- sero_data %>% arrange(country, mid.survey.year, country.survey.year, mid.age.years)
  sero_data <- sero_data %>% mutate(num_ = 1) %>% group_by(survey.ind) %>% mutate(survey_within_age_ind = cumsum(num_)) %>% 
    dplyr::select(-num_) %>% ungroup()
  
  # Rework survey and country numbers for Stan
  sero_data$survey.num <- as.integer(as.factor(sero_data$survey.ind))              # Index id number for each survey 
  sero_data$country.num <- as.integer(as.factor(sero_data$country.ind))            # Index id number for country for each survey
  sero_data$region.num <- as.integer(as.factor(sero_data$region.ind))            # Index id number for country for each survey
  sero_data$survey.age.num <- 1:nrow(sero_data)                                  # Index id number for country for each survey
  
  return(sero_data)
}





##' Functions to reformat things that stan needs after subsetting to specific countries or surveys

reformat_sero_subset <- function(sero_data){
  
  # Rework survey and country numbers for Stan
  sero_data$survey.num  <- as.integer(as.factor(sero_data$survey.num))        # Index id number for each survey 
  sero_data$country.num <- as.integer(as.factor(sero_data$country.num))      # Index id number for country for each survey
  sero_data$region.num  <- as.integer(as.factor(sero_data$region.num))        # Index id number for country for each survey
  sero_data$survey.age.num <- 1:nrow(sero_data) 
  
  return(sero_data)
}


add_foi_years_sero <- function(sero_data){
  
  # First, get calendar years for each survey age group 
  sero_data$survey_mid_year <- floor((sero_data$start.year + sero_data$end.year) / 2)
  sero_data <- sero_data %>% mutate(min_years_surveys = floor(survey_mid_year - mid.age.years),
                                    max_years_surveys = survey_mid_year)
  sero_data <- sero_data %>% group_by(country) %>%
    mutate(country_min_year = min(min_years_surveys, na.rm = TRUE),
           country_max_year = max(max_years_surveys, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(country_year_ind_min = min_years_surveys - country_min_year + 1, 
           country_year_ind_max = max_years_surveys - country_min_year + 1)
  
  #save(sero_data, file="data/sero_data.RData")
  return(sero_data)
}





format_list_sero_for_stan <- function(sero_data){
  
  # # SUMMARY DATA OF SURVEY
  # sero_data_survey <- sero_data[!duplicated(sero_data$survey.ind),] 
  #
  # # MATRIX OF POPULATION BY AGE
  source('source/get_demog_data.R') # needed for scaled mid ages
  
  # N_age <- matrix(NA, nrow = nrow(sero_data_survey), ncol = 100, dimnames = list(NULL, 0:99))
  # for (i in 1:nrow(sero_data_survey)) {
  #   N_age[i,] <- get_N_age(0:99, UNcode=sero_data_survey$UNcode[i], year=sero_data_survey$start.year[i])
  # }
  # N_age[N_age<0] <- 0
  
  # set up stan list data
  sero_data_stan <- list()
  sero_data_stan[["C"]] <- as.integer(length(unique(sero_data$ISO)))          # Number unique countries
  sero_data_stan[["S"]] <- as.integer(length(unique(sero_data$survey.ind)))   # Number unique samples/serosurveys
  sero_data_stan[["R"]] <- as.integer(length(unique(sero_data$region.ind)))   # Number unique regions
  sero_data_stan[["SA"]] <- as.integer(nrow(sero_data))                       # Number unique samples/serosurveys x age groups
  sero_data_stan[["N"]] <- as.integer(sero_data$N)                            # Number indivs in each age grp in each serosurvey/sample
  sero_data_stan[["Y"]] <- as.integer(sero_data$Y)                            # Number seropositive in each age in each serosurvey/sample
  sero_data_stan[["age"]] <- as.numeric(sero_data$scaled.mid.age)             # Mid age value of each age group
  #sero_data_stan[["mu"]] <- as.numeric(sero_data$birth.rate)[!duplicated(sero_data$survey.ind)]     # Birth rate for that country during the year of the survey
  sero_data_stan[["country_name"]] <- sero_data$country2            # Index id number for country for each survey
  #sero_data_stan[["country_ind"]] <- sero_data$country.ind            # Index id number for country for each survey
  sero_data_stan[["country"]] <- sero_data$country.num            # Index id number for country for each survey
  sero_data_stan[["survey"]] <- sero_data$survey.num              # Index id number for each survey 
  sero_data_stan[["survey_age"]] <- sero_data$survey.age.num              # Index id number for each survey/age 
  sero_data_stan[["region_of_survey"]] <- sero_data$region.num[!duplicated(sero_data$survey.num)]     # Index id number for region of each survey 
  # Get which country index belongs to each survey
  sero_data_stan[["country_of_survey"]] <- sero_data_stan[["country"]][!duplicated(sero_data_stan[["survey"]])] 
  sero_data_stan[["region_of_country"]] <- sero_data$region.num[!duplicated(sero_data$country.num)] 
  #sero_data_stan[["age_pred_length"]] <- 100                                      
  #sero_data_stan[["age_pred"]] <- 0:99                                      # Population by age for each surveyPopulation by age for each survey
  sero_data_stan[["mu"]] <- sero_data$birth.rate
  
  
  # Generate design matrices for survey-age-group, survey, and country --- all by year
  
  # First, convert to long
  survey_years_long <- data.frame(surveyage=NA, survey=NA, country=NA, year=NA)#, surveyagemid=NA)
  for (i in 1:nrow(sero_data)){
    survey_years_long <- rbind(survey_years_long, 
                               data.frame(surveyage=sero_data$survey.age.num[i],
                                          survey=sero_data$survey.num[i], country=sero_data$country.num[i], 
                                          year=sero_data$min_years_surveys[i] : sero_data$max_years_surveys[i]))#,
    #surveyagemid=sero_data$scaled.mid.age[i]))
  }
  survey_years_long <- survey_years_long[-1,] %>% mutate(year_present=1) 
  
  # years_country <- survey_years_long %>% filter(year_present==1) %>% group_by(country) %>% filter(!duplicated(year)) %>% ungroup() %>%
  #   arrange(country, year) %>%
  #   mutate(year_ind = as.integer(as.factor(year))) %>% select(country, year, year_ind)
  
  # Next, generate the Design Matrices
  # - for survey-age-group
  survey_years_sa <- survey_years_long  %>% mutate(surveyage_years = paste0(surveyage,"-",year)) %>% filter(!duplicated(surveyage_years)) %>%
    select(-surveyage_years) %>% spread(key=year, value=year_present, fill=0)
  # - for survey
  survey_years_s <- survey_years_long %>% mutate(survey_years = paste0(survey,"-",year)) %>% filter(!duplicated(survey_years)) %>%
    select(-surveyage, -survey_years) %>% spread(key=year, value=year_present, fill=0)
  # - for country
  survey_years_c <- survey_years_long %>% mutate(country_years = paste0(country,"-",year)) %>% filter(!duplicated(country_years)) %>%
    select(-surveyage, -survey, -country_years) %>% spread(key=year, value=year_present, fill=0)
  
  # Add to stan list data
  sero_data_stan[['foi_year_range']] <- range(survey_years_long$year, na.rm=TRUE) # add the range for reference
  sero_data_stan[['years_length']] <- diff(range(survey_years_long$year, na.rm=TRUE)) + 1 # number of unique years
  sero_data_stan[['foi_year_mat_sa']] <- survey_years_sa %>% #arrange(country) %>%
    dplyr::select(-surveyage, -survey, -country) %>% as.matrix()
  sero_data_stan[['foi_year_mat_s']] <- survey_years_s %>% arrange(survey) %>% 
    dplyr::select(-country, -survey) %>% as.matrix()
  sero_data_stan[['foi_year_mat_c']] <- survey_years_c %>% arrange(country) %>% 
    dplyr::select(-country) %>% as.matrix()
  foi_years <- sort(unique(survey_years_long$year))
  
  # Birthrate by Country/year
  mu_tc <- matrix(0, nrow=length(unique(sero_data$country.num)), ncol=sero_data_stan$years_length)
  colnames(mu_tc) <- foi_years
  un_codes <- sero_data$UNcode[!duplicated(sero_data$country.num)]
  row.names(mu_tc) <- sort(unique(sero_data$country.num))
  
  for (unc in 1:length(un_codes)){
    mu_tc[unc,] <- get_birthrate(year=foi_years, UN.code = un_codes[unc])
  }
  sero_data_stan[['mu_tc']] <- mu_tc
  
  if("wpp2017" %in% (.packages())){ detach("package:wpp2017", unload=TRUE) }
  if("reshape2" %in% (.packages())){ detach("package:reshape2", unload=TRUE) }
  if("plyr" %in% (.packages())){ detach("package:plyr", unload=TRUE) }
  if("fitdistrplus" %in% (.packages())){ detach("package:fitdistrplus", unload=TRUE) }
  if("MASS" %in% (.packages())){ detach("package:MASS", unload=TRUE) }
  library(dplyr)
  
  rm(pop, popF, popM, popFprojMed,popFTproj,popFT,popMT,popMprojMed,popMTproj,tfr,tfrprojMed,
     who_regions,percentASFR,popproj,region_data,iso_data,iso_data_full,dhs_countrydata,nationality_data,UNcode)
  
  return(sero_data_stan) 
}






format_list_sero_for_stan_AGE <- function(sero_data){
  
  # SUMMARY DATA OF SURVEY
  sero_data_survey <- sero_data[!duplicated(sero_data$survey.ind),]
  
  # MATRIX OF POPULATION BY AGE
  source('source/get_demog_data.R') # needed for scaled mid ages
  
  N_age <- matrix(NA, nrow = nrow(sero_data_survey), ncol = 100, dimnames = list(NULL, 0:99))
  for (i in 1:nrow(sero_data_survey)) {
    N_age[i,] <- get_N_age(0:99, UNcode=sero_data_survey$UNcode[i], year=sero_data_survey$start.year[i])
  }
  N_age[N_age<0] <- 0
  
  # set up stan list data
  sero_data_stan <- list()
  sero_data_stan[["C"]] <- as.integer(length(unique(sero_data$ISO)))          # Number unique countries
  sero_data_stan[["S"]] <- as.integer(length(unique(sero_data$survey.ind)))   # Number unique samples/serosurveys
  sero_data_stan[["R"]] <- as.integer(length(unique(sero_data$region.ind)))   # Number unique regions
  sero_data_stan[["SA"]] <- as.integer(nrow(sero_data))                       # Number unique samples/serosurveys x age groups
  sero_data_stan[["N"]] <- as.integer(sero_data$N)                            # Number indivs in each age grp in each serosurvey/sample
  sero_data_stan[["Y"]] <- as.integer(sero_data$Y)                            # Number seropositive in each age in each serosurvey/sample
  sero_data_stan[["age"]] <- as.numeric(sero_data$scaled.mid.age)             # Mid age value of each age group
  sero_data_stan[["max_age"]] <- as.integer(round(sero_data$scaled.mid.age,0))             # Mid age value of each age group
  sero_data_stan[["country_name"]] <- sero_data$country2            # Index id number for country for each survey
  sero_data_stan[["country"]] <- sero_data$country.num            # Index id number for country for each survey
  sero_data_stan[["survey"]] <- sero_data$survey.num              # Index id number for each survey 
  sero_data_stan[["survey_age"]] <- sero_data$survey.age.num              # Index id number for each survey/age 
  sero_data_stan[["region_of_survey"]] <- sero_data$region.num[!duplicated(sero_data$survey.num)]     # Index id number for region of each survey 
  # Get which country index belongs to each survey
  sero_data_stan[["country_of_survey"]] <- sero_data_stan[["country"]][!duplicated(sero_data_stan[["survey"]])] 
  sero_data_stan[["region_of_country"]] <- sero_data$region.num[!duplicated(sero_data$country.num)] 
  sero_data_stan[["mu"]] <- sero_data$birth.rate
  sero_data_stan[["N_age"]] <- N_age                                      # Population by age for each survey
  #sero_data_stan[["mu"]] <- sero_data_survey$birth.rate
  sero_data_stan[["M"]] <- 0.5
  
  
  # Generate design matrices for survey-age-group, survey, and country --- all by age
  
  # First, convert to long
  survey_age_long <- data.frame(surveyage=NA, survey=NA, survey_year=NA, country=NA, age=NA, survey_age=NA, age_year=NA)#, surveyagemid=NA)
  for (i in 1:nrow(sero_data)){
    survey_age_long <- rbind(survey_age_long, 
                             data.frame(surveyage=sero_data$survey.age.num[i],
                                        survey=sero_data$survey.num[i], 
                                        survey_year=sero_data$start.year[i],
                                        country=sero_data$country.num[i], 
                                        age= 0 : sero_data$scaled.mid.age[i],
                                        survey_age=sero_data$scaled.mid.age[i],
                                        age_year=NA))  }
  survey_age_long <- survey_age_long[-1,] %>% mutate(age_present=1,
                                                     age_year = survey_year - floor(survey_age - age)) %>%
    dplyr::select(-survey_age, -survey_year)
  
  # Save years for Forces of Infection
  foi_year <- sort(unique(survey_age_long$age_year))
  foi_age <- sort(unique(survey_age_long$age))
  
  
  # Next, generate the Design Matrices
  # - for survey-age-group
  survey_age_sa <- survey_age_long  %>% mutate(surveyage_age = paste0(surveyage,"-",age)) %>% filter(!duplicated(surveyage_age)) %>%
    select(-surveyage_age, -age_year) %>% spread(key=age, value=age_present, fill=0)
  # - for survey
  survey_age_s <- survey_age_long %>% mutate(survey_age = paste0(survey,"-",age)) %>% filter(!duplicated(survey_age)) %>%
    select(-surveyage, -survey_age, -age_year) %>% spread(key=age, value=age_present, fill=0)
  # - for country
  survey_age_c <- survey_age_long %>% mutate(country_age = paste0(country,"-",age)) %>% filter(!duplicated(country_age)) %>%
    select(-surveyage, -survey, -country_age, -age_year) %>%  spread(key=age, value=age_present, fill=0)
  # - year of each age per survey
  survey_ageyear_s <- survey_age_long %>% mutate(survey_age = paste0(survey,"-",age)) %>% filter(!duplicated(survey_age)) %>%
    select(-surveyage, -survey_age, -age_present) %>% spread(key=age, value=age_year, fill=NA)
  
  
  # Add to stan list data
  sero_data_stan[['foi_age_range']] <- range(survey_age_long$age, na.rm=TRUE) # add the range for reference
  sero_data_stan[['age_length']] <- diff(range(survey_age_long$age, na.rm=TRUE)) + 1 # number of unique age
  
  sero_data_stan[['foi_age_mat_sa']] <- survey_age_sa %>% #arrange(country) %>%
    dplyr::select(-surveyage, -survey, -country) %>% as.matrix()
  sero_data_stan[['foi_age_mat_s']] <- survey_age_s %>% arrange(survey) %>% 
    dplyr::select(-country, -survey) %>% as.matrix()
  sero_data_stan[['foi_age_mat_c']] <- survey_age_c %>% arrange(country) %>% 
    dplyr::select(-country) %>% as.matrix()
  
  # Birthrate by Country/year
  mu_ts <- matrix(0, nrow=length(unique(sero_data$survey.num)), ncol=sero_data_stan$age_length)
  colnames(mu_ts) <- foi_age
  row.names(mu_ts) <- 1:sero_data_stan$S
  
  un_codes <- sero_data$UNcode[!duplicated(sero_data$survey.num)]
  
  for (s in 1:sero_data_stan$S){
    mu_ts[s,] <- get_birthrate(year=as.integer(survey_ageyear_s[s,-(1:2)], na.rm=TRUE), UN.code = un_codes[s])
  }
  mu_ts[is.na(mu_ts)] <- 0
  sero_data_stan[['mu_ts']] <- mu_ts
  
  if("wpp2017" %in% (.packages())){ detach("package:wpp2017", unload=TRUE) }
  if("reshape2" %in% (.packages())){ detach("package:reshape2", unload=TRUE) }
  if("plyr" %in% (.packages())){ detach("package:plyr", unload=TRUE) }
  if("fitdistrplus" %in% (.packages())){ detach("package:fitdistrplus", unload=TRUE) }
  if("MASS" %in% (.packages())){ detach("package:MASS", unload=TRUE) }
  library(dplyr)
  
  rm(pop, popF, popM, popFprojMed,popFTproj,popFT,popMT,popMprojMed,popMTproj,tfr,tfrprojMed,
     who_regions,percentASFR,popproj,region_data,iso_data,iso_data_full,dhs_countrydata,nationality_data,UNcode)
  
  return(sero_data_stan) 
}


get_N_age <- function(a, UNcode=UNcode, year=year_tmp){
  
  #needed packages
  require(wpp2017)
  require(zoo)
  
  year <- round(year / 5) * 5 # get closest 5-yr year
  
  data(pop) #populations to 2015
  ages.by5 <- popF[popF$country_code==UNcode,3]
  mid.age.by5 <- seq(2.5, 102.5, 5)
  popF <- popF[popF$country_code==UNcode, match(year, colnames(popF))]
  popM <- popM[popM$country_code==UNcode, match(year, colnames(popM))]
  pop_MF <- popF + popM
  pop_tmp <- data.frame(age_range = ages.by5, age=mid.age.by5, pop=pop_MF*1000)
  pop_tmp[is.na(pop_tmp)] <- 0
  f <- smooth.spline(x=pop_tmp$age, y=pop_tmp$pop)
  predpop <- predict(f,a)$y * (diff(a)[1] / 5) # previously 5yr age groups
  predpop[a>=min(pop_tmp$age[pop_tmp$pop==0])] <- 0
  return(predpop)
}


get_N_age_new <- function(a, UNcode, year=year_tmp){
  
  if (!exists("popB1")){
    require(wpp2022)
    data(popB1) #populations to 2015
  }
  
  pop_age_loc <- popB1 %>% filter(country_code==UNcode, age %in% a) %>% 
    arrange(age) %>%
    dplyr::select(matches(as.character(year))) %>% as.vector() %>% unlist() %>% as.numeric()
  pop_age_loc <- pop_age_loc * 1000
  
  return(pop_age_loc)
}




format_list_sero_for_stan_OLD <- function(sero_data, get_age=TRUE, N_age=NULL){
  
  # SUMMARY DATA OF SURVEY
  sero_data_survey <- sero_data[!duplicated(sero_data$survey.ind),]
  
  if (get_age){
    # MATRIX OF POPULATION BY AGE
    source('source/get_demog_data.R') # needed for scaled mid ages
    
    require(wpp2022)
    data(popB1) #populations to 2015
    
    N_age <- matrix(NA, nrow = nrow(sero_data_survey), ncol = 100, dimnames = list(NULL, 0:99))
    for (i in 1:nrow(sero_data_survey)) {
      N_age[i,] <- get_N_age_new(0:99, UNcode=sero_data_survey$UNcode[i], year=sero_data_survey$start.year[i])
    }
    N_age[N_age<0] <- 0
  }
  
  # set up stan list data
  sero_data_stan <- list()
  sero_data_stan[["survey_ind"]] <- unique(sero_data$survey.ind)          # Number unique countries
  sero_data_stan[["C"]] <- as.integer(length(unique(sero_data$ISO)))          # Number unique countries
  sero_data_stan[["S"]] <- as.integer(length(unique(sero_data$survey.ind)))   # Number unique samples/serosurveys
  sero_data_stan[["R"]] <- as.integer(length(unique(sero_data$region.ind)))   # Number unique regions
  sero_data_stan[["SA"]] <- as.integer(nrow(sero_data))                       # Number unique samples/serosurveys x age groups
  sero_data_stan[["N"]] <- as.integer(sero_data$N)                            # Number indivs in each age grp in each serosurvey/sample
  sero_data_stan[["Y"]] <- as.integer(sero_data$Y)                            # Number seropositive in each age in each serosurvey/sample
  sero_data_stan[["age"]] <- as.numeric(sero_data$scaled.mid.age)             # Mid age value of each age group
  sero_data_stan[["country_name"]] <- sero_data$country2            # Index id number for country for each survey
  sero_data_stan[["country"]] <- sero_data$country.num            # Index id number for country for each survey
  sero_data_stan[["survey"]] <- sero_data$survey.num              # Index id number for each survey 
  sero_data_stan[["survey_age"]] <- sero_data$survey.age.num              # Index id number for each survey/age 
  sero_data_stan[["region_of_survey"]] <- sero_data$region.num[!duplicated(sero_data$survey.num)]     # Index id number for region of each survey 
  # Get which country index belongs to each survey
  sero_data_stan[["country_of_survey"]] <- sero_data_stan[["country"]][!duplicated(sero_data_stan[["survey"]])] 
  sero_data_stan[["region_of_country"]] <- sero_data$region.num[!duplicated(sero_data$country.num)] 
  sero_data_stan[["age_pred_length"]] <- 100                                      # Population by age for each survey
  sero_data_stan[["age_pred"]] <- 0:99                                      # Population by age for each survey
  sero_data_stan[["N_age"]] <- N_age                                      # Population by age for each survey
  sero_data_stan[["mu"]] <- sero_data_survey$birth.rate
  sero_data_stan[["M"]] <- 0.5
  
  # Generate design matrices for survey-age-group, survey, and country --- all by age
  # 
  # # First, convert to long
  # survey_age_long <- data.frame(surveyage=NA, survey=NA, survey_year=NA, country=NA, age=NA, survey_age=NA, age_year=NA)#, surveyagemid=NA)
  # for (i in 1:nrow(sero_data)){
  #   survey_age_long <- rbind(survey_age_long, 
  #                            data.frame(surveyage=sero_data$survey.age.num[i],
  #                                       survey=sero_data$survey.num[i], 
  #                                       survey_year=sero_data$start.year[i],
  #                                       country=sero_data$country.num[i], 
  #                                       age= 0 : sero_data$scaled.mid.age[i],
  #                                       survey_age=sero_data$scaled.mid.age[i],
  #                                       age_year=NA))  }
  # survey_age_long <- survey_age_long[-1,] %>% mutate(age_present=1,
  #                                                    age_year = survey_year - floor(survey_age - age)) %>%
  #   dplyr::select(-survey_age, -survey_year)
  # 
  # # Save years for Forces of Infection
  # foi_year <- sort(unique(survey_age_long$age_year))
  # foi_age <- sort(unique(survey_age_long$age))
  # 
  # 
  # # Next, generate the Design Matrices
  # # - for survey-age-group
  # survey_age_sa <- survey_age_long  %>% mutate(surveyage_age = paste0(surveyage,"-",age)) %>% filter(!duplicated(surveyage_age)) %>%
  #   select(-surveyage_age, -age_year) %>% spread(key=age, value=age_present, fill=0)
  # # - for survey
  # survey_age_s <- survey_age_long %>% mutate(survey_age = paste0(survey,"-",age)) %>% filter(!duplicated(survey_age)) %>%
  #   select(-surveyage, -survey_age, -age_year) %>% spread(key=age, value=age_present, fill=0)
  # # - for country
  # survey_age_c <- survey_age_long %>% mutate(country_age = paste0(country,"-",age)) %>% filter(!duplicated(country_age)) %>%
  #   select(-surveyage, -survey, -country_age, -age_year) %>%  spread(key=age, value=age_present, fill=0)
  # # - year of each age per survey
  # survey_ageyear_s <- survey_age_long %>% mutate(survey_age = paste0(survey,"-",age)) %>% filter(!duplicated(survey_age)) %>%
  #   select(-surveyage, -survey_age, -age_present) %>% spread(key=age, value=age_year, fill=NA)
  # 
  # 
  # # Add to stan list data
  # sero_data_stan[['foi_age_range']] <- range(survey_age_long$age, na.rm=TRUE) # add the range for reference
  # sero_data_stan[['age_length']] <- diff(range(survey_age_long$age, na.rm=TRUE)) + 1 # number of unique age
  # 
  # sero_data_stan[['foi_age_mat_sa']] <- survey_age_sa %>% #arrange(country) %>%
  #   dplyr::select(-surveyage, -survey, -country) %>% as.matrix()
  # sero_data_stan[['foi_age_mat_s']] <- survey_age_s %>% arrange(survey) %>% 
  #   dplyr::select(-country, -survey) %>% as.matrix()
  # sero_data_stan[['foi_age_mat_c']] <- survey_age_c %>% arrange(country) %>% 
  #   dplyr::select(-country) %>% as.matrix()
  # 
  # # Birthrate by Country/year
  # mu_ts <- matrix(0, nrow=length(unique(sero_data$survey.num)), ncol=sero_data_stan$age_length)
  # colnames(mu_ts) <- foi_age
  # row.names(mu_ts) <- 1:sero_data_stan$S
  # 
  # un_codes <- sero_data$UNcode[!duplicated(sero_data$survey.num)]
  # 
  # for (s in 1:sero_data_stan$S){
  #   mu_ts[s,] <- get_birthrate(year=as.integer(survey_ageyear_s[s,-(1:2)], na.rm=TRUE), UN.code = un_codes[s])
  # }
  # mu_ts[is.na(mu_ts)] <- 0
  # sero_data_stan[['mu_ts']] <- mu_ts
  # 
  if("wpp2017" %in% (.packages())){ detach("package:wpp2017", unload=TRUE) }
  if("reshape2" %in% (.packages())){ detach("package:reshape2", unload=TRUE) }
  if("plyr" %in% (.packages())){ detach("package:plyr", unload=TRUE) }
  if("fitdistrplus" %in% (.packages())){ detach("package:fitdistrplus", unload=TRUE) }
  if("MASS" %in% (.packages())){ detach("package:MASS", unload=TRUE) }
  library(dplyr)
  
  rm(pop, popF, popM, popFprojMed,popFTproj,popFT,popMT,popMprojMed,popMTproj,tfr,tfrprojMed,
     who_regions,percentASFR,popproj,region_data,iso_data,iso_data_full,dhs_countrydata,nationality_data,UNcode)
  
  return(sero_data_stan) 
}






# PRINT AND SAVE SUMMARY DATA ABOUT SURVEYS -------------------------------

get_country_summary <- function(save_data=FALSE, version=NULL, sero_data=NULL){
  
  if(is.null(sero_data)){
    if (is.null(version)){
      sero_data <- read_csv("data/sero_data_clean.csv") %>% as.data.frame() 
    } else {
      sero_data <- read_csv(file.path( "data", version, "sero_data_clean.csv")) %>% as.data.frame()
    }
  }
  # sero_data$meet_age_crit <- as.logical(sero_data$meet_age_crit)
  # sero_data$meet_age_crit_partial <- as.logical(sero_data$meet_age_crit_partial)
  # 
  survey_tmp_dat <- sero_data[!duplicated(sero_data$survey.ind),]
  tmp <- data.frame(table(survey_tmp_dat$ISO))
  tmp2 <- data.frame(table(survey_tmp_dat$ISO, survey_tmp_dat$meet_age_crit))
  tmp3 <- data.frame(table(survey_tmp_dat$ISO, survey_tmp_dat$meet_age_crit_partial))
  tmp4 <- data.frame(table(survey_tmp_dat$ISO, (survey_tmp_dat$meet_age_crit | survey_tmp_dat$meet_age_crit_partial)))
  
  country_summary <- left_join(left_join(left_join(tmp %>% rename(ISO=Var1, surveys=Freq), 
                                                   tmp2 %>% filter(Var2==TRUE) %>% select(ISO=Var1, meet_crit=Freq)),
                                         (tmp3 %>% filter(Var2==TRUE) %>% select(ISO=Var1, meet_crit_partial=Freq))),
                               (tmp4 %>% filter(Var2==TRUE) %>% select(ISO=Var1, meet_crit_either=Freq))) %>%
    mutate(country=globaltoolboxlite::get_country_name_ISO3(ISO), 
           region =get.region(ISO)) %>% 
    select(country, ISO, region, everything()) %>% arrange(region, country)
  
  country_summary <- left_join(country_summary, 
                               survey_tmp_dat %>% group_by(ISO) %>% summarize(years = toString(start.year)))
  
  
  # merge summary with priority countries
  
  priority_countries <- read.csv(file="data/vimc_countries_of_interest.csv", header=TRUE, stringsAsFactors = FALSE) %>%
    as.data.frame() %>% mutate(req_country=TRUE)
  
  countries_priority <- full_join(priority_countries %>% rename(ISO=iso3) %>% dplyr::select(-country), 
                                  country_summary, by=c("ISO"="ISO")) %>% 
    mutate(region=get.region(ISO)) %>% 
    mutate(country=globaltoolboxlite::get_country_name_ISO3(ISO),
           surveys=ifelse(is.na(surveys), 0, surveys),
           meet_crit=ifelse(is.na(meet_crit), 0, meet_crit),
           meet_crit_partial=ifelse(is.na(meet_crit_partial), 0, meet_crit_partial),
           meet_crit_either=ifelse(is.na(meet_crit_either), 0, meet_crit_either),
           req_country=ifelse(is.na(req_country), FALSE, req_country))
  
  countries_priority <- countries_priority %>% mutate(priority=ifelse(req_country & surveys==0, 1,  
                                                                      ifelse(req_country & meet_crit_either==0, 2,
                                                                             ifelse(req_country & meet_crit==0, 3, 
                                                                                    ifelse(req_country & meet_crit_either==0, 4, 5))))) %>%
    arrange(priority, region, country)
  
  if (save_data){
    dir.create(file.path("results",version))
    write.csv(countries_priority, file=file.path("results",version,"priority_countries_for_serology.csv"), row.names = FALSE)
  }
  
  return(countries_priority)
}






# # PLOT SEROLOGY BY COUNTRY ------------------------------------------------
# 
# plot_country_sero <- function(sero_data, country){
#   
#   p 
#   
# }




