### TRANSFORM PREM AGE CONTACT TO RATIOS FOR ADJUSTING STAN R0




# SETUP -------------------------------------------------------------------

library(tidyverse)
source("source/get.prem.WAIFW.R")
source('source/get_demog_data.R') 

# analysis_version <- analysis_version_data <- "May2024"#"Sept2021"


# FUNCTIONS ---------------------------------------------------------------





# LOAD DATA ---------------------------------------------------------------

# load country/year data for each survey
countries_priority <- read.csv(file=file.path("results",analysis_version,"priority_countries_for_serology.csv"))

#load sero_data
sero_data_tmp <- read_csv(file.path("data", analysis_version_data, "sero_data_final.csv")) %>% as.data.frame()
country_year <- sero_data_tmp %>% 
    filter(meet_age_crit | meet_age_crit_partial) %>% 
    select(ISO, UNcode, year=start.year) %>%
    distinct() %>%
    arrange(ISO, year)
rm(sero_data_tmp)
# load Prem
countries_ <- unique(country_year$ISO)
country_prem <- list()
for (c in 1:length(countries_)){
    prem_mat <- tryCatch(get.prem2021.WAIFW.RAW(iso3code = countries_[c]), error=function(e) NA)
    country_prem[[c]] <- prem_mat
}
names(country_prem) <- countries_


# load demography

max_year <- lubridate::year(Sys.Date())

pop_1980 <- pop_2020 <- pop_survey <- list()
for (s in 1:nrow(country_year)){
    
    UNcode_ <- country_year$UNcode[s]
    year_   <- country_year$year[s]

    pop_ <- getDemography5yr(uncode = UNcode_)$pop.age.byageclasses.1950.2100 %>% as.data.frame() 
    pop_ <- pop_ %>%
        mutate(age = rownames(pop_)) %>%
        dplyr::select(age, everything())
    
    pop_ <- pop_ %>% select(age, as.character(1950:max_year))
    
    pop_1980[[s]] <- pop_ %>% select(age, pop = as.character(1980)) %>% 
        mutate(pop = ifelse(is.na(pop), 0, pop)) %>%
        mutate(age_l = seq(0,100, 5),
               age_l = ifelse(age_l>=80, 80, age_l),
               age = ifelse(age_l>=80, "80-100", age_l)) %>%
        group_by(age_l) %>%
        summarise(pop = sum(pop, na.rm = TRUE)) 
    pop_2020[[s]] <- pop_ %>% select(age, pop = as.character(2020)) %>% 
        mutate(pop = ifelse(is.na(pop), 0, pop)) %>%
        mutate(age_l = seq(0,100, 5),
               age_l = ifelse(age_l>=80, 80, age_l),
               age = ifelse(age_l>=80, "80-100", age_l)) %>%
        group_by(age_l) %>%
        summarise(pop = sum(pop, na.rm = TRUE))     
    pop_survey[[s]] <- pop_ %>% select(age, pop = as.character(year_)) %>% 
        mutate(pop = ifelse(is.na(pop), 0, pop)) %>%
        mutate(age_l = seq(0,100, 5),
               age_l = ifelse(age_l>=80, 80, age_l),
               age = ifelse(age_l>=80, "80-100", age_l)) %>%
        group_by(age_l) %>%
        summarise(pop = sum(pop, na.rm = TRUE))    
}
names(pop_1980) <- names(pop_2020) <- names(pop_survey) <- country_year$ISO




# Get estimates
get_contactR0 <- function(country = country_year$ISO[1],
                          pop_age = pop_1980,
                          country_prem){
    
    age.structure.vector <- pop_age[[country]] %>% filter(age_l<80) %>% pull(pop) %>% as.numeric()
    denom.scalar <- sum(age.structure.vector, na.rm = TRUE)
    prem.matrix <- country_prem[[country]]
    
    next.gen <- age.structure.vector*(1-exp(-prem.matrix/denom.scalar))
    dominanteigenvalue <- Re(eigen(next.gen)$value[1])
    
    return(dominanteigenvalue)
}

# Get estimates
get_contactR0_2 <- function(country = country_year$ISO[1],
                          pop_age1 = pop_1980,
                          pop_age2 = pop_2020,
                          country_prem){
    
    age_structure1 <- pop_age1[[country]] %>% filter(age_l<80) %>% pull(pop) %>% as.numeric()
    age_structure2 <- pop_age2[[country]] %>% filter(age_l<80) %>% pull(pop) %>% as.numeric()
    
    denom.scalar <- sum(age_structure1, na.rm = TRUE)
    prem.matrix <- country_prem[[country]]
    
    #adj_contact_matrix <- prem.matrix * (age_structure1 / age_structure2)
    adj_contact_matrix <- prem.matrix * ((age_structure1*(sum(age_structure2)/sum(age_structure1))) / age_structure2)
    
    next.gen <- age_structure1*(1-exp(-adj_contact_matrix/denom.scalar))
    dominanteigenvalue <- Re(eigen(next.gen)$value[1])
    
    return(dominanteigenvalue)
}
    

contactR0_1980 <- contactR0_2020 <- contactR0_survey <- as.numeric(rep(NA, nrow(country_year)))
for (y in 1:nrow(country_year)){
    # contactR0_1980[y]   <- tryCatch(get_contactR0(country = country_year$ISO[y], pop_age = pop_1980, country_prem), error=function(e) NA)
    contactR0_2020[y]   <- tryCatch(get_contactR0(country = country_year$ISO[y], pop_age = pop_2020, country_prem), error=function(e) NA)
    # contactR0_survey[y] <- tryCatch(get_contactR0(country = country_year$ISO[y], pop_age = pop_survey, country_prem), error=function(e) NA)
    
    contactR0_1980[y]   <- tryCatch(get_contactR0_2(country = country_year$ISO[y], pop_age1 = pop_1980, pop_age2 = pop_2020, country_prem), error=function(e) NA)
    contactR0_survey[y] <- tryCatch(get_contactR0_2(country = country_year$ISO[y], pop_age1 = pop_survey, pop_age2 = pop_2020, country_prem), error=function(e) NA)
}

r_ratio_1980_survey <- contactR0_1980 / contactR0_survey
r_ratio_1980_survey
country_year$year


r_ratio_1980_survey[is.na(r_ratio_1980_survey)] <- 1
country_year$r_ratio_1980 <- r_ratio_1980_survey

write.csv(country_year, file=file.path("data", analysis_version_data, "r_ratio_1980.csv"))
