
#" Run_R0_est_stan.R
#" Purpose: Setup and Run Stan Estimation of R0 using serological survey data
#" Author: Shaun Truelove, satruelove@gmail.com
#" Date updated: 30 Aug 2018


# SETUP -------------------------------------------------------------------

# source("R0_model/source/R0_functions.R")
# source("R0_model/source/ISO_code_source.R")
#source("R0_model/source/load_and_clean_R0_data.R")
source('R0_model/source/load_and_clean_R0_data_functions.R')

analysis_version_data <- "Jan2020"
analysis_version <- "Jan2020_farrington_method_LOO"
analysis_version_repo <- file.path("R0_model", "results", analysis_version)
analysis_version_data_repo <- file.path("R0_model", "data", analysis_version_data)
dir.create(analysis_version_data_repo, recursive = TRUE)
dir.create(analysis_version_repo, recursive = TRUE)
countries_of_interest_path <- "R0_model/data/vimc_countries_of_interest_201910gavi_v4.csv"
load_saved_data <- FALSE

# ~ Load packages ----
library(rstan)
library(tidyverse)

# ~ Functions -------
x_a_funct <- function(a, b0, b1, b2){
  1 - exp((b0/b1)*a*exp(-b1*a) + (1/b1)*((b0/b1)-b2)*(exp(-b1*a)-1) - b2*a)
}

# ~ Stan options ----
rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores()); print(parallel::detectCores()) # Detect and set to multi cores
options(mc.cores = 4) # Detect and set to multi cores
options(max.print = .Machine$integer.max)
options(scipen=999)
options(digits=5)
set_cppo(mode="fast")

subset_number <- 20
subset_number <- NULL


# DATA -----------------------------------------------------

# Prepared in 
#source("R0_model/source/setup_R0_est_data_stan.R")

# Load needed data

if (load_saved_data){
  load(file="R0_model/data/R0_est_testing_AGE_DATA.RData") # Loads sero_data_stan, sero_data
} else {
  
  #load(file='R0_model/data/sero_data_stan.Rdata') # Loads 'sero_data_stan'
  #load(file="R0_model/data/sero_data.RData")      # Loads 'sero_data'
  sero_data <- read_csv(file.path("R0_model","data", analysis_version_data, "sero_data_final.csv")) %>% as.data.frame()
  sero_data_name <- file.path("R0_model","data", analysis_version_data, "sero_data_final.csv")
  
  # Remove unneeded columns
  sero_data <- sero_data %>% dplyr::select(-further_review_needed, -check_notes, -dupl_othertest, -manual_skip, -entered_by, -URL, -second_reviewer, 
                                           -Checked, -X30, -X31, -author, -unique_ref_entry_withID)
  
  
  # MODIFYING SERO_DATA_STAN
  
  # to remove: R,   country_ind, age_pred_length, foi_years_range, N_age, age_pred, mu
  # to remove (maybe): region_of_country, region_of_survey, 
  
  if (!is.null(subset_number)){
    # Subset to a subset of surveys to test
    surveys_ <- unique(sero_data$survey.num)
    survey_samp <- sample(surveys_, subset_number)
    sero_data <- sero_data %>% filter(survey.num %in% survey_samp)# country=="Bangladesh")
  }
  
  # Reformat indicators
  sero_data <- reformat_sero_subset(sero_data)
  
  # Add number of full survey
  sero_data <- sero_data %>% group_by(survey.num) %>% mutate(N_full_survey = sum(N)) %>% ungroup()
  
  
  # Add Multiyear birth rate
  # -- function located in get_demog_data.R
  source("R0_model/source/get_demog_data.R")
  survey_summary <- sero_data %>% filter(!duplicated(survey.num))
  survey_summary$mu_10yr <- survey_summary$mu_10yr <- survey_summary$mu_1yr <- numeric(nrow(survey_summary))

  for (i in 1:nrow(survey_summary)){
    survey_summary$mu_1yr[i] <- GetSurveyCBR(year=survey_summary$mid.survey.year[i], n_years_mean=1, UN.code=survey_summary$UNcode[i])
    survey_summary$mu_10yr[i] <- GetSurveyCBR(year=survey_summary$mid.survey.year[i], n_years_mean=10, UN.code=survey_summary$UNcode[i])
    survey_summary$mu_20yr[i] <- GetSurveyCBR(year=survey_summary$mid.survey.year[i], n_years_mean=20, UN.code=survey_summary$UNcode[i])
  }  
  
  #merge back to sero_data
  survey_summary <- survey_summary %>% select(country.survey.year, mu_1yr, mu_10yr, mu_20yr)
  sero_data <- sero_data %>% full_join(survey_summary, by=c("country.survey.year"="country.survey.year"))
  
  
  # # Now add FOI years to each survey-age
  # sero_data <- add_foi_years_sero(sero_data) # Not currently needed, only for age-specific estimates
  
  
  # SETUP DATA AS LIST FOR STAN
  sero_data_stan <- format_list_sero_for_stan_OLD(sero_data)
  #sero_data_stan <- format_list_sero_for_stan_AGE(sero_data)
  rm(pop, popF, popM, popFprojMed,popFTproj,popFT,popMT,popMprojMed,popMTproj,tfr,tfrprojMed,
     who_regions,percentASFR,popproj,region_data,iso_data,iso_data_full,dhs_countrydata,nationality_data,UNcode)
  
  countryregion_of_survey  <- sero_data %>% filter(!duplicated(survey.num)) %>% select(survey.num, ISO, country.num, country, sub.region, region.num) %>% as.data.frame()    # Index id number for region of each survey names(region_of_survey) <- sero_data$sub.region[!duplicated(sero_data$survey.num)]     # Index id number for region of each survey 
  region_of_country <- sero_data %>% filter(!duplicated(country.num)) %>% select(country.num, ISO, sub.region, region.num) %>% as.data.frame()
  #mod.mat <- model.matrix(~country.num + region.num, sero_data)

  
  
  # Add mus to stan list
  sero_data_stan$mu <- (sero_data %>% filter(!duplicated(survey.num)) %>% as.data.frame())$mu_10yr
  
  
  # Save copies of the clean data -- we do this for record keeping
  dir.create(file.path(analysis_version_repo,"data_record"), recursive = TRUE)
  write_csv(sero_data, file.path(analysis_version_repo,"data_record", "sero_data.csv"))
  saveRDS(sero_data_stan, file.path(analysis_version_repo,"data_record", "sero_data_stan.rds"))
  write_csv(countryregion_of_survey, file.path(analysis_version_repo,"data_record", "country_of_survey.csv"))
  write_csv(region_of_country, file.path(analysis_version_repo,"data_record", "region_of_country.csv"))
  
}

# Add mus to stan list --> starting with 10yr mean mu
sero_data_stan$mu <- (sero_data %>% filter(!duplicated(survey.num)) %>% as.data.frame())$mu_10yr

# # Add plus1
# plus1_dat <- read_csv(file.path(analysis_version_repo, "pop_struct_data.csv"))
# plus1_dat <- plus1_dat %>% select(-country.survey.year, -mu_1yr, -mu_10yr, -mu_20yr) %>%
#   full_join(sero_data %>% filter(!duplicated(survey.ind)) %>% select(survey.num, country, mid.survey.year))
# 
# sero_data_stan$plus1 <- plus1_dat$plus1
# # sero_data_stan$plus1 <- rep(0, sero_data_stan$S)
# # sero_data_stan$plus1 <- rep(1, sero_data_stan$S)






# ~ Save workspace image  --------------------------------------------------
# -- in case of a crash
# dir.create(file.path("R0_model","data",analysis_version), recursive=TRUE)
# save.image(file.path("R0_model","data",analysis_version,"R0_stan_sero.Rdata"))
# load(file.path("R0_model","data",analysis_version,"R0_stan_sero.Rdata"))



# RUN CHECKS --------------------------------------------------------------

# Plot serological curves 
source("R0_model/source/plot_serosurvey_curves.R")
plot_serosurvey_curves(sero_data, proj_dir=analysis_version_repo)

# Check Data
# -- check that N and N.seropositives make sense
length(unique(sero_data$survey.num))

p <- sero_data$Y / sero_data$N
hist(p, breaks=100)
sum(p>1)
# Check for weird ages
sum(sero_data$max.age.years<sero_data$min.age.years)
# check for missing data
sum(is.na(sero_data$scaled.mid.age))

# Do some more quick checks
mean(unique(sero_data_stan[["country"]])) == (max(sero_data_stan[["country"]])+1)/2  # Checking that the indexes are complete
mean(unique(sero_data_stan[["survey"]])) == (max(sero_data_stan[["survey"]])+1)/2
p <- sero_data_stan[["Y"]] / sero_data_stan[["N"]]
summary(p)
hist(p)
sum(p>1)
rm(p)

summary(unique(sero_data_stan[["mu"]]))
summary(unique(sero_data_stan[["age"]]))

# Unique Surveys info
surveys_unique <- sero_data[!duplicated(sero_data$survey.ind),]
table(surveys_unique$country2)



# look at birth rate by year and country
pop_struct_data <- (sero_data %>% filter(!duplicated(survey.ind) & !duplicated(paste0(country,mid.survey.year))) %>% 
                      select(country.survey.year, country, mid.survey.year, mu_1yr, mu_10yr, mu_20yr))

write_csv(pop_struct_data, file.path(analysis_version_repo, "pop_struct_data.csv"))


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# MODEL
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# SETUP MODEL -------------------------------------------------------------

# Add mus to stan list --> starting with 10yr mean mu
sero_data_stan$mu <- (sero_data %>% filter(!duplicated(survey.num)) %>% as.data.frame())$mu_10yr

#sero_data_stan$plus1 <- plus1_dat$plus1
# sero_data_stan$plus1 <- rep(0, sero_data_stan$S)
# sero_data_stan$plus1 <- rep(1, sero_data_stan$S)


#fileName_R0_est <- "R0_model/source/R0_est_fit_sero.stan"
fileName_R0_est <- "R0_model/source/R0_est_original_betabinom.stan"
ret <- stanc(fileName_R0_est) # Check Stan file
ret_sm <- stan_model(stanc_ret = stanc(fileName_R0_est)) # Compile Stan code
save(fileName_R0_est, ret, ret_sm, file="R0_model/source/R0_est_original_betabinom_compiled.RData")
load(file="R0_model/source/R0_est_original_betabinom_compiled.RData")



# MODEL: TEST ---------------------------------------------------------------
(n_survey <-  sero_data_stan[["S"]])
(n_country <-  sero_data_stan[["C"]])
(n_region <-  sero_data_stan[["R"]])

inits_ <- list(b0=rep(0.15, n_survey),
               b1=rep(0.15, n_survey),
               b2=rep(0.003, n_survey),
               logR0_g = log(5.0), 
               logR0_c = rep(log(5.0),sero_data_stan$C), 
               logR0_s = rep(log(5.0),sero_data_stan$S),
               #logR0_s_raw = rep(0.1,sero_data_stan$S),
               #logR0_c = rep(log(5.0),sero_data_stan$C),
               #logR0_c_raw = rep(0.1,sero_data_stan$C),
               sigma_g = .5, 
               sigma_c = .5,
               gamma = .05,
               p = rep(0.5, sero_data_stan$SA))
inits <- list(inits_, inits_)#, inits_, inits_)
fit_R0_test <- sampling(ret_sm, warmup=100, iter=200, seed=12345, data=sero_data_stan, 
                        chains=2, control=list(adapt_delta=.85, max_treedepth = 10),
                        diagnostic_file=file.path(analysis_version_repo,"diagnostics.csv"),
                        sample_file=file.path(analysis_version_repo, "sample_file.csv"),
                        init=inits)

print(get_elapsed_time(fit_R0_test)) # print elapsed time
# Save STAN results
#save(ret_sm, fit_R0_test, sero_data_stan, sero_data, file="R0_model/results/R0_full_test_500iters.RData")
save(ret_sm, fit_R0_test, sero_data_stan, sero_data, file=file.path(analysis_version_repo, "R0_full_test_200iters.RData"))
load(file=file.path(analysis_version_repo, "R0_full_test_200iters.RData")) # Loads ret_sm, fit_R0_test, sero_data_stan_list


fits <- rstan::extract(fit_R0_test)

# Check Average Age and R0 Estimates
Ave_age <- fits$Ave_age
colMeans(Ave_age)
exp(colMeans(log(Ave_age)))

R0_s <- fits$R0_s
exp(colMeans(log(R0_s)))

logR0_s <- fits$logR0_s
exp(colMeans(logR0_s))



#calc R0_s from average age
R0_s_tmp = sapply(1:20, FUN=function(x) 1 / ((Ave_age[,x] - .5) * sero_data_stan$mu[x]) + 1)  
exp(colMeans(log(R0_s_tmp)))
exp(colMeans(log(R0_s)))


logR0_c <- fits$logR0_c
exp(colMeans(logR0_c))

sigma_c <- fits$sigma_c
(summary(sigma_c))

logR0_g <- fits$logR0_g
exp(mean(logR0_g))

sigma_g <- fits$sigma_g
(summary(sigma_g))





# ~ Print Fitted Serology ------------------------------------------------

b0 <- apply(fits$b0, 2, median)
b1 <- apply(fits$b1, 2, median)
b2 <- apply(fits$b2, 2, median)
b_fits <- data.frame(survey_num = sero_data_stan$survey[!duplicated(sero_data_stan$survey)], b0=b0, b1=b1, b2=b2 )

source("R0_model/source/plot_serosurvey_curves.R")
plot_serosurvey_curves_withfits(sero_data, proj_dir=analysis_version_repo, b_fits)

# calc A and R0 from these median fits

#calc_ave_age <- funciton(age_pred_length=100, b0, b1, b2, mu)

age_pred_length <- 100
age_pred <- 0:99
M <- 0.5
mu <- sero_data_stan$mu
S <- sero_data_stan$S
N_age <- sero_data_stan$N_age


p_inf_given_age <- n_new_inf_age <- numeric(age_pred_length)
Ave_age_calc <- R0_s_calc <- numeric(S)

for (s in 1:S) {
  for (a in 1:age_pred_length){
    p_inf_given_age[a] = x_a_funct(age_pred[a]+1, b0[s], b1[s], b2[s]) - x_a_funct(age_pred[a], b0[s], b1[s], b2[s]); 
    n_new_inf_age[a] = p_inf_given_age[a] * N_age[s, a];
  }
  tot_new_inf = sum(n_new_inf_age);
  Ave_age_calc[s] = sum((n_new_inf_age / tot_new_inf) * age_pred) # // Average age of infection
  R0_s_calc[s] = 1 + 1 / ((Ave_age_calc[s] ) * mu[s])            #               // R0 
}
Ave_age_calc
R0_s_calc




# MODEL: FULL ---------------------------------------------------------------


# Add mus to stan list --> starting with 10yr mean mu
sero_data_stan$mu <- (sero_data %>% filter(!duplicated(survey.num)) %>% as.data.frame())$mu_10yr

#sero_data_stan$plus1 <- plus1_dat$plus1
sero_data_stan$plus1 <- rep(0, sero_data_stan$S) # this is actually likely most correct as we already adjust for the population structure somewhat.
# sero_data_stan$plus1 <- rep(1, sero_data_stan$S)



inits <- list(inits_, inits_, inits_, inits_)
fit_R0 <- sampling(ret_sm, warmup=1500, iter=3000, seed=12345, data=sero_data_stan, 
                        chains=4, control=list(adapt_delta=.85, max_treedepth = 10),
                        diagnostic_file=file.path(analysis_version_repo,"diagnostics.csv"),
                        sample_file=file.path(analysis_version_repo, "sample_file.csv"),
                        init=inits)
print(get_elapsed_time(fit_R0)) # print elapsed time

# Save STAN results
save(ret, ret_sm, fit_R0, sero_data_stan, sero_data, file=file.path(analysis_version_repo, "R0_full_3000iters.RData"))
load(file=file.path(analysis_version_repo, "R0_full_3000iters.RData")) # Loads ret_sm, fit_R0_test, sero_data_stan_list




# ~ Examine Results ---------------------------------------------------------

fits <- rstan::extract(fit_R0)

# Check Average Age and R0 Estimates
Ave_age <- fits$Ave_age
colMeans(Ave_age)
exp(colMeans(log(Ave_age)))

R0_s <- fits$R0_s
exp(colMeans(log(R0_s)))

logR0_s <- fits$logR0_s
exp(colMeans(logR0_s))




#calc R0_s from average age
R0_s_tmp = sapply(1:sero_data_stan$S, FUN=function(x) 1 / ((Ave_age[,x] - .5) * sero_data_stan$mu[x]) + 1)  
exp(colMeans(log(R0_s_tmp)))
exp(colMeans(log(R0_s)))

R0_s_tmp2 = sapply(1:sero_data_stan$S, FUN=function(x) 1 / ((Ave_age[,x] - .5) * mu2[x]) + 1)  
exp(colMeans(log(R0_s_tmp2)))
exp(colMeans(log(R0_s_tmp)))
exp(colMeans(log(R0_s)))

logR0_c <- fits$logR0_c
exp(colMeans(logR0_c))

sigma_c <- fits$sigma_c
(summary(sigma_c))

logR0_g <- fits$logR0_g
exp(mean(logR0_g))

sigma_g <- fits$sigma_g
(summary(sigma_g))






# ~ Print Fitted Serology ------------------------------------------------

b0 <- apply(fits$b0, 2, median)
b1 <- apply(fits$b1, 2, median)
b2 <- apply(fits$b2, 2, median)
b_fits <- data.frame(survey_num = sero_data_stan$survey[!duplicated(sero_data_stan$survey)], b0=b0, b1=b1, b2=b2 )

source("R0_model/source/plot_serosurvey_curves.R")
plot_serosurvey_curves_withfits(sero_data, proj_dir=analysis_version_repo, b_fits)


# calc A and R0 from these median fits

#calc_ave_age <- funciton(age_pred_length=100, b0, b1, b2, mu)

age_pred_length <- 100
age_pred <- 0:99
M <- 0.5
mu <- sero_data_stan$mu
S <- sero_data_stan$S
N_age <- sero_data_stan$N_age


p_inf_given_age <- n_new_inf_age <- numeric(age_pred_length)
Ave_age_calc <- R0_s_calc <- numeric(S)

for (s in 1:S) {
  for (a in 1:age_pred_length){
    p_inf_given_age[a] = x_a_funct(age_pred[a]+1, b0[s], b1[s], b2[s]) - x_a_funct(age_pred[a], b0[s], b1[s], b2[s]); 
    n_new_inf_age[a] = p_inf_given_age[a] * N_age[s, a];
  }
  tot_new_inf = sum(n_new_inf_age);
  Ave_age_calc[s] = sum((n_new_inf_age / tot_new_inf) * age_pred) # // Average age of infection
  R0_s_calc[s] = 1 + 1 / ((Ave_age_calc[s] ) * mu[s])            #               // R0 
}
Ave_age_calc
R0_s_calc







# ~ Survey Summary --------------------------------------------------------
library(globaltoolbox)
survey_summary <- sero_data %>% filter(!duplicated(survey.num))

R0_survey_tmp <- exp(fits$logR0_s)
R0_survey <- data.frame(country = NA,
                        ISO = countryregion_of_survey$ISO[match(sero_data_stan$country_of_survey, countryregion_of_survey$country.num)],
                        region = countryregion_of_survey$sub.region[match(sero_data_stan$country_of_survey, countryregion_of_survey$country.num)],
                        survey = unique((sero_data$survey.ind)), 
                        N = survey_summary$N_full_survey,
                        mean = exp(apply(fits$logR0_s, 2, mean)),
                        median = exp(apply(fits$logR0_s, 2, median)),
                        LL = exp(apply(fits$logR0_s, 2, quantile, probs=0.025)),
                        UL = exp(apply(fits$logR0_s, 2, quantile, probs=0.975)))
R0_survey$country <- get_country_name_ISO3(as.character(R0_survey$ISO))

AveAge_survey <- data.frame(country = NA, 
                            ISO = countryregion_of_survey$ISO[match(sero_data_stan$country_of_survey, countryregion_of_survey$country.num)],
                            region = countryregion_of_survey$sub.region[match(sero_data_stan$country_of_survey, countryregion_of_survey$country.num)],                            survey = unique((sero_data$survey.ind)), 
                            mean = apply(fits$Ave_age, 2, mean),
                            median = apply(fits$Ave_age, 2, median),
                            LL = apply(fits$Ave_age, 2, quantile, probs=0.025),
                            UL = apply(fits$Ave_age, 2, quantile, probs=0.975))
AveAge_survey$country <- get_country_name_ISO3(AveAge_survey$ISO)


# R0 by country
R0_country <- fits$logR0_c
R0_country_summary <- data.frame(country = unique((sero_data_stan[["country"]])), 
                                mean = exp(apply(R0_country, 2, mean)),
                                median = exp(apply(R0_country, 2, median)),
                                LL = exp(apply(R0_country, 2, quantile, probs=0.025)),
                                UL = exp(apply(R0_country, 2, quantile, probs=0.975)),
                                log_mean = apply(R0_country, 2, mean),
                                log_sd = apply(R0_country, 2, sd))
country_codes <- data.frame(country=unique(sero_data$ISO), country.num=unique(sero_data$country.num))
R0_country_summary$country <- country_codes$country[match(R0_country_summary$country, country_codes$country.num)]



# ~ Get regional R0 -------------------------------------------------------

country_res <- fits$logR0_c
country_ <- R0_country_summary$country
region_ <- get.region(country_)
subregion_ <- get.subregion(country_)

region_info <- data.frame(country=country_, region=region_, subregion=subregion_)

colnames(country_res) <- country_
region_res <- country_res %>% as.data.frame() %>% gather(key="country", value="logR0")
region_res <- region_res %>% full_join(region_info)

region_summ <- region_res %>% group_by(region) %>% 
  summarise(mean_logR0 = mean(logR0),
            low_logR0 = quantile(logR0, probs = 0.025),
            high_logR0 = quantile(logR0, probs= 0.975),
            logR0_sd = sd(logR0)) %>%
  mutate(meanR0 = exp(mean_logR0), 
         lowR0 = exp(low_logR0), 
         highR0 = exp(high_logR0))

subregion_summ <- region_res %>% group_by(subregion) %>% 
  summarise(mean_logR0 = mean(logR0),
            low_logR0 = quantile(logR0, probs = 0.025),
            high_logR0 = quantile(logR0, probs= 0.975),
            logR0_sd = sd(logR0)) %>%
  mutate(meanR0 = exp(mean_logR0), 
         lowR0 = exp(low_logR0), 
         highR0 = exp(high_logR0))



# PUT TOGETHER INTO FINAL DATASET -----------------------------------------
library(globaltoolbox)
priority_countries <- read.csv(file=countries_of_interest_path, header=TRUE, stringsAsFactors = FALSE) %>%
  as.data.frame() %>% mutate(req_country=TRUE) %>% rename(ISO=iso3)

priority_countries <- priority_countries %>% 
  mutate(Country = sapply(ISO, get_country_name_ISO3), sub.region=get.subregion(ISO)) %>% 
  mutate(CountryEst=FALSE, RegionalEst=FALSE, GlobalEst=FALSE) %>% select(-country)

priority_countries <- left_join(priority_countries, R0_country_summary, by=c("ISO"="country")) %>% 
  rename(R0_mean=mean, R0_ll=LL, R0_ul=UL, logR0_mean=log_mean, logR0_sd=log_sd) %>% select(-median) %>%
  mutate(CountryEst = ifelse(!is.na(R0_mean), TRUE, FALSE),
         RegionalEst = ifelse(is.na(R0_mean), TRUE, FALSE))

R0_region_summary <- subregion_summ %>% 
  select(sub.region=subregion, mean=meanR0, LL=lowR0, UL=highR0, log_mean=mean_logR0, log_sd=logR0_sd)
priority_countries <- left_join(priority_countries, R0_region_summary, by=c("sub.region"="sub.region")) %>% 
  mutate(RegionalEst = ifelse(is.na(R0_mean) & !is.na(mean), TRUE, FALSE)) %>%
  mutate(R0_mean=ifelse(is.na(R0_mean), mean, R0_mean), 
         R0_ll=ifelse(is.na(R0_ll), LL, R0_ll),
         R0_ul=ifelse(is.na(R0_ul), UL, R0_ul), 
         logR0_mean=ifelse(is.na(logR0_mean), log_mean, logR0_mean),
         logR0_sd=ifelse(is.na(logR0_sd), log_sd, logR0_sd)) %>% 
  select(-mean, -LL, -UL, -log_mean, -log_sd)

# Global Mean
R0g_mean <- exp(mean(fits$logR0_g ))
logR0g_mean <- mean(fits$logR0_g)
logR0g_sd <- sd(fits$logR0_g) + mean(fits$sigma_g)
# R0g_ll <- exp(quantile(fits$logR0_g, 0.025))
# R0g_ul <- exp(quantile(fits$logR0_g, 0.975))
R0g_ll <- exp(logR0g_mean - 1.96*logR0g_sd)
R0g_ul <- exp(logR0g_mean + 1.96*logR0g_sd)


priority_countries <- priority_countries %>% 
  mutate(GlobalEst = ifelse(is.na(R0_mean), TRUE, FALSE)) %>%
  mutate(R0_mean=ifelse(is.na(R0_mean), R0g_mean, R0_mean), 
         R0_ll=ifelse(is.na(R0_ll), R0g_ll, R0_ll),
         R0_ul=ifelse(is.na(R0_ul), R0g_ul, R0_ul), 
         logR0_mean=ifelse(is.na(logR0_mean), logR0g_mean, logR0_mean),
         logR0_sd=ifelse(is.na(logR0_sd), logR0g_sd, logR0_sd))

write_csv(priority_countries, file.path(analysis_version_repo,"R0_distributions.csv"))




