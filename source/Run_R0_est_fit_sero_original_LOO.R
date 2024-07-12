
#" Purpose: 2 Validation Exercises for R0 estimates
#" Author: Shaun Truelove & Amy Winter
#" Date updated: 16 March 2020

setwd("~/Documents/Github/RubellaVIMC/")
# SETUP -------------------------------------------------------------------

source('R0_model/source/load_and_clean_R0_data_functions.R')

analysis_version_data <- "Jan2020"
analysis_version <- "Jan2020_farrington_method" 
analysis_version_repo <- file.path("R0_model", "results", analysis_version)
analysis_version_data_repo <- file.path("R0_model", "data", analysis_version_data)
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

# Load the data
sero_data <- read_csv(file.path(analysis_version_repo,"data_record", "sero_data.csv"))

# Add WHO region for the sero_data
who_regs <- read_csv("transmission_model/data/country_codes.csv")
sero_data <- left_join(sero_data, who_regs, by=c("ISO"="ISO3_code"))
sero_data <- sero_data %>% mutate(Region_Code = ifelse(country=="Taiwan", "SEARO", Region_Code))



#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# MODEL
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# SETUP MODEL -------------------------------------------------------------


# #fileName_R0_est <- "R0_model/source/R0_est_fit_sero.stan"
fileName_R0_est <- "R0_model/source/R0_est_original_betabinom.stan"
ret <- stanc(fileName_R0_est) # Check Stan file
ret_sm <- stan_model(stanc_ret = stanc(fileName_R0_est)) # Compile Stan code

# Can Load if on Shaun's Computer
#load(file="R0_model/source/R0_est_original_betabinom_compiled.RData")


# MODEL: FULL RESULTS ---------------------------------------------------------------


# Two Leave One Out Analyses

# The first calcuates global R0 leaving out one region at a time (6 total)
fits_loo1 <- NULL
region_options <- unique(sero_data$Region_Code)


# get the age matrix
sero_data_stan <- format_list_sero_for_stan_OLD(sero_data) 
N_age <- sero_data_stan$N_age
survey_regions_ <- sero_data$Region_Code[!duplicated(sero_data$survey.num)]


for (r in 1:length(region_options)){
  
  # Get subset of the data
  sero_data_loo <- sero_data %>% filter(Region_Code!=region_options[r]) ### ADD FILTER
  
  # Turn into Stan Data
  sero_data_stan_loo <- format_list_sero_for_stan_OLD(sero_data_loo, get_age = FALSE) 
  sero_data_stan_loo$N_age <- N_age[!(survey_regions_ %in% region_options[r]),]
  
  # Add mu and plus1 to stan list --> starting with 10yr mean mu
  sero_data_stan_loo$mu <- (sero_data_loo %>% filter(!duplicated(survey.num)) %>% as.data.frame())$mu_10yr
  sero_data_stan_loo$plus1 <- rep(0, sero_data_stan_loo$S) # this is actually likely most correct as we already adjust for the population structure somewhat.
  
  # Re-index everything
  sero_data_stan_loo$survey <- as.integer(as.factor(sero_data_stan_loo$survey))
  sero_data_stan_loo$country_of_survey <- as.integer(as.factor(sero_data_stan_loo$country_of_survey))
  sero_data_stan_loo$country <- as.integer(as.factor(sero_data_stan_loo$country))
  
  
  (n_survey <-  sero_data_stan_loo[["S"]])
  (n_country <-  sero_data_stan_loo[["C"]])
  (n_region <-  sero_data_stan_loo[["R"]])
  inits_ <- list(b0=rep(0.15, n_survey),
                 b1=rep(0.15, n_survey),
                 b2=rep(0.003, n_survey),
                 logR0_g = log(5.0), 
                 logR0_c = rep(log(5.0),sero_data_stan_loo$C), 
                 logR0_s = rep(log(5.0),sero_data_stan_loo$S),
                 sigma_g = .5, 
                 sigma_c = .5,
                 gamma = .05,
                 p = rep(0.5, sero_data_stan_loo$SA))
  inits <- list(inits_, inits_, inits_, inits_)
  
  fit_R0 <- sampling(ret_sm, warmup=100, iter=200, seed=12345, data=sero_data_stan_loo, 
                     chains=4, control=list(adapt_delta=.85, max_treedepth = 10),
                     diagnostic_file=file.path(analysis_version_repo,"diagnostics.csv"),
                     sample_file=file.path(analysis_version_repo, "sample_file.csv"),
                     init=inits)
  
  
  # ~ Examine Results ---------------------------------------------------------
  
  fits <- rstan::extract(fit_R0)
  logR0_g <- fits$logR0_g
  sigma_g <- fits$sigma_g
  
  fits_loo1 <- bind_rows(fits_loo1, data.frame(R0g_mean=exp(mean(logR0_g)), 
                                               R0g_sd=mean(sigma_g),
                                               R0g_ub=exp(quantile(logR0_g, 0.975)),
                                               R0g_lb=exp(quantile(logR0_g, 0.025))))
  
}
# Save the output
save(fits_loo1, file="R0_model/results/Jan2020_farrington_method_LOO/LOORegionAnalysisResults.RData")



# The second calcuates global R0 leaving out 1/5 of the serosurveys at a time, 100 times
# Unique Surveys 
unique_surveys <- unique(sero_data$survey.num) 
number_unique_surveys <- length(unique_surveys)

# get the age matrix
sero_data_stan <- format_list_sero_for_stan_OLD(sero_data) 
N_age <- sero_data_stan$N_age

fits_loo2 <- NULL

for (r in 1:100){
  
  # Get random draw of serosurveys
  random_surveys <- sample(1:number_unique_surveys, ceiling(number_unique_surveys*(4/5)), replace=F) 
  
  # Get subset of the data
  sero_data_loo <- sero_data %>% filter(survey.num %in% unique_surveys[random_surveys]) ### ADD FILTER - xxshaun - is survey.ind the right variable? are thene 120 unique serosurveys
  
  # Turn into Stan Data
  sero_data_stan_loo <- format_list_sero_for_stan_OLD(sero_data_loo, get_age = FALSE) 
  sero_data_stan_loo$N_age <- N_age[random_surveys,]
  
  # Add mu and plus1 to stan list --> starting with 10yr mean mu
  sero_data_stan_loo$mu <- (sero_data_loo %>% filter(!duplicated(survey.num)) %>% as.data.frame())$mu_10yr
  sero_data_stan_loo$plus1 <- rep(0, sero_data_stan_loo$S) # this is actually likely most correct as we already adjust for the population structure somewhat.
  
  # Re-index everything
  sero_data_stan_loo$survey <- as.integer(as.factor(sero_data_stan_loo$survey))
  sero_data_stan_loo$country_of_survey <- as.integer(as.factor(sero_data_stan_loo$country_of_survey))
  sero_data_stan_loo$country <- as.integer(as.factor(sero_data_stan_loo$country))
  
  
  (n_survey <-  sero_data_stan_loo[["S"]])
  (n_country <-  sero_data_stan_loo[["C"]])
  (n_region <-  sero_data_stan_loo[["R"]])
  inits_ <- list(b0=rep(0.15, n_survey),
                 b1=rep(0.15, n_survey),
                 b2=rep(0.003, n_survey),
                 logR0_g = log(5.0), 
                 logR0_c = rep(log(5.0),sero_data_stan_loo$C), 
                 logR0_s = rep(log(5.0),sero_data_stan_loo$S),
                 sigma_g = .5, 
                 sigma_c = .5,
                 gamma = .05,
                 p = rep(0.5, sero_data_stan_loo$SA))
  inits <- list(inits_, inits_, inits_, inits_)
  
  fit_R0 <- sampling(ret_sm, warmup=100, iter=200, seed=12345, data=sero_data_stan_loo, 
                     chains=4, control=list(adapt_delta=.85, max_treedepth = 10),
                     diagnostic_file=file.path(analysis_version_repo,"diagnostics.csv"),
                     sample_file=file.path(analysis_version_repo, "sample_file.csv"),
                     init=inits)
  
  # ~ Examine Results ---------------------------------------------------------
  
  fits <- rstan::extract(fit_R0)
  logR0_g <- fits$logR0_g
  sigma_g <- fits$sigma_g
  
  fits_loo2 <- bind_rows(fits_loo2, data.frame(R0g_mean=exp(mean(logR0_g)), 
                                               R0g_sd=mean(sigma_g),
                                               R0g_ub=exp(quantile(logR0_g, 0.975)),
                                               R0g_lb=exp(quantile(logR0_g, 0.025))))
}

# Save the output
save(fits_loo2, file="R0_model/results/Jan2020_farrington_method_LOO/LOOSurveyAnalysisResults.RData")



#### PLOT ####
# Get final results
mean0 <- 5.87
lb0 <- 5.09
ub0 <- 6.74
sd0 <- mean(c((lb0-mean0)/(-1.96),(ub0-mean0)/1.96))

# Region LOO
mean.reg <- c(mean0, fits_loo1$R0g_mean)
lb.reg <- c(lb0,  fits_loo1$R0g_lb)
ub.reg <- c(ub0,  fits_loo1$R0g_ub)
#pdf("~/Google Drive/Rubella.Global.Estimates/VIMC/figs/201910gavi_v4/model_documentation/LOOregion.pdf", width=8, height=5.5)
plot(1:7, mean.reg, bty="n",
     xaxt="n",
     xlab = "Region Left Out", 
     ylab = "global R0",
     main = "Leave Out Region Analysis",
     col = "black",
     ylim = c(min(lb.reg), max(ub.reg)))
axis(1, 1:7, c("none",region_options))
arrows(1:7, lb.reg, 
       1:7, ub.reg, lwd = 1.5, 
       angle = 90, code = 3, length = 0.05) 
#dev.off()

# 1/5 Survey LOO
range(fits_loo2$R0g_mean-mean0)
range(fits_loo2$R0g_sd-sd0)
sum(abs(fits_loo2$R0g_mean-mean0)<sd0)/length(fits_loo2$R0g_mean)
#94% of the 100 analyeses estimated a mean global R0s within one standard deviation of the true global R0
fits_loo2$R0g_mean[which(abs(fits_loo2$R0g_mean-mean0)>=sd0)]
sum(abs(fits_loo2$R0g_mean-mean0)<2*sd0)/length(fits_loo2$R0g_mean)



