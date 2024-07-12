
#' setup_R0_est_data_stan.R
#' Purpose: Setup and Run Stan Estimation of R0 using serological survey data
#' Author: Shaun Truelove, satruelove@gmail.com
#' Date updated: 30 Aug 2018


# SETUP -------------------------------------------------------------------

source("source/ISO_code_source.R")


# ~ Load packages ----
if("wpp2017" %in% (.packages())){ detach("package:wpp2017", unload=TRUE) }
if("plyr" %in% (.packages())){ detach("package:plyr", unload=TRUE) }
if("fitdistrplus" %in% (.packages())){ detach("package:fitdistrplus", unload=TRUE) }
if("MASS" %in% (.packages())){ detach("package:MASS", unload=TRUE) }
#if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('parallel')) install.packages('parallel'); library(parallel)
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
library(rstan)



# LOAD AND CHECK DATA -----------------------------------------------------

# Get stan-formated data
sero_data <- get_sero_data(meets_crit = meets_crit, 
                           adjust_low = TRUE, 
                           run_from_saved = !source_the_data, 
                           analysis_version = analysis_version, 
                           min_age_group = min_age_group, 
                           new_data = sero_update_files,
                           pull_raw_data = pull_raw_data,
                           data_dir = file.path("data", analysis_version))


sero_data <- sero_data %>% filter(seropositive.immune<=1)

# Now add FOI years to each survey-age
sero_data <- add_foi_years_sero(sero_data)

# Save it as "FINAL" data
write.csv(sero_data, file=file.path("data", analysis_version, "sero_data_final.csv"), row.names=FALSE)
sero_data <- read_csv(file.path("data", analysis_version, "sero_data_final.csv")) %>% as.data.frame()


# ~ Data Checks ----
# --->> REDO INTO SINGLE FUNCTION  ####

# check that N and N.seropositives make sense
p <- sero_data$Y / sero_data$N
hist(p, breaks=100)
sum(p>1) # shoule be 0
# Check for weird ages
sum(sero_data$max_age_years<sero_data$min_age_years) # should be 0
# check for missing data
sum(is.na(sero_data$scaled.mid.age)) # Should be 0


# SETUP DATA --------------------------------------------------------------

# Note studies which do not have at least 2 age groups above 15 years and 2 age groups below
# sero_data <- sero_data %>% 
#   group_by(survey.ind) %>% 
#   mutate(meet_age_crit         = ifelse((sum(max_age_years<15)>=2 & sum(min_age_years>=15)>=2) |
#                                           (sum(max_age_years<15)>=2 & sum(min_age_years<=15 & max_age_years>=15)==1 & sum(min_age_years>=15)>=1), TRUE, FALSE), 
#          meet_age_crit_partial = ifelse(sum(max_age_years<15)>=1 & sum(min_age_years<=15 & max_age_years>=15)==1 & sum(min_age_years>=15)>=1, TRUE, FALSE))

# ~ Get country summary data ------------------------------------------------
country_summary <- suppressWarnings(get_country_summary(save_data=save.data, version=analysis_version, sero_data=sero_data))


# ~ Setup data as list for stan ---------------------------------------------

# Format as list
sero_data_stan <- format_list_sero_for_stan(sero_data)

# Do some more quick checks
mean(unique(sero_data_stan[['country']])) == (max(sero_data_stan[['country']])+1)/2  # Checking that the indexes are complete
mean(unique(sero_data_stan[['survey']])) == (max(sero_data_stan[['survey']])+1)/2
p <- sero_data_stan[['Y']] / sero_data_stan[['N']]
summary(p)
hist(p)
sum(p>1)
rm(p)

summary(unique(sero_data_stan[['mu']]))
summary(unique(sero_data_stan[['age']]))

# ~ Unique Surveys info  ####
surveys_unique <- sero_data[!duplicated(sero_data$survey.ind),]
table(surveys_unique$country2)

# Save needed data
save(sero_data_stan, file=file.path("data", analysis_version,"sero_data_stan.Rdata"))

# # Save workspace image in case of a crash
# save.image('data/R0_stan_sero.Rdata')
# load('data/R0_stan_sero.Rdata')

# Clear it all
rm(list = ls())
gc()
