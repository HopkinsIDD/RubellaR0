
#' load_and_clean_R0_data.R
#' Purpose: Load and clean serology data from google sheets source data
#' Author: Shaun Truelove, satruelove@gmail.com
#' Date updated: 3 June 2019


# This script should be run when new data for serology is added. It pulls the data and cleans and formats it correctly.



# SETUP -------------------------------------------------------------------

# Set Version & Preferences
if (!exists("analysis_version")){
  analysis_version <- "Sept2021"
  analysis_version_repo <- file.path( "data", analysis_version)
  source_the_data <- TRUE 
  meets_crit <- "partial" # options include: "full", "partial", "none"
  min_age_group <- 0.6 # restricts to age groups with >7 months max age 
}
save.data <- TRUE

library(globaltoolboxlite)
# library(globaltoolbox)
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)


# Source Code
source('source/R0_functions.R')
source("source/ISO_code_source.R")
source('source/load_and_clean_R0_data_functions.R')


# Load & Unload packages
# if("wpp2017" %in% (.packages())){ detach("package:wpp2017", unload=TRUE) }
# if("plyr" %in% (.packages())){ detach("package:plyr", unload=TRUE) }
# if("fitdistrplus" %in% (.packages())){ detach("package:fitdistrplus", unload=TRUE) }
# if("MASS" %in% (.packages())){ detach("package:MASS", unload=TRUE) }
# if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
# if(!require('ggforce')) install.packages('ggforce'); library(ggforce)






# LOAD AND CHECK DATA -----------------------------------------------------


# ~ Load and do Basic Cleaning --------------------------------------------
getwd()
dir.create(file.path("data", analysis_version), recursive = TRUE)

sero_data_raw <- load_and_clean_sero_data(taiwan = TRUE, new_data = sero_update_files, data_dir=file.path("data", analysis_version)) # use `taiwan = TRUE` to consider taiwan data for China

write_csv(sero_data_raw, file.path("data",analysis_version,"sero_data_raw.csv"))
sero_data_raw <- read_csv(file.path("data",analysis_version,"sero_data_raw.csv")) %>% as.data.frame()


# ~ Format data for Stan --------------------------------------------------
sero_data_clean <- clean_sero_for_analysis(sero_data=sero_data_raw, save.data=TRUE, version=analysis_version, pull_raw_data=FALSE)
sero_data_clean <- sero_data_clean$sero_data_clean %>% dplyr::select(-further_review_needed, -manual_skip)



# Check criteria
tmp <- sero_data_clean %>% filter(reference_id=="abdullah_1984")

# full criteria
sum(tmp$max_age_years<15)>=2 & sum(tmp$min_age_years>=15)>=2
sum(tmp$max_age_years<15)>=2 & sum(tmp$min_age_years<=15 & tmp$max_age_years>=15)==1 & sum(tmp$min_age_years>=15)>=1
# partial criteria
sum(tmp$max_age_years<15)>=1 & sum(tmp$min_age_years<=15 & tmp$max_age_years>=15)==1 & sum(tmp$min_age_years>=15)>=1


sero_data_clean <- sero_data_clean %>% 
  group_by(survey.ind) %>% 
  mutate(meet_age_crit         = ifelse((sum(max_age_years<15)>=2 & sum(min_age_years>=15)>=2) |
                                          (sum(max_age_years<15)>=2 & sum(min_age_years<=15 & max_age_years>=15)==1 & sum(min_age_years>=15)>=1), TRUE, FALSE), 
         meet_age_crit_partial = ifelse(sum(max_age_years<15)>=1 & sum(min_age_years<=15 & max_age_years>=15)==1 & sum(min_age_years>=15)>=1, TRUE, FALSE))
tmp <- sero_data_clean %>% filter(reference_id=="abdullah_1984")
tmp$meet_age_crit
tmp$meet_age_crit_partial

write.csv(sero_data_clean, file=file.path("data", analysis_version, "sero_data_clean.csv"), row.names=FALSE)


# ~ Take a look at China --------------------------------------------------

# Take a look at China
#colnames(data_tmp)
#View(data_tmp %>% filter(country=="China"))
#rm(data_tmp)


# ~ Get Analysis-Ready Data -----------------------------------------------

sero_data <- get_sero_data(meets_crit=meets_crit, adjust_low=TRUE, run_from_saved=TRUE, min_age_group=min_age_group, analysis_version=analysis_version)


# # Change Taiwan to Taiwan (country)
# sero_data_TaiwanAsChina <- sero_data
# sero_data <- sero_data %>% mutate(country = ifelse(!is.na(city.province.district) & city.province.district=="Taiwan", "Taiwan", country),
#                                   ISO = ifelse(!is.na(city.province.district) & city.province.district=="Taiwan", "TWN", ISO))
  

# DATA CHECKS ----
# check that N and N.seropositives make sense
p <- sero_data$Y / sero_data$N
hist(p, breaks=100)
sum(p>1) # should be 0
# Check for weird ages
sum(sero_data$max_age_years<sero_data$min_age_years) # should be 0
# check for missing data
sum(is.na(sero_data$scaled.mid.age)) # should be 0



# Save it as "FINAL" data
write.csv(sero_data, file=file.path("data", analysis_version, "sero_data_final.csv"), row.names=FALSE)
sero_data <- read_csv(file.path("data", analysis_version, "sero_data_final.csv")) %>% as.data.frame()



# GET COUNTRY SUMMARY DATA ------------------------------------------------
country_summary <- suppressWarnings(get_country_summary(save_data=save.data, version=analysis_version, sero_data=sero_data))



# SETUP DATA AS LIST FOR STAN ---------------------------------------------

sero_data_stan <- format_list_sero_for_stan_AGE(sero_data)
countryregion_of_survey  <- sero_data %>% filter(!duplicated(survey.num)) %>% select(survey.num, ISO, country.num, country, sub.region, region.num) %>% as.data.frame()    # Index id number for region of each survey names(region_of_survey) <- sero_data$sub.region[!duplicated(sero_data$survey.num)]     # Index id number for region of each survey 
region_of_country <- sero_data %>% filter(!duplicated(country.num)) %>% select(country.num, ISO, sub.region, region.num) %>% as.data.frame()
#mod.mat <- model.matrix(~country.num + region.num, sero_data)


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

# Unique Surveys info
surveys_unique <- sero_data[!duplicated(sero_data$survey.ind),]
table(surveys_unique$country2)





# PLOT SEROLOGY BY COUNTRY ------------------------------------------------

library(grid)
library(gridExtra)

countries <- unique(sero_data$country)
length(countries)
i=1


# All plotted together, using facet_wrap
(p_all <- ggplot(sero_data, aes(x=mid.age.years, y=seropositive.immune, colour=country.survey.year)) +
              geom_line() + geom_point() +
              ylim(0,1) +
              theme_classic() +
              theme(legend.position = "none") +
              facet_wrap(~country, nrow=9))
            
# Individual country
(p_indiv <- ggplot(sero_data %>% filter(country %in% countries[i]), 
                aes(x=mid.age.years, y=seropositive.immune, colour=country.survey.year)) +
                  geom_line() + geom_point() +
                  ylim(0,1) +
                  xlab(NULL) + ylab(NULL) +
                  ggtitle(countries[i]) +
                  theme_classic() +
                  theme(plot.title = element_text(margin = margin(t = 10, r=1, b = -20, l=20),
                                                  face="bold",
                                                  hjust = 0.05),
                        legend.position = c(.7, .3), legend.title = element_blank()))

# Groups of countries

n_per_page <- 6 #(3x2)
n_pages <- ceiling(length(countries)/n_per_page)

countries_group <- countries[1:n_per_page]

p_names <- paste0("p", 1:n_per_page)
for (i in 1:n_per_page){
  p_ <- ggplot(sero_data %>% filter(country %in% countries_group[i]), 
                  aes(x=mid.age.years, y=seropositive.immune, colour=country.survey.year)) +
                    geom_line() + geom_point() +
                    ylim(0,1) + coord_cartesian(xlim=c(0,40)) +
                    xlab("Mid Age(yr)") + ylab("Seropositive (%)") +
                    ggtitle(countries_group[i]) +
                    theme_classic() +
                    theme(plot.title = element_text(margin = margin(t=10, r=1, b=-20, l=20), face="bold", hjust=0.05),
                          legend.position = c(.6, .3), legend.title = element_blank())
  assign(p_names[i], p_)
}
                  
gridExtra::grid.arrange(p1 + xlab(NULL),
                        p2 + xlab(NULL) + ylab(NULL),
                        p3 + xlab(NULL),
                        p4 + xlab(NULL) + ylab(NULL),
                        p5,
                        p6 + ylab(NULL), 
                        nrow=3)




# ~ Plot Together and Save ------------------------------------------------

dir.create(file.path("results",analysis_version), recursive = TRUE)

# number per page
n_per_page <- 6 #(3x2)
n_pages <- ceiling(length(countries)/n_per_page)
max_N <- max(sero_data$N)

pdf(file=file.path("results",analysis_version,"serology_raw2.pdf"),
    width = 7.5, height=10, onefile = TRUE)
{
  for (p in 1:n_pages){
    
    countries_group <- countries[(((p-1)*n_per_page)+1) : (p*n_per_page)]
    countries_group <- countries_group[!is.na(countries_group)]
    n_country_ <- sum(!is.na(countries_group))
    n_per_page_ <- min(c(n_per_page, n_country_))
    p_names <- paste0("p", 1:n_per_page_)
    
    for (i in 1:n_per_page_){
      p_ <- sero_data %>% filter(country %in% countries_group[i]) %>%
        ggplot(aes(x=mid.age.years, y=seropositive.immune, colour=country.survey.year)) +
        geom_line() + geom_point() +
        geom_point(aes(size=N), alpha=0.3) + scale_size_continuous(guide="none", limits=c(1,6000), range=c(1, 50)) + #scale_size(range = c(0, 10)) + # point size for size of sample
        geom_errorbarh(aes(xmax=max_age_years, xmin=min_age_years, height=0.01), alpha=.2) + # bar for age range
        ylim(0,1) + coord_cartesian(xlim=c(0,40)) +
        xlab("Mid Age(yr)") + ylab("Seropositive (%)") +
        ggtitle(countries_group[i]) +
        theme_classic() +
        theme(plot.title = element_text(margin = margin(t=10, r=1, b=-20, l=20), face="bold", size=9, hjust=0.05),
              legend.position = c(.7, .25), legend.title = element_blank(), legend.text = element_text(size=6)) +
        guides(colour = guide_legend(
              keyheight = 0.1, default.unit = "inch")) 
      assign(p_names[i], p_)
    }
    if(!exists("p2")) p2 <- ggplot(data.frame(NA))
    if(!exists("p3")) p3 <- ggplot(data.frame(NA))
    if(!exists("p4")) p4 <- ggplot(data.frame(NA))
    if(!exists("p5")) p5 <- ggplot(data.frame(NA))
    if(!exists("p6")) p6 <- ggplot(data.frame(NA))
    
    gridExtra::grid.arrange(p1 + xlab(NULL),
                            p2 + xlab(NULL) + ylab(NULL),
                            p3 + xlab(NULL),
                            p4 + xlab(NULL) + ylab(NULL),
                            p5,
                            p6 + ylab(NULL), 
                            nrow=3)
    rm(p1, p2, p3, p4, p5, p6)
  }
}
dev.off()


#

# 
# 
# p=1
# i=4
# countries_group <- countries[(((p-1)*n_per_page)+1) : (p*n_per_page)]
# p_names <- paste0("p", 1:n_per_page)
# 
# #for (i in 1:n_per_page){
# ggplot(data=sero_data %>% filter(country %in% countries_group[i]), 
#        aes(x=mid.age.years, y=seropositive.immune, colour=country.survey.year)) +
#   geom_line() + geom_point() +
#   geom_point(aes(size=N), alpha=0.3) + scale_size(range = c(0, 10))
#   geom_errorbarh(aes(xmax=max.age.years, xmin=min.age.years, height=0.01)) +
#   ylim(0,1) + coord_cartesian(xlim=c(0,40)) +
#   xlab("Mid Age(yr)") + ylab("Seropositive (%)") +
#   ggtitle(countries_group[i]) +
#   theme_classic() +
#   theme(plot.title = element_text(margin = margin(t=10, r=1, b=-20, l=20), face="bold", size=9, hjust=0.05),
#         legend.position = c(.7, .25), legend.title = element_blank(), legend.text = element_text(size=6)) +
#   guides(colour = guide_legend(
#     keyheight = 0.1, default.unit = "inch"))
# 
# 
# 











