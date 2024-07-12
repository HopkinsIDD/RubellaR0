##' Run Diagnostics and Analysis on Results from Stan R0 model
##' Shaun Truelove, June 2019
##' 


# VERSION AND SPECIFICS ---------------------------------------------------
analysis_version <- "Oct_testing"
analysis_version_repo <- file.path("results", "R0", analysis_version)
source_the_data <- FALSE 
meets_crit <- "partial" # options include: "full", "partial", "none"
min_age_group <- 0.6 # restricts to age groups with >7 months max age 

runs <- 1000


# SETUP AND LOAD RESULTS --------------------------------------------------

# Load packages
library(tidyverse)
library(rstan)

# Load Stan results
load(file=file.path(analysis_version_repo, paste0("R0_sero_full_test_",runs,"iters.RData"))) # Loads ret_sm, fit_R0, sero_data_stan_list

fit_sero_R0 <- fit_R0_test
fits <- rstan::extract(fit_sero_R0)





# SUMMARIZE DATA ----------------------------------------------------------
fit_summary <- summary(fit_sero_R0)
print(fit_summary$summary)

# Sampler Diagnostics
sampler_params <- get_sampler_params(fit_sero_R0, inc_warmup = FALSE)
sampler_params_chain1 <- sampler_params[[1]]
colnames(sampler_params_chain1)

mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
print(mean_accept_stat_by_chain)

max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
print(max_treedepth_by_chain)

# Get warmup and sampling times
print(get_elapsed_time(fit_sero_R0))



summary(fit_sero_R0, pars=c('Ave_age'))$summary
summary(fit_sero_R0, pars=c('Ave_age_post'))$summary
summary(fit_sero_R0, pars=c('R0s'))$summary
summary(fit_sero_R0, pars=c('b0', 'b1', 'b2'))$summary
summary(fit_sero_R0, pars=c('sigma', 'sigma_g', 'log_R0mu', 'log_R0g'))$summary


summary(exp(fits$log_R0g))
exp(quantile(fits$log_R0g, c(0.025, 0.975)))
summary(exp(fits$log_R0c))
exp(quantile(fits$log_R0c, c(0.025, 0.975)))

summary(exp(fits$log_R0c[,1]))
exp(quantile(fits$log_R0c[,1], c(0.025, 0.975)))
sd(fits$log_R0c[,1])

summary((fits$log_R0c[,1]))
(quantile(fits$log_R0c[,1], c(0.025, 0.975)))
sd(fits$log_R0c[,1])

summary(exp(fits$log_R0c[,2]))
exp(quantile(fits$log_R0c[,2], c(0.025, 0.975)))
sd(fits$log_R0c[,2])

summary(exp(fits$log_R0s))
exp(quantile(fits$log_R0s[,1], c(0.025, 0.975)))

summary(fits$log_R0s[,1])
sd(fits$log_R0s[,1])


summary(exp(fits$log_R0s))
exp(quantile(fits$log_R0s[,1], c(0.025, 0.975)))
summary(exp(fits$log_R0c))
exp(quantile(fits$log_R0c, c(0.025, 0.975)))

summary(exp(fits$log_R0r))
exp(quantile(fits$log_R0r, c(0.025, 0.975)))

summary(exp(fits$log_R0g + fits$beta_region[,1]))
exp(quantile(fits$log_R0g + fits$beta_region[,1], c(0.025, 0.975)))
summary(exp(fits$log_R0g + fits$beta_region[,2]))
exp(quantile(fits$log_R0g + fits$beta_region[,2], c(0.025, 0.975)))


summary(fits$R0g)
quantile(fits$R0g, c(0.025, 0.975))




R0_survey_tmp <- exp(fits_R0s$log_R0s)
R0_survey <- data.frame(country = NA,
                        ISO = countryregion_of_survey$ISO[match(sero_data_stan$country_of_survey, countryregion_of_survey$country.num)],
                        region = countryregion_of_survey$sub.region[match(sero_data_stan$country_of_survey, countryregion_of_survey$country.num)],
                        survey = unique((sero_data$survey.ind)), 
                        mean = apply(fits_R0s$R0s, 2, mean),
                        median = apply(fits_R0s$R0s, 2, median),
                        LL = apply(fits_R0s$R0s, 2, quantile, probs=0.025),
                        UL = apply(fits_R0s$R0s, 2, quantile, probs=0.975))
R0_survey$country <- get.country.names.ISO3(R0_survey$ISO)


bpars_survey <- data.frame(country = NA,
                           ISO = countryregion_of_survey$ISO[match(sero_data_stan$country_of_survey, countryregion_of_survey$country.num)],
                           region = countryregion_of_survey$sub.region[match(sero_data_stan$country_of_survey, countryregion_of_survey$country.num)],
                           survey = unique((sero_data$survey.num)), 
                           b0_mean = apply(fits_R0s$b0, 2, mean),
                           b0_LL = apply(fits_R0s$b0, 2, quantile, probs=0.025),
                           b0_UL = apply(fits_R0s$b0, 2, quantile, probs=0.975),
                           b1_mean = apply(fits_R0s$b1, 2, mean),
                           b1_LL = apply(fits_R0s$b1, 2, quantile, probs=0.025),
                           b1_UL = apply(fits_R0s$b1, 2, quantile, probs=0.975),
                           b2_mean = apply(fits_R0s$b2, 2, mean),
                           b2_LL = apply(fits_R0s$b2, 2, quantile, probs=0.025),
                           b2_UL = apply(fits_R0s$b2, 2, quantile, probs=0.975))
bpars_survey$country <- get.country.names.ISO3(bpars_survey$ISO)
bpars_survey$survey_id <- sero_data$country.survey.year[match(bpars_survey$survey, sero_data$survey.num)]



plot_sero_fits <- function(data, params=c(b0=0.15, b1=0.15, b2=0.003)){
  fit.dat <- data.frame(x=0:70, 
                        y=1-(x_a_funct(a=0:70, b0=params["b0"], b1=params["b1"], b2=params["b2"])))# + coef(fit)["mu"]))
  (p1 <- ggplot(data, aes(x=x, y=y)) + geom_point() + 
      geom_line(data=fit.dat, aes(x=x, y=y), color='blue') +
      ggtitle(paste0(data$country, "-", data$survey.num)))
  return(p1)
}

for (i in sort(unique(sero_data$survey.num))){
  dat_survey <- sero_data %>% filter(survey.num==i) %>% select(x=scaled.mid.age, y=seropositive.immune, 
                                                               country=country, survey.num=survey.num) 
  print( plot_sero_fits(data=dat_survey, params=c(b0=bpars_survey$b0_mean[i], 
                                                  b1=bpars_survey$b1_mean[i], b2=bpars_survey$b2_mean[i])))
}  



# PLot multiple on one for saving ------------------------

plot_sero_fits_multi <- function(data, page){
  (p1 <- ggplot(data, aes(x=x, y=y)) + 
     geom_point(data=data%>%filter(geom=='pt'), shape=1) + 
     geom_line(data=data%>%filter(geom=='ln'), color='blue') + 
     facet_wrap(~survey, nrow = 3)) 
  return(p1)
}

par.data <- bpars_survey %>% select(survey.num=survey, survey=survey_id, b0=b0_mean, b1=b1_mean, b2=b2_mean)

fit.data <- NULL
for (i in 1:nrow(par.data)){
  fit.data <- rbind(fit.data, 
                    data.frame(x=0:70, 
                               y = 1-(x_a_funct(a=0:70, b0=par.data$b0[i], b1=par.data$b1[i], b2=par.data$b2[i])),
                               survey.num=par.data$survey.num[i],
                               survey=par.data$survey[i]) )
}
dat_survey <- sero_data %>% select(x=scaled.mid.age, y=seropositive.immune, survey.num=survey.num, survey=country.survey.year) 

data <- rbind(dat_survey %>% mutate(geom="pt"), fit.data %>% mutate(geom='ln'))
plot_sero_fits_multi(data=data %>% filter(survey.num <=9))

pages <- ceiling(nrow(bpars_survey)/9)

dir.create(file.path("results", analysis_version), recursive = TRUE)
pdf(file.path("results", analysis_version, "serology_fits.pdf"))
for(p in 1:pages){
  dat.tmp <- data %>% filter(survey.num>=(((p-1)*9)+1) & survey.num<=(p*9))
  print(plot_sero_fits_multi(dat.tmp))
}
dev.off()




# RUN STAN FOR R0 Fitting ONLY ------------------------------------------

fileName_R0fit_est <- "source/R0_estimation/R0_est_fit_sero_fitR0only_noregion.stan"
R0fit_ret <- stanc(fileName_R0fit_est) # Check Stan file
R0fit_ret_sm <- stan_model(stanc_ret = R0fit_ret) # Compile Stan code
save(fileName_R0fit_est, R0fit_ret, R0fit_ret_sm, file="source/R0_estimation/R0_est_R0fit_compiled.RData")
load(file="source/R0_estimation/R0_est_R0fit_compiled.RData")

sero_data_stan_2 <- sero_data_stan
sero_data_stan_2[["log_R0s"]] <- as.numeric(apply(fits_R0s$log_R0s, 2, mean))

# RUN TEST ---------------------------------------------------------------
(n_survey <-  sero_data_stan[["S"]])
(n_country <-  sero_data_stan[["C"]])
(n_region <-  sero_data_stan[["R"]])

inits_ <- list(log_R0g=1.65,
               log_R0c=rep(2, n_country),
               sigma_c=rep(.5, n_country),
               sigma_g=.5)
inits <- list(inits_, inits_, inits_, inits_)
fit_R0_test <- sampling(R0fit_ret_sm, warmup=100, iter=200, seed=123, data=sero_data_stan_2, 
                        chains=4, control=list(adapt_delta=.99, max_treedepth = 15), 
                        init=inits)
print(get_elapsed_time(fit_R0_test)) # print elapsed time


fit_R0 <- sampling(R0fit_ret_sm, warmup=2000, iter=10000, seed=123, data=sero_data_stan_2, 
                   chains=4, control=list(adapt_delta=.999, max_treedepth = 15), 
                   init=inits)
print(get_elapsed_time(fit_R0)) # print elapsed time

fits <- rstan::extract(fit_R0)




# Save STAN results
dir.create(file.path('results',analysis_version), recursive=TRUE)
save(R0fit_ret_sm, fit_R0, sero_data_stan_2, sero_data, file=file.path('results',analysis_version, 'R0_fits_10000iters.RData'))
load(file.path('results',analysis_version, 'R0_fits_10000iters.RData')) # Loads ret_sm, fit_R0_test, sero_data_stan_list


# check pairs
pairs(fit_R0, pars=c('log_R0g', 'sigma_g', 'lp__'), las=1)
pairs(fit_R0, pars=c('log_R0c[1]', 'sigma_c[1]', 'lp__'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[1]', 'log_R0r[1]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[1]', 'log_R0r[1]', 'sigma_r[1]', 'sigma_c[1]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[1]', 'log_R0r[1]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[1]', 'log_R0c[2]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0r[1]', 'log_R0r[2]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0r[1]', 'sigma_r[1]', 'eta_r[1]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0r[2]', 'sigma_r[2]', 'eta_r[2]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[1]', 'sigma_c[1]', 'eta_c[1]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[2]', 'sigma_c[2]', 'eta_c[2]', 'log_R0r[2]'), las=1)


summary(exp(fits$log_R0g))
exp(quantile(fits$log_R0g, c(0.025, 0.975)))
summary(exp(fits$log_R0c))
exp(quantile(fits$log_R0c, c(0.025, 0.975)))

summary(exp(fits$log_R0c[,1]))
exp(quantile(fits$log_R0c[,1], c(0.025, 0.975)))
summary(exp(fits$log_R0c[,2]))
exp(quantile(fits$log_R0c[,2], c(0.025, 0.975)))

summary(exp(fits$log_R0r))
exp(quantile(fits$log_R0r, c(0.025, 0.975)))

summary(exp(fits$log_R0g + fits$beta_region[,1]))
exp(quantile(fits$log_R0g + fits$beta_region[,1], c(0.025, 0.975)))
summary(exp(fits$log_R0g + fits$beta_region[,2]))
exp(quantile(fits$log_R0g + fits$beta_region[,2], c(0.025, 0.975)))


summary(fits$R0g)
quantile(fits$R0g, c(0.025, 0.975))

apply(fits$R0r, 2, summary)
apply(fits$R0r, 2, quantile, c(0.025, 0.975))

mean(fits$log_R0g) + 3*sd(fits$log_R0g)


R0_survey_tmp <- exp(fits$log_R0s)
R0_survey <- data.frame(country = NA,
                        ISO = countryregion_of_survey$ISO[match(sero_data_stan$country_of_survey, countryregion_of_survey$country.num)],
                        region = countryregion_of_survey$sub.region[match(sero_data_stan$country_of_survey, countryregion_of_survey$country.num)],
                        survey = unique((sero_data$survey.ind)), 
                        mean_log = apply(fits$log_R0s, 2, mean),
                        mean = exp(apply(fits$log_R0s, 2, mean)),
                        median = exp(apply(fits$log_R0s, 2, median)),
                        LL = exp(apply(fits$log_R0s, 2, quantile, probs=0.025)),
                        UL = exp(apply(fits$log_R0s, 2, quantile, probs=0.975)))
R0_survey$country <- get.country.names.ISO3(R0_survey$ISO)



mean(sero_data_stan_2$log_R0s)
mean(sero_data_stan_2$log_R0s) + 3*sd(sero_data_stan_2$log_R0s)
sero_data_stan_2$log_R0s
exp(sero_data_stan_2$log_R0s)



# R0 by country

R0_country_summary <- data.frame(country = unique((sero_data_stan[["country"]])), 
                                 mean = exp(apply(fits$log_R0c, 2, mean)),
                                 median = exp(apply(fits$log_R0c, 2, median)),
                                 LL = exp(apply(fits$log_R0c, 2, quantile, probs=0.025)),
                                 UL = exp(apply(fits$log_R0c, 2, quantile, probs=0.975)),
                                 log_mean = apply(fits$log_R0c, 2, mean),
                                 log_sd = apply(fits$log_R0c, 2, sd))
country_codes <- data.frame(country=unique(sero_data$ISO), country.num=unique(sero_data$country.num))
R0_country_summary$country <- country_codes$country[match(R0_country_summary$country, country_codes$country.num)]
# R0_country_summary <- R0_country_summary %>% mutate(
#   LL3 = exp(log_mean - 1.96*log_sd),
#   UL3 = exp(log_mean + 1.96*log_sd)) 




# R0 by Region
R0_region_summary <- data.frame(region = unique((sero_data_stan[["region_of_country"]])), 
                                mean = exp(apply(fits$log_R0r, 2, mean)),
                                median = exp(apply(fits$log_R0r, 2, median)),
                                LL = exp(apply(fits$log_R0r, 2, quantile, probs=0.025)),
                                UL = exp(apply(fits$log_R0r, 2, quantile, probs=0.975)),
                                log_mean = apply(fits$log_R0r, 2, mean),
                                log_sd = apply(fits$log_R0r, 2, sd))
region_codes <- data.frame(region=unique(sero_data$sub.region), region.num=unique(sero_data$region.num))
R0_region_summary$region <- region_codes$region[match(R0_region_summary$region, region_codes$region.num)]
# R0_region_summary <- R0_region_summary %>% mutate(
#   LL3 = exp(log_mean - 1.96*log_sd),
#   UL3 = exp(log_mean + 1.96*log_sd))          






# TAIWAN AND CHINA TOGETHER -----------------------------------------------

# RUN STAN FOR R0 Fitting ONLY ------------------------------------------

fileName_R0fit_est <- "source/R0_estimation/R0_est_fit_sero_fitR0only_new.stan"
R0fit_ret <- stanc(fileName_R0fit_est) # Check Stan file
R0fit_ret_sm <- stan_model(stanc_ret = R0fit_ret) # Compile Stan code
save(fileName_R0fit_est, R0fit_ret, R0fit_ret_sm, file="source/R0_estimation/R0_est_R0fit_compiled.RData")
load(file="source/R0_estimation/R0_est_R0fit_compiled.RData")


# Get stan-formated data
sero_data_allChina <- get_sero_data(taiwan=FALSE, countries=globaltoolbox::get.iso(c("China","Taiwan")),
                                    meets_crit=meets_crit, adjust_low=TRUE, run_from_saved=FALSE, 
                                    min_age_group=min_age_group, version=analysis_version)
sero_data_Taiwan <- get_sero_data(taiwan=TRUE, countries=globaltoolbox::get.iso(c("Taiwan")),
                                  meets_crit=meets_crit, adjust_low=TRUE, run_from_saved=FALSE, 
                                  min_age_group=min_age_group, version=analysis_version)
unique(sero_data_allChina$country)
unique(sero_data_Taiwan$country)
sero_data <- rbind(sero_data_allChina, sero_data_Taiwan)
sero_data_stan_3 <- format_list_sero_for_stan(sero_data)
sero_data_stan_3[["log_R0s"]] <- as.numeric(apply(fits_R0s$log_R0s, 2, mean))


# RUN TEST ---------------------------------------------------------------
(n_survey <-  sero_data_stan_3[["S"]])
(n_country <-  sero_data_stan_3[["C"]])
(n_region <-  sero_data_stan_3[["R"]])

inits_ <- list(log_R0g=1.65,
               log_R0c=rep(2, n_country),
               # beta_region=rep(0, n_region),
               sigma_c=rep(.5, n_country),
               sigma_g=.5)
inits <- list(inits_, inits_, inits_, inits_)
fit_R0_test <- sampling(R0fit_ret_sm, warmup=100, iter=200, seed=123, data=sero_data_stan_3, 
                        chains=4, control=list(adapt_delta=.99, max_treedepth = 15), 
                        init=inits)
print(get_elapsed_time(fit_R0_test)) # print elapsed time



fit_R0 <- sampling(R0fit_ret_sm, warmup=5000, iter=10000, seed=123, data=sero_data_stan_2, 
                   chains=4, control=list(adapt_delta=.999, max_treedepth = 15), 
                   init=inits)
print(get_elapsed_time(fit_R0)) # print elapsed time

fits <- rstan::extract(fit_R0)




# Save STAN results
dir.create(file.path('results',analysis_version), recursive=TRUE)
save(R0fit_ret_sm, fit_R0_test, sero_data_stan_2, sero_data, file=file.path('results',analysis_version, 'R0_fits_10000iters.RData'))
load(file.path('results',analysis_version, 'R0_fits_10000iters.RData')) # Loads ret_sm, fit_R0_test, sero_data_stan_list


# check pairs
pairs(fit_R0, pars=c('log_R0g', 'sigma_g', 'lp__'), las=1)
pairs(fit_R0, pars=c('log_R0c[1]', 'sigma_c[1]', 'lp__'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[1]', 'log_R0r[1]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[1]', 'log_R0r[1]', 'sigma_r[1]', 'sigma_c[1]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[1]', 'log_R0r[1]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[1]', 'log_R0c[2]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0r[1]', 'log_R0r[2]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0r[1]', 'sigma_r[1]', 'eta_r[1]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0r[2]', 'sigma_r[2]', 'eta_r[2]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[1]', 'sigma_c[1]', 'eta_c[1]'), las=1)
pairs(fit_R0, pars=c('log_R0g', 'log_R0c[2]', 'sigma_c[2]', 'eta_c[2]', 'log_R0r[2]'), las=1)


summary(exp(fits$log_R0g))
exp(quantile(fits$log_R0g, c(0.025, 0.975)))
summary(exp(fits$log_R0c))
exp(quantile(fits$log_R0c, c(0.025, 0.975)))

exp(quantile(fits$log_R0c[,1], c(0.025, 0.975)))
exp(quantile(fits$log_R0c[,2], c(0.025, 0.975)))

summary(exp(fits$log_R0r))
exp(quantile(fits$log_R0r, c(0.025, 0.975)))


summary(fits$R0g)
quantile(fits$R0g, c(0.025, 0.975))

apply(fits$R0r, 2, summary)
apply(fits$R0r, 2, quantile, c(0.025, 0.975))

mean(fits$log_R0g) + 3*sd(fits$log_R0g)























# PUT TOGETHER INTO FINAL DATASET -----------------------------------------

priority_countries <- read.csv(file="data/Requested_VIMC_countries_Rubella.csv", header=TRUE, stringsAsFactors = FALSE) %>%
  as.data.frame() %>% mutate(req_country=TRUE)

priority_countries <- priority_countries %>% 
  rename(ISO=ISO3.code) %>%
  mutate(country=sapply(ISO, get.country.names.ISO3), sub.region=get.subregion(ISO)) %>% 
  mutate(CountryEst=FALSE, RegionalEst=FALSE, GlobalEst=FALSE) %>%
  select(-Country)

priority_countries <- left_join(priority_countries, R0_country_summary, by=c("ISO"="country")) %>% 
  rename(R0_mean=mean, R0_ll=LL, R0_ul=UL, logR0_mean=log_mean, logR0_sd=log_sd) %>% select(-median) %>%
  mutate(CountryEst = ifelse(!is.na(R0_mean), TRUE, FALSE),
         RegionalEst = ifelse(is.na(R0_mean), TRUE, FALSE))

priority_countries <- left_join(priority_countries, R0_region_summary, by=c("sub.region"="region")) %>% 
  mutate(RegionalEst = ifelse(is.na(R0_mean) & !is.na(mean), TRUE, FALSE)) %>%
  mutate(R0_mean=ifelse(is.na(R0_mean), mean, R0_mean), 
         R0_ll=ifelse(is.na(R0_ll), LL, R0_ll),
         R0_ul=ifelse(is.na(R0_ul), UL, R0_ul), 
         logR0_mean=ifelse(is.na(logR0_mean), log_mean, logR0_mean),
         logR0_sd=ifelse(is.na(logR0_sd), log_sd, logR0_sd)) %>% select(-median, -mean, -LL, -UL, -log_mean, -log_sd)

# Global Mean
R0g_mean <- exp(mean(fits$log_R0g))
logR0g_mean <- mean(fits$log_R0g)
logR0g_sd <- sd(fits$log_R0g)
R0g_ll <- exp(quantile(fits$log_R0g, 0.025))
R0g_ul <- exp(quantile(fits$log_R0g, 0.975))

priority_countries <- priority_countries %>% 
  mutate(GlobalEst = ifelse(is.na(R0_mean), TRUE, FALSE)) %>%
  mutate(R0_mean=ifelse(is.na(R0_mean), R0g_mean, R0_mean), 
         R0_ll=ifelse(is.na(R0_ll), R0g_ll, R0_ll),
         R0_ul=ifelse(is.na(R0_ul), R0g_ul, R0_ul), 
         logR0_mean=ifelse(is.na(logR0_mean), logR0g_mean, logR0_mean),
         logR0_sd=ifelse(is.na(logR0_sd), logR0g_sd, logR0_sd))

write_csv(priority_countries, "data/R0_distributions.csv")
