

# SETUP -------------------------------------------------------------------

library(tidyverse)

analysis_version <- "Nov_testing_Age"
analysis_version_repo <- file.path("R0_model", "results", analysis_version)
dir.create(file.path("R0_model","results",analysis_version), recursive = TRUE)


# Load Results
load(file=file.path(analysis_version_repo, 'R0_test_500iters.RData')) # Loads ret_sm, fit_R0_test, sero_data_stan




# GET RESULTS -------------------------------------------------------------

fits <- rstan::extract(fit_R0_test)


# summary(fit_R0_test)[[1]]
# summary(fit_R0_test, pars=c('p', 'lambda_ts', 'logR0_c', 'epsilon_s', 'logR0_g', 'sigma_g'))[[1]]
# summary(fit_R0_test, pars=c('logR0_ts'))[[1]]
# summary(fit_R0_test, pars=c('beta_tc'))[[1]]
# summary(fit_R0_test, pars=c('epsilon_s'))[[1]]


# Look at probabilities vs expected
res <- data.frame(SA = sero_data_stan$survey_age,
                  Age = sero_data_stan$age,
                  Year = sero_data$start.year,
                  survey_details = sero_data$country.survey.year,
                  S = sero_data_stan$survey,
                  C = sero_data_stan$country,
                  Country = sero_data_stan$country_name,
                  Y = sero_data_stan$Y,
                  N = sero_data_stan$N,
                  Prior = sero_data_stan$Y / sero_data_stan$N,
                  Posterior = colMeans(fits$p),
                  Post_low = apply(fits$p, 2, quantile, 0.025),
                  Post_high = apply(fits$p, 2, quantile, 0.975))
tmp <- binom::binom.confint(res$Y, res$N, method="wilson")
res <- res %>% mutate(binom_low=round(tmp$lower,3),
                      binom_high=round(tmp$upper,3))
View(res)


summary(fit_R0_test, pars=c('logR0_g', 'sigma_g'))[[1]]
summary(fit_R0_test, pars=c('logR0_c'))[[1]]

# R0
exp(mean(fits$logR0_g)); exp(quantile(fits$logR0_g, c(0.025, 0.975)))

R0c <- exp(apply(fits$logR0_c, 2, mean))
hist(R0c, breaks=10)

# Overdispersion
mean(fits$gamma); quantile(fits$gamma, c(0.025, 0.975))







# PLOT P FITS -------------------------------------------------------------

plot_probs <- function(data, prior_limits=FALSE, age_limits=c(0,30)){
   
    if (length(unique(data$C))>1){
        return(print("Must limit to 1 country"))
    }
    
    data <- data %>% mutate(survey_details = gsub(paste0(Country[1]," - "), "", survey_details))
    data <- data %>% rename(Survey=survey_details)
    
    data_ <- rbind(data %>% dplyr::select(S, Survey, Age, N, p=Prior, p_low=binom_low, p_high=binom_high) %>% mutate(type="Prior"),
                   data %>% filter(!duplicated(Age)) %>% dplyr::select(S, Survey, Age, N, p=Posterior, p_low=Post_low, p_high=Post_high) %>% 
                             mutate(S="Posterior", Survey="Posterior", type="Posterior")) %>% as.data.frame()
    data_ <- data_ %>% mutate(N = ifelse(type=="Posterior", mean(data_$N), N),
                              Survey = as.factor(Survey)) %>% 
                        mutate(Survey = relevel(Survey, ref="Posterior"))
    
    if (!prior_limits){
        data_ <- data_ %>% mutate(p_low = ifelse(type!="Posterior", NA, p_low),
                                 p_high = ifelse(type!="Posterior", NA, p_high))
    }
    pd <- position_dodge(1) # move them .5 to the left and right
    p <- ggplot(data_, aes(Age, p, group=Survey, color=Survey),  na.rm=TRUE) + geom_point(aes(size=N, shape=type)) + geom_line(size=1) +
         geom_errorbar(aes(ymin=p_low, ymax=p_high), width=2, position = pd,  na.rm=TRUE) +
         coord_cartesian(ylim=c(0,1), xlim=age_limits) +
         theme_bw() + guides(size=FALSE, shape=FALSE) +
         theme(legend.title = element_blank(),
               plot.title = element_text(vjust=-10),
               legend.position = c(.75,.1),
               legend.justification = c("right", "bottom"),
               legend.text = element_text(size = 8)) +
         ggtitle(paste0("  ", data$Country[1]))
    return(p)
}

plot_multi <- function(data, X){
    data_ <- data %>% filter(C==X)
    plot_probs(data_, prior_limits=FALSE, age_limits=c(0,40)) 
}


# SAVE SEROLOGY FITS ------------------------------------------------------



plot_multi(data=res, X=4)

plots_per_page <- 4
plot_rows <- 2
pages <- ceiling(length(unique(res$C)) / plots_per_page)
pdf(file=file.path(analysis_version_repo, "seropositivity_byage.pdf"), width = 8.5, height=10, onefile = TRUE)

for(p in 1:pages){
    dat.tmp <- res %>% filter(C>=(((p-1)*plots_per_page)+1) & C<=(p*plots_per_page))
    myGrobs <- lapply(unique(dat.tmp$C), plot_multi, data=dat.tmp)
    gridExtra::grid.arrange(grobs = myGrobs, nrow = plot_rows)
}

dev.off()






# DIVERGENCES -------------------------------------------------------------

#install.packages("shinystan")
library("shinystan")

params_cp <- as.data.frame(extract(fit_R0_test, permuted=FALSE))
names(params_cp) <- gsub("chain:1.", "", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("[", ".", names(params_cp), fixed = TRUE)
names(params_cp) <- gsub("]", "", names(params_cp), fixed = TRUE)
params_cp$iter <- 1:nrow(params_cp)

par(mar = c(4, 4, 0.5, 0.5))
plot(params_cp$iter, exp(params_cp$logR0_g), pch=16, cex=0.8,
     xlab="Iteration", ylab="R0g")




#
















## Looking at first survey, 5th age group

s = 1
sa = 5

sero_data_stan$age
sero_data_stan$survey
sero_data_stan$mu
sero_data_stan$foi_year_mat_sa
#  ---> looks like the matrices are messed up!!!!







