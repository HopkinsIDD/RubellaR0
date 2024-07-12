


# LOAD AND SETUP DATA AS LIST FOR STAN ---------------------------------------------
{
  sero_data_clean <- read.csv('data/sero_data_clean.csv', header=TRUE, stringsAsFactors=FALSE)
  sero_data_clean <- sero_data_clean %>% filter(!is.na(N)) # make sure no surveys with missing data
  sero_data_clean <- sero_data_clean %>% filter(N > 0) # remove age groups with 0 data
  
  # Remove any before 7 months of age
  View(sero_data_clean %>% filter(max.age.years<.6))
  sero_data_clean <- sero_data_clean %>% filter(max.age.years>=.6) # 6.99/12 = 0.58
  
  # Remove studies which do not have at least 2 age groups above 15 years and 2 age groups below
  sero_data_clean <- sero_data_clean %>% group_by(survey.ind) %>% 
    mutate(meet_age_crit=ifelse(sum(max.age.years<15)>=2 & sum(min.age.years>=15)>=2, TRUE, FALSE), 
           meet_age_crit_partial=ifelse(sum(max.age.years<15)>=1 & sum(min.age.years<=15 & max.age.years>=15)==1 & sum(min.age.years>=15)>=1, TRUE, FALSE))
  
  # subtract 6 months from age (maternal immunity)
  sero_data_clean <- sero_data_clean %>% mutate(scaled.mid.age = ifelse(max.age.years>=1, scaled.mid.age-.5, (max.age.years - .5)/2))
  
  # Rework survey and country numbers for Stan
  sero_data_clean$survey.num <- as.integer(as.factor(sero_data_clean$survey.ind))              # Index id number for each survey 
  sero_data_clean$country.num <- as.integer(as.factor(sero_data_clean$country.ind))            # Index id number for country for each survey
  
  # set up stan list data
  sero_data_stan <- list()
  sero_data_stan[["C"]] <- as.integer(length(unique(sero_data_clean$ISO)))          # Number unique countries
  sero_data_stan[["S"]] <- as.integer(length(unique(sero_data_clean$survey.ind)))   # Number unique samples/serosurveys
  sero_data_stan[["SA"]] <- as.integer(nrow(sero_data_clean))                       # Number unique samples/serosurveys x age groups
  sero_data_stan[["N"]] <- as.integer(sero_data_clean$N)                            # Number indivs in each age grp in each serosurvey/sample
  sero_data_stan[["Y"]] <- as.integer(sero_data_clean$N.seropositive)               # Number seropositive in each age in each serosurvey/sample
  sero_data_stan[["age"]] <- as.numeric(sero_data_clean$scaled.mid.age)             # Mid age value of each age group
  sero_data_stan[["mu"]] <- as.numeric(sero_data_clean$birth.rate)                  # Birth rate for that country during the year of the survey
  sero_data_stan[["country_name"]] <- sero_data_clean$country2            # Index id number for country for each survey
  sero_data_stan[["country_ind"]] <- sero_data_clean$country.ind            # Index id number for country for each survey
  sero_data_stan[["country"]] <- sero_data_clean$country.num            # Index id number for country for each survey
  sero_data_stan[["survey"]] <- sero_data_clean$survey.num              # Index id number for each survey 
  # Get which country index belongs to each survey
  sero_data_stan[["country_of_survey"]] <- sero_data_stan[["country"]][!duplicated(sero_data_stan[["survey"]])] 
}

# Put together into df
lengths(sero_data_stan)
res_dat <- data.frame(ISO=sero_data_clean$ISO,   
                      year=(sero_data_clean$start.year + sero_data_clean$end.year) / 2,
                      do.call(cbind.data.frame, 
                         sero_data_stan[c("N", "Y", "age", "mu", "country_name", "country_ind", "country", "survey")]),
                      min.age=sero_data_clean$min.age.years, 
                      max.age=sero_data_clean$max.age.years, 
                      age.span=sero_data_clean$max.age.years - sero_data_clean$min.age.years,
                      meets_age_crit=sero_data_clean$meet_age_crit,
                      meets_age_crit_partial=sero_data_clean$meet_age_crit_partial,
                      meets_age_crit_any=(sero_data_clean$meet_age_crit | sero_data_clean$meet_age_crit_partial))


# Remove some of positives below age 6 mo
#  - to attempt to fix the problem of maternal antibodies, we will try to just remove those children
n <- with(res_dat, ifelse(min.age<0.5, round(N*((.5-min.age)/age.span)), 0))
res_dat <- res_dat %>% mutate(N2 = N - n, Y2 = Y - n)
res_dat$Y2[res_dat$Y2<0] <- 0






# DO CRUDE ANALYSIS OF ESTIMATED FOI AND R0 -------------------------------

# probability of seropositivity
res_dat <- res_dat %>% mutate(p = Y / N, p2 = Y2 / N2) %>% 
                       mutate(p = ifelse(p==1, .99, p), p2 = ifelse(p2==1, .99, p2)) %>%
                       mutate(lambda = -1*(log(1-p) / age), lambda2 = -1*(log(1-p2) / age))

# # Calculate FOI using Cutts method
# lambda = -ln((1 - p(a1)) / ( 1 - p (a))

res_dat <- res_dat %>% mutate(p = Y / N, p2 = Y2 / N2) %>% 
  mutate(p = ifelse(p==1, .99, p), p2 = ifelse(p2==1, .99, p2)) %>%
  mutate(lambda = -1*(log(1-p) / age), lambda2 = -1*(log(1-p2) / age))


res_dat <- res_dat %>% group_by(survey) %>% mutate(survey_tot_age = max(max.age) - min(min.age)) %>%
                                            mutate(lambda_mean=mean(lambda, na.rm=T), lambda_mean2=mean(lambda2, na.rm=T),
                                                   lambda_wt_mean=sum(lambda*((max.age-min.age)/survey_tot_age), na.rm=T),
                                                   lambda_wt_mean2=sum(lambda2*((max.age-min.age)/survey_tot_age), na.rm=T))

res_dat <- res_dat %>% mutate(R0_indiv = (lambda / mu) + 1, 
                              R0_mean  = (lambda_mean / mu) + 1, 
                              R0_wt_mean = (lambda_wt_mean / mu) + 1,
                              R0_indiv2 = (lambda2 / mu) + 1, 
                              R0_mean2  = (lambda_mean2 / mu) + 1, 
                              R0_wt_mean2 = (lambda_wt_mean2 / mu) + 1)

res_dat_mean <- res_dat %>% filter(!duplicated(survey))


# Look at the results

ggplot(res_dat_mean, aes(meets_age_crit, R0_wt_mean)) + geom_boxplot()
ggplot(res_dat_mean, aes(meets_age_crit, R0_mean)) + geom_boxplot()

ggplot(res_dat_mean, aes(meets_age_crit, R0_wt_mean2)) + geom_boxplot()
ggplot(res_dat_mean, aes(meets_age_crit, R0_mean2)) + geom_boxplot()
by(res_dat_mean$R0_wt_mean2, res_dat_mean$meets_age_crit, summary)
by(res_dat_mean$R0_wt_mean2, res_dat_mean$meets_age_crit_any, summary)

ggplot(res_dat, aes(meets_age_crit, R0_indiv)) + geom_boxplot()
ggplot(res_dat, aes(meets_age_crit, R0_indiv2)) + geom_boxplot()


ggplot(res_dat, aes(color=meets_age_crit, x=log(R0_indiv))) + geom_density()
ggplot(res_dat, aes(color=meets_age_crit, x=log(R0_indiv2))) + geom_density()

median(res_dat$R0_indiv, na.rm=TRUE)

ggplot(res_dat %>% filter(lambda2<5), aes(color=meets_age_crit, x=age, y=lambda)) + geom_point()
ggplot(res_dat %>% filter(lambda2<5), aes(color=meets_age_crit, x=age, y=lambda2)) + geom_point()

ggplot(res_dat %>% filter(lambda2<5), aes(color=survey, x=age, y=p2, group=survey)) + 
  geom_point() + geom_line() +
  facet_grid(rows=vars(meets_age_crit))

# India
ggplot(res_dat %>% filter(lambda2<5 & country_name=="India"), aes(color=survey, x=age, y=p2, group=survey)) + 
  geom_point() + geom_line() + 
  facet_grid(rows=vars(meets_age_crit))

# Fit sero data to gamma distribs





# Try diff method
library(bbmle)

# Functions

x_a_funct <- function(a, b0, b1, b2){
  exp((b0/b1)*a*exp(-b1*a) + (1/b1)*((b0/b1)-b2)*(exp(-b1*a)-1)-(b2*a))
}

lambda_a <- function(a, b0, b1, b2){
  (b0*a - b2) * exp(-b1 * a) + b2
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



X_a <- function(a, b0=b0, b1=b1, b2=b2){ 
  x_a(a=a, b0=b0, b1=b1, b2=b2) * get_N_age(a=a, UNcode=UNcode, year=year_tmp)
}


# Find the RSS 
min.RSS <- function(data=dat, par) {
  with(data, sum((y - x_a_funct(a=x, b0=par[1], b1=par[2], b2=par[3]))^2))
}

LL <- function(b0, b1, b2, sigma) {
  # Find the risidual 
  R <- y - x_a_funct(a=a, b0, b1, b2)
  # Calc the likelihood for the residuals
  R <- suppressWarnings(dnorm(R, 0, sigma, log=TRUE))
  
  # Sum the log likelihoods for the data points
  -sum(R)
}




###

survey_num <- 23

dat_tmp <- res_dat %>% filter(survey==survey_num) %>% as.data.frame()
(p1 <- ggplot(dat_tmp, aes(x=age, y=p2)) + geom_point() + geom_line())
fit <- smooth.spline(dat_tmp$age, dat_tmp$p2)
p1 + geom_line(data=as.data.frame(predict(fit, 0:40)), aes(x, y))

dat <- data.frame(x=dat_tmp$age, y=1-dat_tmp$p2)[-1,]
dat <- rbind(c(0, 1), dat)

dat2 <- data.frame(a=dat_tmp$age, n=dat_tmp$N, y=dat_tmp$Y, x=dat_tmp$N-dat_tmp$Y)[-1,]
dat2 <- rbind(c(0,100,0,100), dat2)
dat2$p <- 1 - dat2$x / dat2$n


UNcode <- get.UNcode.from.ISO3(dat_tmp$ISO[1])
year_tmp <- dat_tmp$year[1]
  
y <- dat$y
a <- dat$x

fit <- mle2(LL, start = list(b0=1,b1=1,b2=0,sigma=1))
summary(fit)
coef(fit)
fit.dat <- data.frame(x=0:50, 
                      y= (x_a_funct(a=0:50, b0=coef(fit)["b0"], b1=coef(fit)["b1"], b2=coef(fit)["b2"])))# + coef(fit)["mu"]))
(p1 <- ggplot(dat, aes(x=x, y=1-y)) + geom_point() + 
    geom_line(data=fit.dat, aes(x=x, y=1-y)) +
    geom_line(data=data.frame(x=0:50, y=1-(x_a_funct(a=0:50, b0=.3, b1=.45, b2=0.01))),
              aes(x=x, y=y), color='blue'))


(fit2 <- optim(par = c(.3,.45,.01), fn = min.RSS, data = dat))
(p2 <- ggplot(dat, aes(x=x, y=1-y)) + geom_point() + 
    geom_line(data=data.frame(x=0:50, y=1-x_a_funct(a=0:50, b0=fit2$par[1], b1=fit2$par[2], b2=fit2$par[3])), 
              aes(x=x, y=y)) +
    geom_line(data=data.frame(x=0:50, y=1-x_a_funct(a=0:50, b0=.3, b1=.45, b2=0.01)), 
              aes(x=x, y=y), color='blue'))





fit <- mle(LL2, start = list(b0=.4,b1=.5,b2=.001), fixed=list(data=dat2))
summary(fit)
coef(fit)
fit.dat <- data.frame(x=0:50, 
                      y= (x_a_funct(a=0:50, b0=coef(fit)["b0"], b1=coef(fit)["b1"], b2=coef(fit)["b2"])))# + coef(fit)["mu"]))
(p1 <- ggplot(dat, aes(x=x, y=1-y)) + geom_point() + 
    geom_line(data=fit.dat, aes(x=x, y=1-y)) +
    geom_line(data=data.frame(x=0:50, y=1-(x_a_funct(a=0:50, b0=.3, b1=.45, b2=0.01))),
              aes(x=x, y=y), color='blue'))

## MLE using binomial probability fits the best -- will use this for the final estimation
##  --> Should be easy to implement in a stan model



# CALC AVERAGE AGE OF INFECTION FROM FITTED FOI
b0=coef(fit)["b0"]
b1=coef(fit)["b1"]
b2=coef(fit)["b2"]

#b0=0.203; b1=0.254; b2=0.0207



L <- 90

f_num <- function(a) { a * lambda_a(a=a, b0=b0, b1=b1, b2=b2) * X_a(a=a, b0=b0, b1=b1, b2=b2) }
lambda_mean <- integrate(function(a) {lambda_a(a=a, b0=b0, b1=b1, b2=b2)}, lower=0, upper=L)$value  / L 


f_denom <- function(a) { lambda_mean * X_a(a=a, b0=b0, b1=b1, b2=b2) }
#f_denom <- function(a) { X_a(a=a, b0=b0, b1=b1, b2=b2) }

f_num2 <- function(a) { x_a(a=a, b0=b0, b1=b1, b2=b2) }
f_denom2 <- function(a) { lambda_a(a=a, b0=b0, b1=b1, b2=b2) * x_a(a=a, b0=b0, b1=b1, b2=b2) }

(A <- (integrate(f_num, lower=0, upper=L)$value /  integrate(f_denom, lower=0, upper=L)$value) )
L / A

B <- 1 / dat_tmp$mu[1]
(R0 <- 1 + (B / (A - 0.5)))
(R0 <- 1 + (B / (A)))

### THIS LOOKS WRONG ###





# Basic estimatation of lambda using x_a ------------------------------------
L <- 80

n_new_inf <- c(0, x_a_funct(0:L, b0, b1, b2) - x_a_funct((0:L)+1, b0, b1, b2)) # new infections (proportion per age group)
lambda_manual <- n_new_inf / (c(0, x_a_funct(0:L, b0, b1, b2)) * 1)
lambda_manual
lambda_a(a=0:L, b0=b0, b1=b1, b2=b2)
# fitted lambda function looks pretty good

mean(lambda_manual, na.rm=T)
mean(lambda_a(a=0:L, b0=b0, b1=b1, b2=b2))


# Manual average age of infection -------------------------------------------
N_age_ <- get_N_age(0:L, UNcode=UNcode, year=year_tmp)
N_total_ <- sum(N_age_)

p_inf_given_age <- n_new_inf[-length(n_new_inf)] # pr(new infection | age)

n_new_inf_age <- p_inf_given_age * N_age_
plot(0:L, n_new_inf_age)
tot_new_inf <- sum(n_new_inf_age)

(A <- sum((n_new_inf_age / tot_new_inf) * 0:L))

# R0
1 + ( (1/dat_tmp$mu[1]) / (A) )
1 + ( (1/dat_tmp$mu[1]) / (A-0.5) )
### !!! IT WORKED !!!!  ###################



# Version 2 - using Bayes Theorem
p_age <- N_age_ / N_total_
p_inf <- 1 - (sum(X_a(0:L, b0=b0, b1=b1, b2=b2)) / sum(get_N_age(0:L, UNcode=UNcode, year=year_tmp)))

p_age_given_inf <- ((p_inf_given_age*p_age) / p_inf) 
sum(p_age_given_inf)

A <- sum(p_age_given_inf * (0:L))





library(pracma) # for integrating numerically (trapz function) -- have to do this since we dont have continuous age populations
a <- seq(0, 70,5)
trapz(x = a, y= lambda_a(a, b0=b0, b1=b1, b2=b2))
trapz(x = a, y= 1-x_a(a, b0=b0, b1=b1, b2=b2))



#






# Checking surveys with young groups --------------------------------------

# ind_survey_lessthan1 <- sero_data_clean$survey.ind[sero_data_clean$min.age.years==0]





  
  
  






