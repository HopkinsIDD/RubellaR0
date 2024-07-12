# Functions

# x_a_funct <- function(a, b0, b1, b2){
#   exp((b0/b1)*a*exp(-b1*a) + (1/b1)*((b0/b1)-b2)*(exp(-b1*a)-1)-(b2*a))
# }
x_a_funct <- function(a, b0, b1, b2){
    1 - exp((b0/b1)*a*exp(-b1*a) + (1/b1)*((b0/b1)-b2)*(exp(-b1*a)-1) - b2*a)
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


library(bbmle)


LL2 <- function(data=dat2, b0, b1, b2) {
  # Find the risidual 
  p_i <- 1-x_a_funct(a=data$a, b0, b1, b2)
  # Calc the likelihood given binomial probilities
  lik <- suppressWarnings(dbinom(data$y, data$n, p_i, log=FALSE))
  
  # Sum the log likelihoods for the data points
  -sum(lik)
}