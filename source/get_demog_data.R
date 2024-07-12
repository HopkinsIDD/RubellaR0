

# SETUP PACKAGES AND SOURCE CODE ------------------------------------------

# if(!require('wpp2017')) install.packages('wpp2017');
# if(!require('zoo')) install.packages('zoo');

#source('source/ISO_code_source.R')






# ### Function to link country name to iso3.code and UN.code
# ###
# ##' @param - country - country name that matches the info sheet in the google doc
# ###
# ### Returns UN.code and iso3.code
# GetCountryCodes <- function(country="Madagascar"){
#   setwd("Google Drive/Rubella.Global.Estimates/data/")
#   cc <- read.csv("country_codes.csv")
#   
#   uncode <- as.numeric(country.codes$uncode[country.codes$Report_country_name==country])
#   iso3code <- as.character(country.codes$ISO3_code[country.codes$uncode==uncode & !is.na(country.codes$uncode)])
#   
#   return(list(UN.code=uncode, iso3.code=iso3code))
# }




### Function to get mid-age scaled by the age distribution for each age group and crude birth rate per 1 for each serosurvey
###
##' @param - min.age.year - vector of minimum ages for all age group (in years) for a serosurvey
##' @param - max.age.year - vector of maximum ages for each age group  (in years) for a serosurvey
##' @param - mid.year - mid-year the serosurvey was collected
##' @param - UN.code - uncode for the country where the serosurvey was collected
###
### Returns scaled mid-age based on the population age strucutre and mu (crude birth rate) at the time of the survey
GetSurveyDemo <- function(min.age.year=c(14,18,22,26,30), max.age.year=c(17,21,25,29,45), mid.year=2002, UN.code=586){
  
  age.classes=c(0.5, seq(1, 101, 1))
  out <- getDemography(UN.code, age.classes)
  
  cbr <- out$cbr.1950.2100[mid.year-1950+1]/1000
  age.struct <- out$pop.age.byageclasses.1950.2100[,(mid.year-1950+1)]
  
  mid.age.year <- rep(NA, length(min.age.year))
  for (g in 1:length(min.age.year)){
    index <- which(age.classes>=min.age.year[g] & age.classes<=(max.age.year[g]))
    ages <- age.classes[index]
    age.group.dist <- age.struct[index]/sum(age.struct[index])
    mid.age.year[g] <- round(sum(ages*age.group.dist),2)
  }
  
  cbr.survey <- rep(cbr, length(mid.age.year))
  
  return(list(mid.age.year=mid.age.year,
              mu.survey=cbr.survey))
  
}




### Function to get the mean crude birth rate per 1 for each serosurvey over X years
###
##' @param - mid.year - mid-year the serosurvey was collected
##' @param - UN.code - uncode for the country where the serosurvey was collected
###
### Returns scaled mid-age based on the population age strucutre and mu (crude birth rate) at the time of the survey
GetSurveyCBR <- function(year=2002, n_years_mean=20, UN.code=586){
  
  yr_min <- max(year - n_years_mean + 1, 1950)
  yr_max <- year
  years <- yr_min:yr_max
  years <- years[!duplicated(years)]
  
  age.classes=c(0.5, seq(1, 101, 1))
  out <- getDemography(UN.code, age.classes)
  cbrs <- out$cbr.1950.2100[years-1950+1]/1000
  cbr <- mean(cbrs)
  return(cbr)
}





### Function to get population size, age strucutre, and fertility from 1950 to 2100
###
##' @param - uncode - needs to match country name
###
### Returns demography for each country
getDemography5yr <- function(uncode=NA){
  
  #needed packages
  library(wpp2017)
  library(zoo)
  
  cc <- uncode
  
  time.by5 <- seq(1950,2100,5)
  time <- seq(1950,2100,1)
  
  #Population total in years 1950 to 2100
  data(pop)
  data(popproj)
  if (any(pop$country_code==cc)){
    pop.total.1950.2100.by5 <-  as.numeric(cbind(pop[pop$country_code==cc,3:ncol(pop)],
                                                 popproj[popproj$country_code==cc,3:ncol(popproj)]))
    f <- smooth.spline(time.by5, pop.total.1950.2100.by5)
    pop.total.1950.2100 <- predict(f,time)$y 
    
    #Population Age Structure by Age over time
    data(popFprojMed)
    data(popMprojMed)
    popF.age.1950.2100.by5 <- cbind(popF[popF$country_code==cc,4:ncol(popF)], popFprojMed[popFprojMed$country_code==cc,4:ncol(popFprojMed)])
    popM.age.1950.2100.by5 <- cbind(popM[popM$country_code==cc,4:ncol(popM)], popMprojMed[popMprojMed$country_code==cc,4:ncol(popMprojMed)])
    pop.age.1950.2100.by5 <- popF.age.1950.2100.by5+popM.age.1950.2100.by5
    ages.by5 <- popF[popF$country_code==cc,3]
    mid.age.by5 <- seq(2.5, 102.5, 5)
    popF.age.1950.2100 <- popM.age.1950.2100 <- matrix(NA, length(ages.by5), length(time))
    for (age in 1:length(ages.by5)){
      #new with wpp2017 because there are NA in older ages and oldest cohorts
      time.by5.index <- which(!is.na(popF.age.1950.2100.by5[age,]))
      time.by5.tmp <- time.by5[time.by5.index] 
      time.index <- ((time.by5.index[1]-1)*5+1):length(time)
      time.tmp <- time[time.index]
      f <- smooth.spline(time.by5.tmp, popF.age.1950.2100.by5[age,time.by5.index])
      popF.age.1950.2100[age,time.index] <- predict(f,time.tmp)$y 
      m <- smooth.spline(time.by5.tmp, popM.age.1950.2100.by5[age,time.by5.index])
      popM.age.1950.2100[age,time.index] <- predict(m,time.tmp)$y 
    }
    colnames(popF.age.1950.2100) <- colnames(popM.age.1950.2100) <- time
    rownames(popF.age.1950.2100) <- rownames(popM.age.1950.2100) <- ages.by5
    #forcing 100+ to be 0 if NA - new with wpp2017 - then using zoo::na.spline in the zoo package to fill them in
    popF.age.1950.2100[21,is.na(popF.age.1950.2100[21,])] <- 0
    popM.age.1950.2100[21,is.na(popM.age.1950.2100[21,])] <- 0
    
    pop.tmp <- popF.age.1950.2100+popM.age.1950.2100
    
    #rescale the age distibution to the estimated population each year (pop.total.1950.2100)
    pop.age.byageclasses.1950.2100 <- matrix(NA, nrow(pop.tmp), ncol(pop.tmp))
    
    for (t in 1:ncol(pop.tmp)){
      pop.age.byageclasses.1950.2100[,t] <- pop.total.1950.2100[t]*pop.tmp[,t]/sum(pop.tmp[,t], na.rm=TRUE)
    }
    rownames(pop.age.byageclasses.1950.2100) <- rownames(popF.age.1950.2100)
    colnames(pop.age.byageclasses.1950.2100) <- time
  }
  
  
  if("wpp2017" %in% (.packages())){ detach("package:wpp2017", unload=TRUE) }
  if("wpp2022" %in% (.packages())){ detach("package:wpp2022", unload=TRUE) }
  if("reshape2" %in% (.packages())){ detach("package:reshape2", unload=TRUE) }
  if("plyr" %in% (.packages())){ detach("package:plyr", unload=TRUE) }
  if("fitdistrplus" %in% (.packages())){ detach("package:fitdistrplus", unload=TRUE) }
  if("MASS" %in% (.packages())){ detach("package:MASS", unload=TRUE) }
  
  # Reload dplyr incase it got messed up
  if("dplyr" %in% (.packages())){ library(dplyr) }
  
  rm(pop, popF, popM, popFprojMed,popFTproj,popFT,popMT,popMprojMed,popMTproj,tfr,tfrprojMed,
     who_regions,percentASFR,popproj,region_data,iso_data,iso_data_full,dhs_countrydata,nationality_data, UNcode)
  
  return(list(pop.total.1950.2100=pop.age.byageclasses.1950.2100, 
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100)) 
}








### Function to get population size, age strucutre, and fertility from 1950 to 2100
###
##' @param - uncode - needs to match country name
##' @param - age.classes  - in years
###
### Returns demography for each country
getDemography <- function(uncode=NA, age.classes=c(0.5, seq(1, 101, 1))){
  
  #needed packages
  library(wpp2017)
  library(zoo)
  
  cc <- uncode
  
  time.by5 <- seq(1950,2100,5)
  time <- seq(1950,2100,1)
  
  #Population total in years 1950 to 2100
  data(pop)
  
  if (any(pop$country_code==cc)){  
    
    if (!exists("popF")){
      data(popF, package="wpp2017")
    }
    if (!exists("popM")){
      data(popM, package="wpp2017")
    }
    if (!exists("popFprojMed")){
      data(popFprojMed, package="wpp2017")
    }
    if (!exists("popMprojMed")){
      data(popMprojMed, package="wpp2017")
    }
    if (!exists("popproj")){
      data(popproj, package="wpp2017")
    }
    pop.total.1950.2100.by5 <-  as.numeric(cbind(pop[pop$country_code==cc,3:ncol(pop)],
                                                 popproj[popproj$country_code==cc,3:ncol(popproj)]))
    f <- smooth.spline(time.by5, pop.total.1950.2100.by5)
    pop.total.1950.2100 <- predict(f,time)$y 
    
    #Population Age Structure by Age over time
    data(popFprojMed)
    data(popMprojMed)
    popF.age.1950.2100.by5 <- cbind(popF[popF$country_code==cc,4:ncol(popF)], popFprojMed[popFprojMed$country_code==cc,4:ncol(popFprojMed)])
    popM.age.1950.2100.by5 <- cbind(popM[popM$country_code==cc,4:ncol(popM)], popMprojMed[popMprojMed$country_code==cc,4:ncol(popMprojMed)])
    pop.age.1950.2100.by5 <- popF.age.1950.2100.by5+popM.age.1950.2100.by5
    ages.by5 <- popF[popF$country_code==cc,3]
    mid.age.by5 <- seq(2.5, 102.5, 5)
    popF.age.1950.2100 <- popM.age.1950.2100 <- matrix(NA, length(ages.by5), length(time))
    for (age in 1:length(ages.by5)){
      #new with wpp2017 because there are NA in older ages and oldest cohorts
      time.by5.index <- which(!is.na(popF.age.1950.2100.by5[age,]))
      time.by5.tmp <- time.by5[time.by5.index] 
      time.index <- ((time.by5.index[1]-1)*5+1):length(time)
      time.tmp <- time[time.index]
      f <- smooth.spline(time.by5.tmp, popF.age.1950.2100.by5[age,time.by5.index])
      popF.age.1950.2100[age,time.index] <- predict(f,time.tmp)$y 
      m <- smooth.spline(time.by5.tmp, popM.age.1950.2100.by5[age,time.by5.index])
      popM.age.1950.2100[age,time.index] <- predict(m,time.tmp)$y 
    }
    colnames(popF.age.1950.2100) <- colnames(popM.age.1950.2100) <- time
    rownames(popF.age.1950.2100) <- rownames(popM.age.1950.2100) <- ages.by5
    #forcing 100+ to be 0 if NA - new with wpp2017 - then using zoo::na.spline in the zoo package to fill them in
    popF.age.1950.2100[21,is.na(popF.age.1950.2100[21,])] <- 0
    popM.age.1950.2100[21,is.na(popM.age.1950.2100[21,])] <- 0
    popF.age.byageclasses.1950.2100 <- apply(popF.age.1950.2100, 2, function(x) predict(smooth.spline(mid.age.by5, zoo::na.spline(x)), age.classes)$y)
    popM.age.byageclasses.1950.2100 <- apply(popM.age.1950.2100, 2, function(x) predict(smooth.spline(mid.age.by5, zoo::na.spline(x)), age.classes)$y)
    #adjust for varying bin width!
    popF.age.byageclasses.1950.2100 <- popF.age.byageclasses.1950.2100*diff(c(0,age.classes)) 
    popM.age.byageclasses.1950.2100 <- popM.age.byageclasses.1950.2100*diff(c(0,age.classes))
    pop.tmp <- popF.age.byageclasses.1950.2100+popM.age.byageclasses.1950.2100
    #rescale the age distibution to the estimated population each year (pop.total.1950.2100)
    pop.age.byageclasses.1950.2100 <- matrix(NA, nrow(pop.tmp), ncol(pop.tmp))
    for (t in 1:ncol(pop.tmp)){
      pop.age.byageclasses.1950.2100[,t] <- pop.total.1950.2100[t]*pop.tmp[,t]/sum(pop.tmp[,t])
    }
    rownames(pop.age.byageclasses.1950.2100) <- age.classes
    colnames(pop.age.byageclasses.1950.2100) <- time
    
    #Sex Distribution by age and over time
    pop.tmp <- popF.age.byageclasses.1950.2100+popM.age.byageclasses.1950.2100
    sex.dist <- matrix(NA, nrow(pop.tmp), ncol(pop.tmp))
    for (c in 1:ncol(pop.tmp)){
      sex.dist[,c] <- popF.age.byageclasses.1950.2100[,c]/pop.tmp[,c]
    }
    
    #TFR over time
    data(tfr)
    data(tfrprojMed)
    mid.time.by5 <- seq(1952, 2097, 5) #TFR is given over a range, therefore assume it is the mid-period
    tfr.1950.2100.by5 <-  as.numeric(cbind(tfr[tfr$country_code==cc,3:15],
                                           tfrprojMed[tfrprojMed$country_code==cc,3:ncol(tfrprojMed)]))
    f <- smooth.spline(mid.time.by5, tfr.1950.2100.by5)
    tfr.1950.2100 <- predict(f,time)$y 
    names(tfr.1950.2100) <- time
    
    #Number of women of reproductive age (in five year age groups 15 to 50) over time
    repro.ages <- seq(15,45,5)
    popF.15to50.1950.2100 <- matrix(NA,7,length(time))
    for (c in 1:ncol(pop.age.byageclasses.1950.2100)){
      for (a in 1:length(repro.ages)){
        index <- which(age.classes>=repro.ages[a] & age.classes<(repro.ages[a]+5))
        popF.15to50.1950.2100[a,c] <- sum(pop.age.byageclasses.1950.2100[index,c]*sex.dist[index,c])
      }
    }
    colnames(popF.15to50.1950.2100) <- time
    rownames(popF.15to50.1950.2100) <- repro.ages
    
    #ASFR over time
    data(percentASFR)
    p.asfr <- (percentASFR[percentASFR$country_code==cc,])
    p.asfr.1950.2100 <- matrix(NA, nrow=nrow(p.asfr), ncol=(length(time)))
    for (t in 1:(ncol(p.asfr)-3)){
      p.asfr.1950.2100[,(t*5-4):(t*5)] <- p.asfr[,(t+3)]/100
    }
    p.asfr.1950.2100[,ncol(p.asfr.1950.2100)] <-  p.asfr[,ncol(p.asfr)]/100
    asfr.1950.2100 <- matrix(NA, nrow=nrow(p.asfr), ncol=(length(time)))
    for (t in 1:ncol(asfr.1950.2100)){
      asfr.1950.2100[,t] <- tfr.1950.2100[t]/5*p.asfr.1950.2100[,t]
    }
    colnames(asfr.1950.2100) <- time
    
    #Births over time
    births.1950.2100 <- rep(NA, length(time))
    for (t in 1:length(births.1950.2100)){
      births.1950.2100[t] <- sum(asfr.1950.2100[,t]*popF.15to50.1950.2100[,t])
    }
    names(births.1950.2100) <- time
    
    #Crude Birth Rates over time
    cbr.1950.2100 <- rep(NA, length(time))
    cbr.1950.2100 <- births.1950.2100/pop.total.1950.2100*1000
    names(cbr.1950.2100) <- time
  } else {
    print(paste("No Population Data for ", country, sep=""))
    pop.total.1950.2100=NA
    pop.age.byageclasses.1950.2100=NA
    tfr.1950.2100=NA
    births.1950.2100=NA
    cbr.1950.2100=NA
  }
  
  
  if("wpp2017" %in% (.packages())){ detach("package:wpp2017", unload=TRUE) }
  if("wpp2022" %in% (.packages())){ detach("package:wpp2022", unload=TRUE) }
  if("reshape2" %in% (.packages())){ detach("package:reshape2", unload=TRUE) }
  if("plyr" %in% (.packages())){ detach("package:plyr", unload=TRUE) }
  if("fitdistrplus" %in% (.packages())){ detach("package:fitdistrplus", unload=TRUE) }
  if("MASS" %in% (.packages())){ detach("package:MASS", unload=TRUE) }
  
  # Reload dplyr incase it got messed up
  if("dplyr" %in% (.packages())){ library(dplyr) }
  
  rm(pop, popF, popM, popFprojMed,popFTproj,popFT,popMT,popMprojMed,popMTproj,tfr,tfrprojMed,
     who_regions,percentASFR,popproj,region_data,iso_data,iso_data_full,dhs_countrydata,nationality_data, UNcode)
  
  return(list(pop.total.1950.2100=pop.total.1950.2100, 
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100, 
              tfr.1950.2100=tfr.1950.2100,
              births.1950.2100=births.1950.2100, 
              cbr.1950.2100=cbr.1950.2100)) 
}






### Function to get population size, age strucutre, and fertility from 1950 to 2100
###
##' @param - uncode - needs to match country name
##' @param - age.classes  - in years
###
### Returns demography for each country
getPopAge <- function(uncode=NA, age.classes=c(0.5, seq(1, 101, 1))){
  
  #needed packages
  library(wpp2017)
  library(zoo)
  
  cc <- uncode
  
  time.by5 <- seq(1950,2100,5)
  time <- seq(1950,2100,1)
  
  #Population total in years 1950 to 2100
  data(pop)
  
  if (any(pop$country_code==cc)){
    
    if (!exists("popF")){
      data(popF, package="wpp2017")
    }
    if (!exists("popM")){
      data(popM, package="wpp2017")
    }
    if (!exists("popFprojMed")){
      data(popFprojMed, package="wpp2017")
    }
    if (!exists("popMprojMed")){
      data(popMprojMed, package="wpp2017")
    }
    if (!exists("popproj")){
      data(popproj, package="wpp2017")
    }
    
    pop.total.1950.2100.by5 <-  as.numeric(cbind(pop[pop$country_code==cc, 3:ncol(pop)],
                                                 popproj[popproj$country_code==cc,3:ncol(popproj)]))
    f <- smooth.spline(time.by5, pop.total.1950.2100.by5)
    pop.total.1950.2100 <- predict(f,time)$y 
    
    #Population Age Structure by Age over time
    
    popF.age.1950.2100.by5 <- cbind(popF[popF$country_code==cc,4:ncol(popF)], popFprojMed[popFprojMed$country_code==cc,4:ncol(popFprojMed)])
    popM.age.1950.2100.by5 <- cbind(popM[popM$country_code==cc,4:ncol(popM)], popMprojMed[popMprojMed$country_code==cc,4:ncol(popMprojMed)])
    pop.age.1950.2100.by5 <- popF.age.1950.2100.by5+popM.age.1950.2100.by5
    ages.by5 <- popF[popF$country_code==cc,3]
    mid.age.by5 <- seq(2.5, 102.5, 5)
    popF.age.1950.2100 <- popM.age.1950.2100 <- matrix(NA, length(ages.by5), length(time))
    for (age in 1:length(ages.by5)){
      #new with wpp2017 because there are NA in older ages and oldest cohorts
      time.by5.index <- which(!is.na(popF.age.1950.2100.by5[age,]))
      time.by5.tmp <- time.by5[time.by5.index] 
      time.index <- ((time.by5.index[1]-1)*5+1):length(time)
      time.tmp <- time[time.index]
      f <- smooth.spline(time.by5.tmp, popF.age.1950.2100.by5[age,time.by5.index])
      popF.age.1950.2100[age,time.index] <- predict(f,time.tmp)$y 
      m <- smooth.spline(time.by5.tmp, popM.age.1950.2100.by5[age,time.by5.index])
      popM.age.1950.2100[age,time.index] <- predict(m,time.tmp)$y 
    }
    colnames(popF.age.1950.2100) <- colnames(popM.age.1950.2100) <- time
    rownames(popF.age.1950.2100) <- rownames(popM.age.1950.2100) <- ages.by5
    #forcing 100+ to be 0 if NA - new with wpp2017 - then using zoo::na.spline in the zoo package to fill them in
    popF.age.1950.2100[21,is.na(popF.age.1950.2100[21,])] <- 0
    popM.age.1950.2100[21,is.na(popM.age.1950.2100[21,])] <- 0
    popF.age.byageclasses.1950.2100 <- apply(popF.age.1950.2100, 2, function(x) predict(smooth.spline(mid.age.by5, zoo::na.spline(x)), age.classes)$y)
    popM.age.byageclasses.1950.2100 <- apply(popM.age.1950.2100, 2, function(x) predict(smooth.spline(mid.age.by5, zoo::na.spline(x)), age.classes)$y)
    #adjust for varying bin width!
    popF.age.byageclasses.1950.2100 <- popF.age.byageclasses.1950.2100*diff(c(0,age.classes)) 
    popM.age.byageclasses.1950.2100 <- popM.age.byageclasses.1950.2100*diff(c(0,age.classes))
    pop.tmp <- popF.age.byageclasses.1950.2100+popM.age.byageclasses.1950.2100
    #rescale the age distibution to the estimated population each year (pop.total.1950.2100)
    pop.age.byageclasses.1950.2100 <- matrix(NA, nrow(pop.tmp), ncol(pop.tmp))
    for (t in 1:ncol(pop.tmp)){
      pop.age.byageclasses.1950.2100[,t] <- pop.total.1950.2100[t]*pop.tmp[,t]/sum(pop.tmp[,t])
    }
    rownames(pop.age.byageclasses.1950.2100) <- age.classes
    colnames(pop.age.byageclasses.1950.2100) <- time
    
    
  } else {
    print(paste("No Population Data for ", country, sep=""))
    pop.total.1950.2100=NA
    pop.age.byageclasses.1950.2100=NA
  }
  
  
  if("wpp2017" %in% (.packages())){ detach("package:wpp2017", unload=TRUE) }
  if("wpp2022" %in% (.packages())){ detach("package:wpp2022", unload=TRUE) }
  if("reshape2" %in% (.packages())){ detach("package:reshape2", unload=TRUE) }
  if("plyr" %in% (.packages())){ detach("package:plyr", unload=TRUE) }
  if("fitdistrplus" %in% (.packages())){ detach("package:fitdistrplus", unload=TRUE) }
  if("MASS" %in% (.packages())){ detach("package:MASS", unload=TRUE) }
  
  # Reload dplyr incase it got messed up
  if("dplyr" %in% (.packages())){ library(dplyr) }
  
  rm(pop, popF, popM, popFprojMed,popFTproj,popFT,popMT,popMprojMed,popMTproj,
     who_regions,popproj,region_data,iso_data,iso_data_full,dhs_countrydata,nationality_data, UNcode)
  
  return(list(pop.total.1950.2100=pop.total.1950.2100, 
              pop.age.byageclasses.1950.2100=pop.age.byageclasses.1950.2100)) 
}








### Function to get population size, age strucutre, and fertility from 1950 to 2100
###
##' @param - uncode - needs to match country name
##' @param - age.classes  - in years
###
### Returns demography for each country
getLifeExpectancy <- function(uncode=NA, age.classes=c(0.5, seq(1, 101, 1))){
  
  #needed packages
  library(wpp2017)
  library(zoo)
  
  cc <- uncode
  
  time.by5 <- seq(1950,2100,5)
  time <- seq(1950,2100,1)
  
  #Population total in years 1950 to 2100
  data(pop)
  data(popproj)
  
  if (any(pop$country_code==cc)){
    pop.total.1950.2100.by5 <-  as.numeric(cbind(pop[pop$country_code==cc,3:ncol(pop)],
                                                 popproj[popproj$country_code==cc,3:ncol(popproj)]))
    f <- smooth.spline(time.by5, pop.total.1950.2100.by5)
    pop.total.1950.2100 <- predict(f,time)$y 
    
    #Population Age Structure by Age over time
    data(popFprojMed)
    data(popMprojMed)
    popF.age.1950.2100.by5 <- cbind(popF[popF$country_code==cc,4:ncol(popF)], popFprojMed[popFprojMed$country_code==cc,4:ncol(popFprojMed)])
    popM.age.1950.2100.by5 <- cbind(popM[popM$country_code==cc,4:ncol(popM)], popMprojMed[popMprojMed$country_code==cc,4:ncol(popMprojMed)])
    
    #Life Expency over time
    data(e0M)
    data(e0Mproj)
    data(e0F)
    data(e0Fproj)
    
    mid.time.by5 <- seq(1952, 2097, 5) #Life expectancy is given over a range, therefore assume it is the mid-period
    e0M.1950.2100.by5 <-  as.numeric(cbind(e0M[e0M$country_code==cc,3:15],
                                           e0Mproj[e0Mproj$country_code==cc,3:ncol(e0Mproj)]))
    e0F.1950.2100.by5 <-  as.numeric(cbind(e0F[e0F$country_code==cc,3:15],
                                           e0Fproj[e0Fproj$country_code==cc,3:ncol(e0Fproj)]))
    
    popM.total.1950.2100.by5 <- colSums(popM.age.1950.2100.by5, na.rm = TRUE)[-ncol(popM.age.1950.2100.by5)]
    popF.total.1950.2100.by5 <- colSums(popF.age.1950.2100.by5, na.rm = TRUE)[-ncol(popM.age.1950.2100.by5)]
    sex.prop.M <- popM.total.1950.2100.by5 / (popM.total.1950.2100.by5 + popF.total.1950.2100.by5)
    sex.prop.F <- popF.total.1950.2100.by5 / (popM.total.1950.2100.by5 + popF.total.1950.2100.by5)
    
    e0.1950.2100.by5 <- e0M.1950.2100.by5*sex.prop.M + e0F.1950.2100.by5*sex.prop.F
    f <- smooth.spline(mid.time.by5, e0.1950.2100.by5)
    e0.1950.2100 <- predict(f,time)$y 
    names(e0.1950.2100) <- time
    
    
  } else {
    print(paste("No Population Data for ", country, sep=""))
    e0.1950.2100=NA
  }
  
  if("wpp2017" %in% (.packages())){ detach("package:wpp2017", unload=TRUE) }
  if("reshape2" %in% (.packages())){ detach("package:reshape2", unload=TRUE) }
  if("plyr" %in% (.packages())){ detach("package:plyr", unload=TRUE) }
  if("fitdistrplus" %in% (.packages())){ detach("package:fitdistrplus", unload=TRUE) }
  if("MASS" %in% (.packages())){ detach("package:MASS", unload=TRUE) }
  
  # Reload dplyr incase it got messed up
  if("dplyr" %in% (.packages())){ library(dplyr) }
  
  rm(pop, popF, popM, popFprojMed,popFTproj,popFT,popMT,popMprojMed,popMTproj,tfr,tfrprojMed,
     who_regions,percentASFR,popproj,region_data,iso_data,iso_data_full,dhs_countrydata,nationality_data, UNcode)
  
  return(e0.1950.2100)
}
