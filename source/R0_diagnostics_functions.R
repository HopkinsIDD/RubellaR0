##' R0 Diagnostics Functions
##' 




# SEROPREVALENCE FITS ----------------------------------------------------------------

summ_sero_fits <- function(fits, analysis_version, sero_data_stan, sero_data){
    
    require(dplyr)
    
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
    return(res)
}




plot_probs <- function(data, prior_limits=FALSE, age_limits=c(0,30)){
    
    require(tidyverse)
    
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
    require(tidyverse)
    data_ <- data %>% filter(C==X)
    plot_probs(data_, prior_limits=FALSE, age_limits=c(0,40)) 
}

#plot_multi(data=res, X=4)




plot_sero_fits <- function(fits, analysis_version, sero_data_stan, sero_data){
    
    res <- summ_sero_fits(fits, analysis_version, sero_data_stan, sero_data)
        
    require(tidyverse)
    analysis_version_repo <- file.path("R0_model", "results", analysis_version)

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

    return(print(paste0("All plots save to [ ", file.path(analysis_version_repo, "seropositivity_byage.pdf")," ].")))
}




plot_country_sero_fits <- function(country=1, age_limits=c(0,40),
                                   fits_=fits, analysis_version_=analysis_version, 
                                   sero_data_stan_=sero_data_stan, sero_data_=sero_data){
   
    res <- summ_sero_fits(fits_, analysis_version_, sero_data_stan_, sero_data_)
    
    require(tidyverse)
    if (!is.numeric(country)){
        data_ <- res %>% filter(tolower(country)==tolower(country))
    } else {
        data_ <- res %>% filter(C==country)
    }
    return(plot_probs(data_, prior_limits=FALSE, age_limits))
}






