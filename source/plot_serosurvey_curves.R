##' plot_serosurvey_curves.R 
##' Plot Serosurvey Cumulative Seropositivity Curves
##' 


plot_serosurvey_curves <- function(sero_data, proj_dir){

    sero_data <- sero_data %>% mutate(survey_id = paste0(survey.num, " - ", country.year)) %>% as.data.frame()
    
    # plot 9 per page
    survey_nums <- unique(sero_data$survey.num)
    n_page <- ceiling(length(survey_nums) / 9)
    
    pdf(file.path(proj_dir, "serology_plot_pre.pdf"), width=8.5, height=11)
    
    for (pg in 1:n_page){
        
        surveys_tmp <- survey_nums[(pg*9-8) : (pg*9)] # surveys to plot on this page
        survey_dat <- sero_data %>% filter(survey.num %in% surveys_tmp)
        binom_ci <- binom::binom.confint(survey_dat$Y, survey_dat$N, methods = "wilson")
        
        # gamma=0.1
        # p_ <- survey_dat$Y / survey_dat$N
        # beta = (1/gamma - 1) * (1 - p_);
        # alpha = p_ * (1/gamma - 1);
        # res <- sapply(1, function(x) rmutil::rbetabinom(size=1, survey_dat$N[x], p_[x], gamma))# / survey_dat$N[x] )
        # lapply(res, quantile, prob=0.025)
        # lapply(res, quantile, prob=0.975)

        survey_dat <- survey_dat %>% mutate(binom_low = binom_ci$lower,
                                            binom_high = binom_ci$upper)
        
        max_age <- ceiling(max(survey_dat$scaled.mid.age))
        p1 <- ggplot(survey_dat, aes(x=scaled.mid.age, y=seropositive.immune)) + 
            geom_point(color="navy", size=3) +
            ylab("Seropositive") + xlab("Age (yr)") + 
            coord_cartesian(xlim=c(0, max_age), ylim=c(0,1)) +
            geom_errorbar(aes(ymin=binom_low, ymax=binom_high), width=0, color="navy", size=1.25) +
            theme_classic() +
            facet_wrap(~survey_id)#, labeller = labeller(survey.num = country.year))
        print(p1)
    }    
        
    dev.off()
    print(paste0("Plots saved to ", file.path(proj_dir, "serology_plot_pre.pdf"), "."))
        
}





sero_curve_funct <- function(a, b0, b1, b2){
    1 - exp((b0/b1)*a*exp(-b1*a) + (1/b1)*((b0/b1)-b2)*(exp(-b1*a)-1) - b2*a)
}

#sero_curve_funct(0:70, sero_data$b0[1], sero_data$b1[1], sero_data$b2[1])



plot_serosurvey_curves_withfits_nobands <- function(sero_data, proj_dir, b_fits){
    
    sero_data_tmp <- sero_data %>% mutate(survey_id = paste0(survey.num, " - ", country.year)) %>% as.data.frame()
    
    # add fitted curve b values
    sero_data_tmp <- left_join(sero_data_tmp, b_fits, by=c("survey.num"="survey_num"))
    
    # plot 9 per page
    survey_nums <- unique(sero_data_tmp$survey.num)
    n_page <- ceiling(length(survey_nums) / 9)
    
    pdf(file.path(proj_dir, "serology_plot_post.pdf"), width=8.5, height=11)
    
    for (pg in 1:n_page){
        
        surveys_tmp <- survey_nums[(pg*9-8) : (pg*9)] # surveys to plot on this page
        survey_dat <- sero_data_tmp %>% filter(survey.num %in% surveys_tmp)
        binom_ci <- binom::binom.confint(survey_dat$Y, survey_dat$N, methods = "wilson")
        
        # gamma=0.1
        # p_ <- survey_dat$Y / survey_dat$N
        # beta = (1/gamma - 1) * (1 - p_);
        # alpha = p_ * (1/gamma - 1);
        # res <- sapply(1, function(x) rmutil::rbetabinom(size=1, survey_dat$N[x], p_[x], gamma))# / survey_dat$N[x] )
        # lapply(res, quantile, prob=0.025)
        # lapply(res, quantile, prob=0.975)
        
        survey_dat <- survey_dat %>% mutate(binom_low = binom_ci$lower,
                                            binom_high = binom_ci$upper)
        
        # Add ages for plotting fits
        survey_dat_0 <- survey_dat %>% filter(!duplicated(survey.num)) %>% mutate(scaled.mid.age = 0, seropositive.immune=NA, binom_low=NA, binom_high=NA)
        
        survey_dat <- bind_rows(survey_dat, survey_dat_0, survey_dat_0 %>% mutate(scaled.mid.age = 1), survey_dat_0 %>% mutate(scaled.mid.age = 2),
                                survey_dat_0 %>% mutate(scaled.mid.age = 3),
                                survey_dat_0 %>% mutate(scaled.mid.age = 5),  survey_dat_0 %>% mutate(scaled.mid.age = 10))
        # Add fitted estimates
        survey_dat <- survey_dat %>% mutate(est_p = sero_curve_funct(scaled.mid.age, b0, b1, b2))
            
        max_age <- ceiling(max(survey_dat$scaled.mid.age))
        p1 <- ggplot(survey_dat, aes(x=scaled.mid.age, y=seropositive.immune)) + 
            geom_point(color="navy", size=3) +
            ylab("Seropositive") + xlab("Age (yr)") + 
            geom_line(aes(x=scaled.mid.age, y=est_p), color="maroon", size=1.1) +
            coord_cartesian(xlim=c(0, max_age), ylim=c(0,1)) +
            geom_errorbar(aes(ymin=binom_low, ymax=binom_high), width=0, color="navy", size=1.25) +
            theme_classic() +
            facet_wrap(~survey_id)#, labeller = labeller(survey.num = country.year))
        print(p1)
    }    
    
    dev.off()
    print(paste0("Plots saved to ", file.path(proj_dir, "serology_plot_post_nobands.pdf"), "."))
    
}





##' 
##' Get b_fits and individual b_ trajectories from the stanfit object. 
##' - for plottin confidence bands on seroprev curves
##' 
##' 
##' 
get_b_for_plots <- function(stanfits = file.path(analysis_version_repo, "R0_full_3000iters.RData")){
    
    # Load model fits
    load(file=stanfits) # Loads ret_sm, fit_R0, sero_data_stan_list
    fits <- rstan::extract(fit_R0)
    
    # ~ get median fits
    b0 <- apply(fits$b0, 2, median)
    b1 <- apply(fits$b1, 2, median)
    b2 <- apply(fits$b2, 2, median)
    b_fits <- data.frame(survey_num = sero_data_stan$survey[!duplicated(sero_data_stan$survey)], b0=b0, b1=b1, b2=b2 )
    
    # ~ get full fits
    b_fits_full_b0 <- fits$b0
    b_fits_full_b1 <- fits$b1
    b_fits_full_b2 <- fits$b2
    
    return(list(b_fits=b_fits, 
                b_fits_full_b0=b_fits_full_b0,
                b_fits_full_b1=b_fits_full_b1, 
                b_fits_full_b2=b_fits_full_b2))
    
}







plot_serosurvey_curves_withfits <- function(stanfits = file.path("R0_model","results", results_version, "R0_full_3000iters.RData"), 
                                            proj_dir){
    
    load(file=stanfits) # Loads ret_sm, fit_R0, sero_data_stan_list, sero_data
    
    # Get the b fits
    b_fits_list <- get_b_for_plots(stanfits)
    
    sero_data_tmp <- sero_data %>% mutate(survey_id = paste0(survey.num, " - ", country.year)) %>% as.data.frame() %>% select(survey.num, survey.ind, survey_id, everything())
    
    # add fitted curve b values
    sero_data_tmp <- left_join(sero_data_tmp, b_fits_list$b_fits, by=c("survey.num"="survey_num"))
    
    # plot 9 per page
    survey_nums <- unique(sero_data_tmp$survey.num)
    n_page <- ceiling(length(survey_nums) / 9)
    
    pdf(file.path(proj_dir, "serology_plot_post.pdf"), width=8.5, height=11)
    
    for (pg in 1:n_page){
        
        surveys_tmp <- survey_nums[(pg*9-8) : (pg*9)] # surveys to plot on this page
        surveys_tmp <- surveys_tmp[!is.na(surveys_tmp)]
        survey_dat <- sero_data_tmp %>% filter(survey.num %in% surveys_tmp)
        #binom_ci <- binom::binom.confint(survey_dat$Y, survey_dat$N, methods = "wilson")
        binom_ci <- as_tibble(t(sapply(X=1:nrow(survey_dat), function(i=X) as.numeric(prop.test(x=survey_dat$Y[i], n=survey_dat$N[i], conf.level = .95, correct = TRUE)$conf.int))))
        # gamma=0.1
        # p_ <- survey_dat$Y / survey_dat$N
        # beta = (1/gamma - 1) * (1 - p_);
        # alpha = p_ * (1/gamma - 1);
        # res <- sapply(1, function(x) rmutil::rbetabinom(size=1, survey_dat$N[x], p_[x], gamma))# / survey_dat$N[x] )
        # lapply(res, quantile, prob=0.025)
        # lapply(res, quantile, prob=0.975)
        
        survey_dat <- survey_dat %>% mutate(binom_low = as.numeric(unlist(binom_ci[,1])),
                                            binom_high = as.numeric(unlist(binom_ci[,2])))
        
        # Add ages for plotting fits
        survey_dat_0 <- survey_dat %>% filter(!duplicated(survey.num)) %>% mutate(scaled.mid.age = 0, seropositive.immune=NA, binom_low=NA, binom_high=NA)
        
        survey_dat <- bind_rows(survey_dat, survey_dat_0, survey_dat_0 %>% mutate(scaled.mid.age = 1), survey_dat_0 %>% mutate(scaled.mid.age = 2),
                                survey_dat_0 %>% mutate(scaled.mid.age = 3),
                                survey_dat_0 %>% mutate(scaled.mid.age = 5),  survey_dat_0 %>% mutate(scaled.mid.age = 10))
        # Add fitted estimates
        survey_dat <- survey_dat %>% mutate(est_p = sero_curve_funct(scaled.mid.age, b0, b1, b2))
        tmp <- survey_dat %>% select(survey_id, survey.num, scaled.mid.age, b0, b1, b2)
        
        
        # Add fitted estimates
        survey_id_ <- unique(survey_dat$survey_id)
        survey.num_ <- surveys_tmp
        b_fit_b0 <- b_fits_list$b_fits_full_b0[,surveys_tmp]
        b_fit_b1 <- b_fits_list$b_fits_full_b1[,surveys_tmp]
        b_fit_b2 <- b_fits_list$b_fits_full_b2[,surveys_tmp]
        
        ages_ <- seq(0,60,1)
        tmpcurve <- matrix(NA, nrow=nrow(b_fit_b0), ncol=length(ages_))
        fit_curve_ll <- matrix(NA, nrow=ncol(b_fit_b0), ncol=length(ages_))
        fit_curve_ul <- matrix(NA, nrow=ncol(b_fit_b0), ncol=length(ages_))
        fit_curve_md <- matrix(NA, nrow=ncol(b_fit_b0), ncol=length(ages_))
        
        for (i in 1:ncol(b_fit_b0)){
            for (s in 1:nrow(b_fit_b0)){
                tmpcurve[s,] <- sero_curve_funct(ages_, 
                                                 b0=b_fit_b0[s,i], 
                                                 b1=b_fit_b1[s,i], 
                                                 b2=b_fit_b2[s,i])
            }    
            
            fit_curve_ll[i,] <- apply(tmpcurve, 2, quantile, probs=.025, na.rm=TRUE)       
            fit_curve_ul[i,] <- apply(tmpcurve, 2, quantile, probs=.975, na.rm=TRUE)        
            fit_curve_md[i,] <- apply(tmpcurve, 2, quantile, probs=.5, na.rm=TRUE)        
        }
        colnames(fit_curve_md) <- ages_
        colnames(fit_curve_ll) <- ages_
        colnames(fit_curve_ul) <- ages_
        
        fits_long <- fit_curve_md %>% as.data.frame() %>% mutate(survey_id = survey_id_, survey.num=survey.num_) %>%
            gather(key="age", value = "est_md", -c(survey_id, survey.num))
        fits_longul <- fit_curve_ul %>% as.data.frame() %>% mutate(survey_id = survey_id_, survey.num=survey.num_) %>%
            gather(key="age", value = "est_ul", -c(survey_id, survey.num))
        fits_longll <- fit_curve_ll %>% as.data.frame() %>% mutate(survey_id = survey_id_, survey.num=survey.num_) %>%
            gather(key="age", value = "est_ll", -c(survey_id, survey.num))
        fits_ <- bind_cols(bind_cols(fits_long, fits_longll %>% select(est_ll)), fits_longul %>% select(est_ul))
        fits_$scaled.mid.age <- fits_$age
        
        data <- full_join(survey_dat %>% select(survey_id, survey.num, scaled.mid.age, seropositive.immune, binom_low, binom_high),
                          fits_ %>% mutate(scaled.mid.age = as.integer(scaled.mid.age),
                                           survey_id = as.character(survey_id)), 
                          by=c("scaled.mid.age", "survey_id" ,"survey.num"))
        data <- data %>% arrange(survey.num, scaled.mid.age)
        
        max_age <- ceiling(max(data$scaled.mid.age))
        
        p1 <- ggplot() + 
            geom_ribbon(data=data %>% filter(!is.na(age)), aes(x=scaled.mid.age, ymin=est_ll, ymax=est_ul), fill="maroon", alpha=.4) +
            geom_path(data=data %>% filter(!is.na(age)), aes(x=scaled.mid.age, y=est_md), color="maroon", size=1.1) +
            geom_errorbar(data=data, aes(x=scaled.mid.age, ymin=binom_low, ymax=binom_high), width=0, color="navy", size=1.25) +
            geom_point(data=data , aes(x=scaled.mid.age, y=seropositive.immune), color="navy", size=3) +
            theme_classic() +
            coord_cartesian(xlim=c(0, max_age), ylim=c(0,1)) +
            ylab("Seropositive") + xlab("Age (yr)") + 
            facet_wrap(~survey_id)#, labeller = labeller(survey.num = country.year))
        print(p1)
    }    
    
    dev.off()
    print(paste0("Plots saved to ", file.path(proj_dir, "serology_plot_post.pdf"), "."))
}







plot_serosurvey_curves_withfits_modelsims <- function(modelsims_seroprev, stanfits = file.path("R0_model","results", results_version, "R0_full_3000iters.RData"), 
                                            proj_dir){
    
    load(file=stanfits) # Loads ret_sm, fit_R0, sero_data_stan_list, sero_data
    
    # Get the b fits
    b_fits_list <- get_b_for_plots(stanfits)
    
    sero_data_tmp <- sero_data %>% mutate(survey_id = paste0(survey.num, " - ", country.year)) %>% as.data.frame() %>% select(survey.num, survey.ind, survey_id, everything())
    
    # add fitted curve b values
    sero_data_tmp <- left_join(sero_data_tmp, b_fits_list$b_fits, by=c("survey.num"="survey_num"))
    
    
    # ADD Seroprev from transmission model
    sero_survey_info <- sero_data_tmp %>% select(survey_id, survey.num, survey.ind, country, ISO, start.year ) %>% filter(!duplicated(survey_id))     # get individual survey info 
    sero_survey_info <- sero_survey_info %>% mutate(start.year = ifelse(start.year<=1980, 1981, start.year))
    seroprev_ <- left_join(modelsims_seroprev, sero_survey_info %>% rename(country_name2 = country), by=c("year"="start.year", "country"="ISO")) %>% filter(!is.na(survey_id))  # Reduce seroprev to data wanted
    
    
    # plot 9 per page
    survey_nums <- unique(sero_data_tmp$survey.num)
    n_page <- ceiling(length(survey_nums) / 9)
    
    pdf(file.path(proj_dir, "serology_plot_post_modelsims.pdf"), width=8.5, height=11)
    
    for (pg in 1:n_page){
        
        surveys_tmp <- survey_nums[(pg*9-8) : (pg*9)] # surveys to plot on this page
        surveys_tmp <- surveys_tmp[!is.na(surveys_tmp)]
        survey_dat <- sero_data_tmp %>% filter(survey.num %in% surveys_tmp)
        binom_ci <- binom::binom.confint(survey_dat$Y, survey_dat$N, methods = "wilson")
        
        # gamma=0.1
        # p_ <- survey_dat$Y / survey_dat$N
        # beta = (1/gamma - 1) * (1 - p_);
        # alpha = p_ * (1/gamma - 1);
        # res <- sapply(1, function(x) rmutil::rbetabinom(size=1, survey_dat$N[x], p_[x], gamma))# / survey_dat$N[x] )
        # lapply(res, quantile, prob=0.025)
        # lapply(res, quantile, prob=0.975)
        
        survey_dat <- survey_dat %>% mutate(binom_low = binom_ci$lower,
                                            binom_high = binom_ci$upper)
        
        # Add ages for plotting fits
        survey_dat_0 <- survey_dat %>% filter(!duplicated(survey.num)) %>% mutate(scaled.mid.age = 0, seropositive.immune=NA, binom_low=NA, binom_high=NA)
        
        survey_dat <- bind_rows(survey_dat, survey_dat_0, survey_dat_0 %>% mutate(scaled.mid.age = 1), survey_dat_0 %>% mutate(scaled.mid.age = 2),
                                survey_dat_0 %>% mutate(scaled.mid.age = 3),
                                survey_dat_0 %>% mutate(scaled.mid.age = 5),  survey_dat_0 %>% mutate(scaled.mid.age = 10))
        # Add fitted estimates
        survey_dat <- survey_dat %>% mutate(est_p = sero_curve_funct(scaled.mid.age, b0, b1, b2))
        tmp <- survey_dat %>% select(survey_id, survey.num, scaled.mid.age, b0, b1, b2)
        
        
        # Add fitted estimates
        survey_id_ <- unique(survey_dat$survey_id)
        survey.num_ <- surveys_tmp
        b_fit_b0 <- b_fits_list$b_fits_full_b0[,surveys_tmp]
        b_fit_b1 <- b_fits_list$b_fits_full_b1[,surveys_tmp]
        b_fit_b2 <- b_fits_list$b_fits_full_b2[,surveys_tmp]
        
        ages_ <- seq(0,60,1)
        tmpcurve <- matrix(NA, nrow=nrow(b_fit_b0), ncol=length(ages_))
        fit_curve_ll <- matrix(NA, nrow=ncol(b_fit_b0), ncol=length(ages_))
        fit_curve_ul <- matrix(NA, nrow=ncol(b_fit_b0), ncol=length(ages_))
        fit_curve_md <- matrix(NA, nrow=ncol(b_fit_b0), ncol=length(ages_))
        
        for (i in 1:ncol(b_fit_b0)){
            for (s in 1:nrow(b_fit_b0)){
                tmpcurve[s,] <- sero_curve_funct(ages_, 
                                                 b0=b_fit_b0[s,i], 
                                                 b1=b_fit_b1[s,i], 
                                                 b2=b_fit_b2[s,i])
            }    
            
            fit_curve_ll[i,] <- apply(tmpcurve, 2, quantile, probs=.025, na.rm=TRUE)       
            fit_curve_ul[i,] <- apply(tmpcurve, 2, quantile, probs=.975, na.rm=TRUE)        
            fit_curve_md[i,] <- apply(tmpcurve, 2, quantile, probs=.5, na.rm=TRUE)        
        }
        colnames(fit_curve_md) <- ages_
        colnames(fit_curve_ll) <- ages_
        colnames(fit_curve_ul) <- ages_
        
        fits_long <- fit_curve_md %>% as.data.frame() %>% mutate(survey_id = survey_id_, survey.num=survey.num_) %>%
            gather(key="age", value = "est_md", -c(survey_id, survey.num))
        fits_longul <- fit_curve_ul %>% as.data.frame() %>% mutate(survey_id = survey_id_, survey.num=survey.num_) %>%
            gather(key="age", value = "est_ul", -c(survey_id, survey.num))
        fits_longll <- fit_curve_ll %>% as.data.frame() %>% mutate(survey_id = survey_id_, survey.num=survey.num_) %>%
            gather(key="age", value = "est_ll", -c(survey_id, survey.num))
        fits_ <- bind_cols(bind_cols(fits_long, fits_longll %>% select(est_ll)), fits_longul %>% select(est_ul))
        fits_$scaled.mid.age <- fits_$age
        
        data <- full_join(survey_dat %>% select(country, survey_id, survey.num, scaled.mid.age, seropositive.immune, binom_low, binom_high),
                          fits_ %>% mutate(scaled.mid.age = as.integer(scaled.mid.age),
                                           survey_id = as.character(survey_id)), 
                          by=c("scaled.mid.age", "survey_id" ,"survey.num"))
        data <- data %>% arrange(survey.num, scaled.mid.age)
        data <- left_join(data, seroprev_ %>% dplyr::select(survey_id, age, mean_seroprev, lb_seroprev, ub_seroprev, median_seroprev), by=c("survey_id","scaled.mid.age"="age"))
        
        
        max_age <- ceiling(max(data$scaled.mid.age))
        
        p1 <- ggplot() + 
            geom_ribbon(data=data %>% filter(!is.na(age)), aes(x=scaled.mid.age, ymin=lb_seroprev, ymax=ub_seroprev), fill="orange", alpha=.4) +
            geom_path(data=data %>% filter(!is.na(age)), aes(x=scaled.mid.age, y=mean_seroprev, color="orange"), size=1.1) +
            geom_ribbon(data=data %>% filter(!is.na(age)), aes(x=scaled.mid.age, ymin=est_ll, ymax=est_ul), fill="maroon", alpha=.4) +
            geom_path(data=data %>% filter(!is.na(age)), aes(x=scaled.mid.age, y=est_md, color="maroon"), size=1.1) +
            geom_errorbar(data=data, aes(x=scaled.mid.age, ymin=binom_low, ymax=binom_high, color="navy"), width=0, size=1.25) +
            geom_point(data=data , aes(x=scaled.mid.age, y=seropositive.immune, color="navy"), size=3) +
            theme_classic() +
            # scale_fill_manual(values=c("maroon"="maroon", "orange"="orange"), 
            #                     labels=c("Serology fits","Transmission fits")) +
            scale_colour_manual(name=NULL, values=c("maroon"="maroon","navy"="navy", "orange"="orange"), 
                                labels=c("Serology fits","Survey estimates", "Transmission fits")) +
            coord_cartesian(xlim=c(0, max_age), ylim=c(0,1)) +
            ylab("Seropositive") + xlab("Age (yr)") + 
            theme(legend.position = "bottom") +
            facet_wrap(~survey_id)#, labeller = labeller(survey.num = country.year))
        print(p1)
    }    
    
    dev.off()
    print(paste0("Plots saved to ", file.path(proj_dir, "serology_plot_post_modelsims.pdf"), "."))
}







