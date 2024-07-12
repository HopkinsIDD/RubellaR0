

#Make a WAIFW matrix based on Prem et al. 2017
#using pakistan for afghanistan
#
#Parameters -
#   age class boundaries - the upper age limit for each age class in YEARS
#   uncode - country associated with this code must match the sheet names in the excel spreadhsheet from prem et al.
#   bandwidth - desired smooth bandwidth - default=c(3,3)
#Returns -
#   a WAIFW matrix based on the Polymod results from chosen location with row and col
#   names indicating age classes
get.prem.WAIFW <- function (age.class.boundaries = (1:90),
                            uncode, other.contact.matrix=F,
                            bandwidth=c(3,3), min_age_bound=0) {
    
    #get country name from the uncode
    cc <- read.csv("data/country_codes.csv")
    name <- as.character(cc$Report_country_name[which(cc$uncode==uncode & !is.na(cc$uncode))])
    if (name=="Democratic Republic of the Congo" | name=="Central African Republic") name <- "Congo"
    if (name=="Afghanistan") name <- "Pakistan"
    if (name=="Angola") name <- "Namibia"
    if (name=="Kosovo") name <- "Albania"
    if (name=="Tuvalu") name <- "Tonga"
    if (name=="Turkmenistan") name <- "Iran (Islamic Republic of)"
    if (name=="Togo") name <- "Benin"
    if (name=="Chad") name <- "Niger"
    if (name=="Swaziland") name <- "Lesotho"
    if (name=="Sao Tome and Principe") name <- "Sao Tome and Principe "
    if (name=="South Sudan" | name=="Sudan" | name=="Eritrea" | name=="Djibouti") name <- "Ethiopia"
    if (name=="Somalia") name <- "Ethiopia"
    if (name=="occupied Palestinian territory") name <- "Jordan"
    if (name=="Democratic Peoples Republic of Korea") name <- "Republic of Korea"
    if (name=="Papua New Guinea" | name=="Marshall Islands" | name=="Micronesia (Federated States of )") name <- "Solomon Islands"
    if (name=="Malawi") name <- "Zambia"
    if (name=="Myanmar") name <- "Thailand"
    if (name=="Mali") name <- "Burkina Faso"
    if (name=="Madagascar" | name=="Comoros") name <- "Mozambique"
    if (name=="Republic of Moldova") name <- "Romania"
    if (name=="Lao Peoples Democratic Republic") name <- "Lao People’s Democratic Republi"
    if (name=="Guinea-Bissau" | name=="Gambia") name <- "Guinea"
    if (name=="Cuba") name <- "Jamaica"
    if (name=="Cote d Ivoire") name <-"Ghana"
    if (name=="Bolivia") name <- "Bolivia (Plurinational State of"
    if (name=="Burundi") name <- "Rwanda"
    if (name=="Angola") name <- "Namibia"
    if (name=="Lao People Democratic Republic") name <- "Lao People’s Democratic Republi"
    #bring in the polymod data on contacts for UK
    if(substr(name, 1,3)<"Mos") df.prem <- readxl::read_excel("data/prem_contact_matrices/MUestimates_all_locations_1.xls", sheet=name, col_names=F)
    if(substr(name, 1,3)>"Mos") df.prem <- readxl::read_excel("data/prem_contact_matrices/MUestimates_all_locations_2.xls", sheet=name, col_names=F)
    if (nrow(df.prem)==17) df.prem <- df.prem[-1,]
    df.prem <- t(apply(df.prem, 1, function(x) as.numeric(x)))
    #turn data into long form (2 columns of ages)
    prem <- round(df.prem*1000)
    prem.v <- as.numeric(unlist(c(prem)))
    if (prem.v[length(prem.v)]==0) prem.v[length(prem.v)] <- 1 #hack, code breaks if end on a zero.
    prem.sum <- as.numeric(cumsum(prem.v))
    prem.ages <- seq(3, 78, 5)
    x <- matrix(NA, sum(prem.v), 2)
    for (c in 1:length(prem.v)){
        index <- min(which(is.na(x[,1])))
        index2 <- prem.sum[c]
        x[(index:index2), 1] <- prem.ages[ceiling(c/16)]
        if (!any(c==seq(16,256,16))) x[(index:index2), 2] <- prem.ages[c-(floor(c/16)*16)] #prem.ages[c-(floor(c/16)*16)+1]
        if (c%in%seq(16,256,16)) x[(index:index2), 2] <- prem.ages[(floor(c/16)*16)-(16*(c/16-1))]
    }
    if (other.contact.matrix){
        if (name=="China") name <- "China_Read2014"
        if(substr(name, 1,3)<"Mos") df.prem <- readxl::read_excel("data/prem_contact_matrices/MUestimates_all_locations_1.xls", sheet=name)
        if(substr(name, 1,3)>"Mos") df.prem <- readxl::read_excel("data/prem_contact_matrices/MUestimates_all_locations_2.xls", sheet=name)
        prem <- round(df.prem*100)
        prem.v <- as.numeric(unlist(c(prem)))
        if (prem.v[length(prem.v)]==0) prem.v[length(prem.v)] <- 1 #hack, code breaks if end on a zero.
        prem.sum <- as.numeric(cumsum(prem.v))
        ages.min <- c(1,6,20,65)
        ages.max <- c(5,19,64,80)
        x <- matrix(NA, sum(prem.v), 2)
        for (c in 1:length(prem.v)){
            index <- min(which(is.na(x[,1])))
            index2 <- prem.sum[c]
            x[(index:index2), 1] <- ages.min[c]
            x[(index:index2), 2] <- ages.max[c]
        }
    }
    #smooth monthly from 0 to 91 yrs of age and get the corresponding ages
    est <- KernSmooth::bkde2D(x, bandwidth=bandwidth, gridsize=c(12*101, 12*101), range.x=list(c((1/24),100),c((1/24),100)))
    ages.polymod.smooth <- est$x1
    #image(ages.polymod.smooth,ages.polymod.smooth,est$fhat,xlim=c(0,10),ylim=c(0,10))
    #find which fitted ages (ages.polymod.smooth) each of the desired lower boundaries is between
    n.age.cats <- length(age.class.boundaries)
    lowpoints <-  c(min_age_bound,age.class.boundaries[2:n.age.cats-1])
    index <- (findInterval(lowpoints,ages.polymod.smooth));
    #set very large ages to all be the same as the largest age
    index[index>=length(ages.polymod.smooth)]=length(ages.polymod.smooth)-1
    #extract appropriate matrix, taking smoothed estimates for the ranges
    foi.matrix <- est$fhat[index+1,index+1]
    #levelplot(foi.matrix*100)
    #adjust for fact that the width of your age class should not affect the number of contacts you make
    #foi.matrix <- foi.matrix/diff(c(0,age.class.boundaries))
    colnames(foi.matrix) <- age.class.boundaries
    rownames(foi.matrix) <- age.class.boundaries
    return(foi.matrix)
}






#Make a WAIFW matrix based on Prem et al. 2017 -- keep it scaled to original scale (daily average number of contacts)
#using pakistan for afghanistan
#
#Parameters -
#   age class boundaries - the upper age limit for each age class in YEARS
#   uncode - country associated with this code must match the sheet names in the excel spreadhsheet from prem et al.
#   bandwidth - desired smooth bandwidth - default=c(3,3)
#Returns -
#   a WAIFW matrix based on the Polymod results from chosen location with row and col
#   names indicating age classes
get.prem.WAIFW.daily <- function (age.class.boundaries = (1:90),
                            uncode, other.contact.matrix=F,
                            bandwidth=c(3,3), min_age_bound=0) {
    
    #get country name from the uncode
    cc <- read.csv("data/country_codes.csv")
    name <- as.character(cc$Report_country_name[which(cc$uncode==uncode & !is.na(cc$uncode))])
    if (name=="Democratic Republic of the Congo" | name=="Central African Republic") name <- "Congo"
    if (name=="Afghanistan") name <- "Pakistan"
    if (name=="Angola") name <- "Namibia"
    if (name=="Kosovo") name <- "Albania"
    if (name=="Tuvalu") name <- "Tonga"
    if (name=="Turkmenistan") name <- "Iran (Islamic Republic of)"
    if (name=="Togo") name <- "Benin"
    if (name=="Chad") name <- "Niger"
    if (name=="Swaziland") name <- "Lesotho"
    if (name=="Sao Tome and Principe") name <- "Sao Tome and Principe "
    if (name=="South Sudan" | name=="Sudan" | name=="Eritrea" | name=="Djibouti") name <- "Ethiopia"
    if (name=="Somalia") name <- "Ethiopia"
    if (name=="occupied Palestinian territory") name <- "Jordan"
    if (name=="Democratic Peoples Republic of Korea") name <- "Republic of Korea"
    if (name=="Papua New Guinea" | name=="Marshall Islands" | name=="Micronesia (Federated States of )") name <- "Solomon Islands"
    if (name=="Malawi") name <- "Zambia"
    if (name=="Myanmar") name <- "Thailand"
    if (name=="Mali") name <- "Burkina Faso"
    if (name=="Madagascar" | name=="Comoros") name <- "Mozambique"
    if (name=="Republic of Moldova") name <- "Romania"
    if (name=="Lao Peoples Democratic Republic") name <- "Lao People’s Democratic Republi"
    if (name=="Guinea-Bissau" | name=="Gambia") name <- "Guinea"
    if (name=="Cuba") name <- "Jamaica"
    if (name=="Cote d Ivoire") name <-"Ghana"
    if (name=="Bolivia") name <- "Bolivia (Plurinational State of"
    if (name=="Burundi") name <- "Rwanda"
    if (name=="Angola") name <- "Namibia"
    if (name=="Lao People Democratic Republic") name <- "Lao People’s Democratic Republi"
    #bring in the polymod data on contacts for UK
    if(substr(name, 1,3)<"Mos") df.prem <- readxl::read_excel("data/prem_contact_matrices/MUestimates_all_locations_1.xls", sheet=name, col_names=F)
    if(substr(name, 1,3)>"Mos") df.prem <- readxl::read_excel("data/prem_contact_matrices/MUestimates_all_locations_2.xls", sheet=name, col_names=F)
    if (nrow(df.prem)==17) df.prem <- df.prem[-1,]
    df.prem <- t(apply(df.prem, 1, function(x) as.numeric(x)))
    prem.ages <- seq(3, 78, 5)

    
    df.prem <- df.prem[1:4,1:4]
    prem.ages <- prem.ages[1:4]
    
    #turn data into long form (2 columns of ages)
    prem <- round(df.prem*1000)
    prem.age.num <- length(prem.ages)
    
    
    
    # for scaling back to absolute contacts per day, get total number of contacts by age
    n_contacts_age <- rowSums(df.prem) / 5    # prem age contacts are per 5 years of age
    n_contacts_age <- data.frame(age=prem.ages, n_contacts=n_contacts_age)
    contacts_spline <- smooth.spline(x=n_contacts_age$age, y=n_contacts_age$n_contacts , spar=.3) 
    est_n_conts <- predict(contacts_spline, x=0.5:max(age.class.boundaries))
    # plot(n_contacts_age, type="p")
    # lines(contacts_spline)
    # lines(est_n_conts, col="red")
    
    
    prem.v <- as.numeric(unlist(c(prem)))
    if (prem.v[length(prem.v)]==0) prem.v[length(prem.v)] <- 1 #hack, code breaks if end on a zero.
    #if (prem.v[length(prem.v)]==0) prem.v[length(prem.v)] <- 1000 #hack, code breaks if end on a zero.
    
    # this turns each contact between ages i,j into a row 
    prem.sum <- as.numeric(cumsum(prem.v))
    x <- matrix(NA, sum(prem.v), 2)
    for (c in 1:length(prem.v)){
        index <- min(which(is.na(x[,1])))
        index2 <- prem.sum[c]
        x[(index:index2), 1] <- prem.ages[ceiling(c/prem.age.num)]
        x[(index:index2), 2] <- prem.ages[c-(floor(c/prem.age.num)*prem.age.num)+1]
    }
    
    
    if (other.contact.matrix){
        if (name=="China") name <- "China_Read2014"
        if(substr(name, 1,3)<"Mos") df.prem <- readxl::read_excel("data/prem_contact_matrices/MUestimates_all_locations_1.xls", sheet=name)
        if(substr(name, 1,3)>"Mos") df.prem <- readxl::read_excel("data/prem_contact_matrices/MUestimates_all_locations_2.xls", sheet=name)
        prem <- round(df.prem*100)
        prem.v <- as.numeric(unlist(c(prem)))
        if (prem.v[length(prem.v)]==0) prem.v[length(prem.v)] <- 1 #hack, code breaks if end on a zero.
        prem.sum <- as.numeric(cumsum(prem.v))
        ages.min <- c(1,6,20,65)
        ages.max <- c(5,19,64,80)
        x <- matrix(NA, sum(prem.v), 2)
        for (c in 1:length(prem.v)){
            index <- min(which(is.na(x[,1])))
            index2 <- prem.sum[c]
            x[(index:index2), 1] <- ages.min[c]
            x[(index:index2), 2] <- ages.max[c]
        }
    }
    
    
    #smooth monthly from 0 to 91 yrs of age and get the corresponding ages
    est <- KernSmooth::bkde2D(x, bandwidth=bandwidth, gridsize=c(12*101, 12*101), range.x=list(c((1/24),100),c((1/24),100)))
    
    est <- KernSmooth::bkde2D(x, bandwidth=bandwidth, gridsize=c(101, 101), range.x=list(c((1/2),100),c((1/2),100)))
    
    ages.polymod.smooth <- est$x1
    # contour(est$x1, est$x2, est$fhat) #Plot contour plots
    # persp(est$fhat) # takes forever, dont run
    #image(ages.polymod.smooth,ages.polymod.smooth,est$fhat,xlim=c(0,10),ylim=c(0,10))
    
    #find which fitted ages (ages.polymod.smooth) each of the desired lower boundaries is between
    n.age.cats <- length(age.class.boundaries)
    #lowpoints <-  c(min_age_bound,age.class.boundaries[2:n.age.cats-1])
    lowpoints <-  c(min_age_bound,age.class.boundaries[2:n.age.cats-1])
    max_age <- max(age.class.boundaries)
    index <- (findInterval(c(lowpoints,max_age),ages.polymod.smooth));

    #set very large ages to all be the same as the largest age
    index[index>=length(ages.polymod.smooth)]=length(ages.polymod.smooth)-1
    
    
    #extract appropriate matrix, taking smoothed estimates for the ranges
    #foi.matrix <- est$fhat[index+1,index+1] # old method --> does not appropriately sum the age groups (particularly problem when large)
    foi.matrix <- matrix(NA, nrow=n.age.cats, ncol = n.age.cats)
    for(i in 1:n.age.cats){
        for (j in 1:n.age.cats){
            foi.matrix[i,j] <- sum(est$fhat[(index[i]+1):(index[i+1]+1),(index[j]+1):(index[j+1]+1)])
        }
    }
    #levelplot(foi.matrix*100)
    #adjust for fact that the width of your age class should not affect the number of contacts you make
    #foi.matrix <- foi.matrix/diff(c(0,age.class.boundaries))
    colnames(foi.matrix) <- age.class.boundaries
    rownames(foi.matrix) <- age.class.boundaries
    
    tot_contacts <- sum(est_n_conts$y)
    foi.matrix <- foi.matrix * tot_contacts
    
    return(foi.matrix)
}




#Make a WAIFW matrix based on Prem et al. 2017 -- RAW DATA AND AGE GROUPS
#   uncode - country associated with this code must match the sheet names in the excel spreadhsheet from prem et al.
#   bandwidth - desired smooth bandwidth - default=c(3,3)
#Returns -
#   Prem Data for 5 year age groups 0-80 years
get.prem.WAIFW.RAW <- function (uncode, other.contact.matrix=F) {
    
    #get country name from the uncode
    cc <- read.csv("data/country_codes.csv")
    name <- as.character(cc$Report_country_name[which(cc$uncode==uncode & !is.na(cc$uncode))])
    if (name=="Democratic Republic of the Congo" | name=="Central African Republic") name <- "Congo"
    if (name=="Afghanistan") name <- "Pakistan"
    if (name=="Angola") name <- "Namibia"
    if (name=="Kosovo") name <- "Albania"
    if (name=="Tuvalu") name <- "Tonga"
    if (name=="Turkmenistan") name <- "Iran (Islamic Republic of)"
    if (name=="Togo") name <- "Benin"
    if (name=="Chad") name <- "Niger"
    if (name=="Swaziland") name <- "Lesotho"
    if (name=="Sao Tome and Principe") name <- "Sao Tome and Principe "
    if (name=="South Sudan" | name=="Sudan" | name=="Eritrea" | name=="Djibouti") name <- "Ethiopia"
    if (name=="Somalia") name <- "Ethiopia"
    if (name=="occupied Palestinian territory") name <- "Jordan"
    if (name=="Democratic Peoples Republic of Korea") name <- "Republic of Korea"
    if (name=="Papua New Guinea" | name=="Marshall Islands" | name=="Micronesia (Federated States of )") name <- "Solomon Islands"
    if (name=="Malawi") name <- "Zambia"
    if (name=="Myanmar") name <- "Thailand"
    if (name=="Mali") name <- "Burkina Faso"
    if (name=="Madagascar" | name=="Comoros") name <- "Mozambique"
    if (name=="Republic of Moldova") name <- "Romania"
    if (name=="Lao Peoples Democratic Republic") name <- "Lao People’s Democratic Republi"
    if (name=="Guinea-Bissau" | name=="Gambia") name <- "Guinea"
    if (name=="Cuba") name <- "Jamaica"
    if (name=="Cote d Ivoire") name <-"Ghana"
    if (name=="Bolivia") name <- "Bolivia (Plurinational State of"
    if (name=="Burundi") name <- "Rwanda"
    if (name=="Angola") name <- "Namibia"
    if (name=="Lao People Democratic Republic") name <- "Lao People’s Democratic Republi"
    #bring in the polymod data on contacts for UK
    if(substr(name, 1,3)<"Mos") df.prem <- readxl::read_excel("data/prem_contact_matrices/MUestimates_all_locations_1.xls", sheet=name, col_names=F)
    if(substr(name, 1,3)>"Mos") df.prem <- readxl::read_excel("data/prem_contact_matrices/MUestimates_all_locations_2.xls", sheet=name, col_names=F)
    if (nrow(df.prem)==17) df.prem <- df.prem[-1,]
    df.prem <- t(apply(df.prem, 1, function(x) as.numeric(x)))
    
    return(df.prem)
}


#Make a WAIFW matrix based on Prem et al. 2021 -- RAW DATA AND AGE GROUPS
#   iso3code - iso3code
#Returns -
#   Prem Data for 5 year age groups 0-80 years
get.prem2021.WAIFW.RAW <- function (iso3code) {
    
    load("data/prem_age_contacts_2021/contact_all.rdata")
    ## https://github.com/kieshaprem/synthetic-contact-matrices/tree/master/output/syntheticcontactmatrices2020/overall
    
    index <- which(names(contact_all)==iso3code)
    if (length(index)==0) {
        print(paste("missing contact matrix for iso3code:", iso3code))
    } else {
        df.prem <- contact_all[[index]]
    }
    
    return(df.prem)
}






