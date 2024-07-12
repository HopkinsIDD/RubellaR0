# Get and Use ISO COdes


# Install and load needed packages
if(!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if(!require('RCurl')) install.packages('RCurl'); library(RCurl)
if(!require('tibble')) install.packages('tibble'); library(tibble)
if(!require('countrycode')) install.packages('countrycode'); library(countrycode)
if(!require('rstudioapi')) install.packages('rstudioapi'); library(rstudioapi)

# Get directory of current file (of ISO_code_source.R). We will use this as the directory to save files downloaded locally.
current.wd <- getwd()
print(paste0('Current working directory: ', getwd()))
# source.path <- file.path(rstudioapi::getActiveDocumentContext()$path)
# tmp <- as.integer(gregexpr("/", source.path)[[1]])
# source.path <- substr(source.path,1,tmp[length(tmp)]-1)

# If sourced from other file, this will give this file's path
source.path <- NULL
tryCatch({
  source.path <- dirname(sys.frame(1)$ofile)
}, error= function(x) print('Local Run'))


# If run directly from file, this will give this file's path
if (is.null(source.path)) {
  source.path <- dirname(rstudioapi::getSourceEditorContext()$path)
}

print(source.path)




#####################################################
# Import and Clean the Data

# First Check that you have internet connection and can access github. If not, import local data source
has_internet <- TRUE
tryCatch({
  is.character(getURL("https://raw.githubusercontent.com/datasets/country-codes/master/data/"))==TRUE &
    !is.null(curl::nslookup("www.github.com")) }, 
  error= function(x){
    has_internet <- FALSE
    print('No Internet Connection')
  }
)


# If you have internet and github is accessible, update the data
if (has_internet) {
  iso_data <- region_data <- NULL
  
  # Get ISO Country Codes
  tryCatch({
    iso_data <- read.csv("https://raw.githubusercontent.com/datasets/country-codes/master/data/country-codes.csv", stringsAsFactors=F)
    region_data <- read.csv("https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv", stringsAsFactors=F)
    who_regions <- read.csv('https://raw.githubusercontent.com/shauntruelove/RandomDataAndCode/master/countries_and_regions/data/who_regions.csv', stringsAsFactors=FALSE)
    dhs_countrydata <- read.csv('https://raw.githubusercontent.com/shauntruelove/RandomDataAndCode/master/countries_and_regions/data/DHS_countrydata.csv', stringsAsFactors=FALSE)
    nationality_data <- read.csv('https://raw.githubusercontent.com/shauntruelove/RandomDataAndCode/master/countries_and_regions/data/nationalities_and_languages.csv', stringsAsFactors=FALSE)
  },
  error= function(x) print('No Internet Connection')
  )
  
  dir.create(file.path(source.path, 'iso_source'), showWarnings = FALSE)
  source.path.folder <- file.path(source.path, 'iso_source')
  print(paste0("Data will be save in ", file.path(source.path.folder, ".")))
  
  # Write these data to directory if they are available so they are updated
  write.csv(iso_data, file.path(source.path.folder, 'iso_data.csv'), row.names = FALSE)
  write.csv(region_data, file.path(source.path.folder, 'region_data.csv'), row.names = FALSE)    
  write.csv(who_regions, file.path(source.path.folder, 'who_regions.csv'), row.names = FALSE)
  write.csv(dhs_countrydata, file.path(source.path.folder, 'dhs_countrydata.csv'), row.names = FALSE)
  write.csv(nationality_data, file.path(source.path.folder, 'nationalities_and_languages.csv'), row.names = FALSE)
  
  
  # If not, use the local data
} else {
  source.path.folder <- file.path(source.path, 'iso_source')
  
  if (file.exists(file.path(source.path.folder, 'iso_data.csv'))) {
    iso_data <- read.csv(file.path(source.path.folder, 'iso_data.csv'), header=TRUE, stringsAsFactors=FALSE)
    region_data <- read.csv(file.path(source.path.folder, 'region_data.csv'),  header=TRUE, stringsAsFactors=FALSE)
    who_regions <- read.csv(file.path(source.path.folder, 'who_regions.csv'),  header=TRUE, stringsAsFactors=FALSE)
    dhs_countrydata <- read.csv(file.path(source.path.folder, 'dhs_countrydata.csv'),  header=TRUE, stringsAsFactors=FALSE)
    nationality_data <- read.csv(file.path(source.path.folder, 'nationalities_and_languages.csv'),  header=TRUE, stringsAsFactors=FALSE)
    
  } else {
    print('No access to internet/github and local files could not be located')
    stop()
  }
} 





# Extract needed vectors
ISO2=as.vector(iso_data$ISO3166.1.Alpha.2)
ISO3=as.vector(iso_data$ISO3166.1.Alpha.3)
UNcode=as.vector(iso_data$ISO3166.1.numeric)

Country=as.character(iso_data$official_name_en)
Country2=as.character(iso_data$CLDR.display.name)

Country[grep("Ivoire", Country)] <- "Cote d'Ivoire"
Country2[grep("Ivoire", Country2)] <- "Cote d'Ivoire"
Country[grep("Micronesia (Federated States of)", Country)] <- "Micronesia, Federated States of"

region.row <- as.integer(sapply(X=ISO3, FUN=function(X) which(toupper(region_data$alpha.3)==toupper(X))))
who.row <-as.integer(sapply(X=ISO3, FUN=function(X) which(toupper(who_regions$iso)==toupper(X))))
iso_data_full <- iso_data
iso_data <- tibble(ISO2, ISO3, UNcode, Country, Country2,
                   region_data$region[region.row], region_data$sub.region[region.row], 
                   who_regions$who.region[who.row], who_regions$region.mortality[who.row])
colnames(iso_data) <- c('ISO2', 'ISO3', 'UNcode', 'Country', 'Country2', 'Region', 'Sub.Region', 'WHO.Region', 'WHO.Region.Mortality')

# Fill in missing data
iso_data$Sub.Region[is.na(iso_data$Sub.Region)] <- iso_data_full$Intermediate.Region.Name[is.na(iso_data$Sub.Region)]
rm(who.row, region.row, Country2, Country, ISO3, ISO2)

iso_data[iso_data$Country=='Namibia', 'ISO2'] <- 'NA'
iso_data[iso_data$Country=='Georgia', 'Country2'] <- 'Republic of Georgia'


###################################################
# Functions For Calling the codes

get.country.names.ISO3 <- function(ISO){
  return(as.character(iso_data$Country2[match(toupper(ISO), iso_data$ISO3)]))
}

get.country.names.ISO2 <- function(ISO2){
  return(as.character(iso_data$Country2[match(toupper(ISO2), iso_data$ISO2)]))
}

get.ISO3.from.ISO2 <- function(ISO2){
  return(as.character(iso_data$ISO3[match(toupper(ISO2), iso_data$ISO2)]))
}

get.ISO2.from.ISO3 <- function(ISO3){
  return(as.character(iso_data$ISO2[match(toupper(ISO3), iso_data$ISO3)]))
}

get.country.DHScode <- function(DHS_code){
  return(as.character(dhs_countrydata$CountryName[match(toupper(DHS_code), toupper(dhs_countrydata$DHS_CountryCode))]))
}

get.UNcode <- function(country){
  iso <- get.iso(country, ISO.only=T)
  return(iso_data$UNcode[match(toupper(iso),iso_data$ISO3)])
}

get.UNcode.from.ISO3 <- function(ISO3){
  return(iso_data$UNcode[match(toupper(ISO3), iso_data$ISO3)])
}

# Get ISO3 from Country Name
get.iso <- function(country, ISO.only=T){
  
  # First try "countrycode" package
  ISO <- countrycode(country, 'country.name', 'iso3c')
  
  if (!is.na(ISO)){
    return(ISO)
    
  } else {
    
    country <- tolower(country)
    
    iso.row.tmp <- which(country==tolower(iso_data$Country) | country==tolower(iso_data$Country2))
    if (length(iso.row.tmp)==0){
      iso.row.tmp <- unique(c(grep(country, tolower(iso_data$Country)), grep(country, tolower(iso_data$Country2))))
    }
    
    # Try matching individual words
    country2 <- gsub(' {2,}',' ',country)
    length.country <- length(strsplit(country2,' ')[[1]])
    country2 <- gsub('\\(', '', country2) # Get rid of parentheses
    country2 <- gsub('\\)', '', country2) # Get rid of parentheses
    
    if (length(iso.row.tmp)==0 & length.country>1){
      country.words <- strsplit(country2,' ')[[1]]
      matches <- lapply(country.words, grep, tolower(iso_data$Country)) # match each word of the country name with words in the iso data
      match.row <- Reduce(intersect, matches)  # The match row is identified as the one for which matches of multiple word occurs (intersect)
      
      if (length(match.row)>0){
        iso.row.tmp <- match.row
      } else if (length(match.row)==0 & sum(unlist(matches))>0) {
        iso.row.tmp <- as.integer(matches[which.min(lengths(matches))])  # Match row is the one with the least matches
      }
    }
    
    if (length(iso.row.tmp)==1){
      if (ISO.only){
        return(as.character(iso_data$ISO3[iso.row.tmp]))
      } else{
        return(as.vector(iso_data$ISO3[iso.row.tmp]))
      }
    } else if (length(iso.row.tmp==0)){
      #print('ISO Not Found')
      return('ISO Not Found')
    } else if (length(iso.row.tmp>1)){
      #print(paste0(length(iso.row.tmp),' ISO matches found for ', country))
      return(paste0(length(iso.row.tmp),' ISO matches found for ', country))
    }
  }
}

# Get ISO2 from Country Name
get.iso2 <- function(country, ISO.only=T){
  
  # First try "countrycode" package
  ISO <- countrycode(country, 'country.name', 'iso2c')
  
  if (!is.na(ISO)){
    return(ISO)
    
  } else {
    
    country <- tolower(country)
    
    iso.row.tmp <- which(country==tolower(iso_data$Country) | country==tolower(iso_data$Country2))
    if (length(iso.row.tmp)==0){
      iso.row.tmp <- unique(c(grep(country, tolower(iso_data$Country)), grep(country, tolower(iso_data$Country2))))
    }
    
    # Try matching individual words
    country2 <- gsub(' {2,}',' ',country)
    length.country <- length(strsplit(country2,' ')[[1]])
    country2 <- gsub('\\(', '', country2) # Get rid of parentheses
    country2 <- gsub('\\)', '', country2) # Get rid of parentheses
    
    if (length(iso.row.tmp)==0 & length.country>1){
      country.words <- strsplit(country2,' ')[[1]]
      matches <- lapply(country.words, grep, tolower(iso_data$Country)) # match each word of the country name with words in the iso data
      match.row <- Reduce(intersect, matches)  # The match row is identified as the one for which matches of multiple word occurs (intersect)
      
      if (length(match.row)>0){
        iso.row.tmp <- match.row
      } else if (length(match.row)==0 & sum(unlist(matches))>0) {
        iso.row.tmp <- as.integer(matches[which.min(lengths(matches))])  # Match row is the one with the least matches
      }
    }
    
    if (length(iso.row.tmp)==1){
      if (ISO.only){
        return(as.character(iso_data$ISO2[iso.row.tmp]))
      } else{
        return(as.vector(iso_data$ISO2[iso.row.tmp]))
      }
    } else if (length(iso.row.tmp==0)){
      #print('ISO Not Found')
      return('ISO Not Found')
    } else if (length(iso.row.tmp>1)){
      #print(paste0(length(iso.row.tmp),' ISO matches found for ', country))
      return(paste0(length(iso.row.tmp),' ISO matches found for ', country))
    }
  }
  
}


# Get WHO Regions using ISO
get.who.region <- function(ISO){
  # Get WHO Regions using ISO
  return(as.character(who_regions$who.region[match(toupper(ISO), who_regions$iso)]))
}


# Get Other Regions using ISO
get.region <- function(ISO){
  # Get Regions using ISO
  return(as.character(iso_data$Region[match(toupper(ISO), iso_data$ISO3)]))
}
# Get Other Regions using ISO
get.subregion <- function(ISO){
  # Get Sub Regions using ISO
  return(as.character(iso_data$Sub.Region[match(toupper(ISO), iso_data$ISO3)]))
}






# Get ISO from nationality

# Get ISO3 from Country Name
get.iso.from.nationality <- function(nationality, ISO.only=T){
    
    nationality <- tolower(nationality)
    
    iso.row.tmp <- which(grep(nationality, tolower(nationality_data$Nationalities)) | nationality==tolower(nationality_data$Country))
    if (length(iso.row.tmp)==0){
        iso.row.tmp <- unique(c(grep(country, tolower(iso_data$Country)), grep(country, tolower(iso_data$Country2))))
    }
    
    # Try matching individual words
    country2 <- gsub(' {2,}',' ',country)
    length.country <- length(strsplit(country2,' ')[[1]])
    country2 <- gsub('\\(', '', country2) # Get rid of parentheses
    country2 <- gsub('\\)', '', country2) # Get rid of parentheses
    
    if (length(iso.row.tmp)==0 & length.country>1){
        country.words <- strsplit(country2,' ')[[1]]
        matches <- lapply(country.words, grep, tolower(iso_data$Country)) # match each word of the country name with words in the iso data
        match.row <- Reduce(intersect, matches)  # The match row is identified as the one for which matches of multiple word occurs (intersect)
        
        if (length(match.row)>0){
            iso.row.tmp <- match.row
        } else if (length(match.row)==0 & sum(unlist(matches))>0) {
            iso.row.tmp <- as.integer(matches[which.min(lengths(matches))])  # Match row is the one with the least matches
        }
    }
    
    if (length(iso.row.tmp)==1){
        if (ISO.only){
            return(as.character(iso_data$ISO3[iso.row.tmp]))
        } else{
            return(as.vector(iso_data$ISO3[iso.row.tmp]))
        }
    } else if (length(iso.row.tmp==0)){
        #print('ISO Not Found')
        return('ISO Not Found')
    } else if (length(iso.row.tmp>1)){
        #print(paste0(length(iso.row.tmp),' ISO matches found for ', country))
        return(paste0(length(iso.row.tmp),' ISO matches found for ', country))
    }
}

