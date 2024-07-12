# merge R0 results and serosurvey priority data

library(tidyverse)
library(globaltoolbox)

sero_priorities <-read_csv("results/serosurvey_priorities.csv")
R0_and_sero <- read_csv("results/R0_and_serosurveys.csv")

sero_priorities <- sero_priorities %>% mutate(ISO = get.iso(country)) %>% 
  mutate(ISO = ifelse(country=="Kosovo", "RKS", ISO))
mergeddata <- full_join(sero_priorities, R0_and_sero, by=c("ISO"="ISO")) %>% as.data.frame() 
mergeddata <- mergeddata %>% select(country, Country, ISO, sub.region, Region, VIMC_98, Priority, everything())

mergeddata <- mergeddata %>% select(-sub.region, -country) %>%
  mutate(Priority = ifelse(VIMC_98==FALSE, "NONE", Priority))
mergeddata <- mergeddata %>% mutate(WHO_region = get.who.region(ISO)) %>%
  select(Country, ISO, WHO_region, everything())

# Save merged data
write_csv(mergeddata, "results/merged_R0_and_surveyinfo.csv")
