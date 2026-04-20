#-----------------------------------------------------------------#
##################    EXPLORATIVE DATA ANALYSIS   #################
#-----------------------------------------------------------------#


# R packages --------------------------------------------------------------
library(tidyverse)
library(data.table)
library(ggplot2)
library(readr)
library(gridExtra)
library(xtable)
library(latex2exp)

# Pollutant data ----------------------------------------------------------
data = read.csv("data/daily_SPEC_2021.csv")

data = data[data$State.Name == "California", ]
unique(data$Local.Site.Name)
unique(data$Parameter.Name)


# Data polishing ----------------------------------------------------------
# Some stations include only meteorological data
# Keep only stations with air pollutants:
unique(data$Parameter.Name)
AQdata = data[data$Parameter.Name %like% "PM2.5", ]
length(unique(AQdata$Local.Site.Name)) 

unique(AQdata$Parameter.Name)
# 71 different pollutants, but we focus only on a subset 


graph_data = AQdata[ ,c(12, 23, 9, 17)]
graph_data$Date.Local = as.Date(graph_data$Date.Local)

graph_data %>% pivot_wider(id_cols=c(Date.Local, Local.Site.Name), 
                           names_from = Parameter.Name, 
                           values_from = Arithmetic.Mean)
graph_data %>% dplyr::group_by(Date.Local, Local.Site.Name, Parameter.Name) %>% 
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>% 
  dplyr::filter(n > 1L) 

weird = AQdata[AQdata$Date.Local=="2021-01-01" & AQdata$Local.Site.Name == "", ]
# some double obs: one with "event included" and the other with "none" 
table(AQdata$Event.Type)
# We cannot exclude either all "included" or "none" events 
# --> we pick an average of the two values recorded
summary_data = graph_data %>% pivot_wider(id_cols=c(Date.Local, Local.Site.Name), 
                                          names_from = Parameter.Name, 
                                          values_from = Arithmetic.Mean,
                                          values_fn = mean)

#We only focus on: sulfur, aluminium, lead, total nitrate, sulfate
new_data = summary_data[ , c("Date.Local", "Local.Site.Name", "Aluminum PM2.5 LC", "Sulfur PM2.5 LC", "Total Nitrate PM2.5 LC", "Sulfate PM2.5 LC",
                         "OC PM2.5 LC TOR", "EC PM2.5 LC TOR")]




# New data for first approach
missing_per_staz = tapply(rowSums(is.na(new_data[ ,-c(1,2)])), new_data$Local.Site.Name, sum)
sort(missing_per_staz)
# keeping only stations with less than 100 NAs
new_data = new_data[which(new_data$Local.Site.Name %in% names(missing_per_staz)[missing_per_staz<100]), ]
new_data = new_data[-which(new_data$Local.Site.Name==""), ]

selected_staz = names(tapply(new_data$Date.Local, new_data$Local.Site.Name, length)[tapply(new_data$Date.Local, new_data$Local.Site.Name, length)>100])
new_data = new_data[new_data$Local.Site.Name %in% selected_staz, ]
new_data$t = yday(new_data$Date.Local)

write_csv(new_data, "realData_subset.csv")


# Stations EDA ------------------------------------------------------------
site_data = read_csv("data/aqs_sites.csv")
site_data = site_data[site_data$`State Name`=="California", ]
new_sites = site_data[site_data$`Local Site Name`%in% unique(new_data$Local.Site.Name), c(4:5,7:9,20,24:25)]
names(new_sites)[4:5] = c("LandUse","LandSetting")

write_csv(new_sites, "FinalStations.csv")