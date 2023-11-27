#code for Kamazima

# Mon Nov 27 11:41:16 2023 ------------------------------

#load packages
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

#read in data
fish<-read.csv(here("data/Barcoding_Results_Full.csv"))
stations<-read.csv(here("data/Station_Info.csv"))

#rename columns, change station ID to factor
fish<- fish %>% rename(Stage = X, Cruise_ID = Cruise.ID, Station_ID = Station.ID) %>% 
  mutate(Station_ID = as.factor(Station_ID))

stations<-stations %>% mutate(Station_ID = as.factor(Station_ID))

#filter results from May 2021, eggs only, remove ones with no data
egg_2105<- fish %>% dplyr::select("Genus.species","Common_Name","Cruise_ID","Station_ID","Stage") %>% 
  filter(Cruise_ID == 2105 & Stage == "Egg" & Common_Name != "")

#tally fish by common name
egg_2105_tally<-egg_2105 %>% group_by(Common_Name, Cruise_ID,Station_ID) %>% tally()

#merge with station data, calculate density
egg_2105_stations<-left_join(egg_2105_tally,stations,by=c("Cruise_ID","Station_ID")) %>% 
  mutate(Density = n/Tow_Depth*0.3^2*pi)

#save
write.csv(egg_2105_stations,"2105_eggs.csv")
