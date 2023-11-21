# Sat Nov 18 15:27:11 2023 ------------------------------

#load packages
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)

#read in data
fish<-read.csv(here("data/Barcoding_Results_Full.csv"))
stations<-read.csv(here("data/Station_Info.csv"))

#convert Cruise_ID, Station_IN to factor
stations$Cruise_ID<-as.factor(stations$Cruise_ID)
fish$Cruise.ID<-as.factor(fish$Cruise.ID)
stations$Station_ID<-as.factor(stations$Station_ID)
fish$Station.ID<-as.factor(fish$Station.ID)

#remove all columns but genus, species, common name, station, type
fish_filter<- fish %>% dplyr::select(c("Genus.species","Common_Name","Cruise.ID","Station.ID","X", "SST")) %>% 
  na.omit() %>% rename("Cruise_ID" = "Cruise.ID") %>% rename("Station_ID"="Station.ID") %>%
  rename("Stage" = "X")
  

#filter out 2103, juveniles
fish_filter<-fish_filter %>% filter(Cruise_ID != "2103") %>% filter(Stage!="Juvenile")

#silver hake eggs
silver<-fish_filter %>% filter(Common_Name == "Silver Hake" & Stage == "Egg") 
#only in 2105 to date

#plot temp data - all
stations %>%
  ggplot(aes(x=SST, color=Cruise_ID, fill=Cruise_ID)) +
  geom_histogram(binwidth = 0.5, alpha=0.5)

#filter stations already processed
stations_filter <- stations %>% filter (Cruise_ID %in% c(2105,2202,2205,2207)) 
stations_filter <- stations_filter[1:76,]
#plot
stations_filter %>%
  ggplot(aes(x=SST, color=Cruise_ID, fill=Cruise_ID)) +
  geom_histogram(binwidth = 0.5, alpha=0.5)

silver_tally<-silver %>% group_by(Cruise_ID,Station_ID) %>% tally()

silver_tally <- left_join(stations_filter,silver_tally,by=c("Cruise_ID","Station_ID")) %>%
  rename("SilverHake" = "n") %>% mutate(SilverHake = ifelse(is.na(SilverHake),0,SilverHake))

ggplot() +
  geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
  geom_histogram(data=silver,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")


#what about GS Flounder?
gs_egg<-fish_filter %>% filter(Common_Name == "Gulfstream Flounder" & Stage == "Egg")
gs_larva<-fish_filter %>% filter(Common_Name == "Gulfstream Flounder" & Stage == "Larva")
ggplot() +
  geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
  geom_histogram(data=gs_egg,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")

#what about Mackerel?
mackerel<-fish_filter %>% filter(Common_Name == "Atlantic Mackerel" & Stage == "Egg")
ggplot() +
  geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
  geom_histogram(data=mackerel,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")


