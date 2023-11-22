#Code to calculate thermal niches of fish eggs 
#from EcoMon survey program & barcoding results reported in Lewis et al (2016)

# Tue Nov 21 09:55:08 2023 ------------------------------

#load packages
library(readxl)
library(dplyr)
library(tidyr)
library(MASS)
library(here)
library(ggplot2)

#Load data
#Load station data from subset of stations with ethanol-preserved samples
db <- read_excel("data/FishEggDatabase_GenID.accdb.xls",sheet = "StationData")
#remove rows with NA egg count
db <- db %>% filter(is.na(Egg_Count) == F)
#convert egg counts to numeric, station to factor
db$Egg_Count <- as.numeric(db$Egg_Count)
db$Station <- as.double(db$Station)
#rename columns to match ecomon
db <- db %>% rename(cruise_name = Cruise_Name, station = Station)

#Load all EcoMon station data
ecomon <- read.csv("data/EcoMon_Plankton_Data.csv")
ecomon$station<-as.double(ecomon$station)
#join datasets, select for relevant columns
db_trim<- left_join(db,ecomon,by=c("cruise_name","station")) %>% dplyr::select(cruise_name,station,sfc_temp,date,Egg_Count)

#histogram of SST of all stations
db_trim %>% 
  ggplot(aes(x=sfc_temp, color=cruise_name, fill=cruise_name,binwidth = 0.5))+
  geom_histogram(alpha=0.5)

#load in DNA barcoding results
eggs<- read_excel("data/FishEggDatabase_GenID.accdb.xls", sheet = "EggIDData")
#rename eggs columns
eggs<- eggs %>% rename(Sample_ID=`Sample ID`, Common_Name = `Common Name`,cruise_name=Cruise, station=Station)
eggs_trim<- eggs %>% dplyr::select(Sample_ID,cruise_name,station,Common_Name)
eggs_trim$station<-as.double(eggs_trim$station)

#join egg data with temperature data
egg_temp<-left_join(eggs_trim,db_trim, by = c("cruise_name","station"))

stations_ecomon<-egg_cruise %>% dplyr::select(cruise_name,station,sfc_temp) %>% distinct()
stations_all<-bind_rows(stations,stations_ecomon) 
ggplot(data = stations_all,aes(x=sfc_temp))+
  geom_histogram(alpha=0.5,fill="magenta",color="magenta")

#Look at one species = silver hake
silver<- egg_temp %>% filter(Common_Name == "Silver hake")

ggplot() +
  geom_histogram(data=stations_ecomon,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=silver,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,color="black")

#what might niche breadth be?
#establish a column that is rounded to nearest 0.5degC
db_trim <- db_trim %>% mutate(temp_bin = round(sfc_temp*2)/2)
#tally, calculate proportion
temp_dist <- db_trim %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_samples=n/total) %>% dplyr::select(temp_bin,prop_samples)

#apply to silver
silver <- silver %>% mutate(temp_bin = round(sfc_temp*2)/2)
silver_dist <- silver %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_fish=n/total) %>% dplyr::select(temp_bin,prop_fish)

silver_niche<-left_join(temp_dist,silver_dist,by="temp_bin") %>% 
  mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
  mutate (niche = sqrt(prop_samples*prop_fish))

silver_NB <- sum(silver_niche$niche)

#apply to Gulf stream flounder
gs_flounder<-egg_temp %>% filter(Common_Name == "Gulf Stream flounder")
ggplot() +
  geom_histogram(data=stations_ecomon,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=gs_flounder,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,color="black")


gs_flounder <- gs_flounder %>% mutate(temp_bin = round(sfc_temp*2)/2)
gs_flounder_dist <- gs_flounder %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_fish=n/total) %>% dplyr::select(temp_bin,prop_fish)

gs_flounder_niche<-left_join(temp_dist,gs_flounder_dist,by="temp_bin") %>% 
  mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
  mutate (niche = sqrt(prop_samples*prop_fish))

gs_flounder_NB <- sum(gs_flounder_niche$niche)

#cod
cod<- egg_cruise %>% filter(Common_Name == "Atlantic cod")
ggplot() +
  geom_histogram(data=stations_ecomon,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=cod,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,color="black")

#Haddock
haddock<-subset(eggs.trim,eggs.trim$`Common Name` == "Haddock")
hist(haddock$temp,main="Haddock (n=44)", xlab = "Temp",prob=T,xlim = c(0,30))

#Red hake
red<-subset(eggs.trim,eggs.trim$`Common Name` == "Red hake")
hist(red$temp,main="Red Hake (n=139)", xlab = "Temp",prob=T,xlim = c(0,30))

#yellowtail
yellowtail<-subset(eggs.trim,eggs.trim$`Common Name` == "Yellowtail flounder")
hist(yellowtail$temp, main="Yellowtail Flounder (n=74)", xlab = "Temp",prob=T,xlim = c(0,30))



#Fourspot flounder
FS.flounder<-subset(eggs.trim,eggs.trim$`Common Name` == "American fourspot flounder")
hist(FS.flounder$temp,main="Fourspot Flounder (n=182)", xlab = "Temp",prob=T,xlim = c(0,30))

#Windowpane
windowpane<-subset(eggs.trim,eggs.trim$`Common Name` == "Windowpane")
hist(windowpane$temp,main="Windowpane (n=79)", xlab = "Temp",prob=T,xlim = c(0,30))

#Butterfish
butterfish<-subset(eggs.trim,eggs.trim$`Common Name` == "Atlantic butterfish")
hist(butterfish$temp,main="Butterfish (n=52)", xlab = "Temp",prob=T,xlim = c(0,30))

#Summer flounder
summer<-subset(eggs.trim,eggs.trim$`Common Name` == "Summer flounder")
hist(summer$temp,main="Summer Flounder (n=52)", xlab = "Temp",prob=T,xlim = c(0,30))

#Bay anchovy
bay.anchovy<-subset(eggs.trim,eggs.trim$`Common Name` == "Bay anchovy")
hist(bay.anchovy$temp,main="Bay Anchovy (n=34)", xlab = "Temp",prob=T,xlim = c(0,30))

#Searobin
searobin<-subset(eggs.trim,eggs.trim$`Common Name` == "Northern searobin")
hist(searobin$temp,main="Searobin (n=49)", xlab = "Temp",prob=T,xlim = c(0,30))

silver.count<-na.exclude(count(silver,temp.round))
fit<-fitdistr(silver.count$temp.round,"Poisson")
hist(silver$temp,prob=T)
lines(dpois(0:max(silver.count$temp.round),lambda = fit$estimate[1]))

#Now tally cod counts
cod.count<-count(cod,temp.round)
fit<-fitdistr(cod.count$temp.round,"Poisson")
hist(cod$temp)
lines(dpois(cod$temp.round,lambda = fit$estimate[1]))

#data exploration
bsb.1<-subset(db,db$Cruise_Name == "AL0605")
bsb.2<-subset(db,db$Cruise_Name == "AL0408")

#Try fitting different distributions to silver hake data
silver.count<-na.exclude(count(silver,temp.round))
fit<-fitdistr(silver.count$temp.round,"Poisson")
fit2<-fitdistr(silver.count$temp.round,"normal")
fit3<-fitdistr(silver.count$temp.round,"negative binomial")
hist(silver$temp,prob=T)
lines(dpois(0:max(silver.count$temp.round),lambda = fit$estimate[1]),col="red")
lines(dnorm(0:max(silver.count$temp.round),mean = fit2$estimate[1],sd=fit2$estimate[2]),col="magenta")
lines(dnbinom(0:max(silver.count$temp.round),size = fit3$estimate[1],mu=fit3$estimate[2]),col="blue")

