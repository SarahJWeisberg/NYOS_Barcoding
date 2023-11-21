# Tue Nov 21 09:55:08 2023 ------------------------------


library(readxl)
library(dplyr)
library(tidyr)
library(MASS)
library(here)
library(ggplot2)

#Load data
#db <- read_excel("FishEggDatabase_GenID.accdb.xls",sheet = "StationData")
db <- read.csv("data/EcoMon_Plankton_Data.csv")
eggs<- read_excel("FishEggDatabase_GenID.accdb.xls", sheet = "EggIDData")

#rename eggs columns
eggs<- eggs %>% rename(Sample_ID=`Sample ID`, Common_Name = `Common Name`,cruise_name=Cruise, station=Station)

#Trim data for relevant columns
db_trim<- db %>% dplyr::select(cruise_name,station,sfc_temp,date)
db_trim$station<-as.factor(db_trim$station)
eggs_trim<- eggs %>% dplyr::select(Sample_ID,cruise_name,station,Common_Name)
eggs_trim$station<-as.factor(eggs_trim$station)

#join egg data with temperature data
egg_cruise<-left_join(eggs_trim,db_trim, by = c("cruise_name","station"))

#histogram of all stations
egg_cruise %>% dplyr::select(cruise_name,station,sfc_temp) %>% distinct() %>%
  ggplot(aes(x=sfc_temp, color=cruise_name, fill=cruise_name,binwidth = 0.5))+
  geom_histogram(alpha=0.5)

stations_ecomon<-egg_cruise %>% dplyr::select(cruise_name,station,sfc_temp) %>% distinct()
stations_all<-bind_rows(stations,stations_ecomon) 
ggplot(data = stations_all,aes(x=sfc_temp))+
  geom_histogram(alpha=0.5,fill="magenta",color="magenta")

#Look at one species = silver hake
silver<- egg_cruise %>% filter(Common_Name == "Silver hake")

ggplot() +
  geom_histogram(data=stations_ecomon,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=silver,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,color="black")

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


#Gulf stream flounder
gs_flounder<-egg_cruise %>% filter(Common_Name == "Gulf Stream flounder")
ggplot() +
  geom_histogram(data=stations_ecomon,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=gs_flounder,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,color="black")


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

