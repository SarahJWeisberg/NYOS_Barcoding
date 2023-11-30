#Code to calculate thermal niches of fish eggs 
#from EcoMon survey program & barcoding results reported in Lewis et al (2016)

# Mon Nov 27 14:18:41 2023 ------------------------------

# prep --------------------------------------------------------------------
#load packages
library(readxl)
library(dplyr)
library(tidyr)
library(MASS)
library(here)
library(ggplot2)

# Load and tidy data ------------------------------------------------------
#Load station data from subset of stations with ethanol-preserved samples
db <- read_excel("data/FishEggDatabase_GenID.accdb.xls",sheet = "StationData")
#remove rows with NA egg count
#convert egg counts to numeric, station to double
#rename columns to match ecomon
db <- db %>% filter(is.na(Egg_Count) == F) %>% 
  mutate(Egg_Count = as.numeric(Egg_Count),Station=as.double(Station),Ich_Haul_Fctr = as.numeric(Ich_Haul_Fctr)) %>% 
  rename(cruise_name = Cruise_Name, station = Station)

#Load all EcoMon station data
#convert station to double
ecomon <- read.csv("data/EcoMon_Plankton_Data.csv") %>% mutate(station = as.double(station))

#load in DNA barcoding results
eggs<- read_excel("data/FishEggDatabase_GenID.accdb.xls", sheet = "EggIDData")
#rename eggs columns
eggs<- eggs %>% rename(Sample_ID=`Sample ID`, Common_Name = `Common Name`,cruise_name=Cruise, station=Station)
#fix errors in common names
eggs <- eggs %>% mutate(Common_Name = replace(Common_Name,Common_Name == "American","American fourspot flounder")) %>%
  mutate(Common_Name = replace(Common_Name,Common_Name == "Witch Flounder","Witch flounder")) %>%
  mutate(Common_Name = replace(Common_Name,Common_Name == "Yellowtail Flounder","Yellowtail flounder"))

#join datasets, select for relevant columns
db_trim<- left_join(db,ecomon,by=c("cruise_name","station")) %>% 
  dplyr::select(cruise_name,station,sfc_temp,date,Egg_Count,Ich_Haul_Fctr)

#replace NA haul factors with average
db_trim <- db_trim %>% 
  mutate(Ich_Haul_Fctr = ifelse(is.na(Ich_Haul_Fctr),mean(Ich_Haul_Fctr,na.rm = T),Ich_Haul_Fctr))

#histogram of SST of all stations
db_trim %>% 
  ggplot(aes(x=sfc_temp, color=cruise_name, fill=cruise_name,binwidth = 0.5))+
  geom_histogram(alpha=0.5)

#select relevant columns of egg data
eggs_trim<- eggs %>% dplyr::select(Sample_ID,cruise_name,station,Common_Name) %>% 
  mutate(station = as.double(station))

#join egg data with temperature data
egg_temp<-left_join(eggs_trim,db_trim, by = c("cruise_name","station"))

#What is tally of eggs barcoded per station, compared with egg count?
station_tally<-eggs_trim %>% group_by(cruise_name,station) %>% tally()
station_tally <- left_join(db_trim,station_tally,by = c("cruise_name","station")) %>%
  rename(barcoded_count = n) %>% mutate(barcoded_count = ifelse(is.na(barcoded_count),0,barcoded_count))
#there are some whole cruises with non-zero egg counts but none barcoded
#six other stations with this
#let's ignore those for now
station_tally <- station_tally %>% filter(!cruise_name %in% c("GU1302","PC1301","HB1202","DE1105","DE1102")) %>% 
  filter(!(Egg_Count > 0 & barcoded_count == 0)) %>% mutate(diff=Egg_Count-barcoded_count) %>% 
  mutate(cruise_station = paste(cruise_name,station,sep="_"))

#Dealing with subsampling issue
#make new column with cruise and station merged
subsampled<-station_tally %>% filter(diff>0)
#need to bin into 0.05mm diameter bins

eggs_sub <-egg_temp %>% mutate(cruise_station = paste(cruise_name,station,sep="_")) %>% 
  filter(cruise_station %in% subsampled$cruise_station)

eggs_sub <- eggs_sub %>% group_by(cruise_station,Common_Name) %>% tally() %>% rename(abd = n)
eggs_sub <- left_join(eggs_sub,subsampled,by="cruise_station") %>% 
  mutate(est = abd/barcoded_count * Egg_Count * Ich_Haul_Fctr) %>% 
  dplyr::select(Common_Name,cruise_station,sfc_temp,date,est)

eggs_no_sub <- egg_temp %>% mutate(cruise_station = paste(cruise_name,station,sep="_")) %>% 
  filter(!cruise_station %in% subsampled$cruise_station) %>% group_by(Common_Name,cruise_station) %>% 
  tally() %>% rename(abd = n)

eggs_no_sub <- left_join(eggs_no_sub,station_tally,by = "cruise_station") %>% 
  mutate(est = abd * Ich_Haul_Fctr) %>% dplyr::select(Common_Name,cruise_station,sfc_temp,date,est)

egg_estimates<- bind_rows(eggs_no_sub,eggs_sub)

# Analysis ----------------------------------------------------------------
#what might niche breadth be?
#establish a column that is rounded to nearest 0.5degC
#based on Asch & Erisman, 2018
db_trim <- db_trim %>% mutate(temp_bin = round(sfc_temp*2)/2)
#tally, calculate proportion
temp_dist <- db_trim %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_samples=n/total) %>% dplyr::select(temp_bin,prop_samples)


# species-specific analyses -----------------------------------------------
#Look at one species = silver hake
silver<- egg_estimates %>% filter(Common_Name == "Silver hake")
max(silver$sfc_temp,na.rm = T) - min(silver$sfc_temp,na.rm = T)
quantile(silver$sfc_temp,0.95,na.rm = T) - quantile(silver$sfc_temp,0.05,na.rm = T)

ggplot() +
  geom_histogram(data=db_trim,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_line(data=silver,aes(x=sfc_temp,y=est),color="black")

#apply to silver
silver <- silver %>% mutate(temp_bin = round(sfc_temp*2)/2)
silver_dist <- silver %>% mutate(total = sum(est,na.rm=T)) %>% 
  group_by(temp_bin) %>% filter(is.na(temp_bin) == F) %>% 
  mutate(prop_fish=sum(est)/total) %>% dplyr::select(temp_bin,prop_fish) %>% distinct()

ggplot() +
  geom_line(data=silver_dist,aes(x=temp_bin,y=prop_fish))

silver_niche<-left_join(temp_dist,silver_dist,by="temp_bin") %>% 
  mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
  mutate (niche = sqrt(prop_samples*prop_fish)) %>%
  mutate(q = prop_fish/prop_samples) %>%
  mutate(q_norm = q/sum(q)) %>%
  mutate(cdf_q_norm = cumsum(q_norm))

ggplot(silver_niche) +
  geom_line(aes(x=temp_bin,y=cdf_q_norm),color="magenta")
#geom_line(data=cod_niche,aes(x=temp_bin,y=cdf_q_norm),color="navy")

silver_NB <- sum(silver_niche$niche)


#Look at one species = silver hake (old way)
silver<- egg_temp %>% filter(Common_Name == "Silver hake")
max(silver$sfc_temp,na.rm = T) - min(silver$sfc_temp,na.rm = T)
quantile(silver$sfc_temp,0.95,na.rm = T) - quantile(silver$sfc_temp,0.05,na.rm = T)

ggplot() +
  geom_histogram(data=db_trim,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=silver,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,color="black")

#apply to silver
silver <- silver %>% mutate(temp_bin = round(sfc_temp*2)/2)
silver_dist <- silver %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_fish=n/total) %>% dplyr::select(temp_bin,prop_fish)

silver_niche<-left_join(temp_dist,silver_dist,by="temp_bin") %>% 
  mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
  mutate (niche = sqrt(prop_samples*prop_fish)) %>%
  mutate(q = prop_fish/prop_samples) %>%
  mutate(q_norm = q/sum(q)) %>%
  mutate(cdf_q_norm = cumsum(q_norm))

ggplot(silver_niche) +
  geom_line(aes(x=temp_bin,y=cdf_q_norm),color="magenta")+
  geom_line(data=gs_flounder_niche,aes(x=temp_bin,y=cdf_q_norm),color="navy")

silver_NB <- sum(silver_niche$niche)

#apply to Gulf stream flounder
gs_flounder<- egg_estimates %>% filter(Common_Name == "Gulf Stream flounder")
max(gs_flounder$sfc_temp,na.rm = T) - min(gs_flounder$sfc_temp,na.rm = T)
quantile(gs_flounder$sfc_temp,0.95,na.rm = T) - quantile(gs_flounder$sfc_temp,0.05,na.rm = T)

ggplot() +
  geom_histogram(data=db_trim,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_line(data=gs_flounder,aes(x=sfc_temp,y=est),color="black")

#apply to gs_flounder
gs_flounder <- gs_flounder %>% mutate(temp_bin = round(sfc_temp*2)/2)
gs_flounder_dist <- gs_flounder %>% mutate(total = sum(est,na.rm=T)) %>% 
  group_by(temp_bin) %>% filter(is.na(temp_bin) == F) %>% 
  mutate(prop_fish=sum(est)/total) %>% dplyr::select(temp_bin,prop_fish) %>% distinct()

ggplot() +
  geom_line(data=gs_flounder_dist,aes(x=temp_bin,y=prop_fish))

gs_flounder_niche<-left_join(temp_dist,gs_flounder_dist,by="temp_bin") %>% 
  mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
  mutate (niche = sqrt(prop_samples*prop_fish)) %>%
  mutate(q = prop_fish/prop_samples) %>%
  mutate(q_norm = q/sum(q)) %>%
  mutate(cdf_q_norm = cumsum(q_norm))

#old way
gs_flounder<-egg_temp %>% filter(Common_Name == "Gulf Stream flounder")
ggplot() +
  geom_histogram(data=db_trim,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=gs_flounder,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,color="black")

#max-min
max(gs_flounder$sfc_temp,na.rm = T) - min(gs_flounder$sfc_temp,na.rm = T)
#90 CI
quantile(gs_flounder$sfc_temp,0.95,na.rm = T) - quantile(gs_flounder$sfc_temp,0.05,na.rm = T)

gs_flounder <- gs_flounder %>% mutate(temp_bin = round(sfc_temp*2)/2)
gs_flounder_dist <- gs_flounder %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_fish=n/total) %>% dplyr::select(temp_bin,prop_fish)

gs_flounder_niche<-left_join(temp_dist,gs_flounder_dist,by="temp_bin") %>% 
  mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
  mutate (niche = sqrt(prop_samples*prop_fish))%>%
  mutate(q = prop_fish/prop_samples) %>%
  mutate(q_norm = q/sum(q)) %>%
  mutate(cdf_q_norm = cumsum(q_norm))

gs_flounder_NB <- sum(gs_flounder_niche$niche)

#apply to yellowtail flounder
yt_flounder<-egg_temp %>% filter(Common_Name == "Yellowtail flounder")
ggplot() +
  geom_histogram(data=db_trim,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=yt_flounder,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,color="black")

#max-min
max(yt_flounder$sfc_temp,na.rm = T) - min(yt_flounder$sfc_temp,na.rm = T)
#90 CI
quantile(yt_flounder$sfc_temp,0.95,na.rm = T) - quantile(yt_flounder$sfc_temp,0.05,na.rm = T)

yt_flounder <- yt_flounder %>% mutate(temp_bin = round(sfc_temp*2)/2)
yt_flounder_dist <- yt_flounder %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_fish=n/total) %>% dplyr::select(temp_bin,prop_fish)

yt_flounder_niche<-left_join(temp_dist,yt_flounder_dist,by="temp_bin") %>% 
  mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
  mutate (niche = sqrt(prop_samples*prop_fish))%>%
  mutate(q = prop_fish/prop_samples) %>%
  mutate(q_norm = q/sum(q)) %>%
  mutate(cdf_q_norm = cumsum(q_norm))

yt_flounder_NB <- sum(yt_flounder_niche$niche)

#cod
cod<- egg_temp %>% filter(Common_Name == "Atlantic cod")
max(cod$sfc_temp,na.rm = T) - min(cod$sfc_temp,na.rm = T)
quantile(cod$sfc_temp,0.95,na.rm = T) - quantile(cod$sfc_temp,0.05,na.rm = T)

ggplot() +
  geom_histogram(data=db_trim,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=cod,aes(x=sfc_temp),binwidth = 0.5, alpha=0.5,color="black")

#apply to cod
cod <- cod %>% mutate(temp_bin = round(sfc_temp*2)/2)
cod_dist <- cod %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_fish=n/total) %>% dplyr::select(temp_bin,prop_fish)

cod_niche<-left_join(temp_dist,cod_dist,by="temp_bin") %>% 
  mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
  mutate (niche = sqrt(prop_samples*prop_fish)) %>%
  mutate(q = prop_fish/prop_samples) %>%
  mutate(q_norm = q/sum(q)) %>%
  mutate(cdf_q_norm = cumsum(q_norm))

cod_NB <- sum(cod_niche$niche)



#Try fitting different distributions to silver hake data
silver.count<-na.exclude(count(silver,temp.round))
fit<-fitdistr(silver.count$temp.round,"Poisson")
fit2<-fitdistr(silver.count$temp.round,"normal")
fit3<-fitdistr(silver.count$temp.round,"negative binomial")
hist(silver$temp,prob=T)
lines(dpois(0:max(silver.count$temp.round),lambda = fit$estimate[1]),col="red")
lines(dnorm(0:max(silver.count$temp.round),mean = fit2$estimate[1],sd=fit2$estimate[2]),col="magenta")
lines(dnbinom(0:max(silver.count$temp.round),size = fit3$estimate[1],mu=fit3$estimate[2]),col="blue")

