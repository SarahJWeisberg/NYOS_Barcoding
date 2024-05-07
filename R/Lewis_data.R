#Code to calculate thermal niches of fish eggs 
#from EcoMon survey program & barcoding results reported in Lewis et al (2016)

# Mon Jan 15 10:36:51 2024 ------------------------------

# prep --------------------------------------------------------------------
#load packages
library(readxl)
library(dplyr)
library(tidyr)
library(MASS)
library(here)
library(ggplot2)
library(rio)

# Load and tidy data ------------------------------------------------------
#Load station data from subset of stations with ethanol-preserved samples
db_orig <- read_excel("data/FishEggDatabase_GenID.accdb.xls",sheet = "StationData")
#remove rows with NA egg count
#convert egg counts to numeric, create new column with cruise/station unique ID
db <- db_orig %>% filter(is.na(Egg_Count) == F) %>% 
  mutate(Egg_Count = as.numeric(Egg_Count),Station=as.double(Station),Ich_Haul_Fctr = as.numeric(Ich_Haul_Fctr)) %>%
  mutate(cruise_station = paste(Cruise_Name,Station,sep="_"))

#read in database with haul factor info
stations <- read.csv("data/DensitiesBySpecies.csv")
haul <- stations %>% rename(IHF = Ich_Haul_Fctr) %>% mutate(cruise_station = paste(Cruise_Name,Station,sep="_")) %>%
  dplyr::select(cruise_station,IHF)

#merge back with egg db, remove rows with haul factor = NA
#these are from 4 cruises without barcoding results
db <- left_join(db,haul, by = "cruise_station") %>% dplyr::select(-Ich_Haul_Fctr) %>% filter(is.na(IHF) ==F)

#load in DNA barcoding results
eggs_orig<- read_excel("data/FishEggDatabase_GenID.accdb.xls", sheet = "EggIDData")
#rename eggs columns
eggs<- eggs_orig %>% rename(Sample_ID=`Sample ID`, Common_Name = `Common Name`,Cruise_Name=Cruise) %>%
  mutate(cruise_station = paste(Cruise_Name,Station,sep="_"))
#fix errors in common names
eggs <- eggs %>% mutate(Common_Name = replace(Common_Name,Common_Name == "American","American fourspot flounder")) %>%
  mutate(Common_Name = replace(Common_Name,Common_Name == "Witch Flounder","Witch flounder")) %>%
  mutate(Common_Name = replace(Common_Name,Common_Name == "Yellowtail Flounder","Yellowtail flounder"))

#Load station CTD data
ecomon <- read.csv("data/station_CTD.csv")
#rename columns, select for stations of interest, filter for columns of interest
sst <- ecomon %>% rename(Cruise_Name = CRUISE_NAME, Station = STATION, SST = OCTEMPS_SFC_TEMP) %>% 
  mutate(cruise_station = paste(Cruise_Name,Station,sep="_")) %>% 
  filter(cruise_station %in% eggs$cruise_station | cruise_station %in% db$cruise_station) %>%
  dplyr::select(cruise_station,SST)

#average SST for stations with more than 1
sst <- sst %>% group_by(cruise_station) %>% mutate(SST = mean(SST)) %>% distinct()

#join datasets, select for relevant columns, remove NAs
db_trim<- left_join(db,sst,by="cruise_station") %>% 
  dplyr::select(cruise_station,SST,Egg_Count,IHF,Cruise_Name,Station) %>% filter(is.na(SST) == F)

#histogram of SST of all stations
db_trim %>%
  ggplot(aes(x=SST, ,binwidth = 0.5))+
  geom_histogram(alpha=0.5,color="orange",fill="orange")

#join egg data with temperature data
egg_temp<-left_join(eggs,sst, by = "cruise_station") 

#Dealing with biased subsampling
#see Lewis et al., 2016 for details
#What is tally of eggs barcoded per station, compared with egg count?
station_tally<-eggs %>% group_by(cruise_station) %>% tally()
station_tally <- left_join(station_tally,db_trim,by = "cruise_station") %>%
  rename(barcoded_count = n) %>% mutate(barcoded_count = ifelse(is.na(barcoded_count),0,barcoded_count)) %>% 
  mutate(diff=Egg_Count-barcoded_count) %>% dplyr::select(-SST)
station_tally <- left_join(station_tally, sst, by="cruise_station") %>% filter(is.na(SST) == F)

#pull out eggs from stations that were subsampled
subsampled<-station_tally %>% filter(diff>0)
eggs_sub <-egg_temp %>% 
  filter(cruise_station %in% subsampled$cruise_station)

#pull in information about diameter-based subsampling
#add number of unid'd eggs + id'd per size bin to get total - can then calculate proportionality
sizes <- dplyr::bind_rows(import_list("data/FishEgg_MsmtData.xlsx",setclass = "tbl"),.id = "sheet")
sizes <- sizes %>% rename(cruise_station = sheet) %>% mutate(Tot_Frequency = UnID_Frequency + ID_Frequency) %>%
  dplyr::select(-c(UnID_Frequency,ID_Frequency))

#separate which eggs we have size info for, which not
eggs_sub_size <- eggs_sub %>% filter(cruise_station %in% sizes$cruise_station)
eggs_sub_no_size <- eggs_sub %>% filter(!cruise_station %in% sizes$cruise_station)

#estimate total number of eggs per species per subsampled station
eggs_sub_final<-c()
for (i in 1:length(unique(sizes$cruise_station))){
  cs<-unique(sizes$cruise_station)[i] #pull out unique station
  #get barcorded eggs from that station
  #bin species IDs into size categories
  eggs_to_bin <- eggs_sub_size %>% filter(cruise_station == cs) %>% mutate(diameter = as.numeric(`Egg Diameter (mm)`)) %>% 
    mutate(Bin = ceiling(diameter*20)/20) %>% group_by(Common_Name,Bin,cruise_station) %>% tally()
  #get size bin info from that station
  sizes_bin <- sizes %>% filter(cruise_station == cs)
  #join
  bin_est <- right_join(eggs_to_bin,sizes_bin,by=c("cruise_station","Bin"))
  bin_est <- bin_est %>% group_by(Bin) %>% mutate(tot_barcoded = sum(n)) %>% ungroup()
  #deal with nas, have to join with closest bin with barcoding results
  bin_bad<- bin_est %>% filter(Tot_Frequency >0 & is.na(tot_barcoded)==T)
  if (length(bin_bad$Bin) >0){
    for (j in 1:length(bin_bad$Bin)){
    bin_good<-bin_est %>% filter(is.na(tot_barcoded)==F)
    bin_est <- bin_good %>% 
      mutate(Tot_Frequency= ifelse(Bin == bin_good$Bin[which.min(abs(bin_good$Bin-bin_bad$Bin[j]))],
                                   bin_bad$Tot_Frequency[j]+Tot_Frequency,Tot_Frequency))
    }
  }
  bin_est <- bin_est %>% mutate(est = n*Tot_Frequency/tot_barcoded)
  eggs_sub_size_est <- bin_est %>% group_by(Common_Name,cruise_station) %>% mutate(tot=sum(est)) %>% 
    dplyr::select(Common_Name,cruise_station,tot) %>% distinct()
  eggs_sub_final <- rbind(eggs_sub_final,eggs_sub_size_est)
}

# 11 stations without size info - am going to assume proportionality
eggs_sub_no_size <- eggs_sub_no_size %>% group_by(cruise_station,Common_Name) %>% tally() %>% rename(abd = n)
eggs_sub_no_size <- left_join(eggs_sub_no_size,subsampled,by="cruise_station") %>% 
  mutate(tot = abd/barcoded_count * Egg_Count) %>% 
  dplyr::select(Common_Name,cruise_station,tot)

eggs_no_sub <- egg_temp %>% filter(!cruise_station %in% subsampled$cruise_station) %>% group_by(Common_Name,cruise_station) %>% 
  tally() %>% rename(tot = n)

#put it all together
#remove NAs, replace 0's with 1's
egg_estimates <- bind_rows(eggs_no_sub,eggs_sub_final,eggs_sub_no_size) %>% filter(is.na(Common_Name)==F) %>%
  mutate(tot=ifelse(tot==0,1,tot))

#add in SST, haul factor info
egg_estimates <- left_join(egg_estimates,sst,by='cruise_station') %>% filter(is.na(SST)==F)
egg_estimates <- left_join(egg_estimates,haul,by="cruise_station") %>% filter(is.na(IHF)==F)

#multiply count by haul factor to get estimates normalized for effort
egg_estimates <- egg_estimates %>% mutate(est=tot*IHF)

# Analysis ----------------------------------------------------------------
#what might niche breadth be?
#establish a column that is rounded to nearest 0.5degC
#based on Asch & Erisman, 2018
db_trim <- db_trim %>% mutate(temp_bin = round(SST*2)/2)
#tally, calculate proportion
temp_dist <- db_trim %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_samples=n/total) %>% dplyr::select(temp_bin,prop_samples)

# species-specific analyses -----------------------------------------------
#get cdf plots for several species of interest
spp<-egg_estimates %>% group_by(Common_Name) %>% mutate(total=sum(est)) %>% 
  dplyr::select(Common_Name,total) %>% distinct()
high_abd_spp<-spp %>% ungroup() %>% slice_max(total,n=8)

#or select 3 species of interest
spp_interest<-c("Silver hake","American fourspot flounder","Gulf Stream flounder")

cdfs_lewis<-c()
for(i in 1:length(high_abd_spp$Common_Name)){
  sp <- high_abd_spp$Common_Name[i]
  eggs_sp <- egg_estimates %>% filter(Common_Name == sp)
  #round temp to nearest 0.5deg bin
  eggs_sp <- eggs_sp %>% mutate(temp_bin = round(SST*2)/2)
  eggs_sp_dist <- eggs_sp %>% mutate(total = sum(est,na.rm=T)) %>% 
    group_by(temp_bin) %>% filter(is.na(temp_bin) == F) %>% 
    mutate(prop_fish=sum(est)/total) %>% dplyr::select(temp_bin,prop_fish) %>% distinct()
  eggs_niche<-left_join(temp_dist,eggs_sp_dist,by="temp_bin") %>% 
    mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
    mutate (niche = sqrt(prop_samples*prop_fish)) %>%
    mutate(q = prop_fish/prop_samples) %>%
    mutate(q_norm = q/sum(q)) %>%
    mutate(cdf_q_norm = cumsum(q_norm)) %>%
    mutate(spp = sp)
  cdfs_lewis <- rbind(cdfs_lewis,eggs_niche)
}

#plot results
ggplot(cdfs_lewis) +
  geom_line(aes(x=temp_bin,y=cdf_q_norm,color=spp))

#fourspot
fourspot<-cdfs_lewis %>% filter(spp == "American fourspot flounder")
#Look at one species = silver hake
silver<- egg_estimates %>% filter(Common_Name == "Silver hake")
max(silver$SST,na.rm = T) - min(silver$SST,na.rm = T)
quantile(silver$SST,0.9,na.rm = T) - quantile(silver$SST,0.1,na.rm = T)

ggplot() +
  geom_histogram(data=db_trim,aes(x=SST),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_line(data=silver,aes(x=SST,y=est),color="black")

#apply to silver
silver <- silver %>% mutate(temp_bin = round(SST*2)/2)
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
  # geom_vline(xintercept = 13.75, color= "black")+
  # geom_vline(xintercept = 18.25, color= "black")
#geom_line(data=cod_niche,aes(x=temp_bin,y=cdf_q_norm),color="navy")

ggplot(silver_niche) +
  geom_line(aes(x=temp_bin,y=q_norm),color="magenta")

  

silver_NB <- sum(silver_niche$niche)


#Look at one species = silver hake (old way)
silver<- egg_temp %>% filter(Common_Name == "Silver hake")
max(silver$SST,na.rm = T) - min(silver$SST,na.rm = T)
quantile(silver$SST,0.95,na.rm = T) - quantile(silver$SST,0.05,na.rm = T)

ggplot() +
  geom_histogram(data=db_trim,aes(x=SST),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=silver,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")

#apply to silver
silver <- silver %>% mutate(temp_bin = round(SST*2)/2)
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
  
ggplot()+
  geom_line(data=gs_flounder_niche,aes(x=temp_bin,y=cdf_q_norm),color="navy")

silver_NB <- sum(silver_niche$niche)

#apply to Gulf stream flounder
gs_flounder<- egg_estimates %>% filter(Common_Name == "Gulf Stream flounder")
max(gs_flounder$SST,na.rm = T) - min(gs_flounder$SST,na.rm = T)
quantile(gs_flounder$SST,0.95,na.rm = T) - quantile(gs_flounder$SST,0.05,na.rm = T)

ggplot() +
  geom_histogram(data=db_trim,aes(x=SST),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_line(data=gs_flounder,aes(x=SST,y=est),color="black")

#apply to gs_flounder
gs_flounder <- gs_flounder %>% mutate(temp_bin = round(SST*2)/2)
gs_flounder_dist <- gs_flounder %>% mutate(total = sum(est,na.rm=T)) %>% 
  group_by(temp_bin) %>% filter(is.na(temp_bin) == F) %>% 
  mutate(prop_fish=sum(est)/total) %>% dplyr::select(temp_bin,prop_fish) %>% distinct()

ggplot(data=gs_flounder_dist, aes(x=temp_bin,y=prop_fish))+
  geom_line()

gs_flounder_niche<-left_join(temp_dist,gs_flounder_dist,by="temp_bin") %>% 
  mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
  mutate (niche = sqrt(prop_samples*prop_fish)) %>%
  mutate(q = prop_fish/prop_samples) %>%
  mutate(q_norm = q/sum(q)) %>%
  mutate(cdf_q_norm = cumsum(q_norm))

ggplot() +
  geom_line(data=gs_flounder_niche,aes(x=temp_bin,y=q_norm),color= "navy")

ggplot() +
  geom_line(data=gs_flounder_niche,aes(x=temp_bin,y=cdf_q_norm),color= "navy")+
  geom_vline(xintercept = 16.5) +
  geom_vline(xintercept = 25.5)

#old way
gs_flounder<-egg_temp %>% filter(Common_Name == "Gulf Stream flounder")
ggplot() +
  geom_histogram(data=db_trim,aes(x=SST),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=gs_flounder,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")

#max-min
max(gs_flounder$SST,na.rm = T) - min(gs_flounder$SST,na.rm = T)
#90 CI
quantile(gs_flounder$SST,0.95,na.rm = T) - quantile(gs_flounder$SST,0.05,na.rm = T)

gs_flounder <- gs_flounder %>% mutate(temp_bin = round(SST*2)/2)
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
  geom_histogram(data=db_trim,aes(x=SST),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=yt_flounder,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")

#max-min
max(yt_flounder$SST,na.rm = T) - min(yt_flounder$SST,na.rm = T)
#90 CI
quantile(yt_flounder$SST,0.95,na.rm = T) - quantile(yt_flounder$SST,0.05,na.rm = T)

yt_flounder <- yt_flounder %>% mutate(temp_bin = round(SST*2)/2)
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
max(cod$SST,na.rm = T) - min(cod$SST,na.rm = T)
quantile(cod$SST,0.95,na.rm = T) - quantile(cod$SST,0.05,na.rm = T)

ggplot() +
  geom_histogram(data=db_trim,aes(x=SST),binwidth = 0.5, alpha=0.5,fill="magenta")+
  geom_histogram(data=cod,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")

#apply to cod
cod <- cod %>% mutate(temp_bin = round(SST*2)/2)
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

