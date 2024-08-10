#data prep for SPQ analysis
#data = fish egg species ID from DNA barcoding & temp from CTD casts
#samples from EcoMon survey program
#barcoding results as reported in Lewis et al (2016)

# Fri Aug  9 13:52:27 2024 ------------------------------


# prep --------------------------------------------------------------------
#load packages
library(readxl)
library(dplyr)
library(tidyr)
library(here)
library(ggplot2)

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

#bin temperature into 0.5degC bins
db_trim <- db_trim %>% mutate(temp_bin = round(SST*2)/2)
#tally, calculate proportion
temp_dist_lewis <- db_trim %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_samples=n/total) %>% dplyr::select(temp_bin,prop_samples)

#clear environment
rm(list = setdiff(ls(),c("egg_estimates","temp_dist_lewis","db_trim")))





