# Mon Nov 27 14:17:51 2023 ------------------------------


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

#use CTD SST data
stations <- stations %>% rename(SST = SST_CTD) %>% mutate(SST = as.numeric(SST))
#filter just SST data
stations_filter <- stations %>% dplyr::select(Cruise_ID,Station_ID,SST)

#remove all columns but genus, species, common name, station, type
fish_filter<- fish %>% 
  dplyr::select(c("Genus.species","Common_Name","Cruise.ID","Station.ID","Stage")) %>% 
  rename("Cruise_ID" = "Cruise.ID") %>% rename("Station_ID"="Station.ID")
  
#filter out 2103, juveniles
fish_filter<-fish_filter %>% filter(Cruise_ID != "2103") %>% filter(Stage!="Juvenile")

#merge with SST data
fish_filter <- left_join(fish_filter,stations_filter,by=c("Cruise_ID","Station_ID"))

#silver hake eggs
silver<-fish_filter %>% filter(Common_Name == "Silver Hake" & Stage == "Egg") 

#plot temp data - all
stations %>%
  ggplot(aes(x=SST, color=Cruise_ID, fill=Cruise_ID)) +
  geom_histogram(binwidth = 0.5, alpha=0.5)

#tally all
fish_tally <- fish_filter %>% group_by(Common_Name,Stage) %>% tally()

silver_tally<-silver %>% group_by(Cruise_ID,Station_ID) %>% tally() 

silver_tally <- left_join(stations_filter,silver_tally,by=c("Cruise_ID","Station_ID")) %>%
  rename("SilverHake" = "n") %>% mutate(SilverHake = ifelse(is.na(SilverHake),0,SilverHake))

ggplot() +
  geom_histogram(data=stations,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
  geom_histogram(data=silver,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")


#what about GS Flounder?
gs_egg<-fish_filter %>% filter(Common_Name == "Gulfstream Flounder" & Stage == "Egg")
gs_larva<-fish_filter %>% filter(Common_Name == "Gulfstream Flounder" & Stage == "Larva")
ggplot() +
  geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
  geom_histogram(data=gs_egg,aes(x=SST),binwidth = 0.5, alpha=0.5,fill="black")+
  geom_histogram(data=gs_larva,aes(x=SST),binwidth = 0.5, alpha=0.5,fill="navy")

#what about Mackerel?
mackerel<-fish_filter %>% filter(Common_Name == "Atlantic Mackerel" & Stage == "Egg")
ggplot() +
  geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
  geom_histogram(data=mackerel,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")

#what about Red Hake?
red_hake<-fish_filter %>% filter(Common_Name == "Red Hake" & Stage == "Egg")
ggplot() +
  geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
  geom_histogram(data=red_hake,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")

#what about Northern Searobin?
n_searobin<-fish_filter %>% filter(Common_Name == "Northern Searobin" & Stage == "Larva")
ggplot() +
  geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
  geom_histogram(data=n_searobin,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")

# niche estimates -----------------------------------------------
#what might niche breadth be?
#establish a column that is rounded to nearest 0.5degC
#based on Asch & Erisman, 2018
temp_dist_NYOS <- stations_filter %>% mutate(temp_bin = round(SST*2)/2)
#tally, calculate proportion
temp_dist_NYOS <- temp_dist_NYOS %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_samples=n/total) %>% dplyr::select(temp_bin,prop_samples)

#get cdf plots for several species of interest
spp_interest<- c("Silver Hake","Gulfstream Flounder", "Northern Searobin","Butterfish","Fourspot Flounder","Bullet Tuna")

cdfs<-c()
for(i in 1:length(spp_interest)){
  sp <- spp_interest[i]
  n_sp <- fish_tally %>% filter(Common_Name == sp) %>% rename(total=n)
  eggs_larv_sp <- fish_filter %>% filter(Common_Name == sp)
  #round temp to nearest 0.5deg bin
  eggs_larv_sp <- eggs_larv_sp %>% mutate(temp_bin = round(SST*2)/2)
  sp_dist <- eggs_larv_sp %>% group_by(temp_bin,Stage) %>% tally()
  #add column with totals
  sp_dist <- left_join(sp_dist,n_sp,by="Stage") %>% group_by(Stage) %>%
    mutate(prop_fish=n/total) %>% dplyr::select(temp_bin,prop_fish) %>% distinct()
  egg_niche<- sp_dist %>% 
    filter(Stage == "Egg")%>% ungroup %>%
    dplyr::select(-Stage)
  egg_niche<- left_join(temp_dist_NYOS,egg_niche,by="temp_bin") %>%
    mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
    mutate (niche = sqrt(prop_samples*prop_fish)) %>%
    mutate(q = prop_fish/prop_samples) %>%
    mutate(q_norm = q/sum(q)) %>%
    mutate(cdf_q_norm = cumsum(q_norm)) %>%
    mutate(spp = sp) %>%
    mutate(Stage = "Egg")
  larv_niche<- sp_dist %>% 
    filter(Stage == "Larva")%>% ungroup %>%
    dplyr::select(-Stage)
  larv_niche<- left_join(temp_dist_NYOS,larv_niche,by="temp_bin") %>%
    mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
    mutate (niche = sqrt(prop_samples*prop_fish)) %>%
    mutate(q = prop_fish/prop_samples) %>%
    mutate(q_norm = q/sum(q)) %>%
    mutate(cdf_q_norm = cumsum(q_norm)) %>%
    mutate(spp = sp) %>%
    mutate(Stage = "Larva")
  cdfs <- rbind(cdfs,egg_niche,larv_niche)
}

#plot
ggplot(data=cdfs,aes(x=temp_bin,y=cdf_q_norm,color=Stage))+
  geom_line()+
  facet_wrap(~spp)

cdfs %>% filter(spp == "Silver Hake") %>%
  ggplot()+
  geom_line(aes(x=temp_bin,y=cdf_q_norm,color=Stage))+
  geom_line(data=silver_niche,aes(x=temp_bin,y=cdf_q_norm,color="Lewis"))+
  labs(title="Silver Hake")

cdfs %>% filter(spp == "Fourspot Flounder") %>%
  ggplot()+
  geom_line(aes(x=temp_bin,y=cdf_q_norm,color=Stage))+
  geom_line(data=fourspot,aes(x=temp_bin,y=cdf_q_norm,color="Lewis"))+
  labs(title="Fourspot Flounder")

cdfs %>% filter(spp == "Gulfstream Flounder") %>%
  ggplot()+
  geom_line(aes(x=temp_bin,y=cdf_q_norm,color=Stage))+
  geom_line(data=gs_flounder_niche,aes(x=temp_bin,y=cdf_q_norm,color="Lewis"))+
  labs(title="Gulfstream Flounder")
