#single parameter quotient calculations for thermal habitats of fish eggs and larvae
#samples from Northeast US continental shelf

#larval samples from NOAA's EcoMon program
#eggs ID'd via DNA barcoding - two sources
#1 - samples from Lewis et al (2016): subset of EcoMon samples
#2 - samples from NY Offshore Project (NYOS)

# Fri Aug  9 13:47:58 2024 ------------------------------


# prep --------------------------------------------------------------------

#load packages
library(here)

#plotting
library(ggplot2) 
library(ggthemes) 
library(ggpubr)
library(gridExtra)
library(RColorBrewer)

#data wrangling
library(dplyr) #tidy 
library(tidyverse) #tidy
library(zoo) #ordering
library(broom) #converts to tidy 
library(readxl) #reading 
library(here)

#modeling
library(rstatix) #modeling with tidy
library(glmmTMB)#mixed models 

#set overall theme 
theme_Publication <- function(base_size=10, base_family="Arial") {
  
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.spacing  = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0")
   ))
  
}

select <- dplyr::select

#load prepped data sets & merge
source(here("R/Lewis_data.R")) #egg data from Lewis_et_al_2016
source(here("R/EcoMon_larval_data.R"))#larval data from EcoMon
source(here("R/NYOS_egg_data.R"))

#rename so it's easier
#tweak dataframes so they match, merge
#need cruise_station, ID, temp, abundance estimate for each
eggs_lewis <- egg_estimates %>% select(Common_Name,cruise_station,SST,est) %>% relocate(Common_Name,.after = SST) %>%
  rename(spp = Common_Name, abd = est) %>% mutate(data="Lewis")

eggs_NYOS <- fish_filter %>% relocate(SST,.after = cruise_station) %>%
  rename(spp = Common_Name,abd=n) %>% mutate(data="NYOS")

larvae_ecomon <- ecomon_long %>% rename(SST = sfc_temp,abd=value) %>% mutate(data="EcoMon")

#rename species of interest so they match
#note that I've already filtered EcoMon larval dataset for spp of interest
larvae_ecomon <- larvae_ecomon %>% mutate(spp = ifelse(spp=="hipobl_10m2","Fourspot Flounder",
                                                       ifelse(spp=="merbil_10m2","Silver Hake",
                                                              ifelse(spp=="citarc_10m2","Gulfstream Flounder",
                                                                     ifelse(spp=="pep_tri","Butterfish",spp)))))

eggs_lewis <- eggs_lewis %>% mutate(spp = ifelse(spp=="American fourspot flounder","Fourspot Flounder",
                                                 ifelse(spp=="Silver hake","Silver Hake",
                                                        ifelse(spp=="Gulf Stream flounder","Gulfstream Flounder",
                                                               ifelse(spp=="Atlantic butterfish","Butterfish",spp)))))


ichthyo_all <- rbind(eggs_lewis,eggs_NYOS,larvae_ecomon)
ichthyo_all <- ichthyo_all %>% mutate(temp_bin=round(SST*2)/2)

#merge temperature proportion data into one dataframe
temp_dist_lewis <- temp_dist_lewis %>% select(temp_bin,prop_samples) %>% mutate(data="Lewis")
temp_dist_NYOS <- temp_dist_NYOS %>% select(temp_bin,prop_samples) %>% mutate(data="NYOS")
temp_dist_ecomon <- temp_dist_ecomon %>% select(temp_bin,prop_samples) %>% mutate(data="EcoMon")

temp_dist_all <- rbind(temp_dist_lewis,temp_dist_NYOS,temp_dist_ecomon)

# Single parameter quotient (SPQ) calculations --------------------------------------------------------
#select 4 species of interest
spp_interest<-c("Silver Hake","Fourspot Flounder","Gulfstream Flounder", "Butterfish")
datasets <- c("Lewis","NYOS","EcoMon")

ichthyo_spp <- ichthyo_all %>% filter(spp %in% spp_interest)

#calculate spq, normalized q, cdf of normalized q for each species x dataset combo
cdfs<-c()
for(j in 1: length(datasets)){
  temp_dist <- temp_dist_all %>% filter(data == datasets[j])
  ichthyo <- ichthyo_spp %>% filter(data == datasets[j])
  for(i in 1:length(spp_interest)){
    fish <- ichthyo %>% filter(spp == spp_interest[i])
    #round temp to nearest 0.5deg bin
    fish <- fish %>% mutate(temp_bin = round(SST*2)/2)
    fish_dist <- fish %>% mutate(total = sum(abd,na.rm=T)) %>% 
      group_by(temp_bin) %>% filter(is.na(temp_bin) == F) %>% 
      mutate(prop_fish=sum(abd)/total) %>% dplyr::select(temp_bin,prop_fish) %>% distinct()
    fish_spq<-left_join(temp_dist,fish_dist,by="temp_bin") %>% 
      mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
      mutate(q = prop_fish/prop_samples) %>%
      mutate(q_norm = q/sum(q)) %>%
      mutate(cdf_q_norm = cumsum(q_norm)) %>%
      mutate(spp = spp_interest[i],data = datasets[j])
    cdfs <- rbind(cdfs,fish_spq)
  }
}


# bootstrapping -----------------------------------------------------------
#need stations with true 0s for bootstrapping
stations_NYOS <- stations_filter %>% select(cruise_station,SST) %>% mutate(data="NYOS")
stations_lewis <- db_trim %>% select(cruise_station,SST) %>% mutate(data="Lewis")
stations_ecomon <- ecomon_long %>% select(cruise_station,sfc_temp) %>% unique() %>% 
 group_by(cruise_station) %>% summarise(SST= mean(sfc_temp)) %>% mutate(data="EcoMon")

stations_all <- rbind(stations_NYOS,stations_lewis,stations_ecomon)

#bootstrapping runs
boot<-1000
boot_out<-c()
for (s in 1:length(datasets)){
  stations_boot <- stations_all %>% filter(data == datasets[s])
  ichthyo_boot <- ichthyo_spp %>% filter(data == datasets[s])
  for (k in 1:length(spp_interest)) {
    spp_boot <- ichthyo_boot %>% filter(spp == spp_interest[k]) %>% ungroup %>% select(cruise_station,abd)
    boot_data <- left_join(stations_boot,spp_boot,by="cruise_station") %>% 
      mutate(abd = replace(abd,is.na(abd),0)) %>% select(-data)
    for (j in 1:boot){
      station_n <- nrow(boot_data)
      newdata<-as.data.frame(matrix(nrow = station_n,ncol=ncol(boot_data)))
      for(i in 1:station_n){
        rownums<-sample(1:station_n,length(station_n),replace = T) 
        newdata[i,]<-boot_data[rownums,]
      }
      colnames(newdata)<-c("cruise_station","SST","abd")
      newdata <- newdata %>% mutate(temp_bin = round(SST*2)/2)
      newdata_temp <- newdata %>% group_by(temp_bin) %>% tally() %>% rename(n_samples = n)
      newdata <- newdata %>% group_by(temp_bin) %>% mutate(n_fish = sum(abd)) %>%
        dplyr::select(temp_bin,n_fish) %>% unique()
      newdata<-left_join(newdata,newdata_temp,by="temp_bin")
      new_q <- newdata %>% ungroup() %>% mutate(q = (n_fish/sum(n_fish)/(n_samples/sum(n_samples)))) %>% 
        dplyr::select(temp_bin,q) %>% mutate(data=datasets[s],spp=spp_interest[k])
      boot_out<-rbind(boot_out,new_q)
    }
  }
}

#get 5,95 CI
boot_summary <- boot_out %>% group_by(temp_bin,spp,data) %>% 
  summarise(lower_CI = quantile(q,0.05),upper_CI = quantile(q,0.95)) 

#think about estimating thermal breadth from each run instead
#remove any q's less than 1
boot_range <- boot_out %>% filter(q >=1) %>%
  group_by(run_num) %>% 
  summarise(lower = min(temp_bin),upper = max(temp_bin)) %>%
  mutate(range = upper-lower)



# modeling ----------------------------------------------------------------
test_glm<-glm(data=butter_cdf,formula=cdf_q_norm~temp_bin, family = "quasibinomial")
test_predictors<-list(temp_bin=seq(min(temp_dist$temp_bin), max(temp_dist$temp_bin),0.5))
glm_vals <- predict(test_glm, newdata=test_predictors, type = "response")
glm_predict <- data.frame(temp = test_predictors[[1]], prop = glm_vals)
ggplot(data=glm_predict, aes(x=temp,y=prop))+
  geom_line()+
  geom_line(data=butter_cdf,aes(x=temp_bin,y=cdf_q_norm),color="magenta")


# plotting ----------------------------------------------------------------
#plot results
cdfs_lewis %>% filter(spp %in% spp_interest) %>%
  ggplot() +
  geom_line(aes(x=temp_bin,y=cdf_q_norm,color=spp))+
  #geom_hline(yintercept = 1)+
  facet_wrap(~spp)+
  theme(legend.position = "none")
boot_summary %>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI))+
  facet_wrap(~spp)


