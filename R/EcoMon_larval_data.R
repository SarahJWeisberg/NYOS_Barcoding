#data prep for SPQ analysis
#data = fish larvae species ID & temp from CTD casts
#samples from EcoMon survey program

# Fri Aug  9 14:37:03 2024 ------------------------------

#load in main EcoMon dataset
#accession date: May 6, 2024
ecomon <- read.csv(here("data/EcoMon_Data/EcoMon_Plankton_Data_v3_8.csv"))

#load in butterfish dataset from Dave Richardson
#has species-level classification of butterfish (genus-level reported in public EcoMon dataset)
butterfish_raw <- read_xlsx(here("data/DatasetExplanation_Peprilus_May2024.xlsx"))

#get rid of leading zeros in station column
butterfish <- butterfish_raw %>% mutate(station = round(as.numeric(STATION))) %>% 
  mutate(cruise_station = paste(CRUISE_NAME,station, sep = "_"))
#get rid of unneeded columns
butterfish <- butterfish %>% select(cruise_station,TAXA_ICHTHYO,ABUNDANCE)
#calculate abundance by butterfish species
#170511100 = Peprilus sp., 170511101 = P. paru, 170511103 = P. burti (presumed mis-ID), 170511104 = P. triacanthus
butterfish <- butterfish %>% group_by(cruise_station) %>% 
  mutate(pep_sp = sum(ABUNDANCE[which(TAXA_ICHTHYO == 170511100)])) %>%
  mutate(pep_par = sum(ABUNDANCE[which(TAXA_ICHTHYO == 170511101)])) %>%
  mutate(pep_tri = sum(ABUNDANCE[which(TAXA_ICHTHYO >= 170511102)])) %>%
  mutate(pep_tri = replace(pep_tri, is.na(pep_tri), 0)) %>%
  select(cruise_station,pep_par,pep_sp,pep_tri) %>%
  unique

butterfish <- butterfish %>% mutate(pep_all = pep_tri+pep_sp+pep_par)

#check prorportion of abundance that is pep_tri = 99.9%
sum(butterfish$pep_tri)/(sum(butterfish$pep_tri)+sum(butterfish$pep_par))

#merge butterfish with rest
ecomon <- ecomon %>% mutate(cruise_station = paste(cruise_name,station, sep = "_"))

ecomon <- left_join(ecomon,butterfish,by="cruise_station")

## --------------------------------------------------------------------------------------------------------------------------------------------------
ecomon[ecomon == 9999] <- NA

#first split out x (environment) and y (species) to convert NAs to true 0s (because it is a survey)
#select only species of interest
ecomon_x <- ecomon %>% dplyr::select(cruise_station,sfc_temp)
ecomon_y <- ecomon %>% select("hipobl_10m2","merbil_10m2","citarc_10m2","scoaqu_10m2","pep_tri")
ecomon_y[is.na(ecomon_y)] <- 0 #NAs should be true 0s

ecomon <- cbind(ecomon_x, ecomon_y)
ecomon <- ecomon %>% filter(!is.na(sfc_temp))
ecomon_ynames <- colnames(ecomon_y)

#histogram of SSTs
ecomon %>% 
  ggplot(aes(x=sfc_temp,binwidth = 0.5))+
  geom_histogram(alpha=0.5,color="#F8766D",fill="#F8766D")+
  theme(legend.position = "none")

ecomon <- ecomon %>% mutate(temp_bin = round(sfc_temp*2)/2)

temp_dist_ecomon <- ecomon %>% group_by(temp_bin) %>% tally %>% 
  mutate(total=sum(n),prop_samples=n/total) %>% dplyr::select(temp_bin,prop_samples,n)

#create long dataset of larval abundances for tidy manipulation 
#around temp to nearest 0.5deg bin
ecomon_long <- ecomon %>% pivot_longer(all_of(ecomon_ynames), names_to = "spp") %>% select(-temp_bin)
 
rm(butterfish,butterfish_raw,ecomon_y,ecomon_x,ecomon_ynames)



