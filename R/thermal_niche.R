# Wed Aug  7 13:52:00 2024 ------------------------------



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


#plot temp data - all
stations %>%
  ggplot(aes(x=SST)) +
  geom_histogram(binwidth = 0.5, alpha=0.5, color="#A58AFF",fill="#A58AFF")

#tally all
fish_tally <- fish_filter %>% group_by(Common_Name,Stage) %>% tally()

#single species example: silver hake
#silver hake eggs
silver<-fish_filter %>% filter(Common_Name == "Silver Hake" & Stage == "Egg") 
silver_tally<-silver %>% group_by(Cruise_ID,Station_ID) %>% tally()

silver_tally <- left_join(stations_filter,silver_tally,by=c("Cruise_ID","Station_ID")) %>%
  rename("SilverHake" = "n") %>% mutate(SilverHake = ifelse(is.na(SilverHake),0,SilverHake))

ggplot() +
  geom_histogram(data=stations,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
  geom_histogram(data=silver,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")


# #what about GS Flounder?
# gs_egg<-fish_filter %>% filter(Common_Name == "Gulfstream Flounder" & Stage == "Egg")
# gs_larva<-fish_filter %>% filter(Common_Name == "Gulfstream Flounder" & Stage == "Larva")
# ggplot() +
#   geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
#   geom_histogram(data=gs_egg,aes(x=SST),binwidth = 0.5, alpha=0.5,fill="black")+
#   geom_histogram(data=gs_larva,aes(x=SST),binwidth = 0.5, alpha=0.5,fill="navy")
# 
# #what about Mackerel?
# mackerel<-fish_filter %>% filter(Common_Name == "Atlantic Mackerel" & Stage == "Egg")
# ggplot() +
#   geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
#   geom_histogram(data=mackerel,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")
# 
# #what about Red Hake?
# red_hake<-fish_filter %>% filter(Common_Name == "Red Hake" & Stage == "Egg")
# ggplot() +
#   geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
#   geom_histogram(data=red_hake,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")
# 
# #what about Northern Searobin?
# n_searobin<-fish_filter %>% filter(Common_Name == "Northern Searobin" & Stage == "Larva")
# ggplot() +
#   geom_histogram(data=stations_filter,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
#   geom_histogram(data=n_searobin,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")

# niche estimates -----------------------------------------------
#what might niche breadth be?
#establish a column that is rounded to nearest 0.5degC
#based on Asch & Erisman, 2018
temp_dist_NYOS <- stations_filter %>% mutate(temp_bin = round(SST*2)/2)

bins<- temp_dist_NYOS %>% dplyr::select(temp_bin) %>% filter(!is.na(temp_bin)) %>% unique()

#plotstations
temp_dist_NYOS%>%
ggplot(aes(x=SST, color=Cruise_ID, fill=Cruise_ID)) +
  geom_histogram(binwidth = 0.5, alpha=0.5)

#tally, calculate proportion
temp_dist_NYOS <- temp_dist_NYOS %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_samples=n/total) %>% dplyr::select(temp_bin,prop_samples,n) %>% rename(n_temp = n)

#get cdf plots for several species of interest
spp_interest<- c("Silver Hake","Gulfstream Flounder","Butterfish","Fourspot Flounder")

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
    mutate(prop_fish=n/total) %>% dplyr::select(temp_bin,prop_fish,n) %>% distinct()
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
    mutate(Stage = "Egg") %>%
    mutate(n = replace(n,is.na(n),0))
  # larv_niche<- sp_dist %>% 
  #   filter(Stage == "Larva")%>% ungroup %>%
  #   dplyr::select(-Stage)
  # larv_niche<- left_join(temp_dist_NYOS,larv_niche,by="temp_bin") %>%
  #   mutate(prop_fish=ifelse(is.na(prop_fish),0,prop_fish)) %>%
  #   mutate (niche = sqrt(prop_samples*prop_fish)) %>%
  #   mutate(q = prop_fish/prop_samples) %>%
  #   mutate(q_norm = q/sum(q)) %>%
  #   mutate(cdf_q_norm = cumsum(q_norm)) %>%
  #   mutate(spp = sp) %>%
  #   mutate(Stage = "Larva")
  cdfs <- rbind(cdfs,egg_niche)
}

#plot
#test_avg %>% filter(spp %in% c("Silver Hake", "Gulfstream Flounder","Fourspot Flounder")) %>%
cdfs %>%
ggplot(aes(x=temp_bin,y=cdf_q_norm,color=spp))+
  geom_line(linewidth=1)+
  #geom_hline(yintercept = 1)+
  #scale_color_manual(values = c("#00BA38","#00BFC4","#F564E3"))+
  facet_wrap(~spp)+
  theme(legend.position = "none")

cdfs %>% filter(spp == "Silver Hake") %>%
  ggplot()+
  geom_line(aes(x=temp_bin,y=cdf_q_norm))+
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


cdfs_all %>% filter(!data=="EcoMon") %>%
  ggplot()+
  geom_line(aes(x=temp_bin,y=cdf_q_norm,color=spp,linetype = data),linewidth=1)+
  scale_color_manual(values = c("#00BA38","#00BFC4","#F564E3"))+
  facet_wrap(~spp)+
  guides(color="none")




# modeling ----------------------------------------------------------------

#trying out using glms to model these 
test_silver<-cdfs %>% filter(spp=="Silver Hake")
test_glm<-glm(data=test_avg_silver,formula=cdf_q_norm~temp_bin, family = "quasibinomial")
test_predictors<-list(temp_bin=seq(min(test_silver$temp_bin), max(test_silver$temp_bin),0.5))
glm_vals <- predict(test_glm, newdata=test_predictors, type = "response")
glm_predict <- data.frame(temp = test_predictors[[1]], prop = glm_vals)
ggplot(data=glm_predict, aes(x=temp,y=prop))+
  geom_line()+
  geom_line(data=test_silver,aes(x=temp_bin,y=cdf_q_norm),color="magenta")+
  geom_line(data=silver_q,aes(x=temp_bin,y=cdf_q_norm),color="blue")
#ok now mixed model?
# NYOS_silver<-cdfs %>% filter(spp=="Silver Hake") %>% mutate(survey = "NYOS") %>% dplyr::select(-Stage,-spp)
# Lewis_silver <- cdfs_lewis %>% filter(spp == "Silver hake") %>% mutate(survey = "EcoMon") %>% dplyr::select(-spp)
# test_silver <- rbind(NYOS_silver,Lewis_silver)
# test_glm<-glmer(data=test_silver,formula=cdf_q_norm~temp_bin|survey, family = "Gamma")
# test_predictors<-list(temp_bin=seq(min(test_silver$temp_bin), max(test_silver$temp_bin),0.5))
# glm_vals <- predict(test_glm, newdata=test_predictors, type = "response")
# glm_predict <- data.frame(temp = test_predictors[[1]], prop = glm_vals)
# ggplot(data=glm_predict, aes(x=temp,y=prop))+
#   geom_line()+
#   geom_line(data=test_silver,aes(x=temp_bin,y=cdf_q_norm),color="magenta")

#double logistic fitting
library(minpack.lm)
double_logistic <- function(x, A1, k1, x01, A2, k2, x02) {
  A1 / (1 + exp(-k1 * (x - x01))) + A2 / (1 + exp(-k2 * (x - x02)))
}

A1 <- 0.1
k1 <- 60
x01 <- 14.5
A2 <- 0.1
k2 <- (-60) 
x02<- 18

plot(x=x,y=double_logistic(x, A1, k1, x01, A2, k2, x02))
plot(x=Lewis_silver$temp_bin,y=Lewis_silver$q_norm)

x<-seq(5,27,by=0.5)

# Example data
df <- data.frame(x = seq(1, 100, length.out = 100), y = double_logistic(seq(1, 100, length.out = 100), 1, 0.1, 20, 0.5, 0.2, 60) + rnorm(100, sd = 0.05))

# Fit the model
fit <- nls(y ~ double_logistic(x, A1, k1, x01, A2, k2, x02), data = df,
           start = list(A1 = 1, k1 = 0.1, x01 = 20, A2 = 0.5, k2 = 0.2, x02 = 60))


fit<-nls(data=Lewis_silver,q_norm~double_logistic(temp_bin, A1, k1, x01, A2, k2, x02),
       start = list(A1 = 0.1, k1 = 20, x01 = 12.5, A2 = 0.1, k2 = (-20), x02 = 18))

# Summary of the fit
summary(fit)

# Plot
plot(df$x, df$y, main = "Double Logistic Fit", xlab = "x", ylab = "y")
lines(df$x, predict(fit, newdata = df), col = "red")

# averaging ---------------------------------------------------------------
#what happens if you just simply average
cdf<-cdfs %>% dplyr::select(-c(Stage,n_temp,n))
cdfs_all <- rbind(cdf,cdfs_lewis)

#rename Lewis spp to match
cdfs_all<-cdfs_all %>%
  mutate(spp = replace(spp,spp == "Silver hake","Silver Hake")) %>%
  mutate(spp = replace(spp,spp == "American fourspot flounder","Fourspot Flounder")) %>%
  mutate(spp = replace(spp,spp == "Gulf Stream flounder","Gulfstream Flounder"))

spp_interest<-c("Silver Hake", "Fourspot Flounder", "Gulfstream Flounder")

test_avg <- cdfs_all %>% filter(spp %in% spp_interest) %>% group_by(spp,temp_bin) %>% mutate(mean_q=mean(q))

test_avg_silver<-test_avg %>% filter(spp == "Silver Hake")


# bootstrapping -----------------------------------------------------------
stations_filter <-stations_filter %>% filter(!Station_ID =="7.1.1", !is.na(SST)) %>% mutate(cruise_station = paste(Cruise_ID,Station_ID,sep="_")) %>%
  dplyr::select(SST,cruise_station)
silver_boot <- silver_tally %>% mutate(cruise_station = paste(Cruise_ID,Station_ID,sep="_")) %>% ungroup %>%
  dplyr::select(SilverHake,cruise_station) 
silver_boot<-left_join(stations_filter,silver_boot,by="cruise_station") 
silver_boot<-silver_boot %>% mutate(SilverHake = replace(SilverHake,is.na(SilverHake),0))

#bootstrapping runs
boot<-100
station_n<-length(stations_filter$cruise_station)
boot_out<-c()
for (j in 1:boot){
  newdata<-as.data.frame(matrix(nrow = station_n,ncol=ncol(silver_boot)))
  for(i in 1:station_n){
   rownums<-sample(1:station_n,length(station_n),replace = T) 
   newdata[i,]<-silver_boot[rownums,]
  }
  colnames(newdata)<-c("SST","cruise_station","n")
  newdata <- newdata %>% mutate(temp_bin = round(SST*2)/2)
  newdata_temp <- newdata %>% group_by(temp_bin) %>% tally() %>% rename(n_samples = n)
  newdata <- newdata %>% group_by(temp_bin) %>% mutate(n_fish = sum(n)) %>%
    dplyr::select(temp_bin,n_fish) %>% unique()
  newdata<-left_join(newdata,newdata_temp,by="temp_bin")
  new_q <- newdata %>% ungroup() %>% mutate(q = (n_fish/sum(n_fish)/(n_samples/sum(n_samples)))) %>% 
    dplyr::select(temp_bin,q)
  new_q$run_num<-j
  boot_out<-rbind(boot_out,new_q)
}

#get 5,95 CI
boot_summary <- boot_out %>% group_by(temp_bin) %>% 
  summarise(lower_CI = quantile(q,0.05),upper_CI = quantile(q,0.95)) 

boot_summary %>%
  ggplot(aes(x=temp_bin))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_point(data=silver_q,aes(x=temp_bin,y=q),color="magenta")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI))


boot_range <- boot_out %>% filter(q >=1) %>%
  group_by(run_num) %>% 
  summarise(lower = min(temp_bin),upper = max(temp_bin)) %>%
  mutate(range = upper-lower)

# bootstrapping q=1 line --------------------------------------------------

#modified from Hare, Richardson code
#assumming
silver_q<-cdfs %>% filter(spp == "Silver Hake")
bins<-temp_dist_NYOS$temp_bin

#make vectors to match matlab code
silver_tally <- silver_tally %>% filter(!is.na(SST))
abund_sub<-silver_tally$SilverHake
stations_filter <- stations_filter %>% filter(!is.na(SST))
env_sub<-stations_filter$SST

# Parameters# Parameterstemp_bin
boot <- 10000  # Number of bootstrap iterations (replace with your value)
ci_l <- 0.05 # Lower confidence interval level (replace with your value)
ci_u <- 0.95 # Upper confidence interval level (replace with your value)

# Initialize matrices to store results
q_rand <- matrix(0, nrow = boot, ncol = length(bins))
sta_bins_rand <- numeric(length(bins))
sta_bins_pro_rand <- numeric(length(bins))
abund_sub_bins_rand <- numeric(length(bins))
abund_sub_bins_pro_rand <- numeric(length(bins))

set.seed(19)

for (jj in seq_len(boot)) {
  # Generate random bootstrap samples
  bootstrap <- runif(length(abund_sub))
  
  # Sort bootstrap values and reorder abund_sub accordingly
  IX <- order(bootstrap)
  abund_sub_boot <- abund_sub[IX]
  
  for (yy in seq_len(length(bins))) {
    # Find indices within the bin range
    a <- which(env_sub > bins[yy] & env_sub <= bins[yy + 1])
    
    # Bin statistics
    sta_bins_rand[yy] <- length(a)
    sta_bins_pro_rand[yy] <- length(a) / length(env_sub)
    
    # Calculate bootstrap statistics based on 'type'
      abund_sub_bins_rand[yy] <- sum(abund_sub_boot[a])
      abund_sub_bins_pro_rand[yy] <- sum(abund_sub_boot[a]) / sum(abund_sub_boot)
  }
  
  # Calculate quotient
  q_rand[jj, ] <- abund_sub_bins_pro_rand / sta_bins_pro_rand
}

# Compute confidence intervals for the quotients
q_rand_ci_lower <- numeric(length(bins))
q_rand_ci_upper <- numeric(length(bins))

for (yy in seq_len(length(bins))) {
  if (any(is.nan(q_rand[, yy]))) {
    q_rand_ci_lower[yy] <- NA
    q_rand_ci_upper[yy] <- NA
  } else {
    # Calculate empirical cumulative distribution function (ECDF)
    ecdf_result <- ecdf(q_rand[, yy])
    
    # Interpolation for confidence intervals
    q_rand_ci_lower[yy] <- quantile(q_rand[, yy], ci_l)
    q_rand_ci_upper[yy] <- quantile(q_rand[, yy], ci_u)
  }
}

# Output results
results <- list(
  q_rand = q_rand,
  q_rand_ci_lower = q_rand_ci_lower,
  q_rand_ci_upper = q_rand_ci_upper
)

plot(x=temp_dist_NYOS$temp_bin,y=results[["q_rand_ci_lower"]])

print(results)

silver_test<-as.data.frame(cbind(results[["q_rand_ci_upper"]],silver_q$temp_bin,silver_q$q))
colnames(silver_test) <- c("upper","temp","q")
ggplot(silver_test)+
  geom_point(aes(x=temp,y=q))+
  geom_point(aes(x=temp,y=upper),col="blue")
