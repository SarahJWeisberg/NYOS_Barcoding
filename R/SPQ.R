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
library(patchwork)
library(paletteer)

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
           panel.grid.minor = element_line(colour="#f0f0f0"),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.spacing  = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0")
   ))
  
}

select <- dplyr::select
#choose color palette
my_pal <- paletteer_d("PrettyCols::Rainbow")

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
                                                                     ifelse(spp=="scoaqu_10m2","Windowpane Flounder",
                                                                     ifelse(spp=="pep_tri","Butterfish",spp))))))

eggs_lewis <- eggs_lewis %>% mutate(spp = ifelse(spp=="American fourspot flounder","Fourspot Flounder",
                                                 ifelse(spp=="Silver hake","Silver Hake",
                                                        ifelse(spp=="Gulf Stream flounder","Gulfstream Flounder",
                                                               ifelse(spp=="Windowpane", "Windowpane Flounder",
                                                               ifelse(spp=="Atlantic butterfish","Butterfish",spp))))))


ichthyo_all <- rbind(eggs_lewis,eggs_NYOS,larvae_ecomon)
ichthyo_all <- ichthyo_all %>% mutate(temp_bin=round(SST*2)/2)

#merge temperature proportion data into one dataframe
temp_dist_lewis <- temp_dist_lewis %>% select(temp_bin,prop_samples) %>% mutate(data="Lewis")
temp_dist_NYOS <- temp_dist_NYOS %>% select(temp_bin,prop_samples) %>% mutate(data="NYOS")
temp_dist_ecomon <- temp_dist_ecomon %>% select(temp_bin,prop_samples) %>% mutate(data="EcoMon")

temp_dist_all <- rbind(temp_dist_lewis,temp_dist_NYOS,temp_dist_ecomon)

# Single parameter quotient (SPQ) calculations --------------------------------------------------------
#select 4 species of interest
spp_interest<-c("Silver Hake","Fourspot Flounder","Gulfstream Flounder", "Butterfish","Windowpane Flounder")
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

#bootstrapping runs -- this take a long time
# set.seed(19)
# boot<-1000
# boot_out<-c()
# for (s in 1:length(datasets)){
#   stations_boot <- stations_all %>% filter(data == datasets[s])
#   ichthyo_boot <- ichthyo_spp %>% filter(data == datasets[s])
#   for (k in 1:length(spp_interest)) {
#     spp_boot <- ichthyo_boot %>% filter(spp == spp_interest[k]) %>% ungroup %>% select(cruise_station,abd)
#     boot_data <- left_join(stations_boot,spp_boot,by="cruise_station") %>%
#       mutate(abd = replace(abd,is.na(abd),0)) %>% select(-data)
#     for (j in 1:boot){
#       station_n <- nrow(boot_data)
#       newdata<-as.data.frame(matrix(nrow = station_n,ncol=ncol(boot_data)))
#       for(i in 1:station_n){
#         rownums<-sample(1:station_n,length(station_n),replace = T)
#         newdata[i,]<-boot_data[rownums,]
#       }
#       colnames(newdata)<-c("cruise_station","SST","abd")
#       newdata <- newdata %>% mutate(temp_bin = round(SST*2)/2)
#       newdata_temp <- newdata %>% group_by(temp_bin) %>% tally() %>% rename(n_samples = n)
#       newdata <- newdata %>% group_by(temp_bin) %>% mutate(n_fish = sum(abd)) %>%
#         dplyr::select(temp_bin,n_fish) %>% unique()
#       newdata<-left_join(newdata,newdata_temp,by="temp_bin")
#       new_q <- newdata %>% ungroup() %>% mutate(q = (n_fish/sum(n_fish)/(n_samples/sum(n_samples)))) %>%
#         dplyr::select(temp_bin,q) %>% mutate(data=datasets[s],spp=spp_interest[k])
#       boot_out<-rbind(boot_out,new_q)
#     }
#   }
# }
# 
# #save
# write.csv(boot_out,"boot_windowpane_1000x.csv")

#read in bootstrapping results 
boot_out <- read_csv("boot_1000x.csv")
window_out <- read_csv("boot_windowpane_1000x.csv")
boot_out <- boot_out[,-1]
window_out <- window_out[,-1]

boot_out <- rbind(boot_out,window_out)

#get 5,95 CI
boot_summary <- boot_out %>% group_by(temp_bin,spp,data) %>% 
  summarise(lower_CI = quantile(q,0.05,na.rm=T),upper_CI = quantile(q,0.95,na.rm=T)) 

#think about estimating thermal breadth from each run instead
#remove any q's less than 1
# boot_range <- boot_out %>% filter(q >=1) %>%
#   group_by(run_num) %>% 
#   summarise(lower = min(temp_bin),upper = max(temp_bin)) %>%
#   mutate(range = upper-lower)



# modeling ----------------------------------------------------------------
spp_cdf<-cdfs %>% filter(data == "EcoMon" & spp == "Silver Hake")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
test_predictors<-list(temp_bin=seq(min(temp_dist$temp_bin), max(temp_dist$temp_bin),0.5))
glm_vals <- predict(test_glm, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
conf_interval <- qnorm(0.975) * glm_vals$se.fit
glm_vals$lower <- glm_vals$fit - conf_interval
glm_vals$upper <- glm_vals$fit + conf_interval

#prep for plotting
glm_predict <- data.frame(temp = test_predictors[[1]], prop = glm_vals$fit,lower=glm_vals$lower,upper=glm_vals$upper)
glm_predict <- glm_predict %>% mutate(stage="Larvae")

#plot
ggplot(data=glm_predict, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm),color=my_pal[1],linewidth=1)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(title = "Windowpane Larvae")

spp_cdf<-cdfs %>% filter(!data == "EcoMon" & spp == "Silver Hake")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm2<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
test_predictors<-list(temp_bin=seq(min(spp_cdf$temp_bin), max(spp_cdf$temp_bin),0.5))
glm_vals2 <- predict(test_glm2, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
lower <- glm_vals2$fit - qnorm(0.975) * glm_vals2$se.fit
upper <- glm_vals2$fit + qnorm(0.975) * glm_vals2$se.fit

#prep for plotting
glm_predict2 <- data.frame(temp = test_predictors[[1]], prop = glm_vals2$fit,lower=lower,upper=upper)
glm_predict2 <- glm_predict2 %>% mutate(stage="Egg")

ggplot(data=glm_predict2, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm,color=data),linewidth=1)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(title = "Windowpane Eggs")

glm_out <- rbind(glm_predict,glm_predict2)

#plot
ggplot(data=glm_out, aes(x=temp,y=prop,color=stage))+
  geom_line(linewidth=2)+
  theme_Publication()+
  scale_color_manual(values = c(my_pal[8],my_pal[12]))+
  labs(title = "Windowpane")
Fig_3_top <- (Fig_3a + Fig_3b)
Fig_3_bot <- (Fig_3a / Fig_3b)| Fig_3c
Fig_3_top / Fig_3_bot  + plot_annotation(tag_levels = "A")

# plotting ----------------------------------------------------------------

# Figure 2 - SST histograms ----------------------------------------------------------------

Fig_2a <- stations_all %>% filter(!data == "EcoMon") %>%
  ggplot()+
  geom_histogram(aes(x=SST,fill=data),alpha=0.8,binwidth = 0.5)+
  scale_fill_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  geom_histogram(aes(x=SST,color=data,fill=NULL),alpha=0,binwidth = 0.5,show.legend = F)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]))+
  coord_cartesian(expand=F)+
  theme_Publication()+
  theme(legend.position = "inside",legend.position.inside = c(0.1,0.8))+
  xlim(c(-0.5,29))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Station Count")
Fig_2a

Fig_2b <- stations_all %>% filter(data == "EcoMon") %>%
  ggplot()+
  geom_histogram(aes(x=SST,fill=data),alpha=0.8,binwidth = 0.5)+
  scale_fill_manual(values = my_pal[1],labels="EcoMon Larvae")+
  geom_histogram(aes(x=SST,color=data,fill=NULL),alpha=0,binwidth = 0.5,show.legend = F)+
  scale_color_manual(values = my_pal[1])+
  coord_cartesian(expand=F)+
  theme_Publication()+
  theme(legend.position = "inside",legend.position.inside = c(0.1,0.875))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Station Count")
Fig_2b

Fig_2b / Fig_2a + plot_annotation(tag_levels = "A")


# Figure 3 - Silver Hake --------------------------------------------------
Fig_3a <- boot_summary %>%
  filter(!data == "EcoMon" & spp == "Silver Hake")%>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI),linewidth=1.5,alpha=0.5)+
  geom_point(data = cdfs %>% filter(!data == "EcoMon" & spp == "Silver Hake"),aes(x=temp_bin,y=q,color=data),size=2.25)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Single Parameter Quotient (SPQ)",
       title="Silver Hake Eggs")+
  theme_Publication()+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))

Fig_3b <- boot_summary %>%
  filter(data == "EcoMon" & spp == "Silver Hake")%>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI),linewidth=1.5,alpha=0.5,show.legend=F)+
  geom_point(data = cdfs %>% filter(data == "EcoMon" & spp == "Silver Hake"),aes(x=temp_bin,y=q,color=data),size=2.25,
             show.legend=F)+
  scale_color_manual(values = c(my_pal[1]))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Single Parameter Quotient (SPQ)",title="Silver Hake Larvae")+
  theme_Publication()

Fig_3_top <- Fig_3a+Fig_3b

#modeling larvae
spp_cdf<-cdfs %>% filter(data == "EcoMon" & spp == "Silver Hake")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
test_predictors<-list(temp_bin=seq(min(temp_dist$temp_bin), max(temp_dist$temp_bin),0.5))
glm_vals <- predict(test_glm, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
conf_interval <- qnorm(0.975) * glm_vals$se.fit
glm_vals$lower <- glm_vals$fit - conf_interval
glm_vals$upper <- glm_vals$fit + conf_interval

#prep for plotting
glm_predict <- data.frame(temp = test_predictors[[1]], prop = glm_vals$fit,lower=glm_vals$lower,upper=glm_vals$upper)
glm_predict <- glm_predict %>% mutate(stage="Larvae")

#plot
Fig_3d <- ggplot(data=glm_predict, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm),color=my_pal[1],linewidth=1)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Silver Hake Larvae")

#modeling eggs
spp_cdf<-cdfs %>% filter(!data == "EcoMon" & spp == "Silver Hake")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm2<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
test_predictors<-list(temp_bin=seq(min(spp_cdf$temp_bin), max(spp_cdf$temp_bin),0.5))
glm_vals2 <- predict(test_glm2, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
lower <- glm_vals2$fit - qnorm(0.975) * glm_vals2$se.fit
upper <- glm_vals2$fit + qnorm(0.975) * glm_vals2$se.fit

#prep for plotting
glm_predict2 <- data.frame(temp = test_predictors[[1]], prop = glm_vals2$fit,lower=lower,upper=upper)
glm_predict2 <- glm_predict2 %>% mutate(stage="Eggs")


Fig_3c <- ggplot(data=glm_predict2, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm,color=data),linewidth=1)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Silver Hake Eggs")+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))

#compare both models
glm_out <- rbind(glm_predict,glm_predict2)

#plot
Fig_3e <- ggplot(data=glm_out, aes(x=temp,y=prop,color=stage))+
  geom_line(linewidth=1.5)+
  geom_ribbon(aes(ymin = lower, ymax = upper,fill=stage), alpha = 0.4, linetype = 0,show.legend = F) +
  theme_Publication()+
  scale_color_manual(values = c(my_pal[8],my_pal[12]))+
  scale_fill_manual(values = c(my_pal[8],my_pal[12]))+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Silver Hake: Eggs & Larvae")

Fig_3_top <- (Fig_3a + Fig_3b)
Fig_3_bot <- (Fig_3c / Fig_3d)| Fig_3e
Fig_3 <- Fig_3_top / Fig_3_bot  + plot_annotation(tag_levels = "A")

ggsave("Figure_3.png", Fig_3, width=11, height=8.5,units=c("in"))

# Figure 4 - Windowpane--------------------------------------------------
Fig_4a <- boot_summary %>%
  filter(!data == "EcoMon" & spp == "Windowpane Flounder")%>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI),linewidth=1.5,alpha=0.5)+
  geom_point(data = cdfs %>% filter(!data == "EcoMon" & spp == "Windowpane Flounder"),aes(x=temp_bin,y=q,color=data),size=2.25)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Single Parameter Quotient (SPQ)",
       title="Windowpane Flounder Eggs")+
  theme_Publication()+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))

Fig_4b <- boot_summary %>%
  filter(data == "EcoMon" & spp == "Windowpane Flounder")%>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI),linewidth=1.5,alpha=0.5,show.legend=F)+
  geom_point(data = cdfs %>% filter(data == "EcoMon" & spp == "Windowpane Flounder"),aes(x=temp_bin,y=q,color=data),size=2.25,
             show.legend=F)+
  scale_color_manual(values = c(my_pal[1]))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Single Parameter Quotient (SPQ)",title="Windowpane Flounder Larvae")+
  theme_Publication()

Fig_4_top <- Fig_4a+Fig_4b

#modeling larvae
spp_cdf<-cdfs %>% filter(data == "EcoMon" & spp == "Windowpane Flounder")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
test_predictors<-list(temp_bin=seq(min(temp_dist$temp_bin), max(temp_dist$temp_bin),0.5))
glm_vals <- predict(test_glm, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
conf_interval <- qnorm(0.975) * glm_vals$se.fit
glm_vals$lower <- glm_vals$fit - conf_interval
glm_vals$upper <- glm_vals$fit + conf_interval

#prep for plotting
glm_predict <- data.frame(temp = test_predictors[[1]], prop = glm_vals$fit,lower=glm_vals$lower,upper=glm_vals$upper)
glm_predict <- glm_predict %>% mutate(stage="Larvae")

#plot
Fig_4d <- ggplot(data=glm_predict, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm),color=my_pal[1],linewidth=1)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Windowpane Flounder Larvae")

#modeling eggs
spp_cdf<-cdfs %>% filter(!data == "EcoMon" & spp == "Windowpane Flounder")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm2<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
test_predictors<-list(temp_bin=seq(min(spp_cdf$temp_bin), max(spp_cdf$temp_bin),0.5))
glm_vals2 <- predict(test_glm2, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
lower <- glm_vals2$fit - qnorm(0.975) * glm_vals2$se.fit
upper <- glm_vals2$fit + qnorm(0.975) * glm_vals2$se.fit

#prep for plotting
glm_predict2 <- data.frame(temp = test_predictors[[1]], prop = glm_vals2$fit,lower=lower,upper=upper)
glm_predict2 <- glm_predict2 %>% mutate(stage="Egg")


Fig_4c <- ggplot(data=glm_predict2, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm,color=data),linewidth=1)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Windowpane Flounder Eggs")+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))

#compare both models
glm_out <- rbind(glm_predict,glm_predict2)

#plot
Fig_4e <- ggplot(data=glm_out, aes(x=temp,y=prop,color=stage))+
  geom_line(linewidth=1.5)+
  geom_ribbon(aes(ymin = lower, ymax = upper,fill=stage), alpha = 0.4, linetype = 0,show.legend = F) +
  theme_Publication()+
  scale_color_manual(values = c(my_pal[8],my_pal[12]))+
  scale_fill_manual(values = c(my_pal[8],my_pal[12]))+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Windowpane Flounder: Eggs & Larvae")

Fig_4_top <- (Fig_4a + Fig_4b)
Fig_4_bot <- (Fig_4c / Fig_4d)| Fig_4e
Fig_4 <- Fig_4_top / Fig_4_bot  + plot_annotation(tag_levels = "A")

ggsave("Figure_4.png", Fig_4, width=11, height=8.5,units=c("in"))

# Figure 5 - Gulfstream--------------------------------------------------
Fig_5a <- boot_summary %>%
  filter(!data == "EcoMon" & spp == "Gulfstream Flounder")%>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI),linewidth=1.5,alpha=0.5)+
  geom_point(data = cdfs %>% filter(!data == "EcoMon" & spp == "Gulfstream Flounder"),aes(x=temp_bin,y=q,color=data),size=2.25)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Single Parameter Quotient (SPQ)",
       title="Gulf Stream Flounder Eggs")+
  theme_Publication()+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))

Fig_5b <- boot_summary %>%
  filter(data == "EcoMon" & spp == "Gulfstream Flounder")%>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI),linewidth=1.5,alpha=0.5,show.legend=F)+
  geom_point(data = cdfs %>% filter(data == "EcoMon" & spp == "Gulfstream Flounder"),aes(x=temp_bin,y=q,color=data),size=2.25,
             show.legend=F)+
  scale_color_manual(values = c(my_pal[1]))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Single Parameter Quotient (SPQ)",title="Gulf Stream Flounder Larvae")+
  theme_Publication()

Fig_5_top <- Fig_5a+Fig_5b

#modeling larvae
test_predictors<-list(temp_bin=seq(-0.5,30,0.5))
spp_cdf<-cdfs %>% filter(data == "EcoMon" & spp == "Gulfstream Flounder")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
glm_vals <- predict(test_glm, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
conf_interval <- qnorm(0.975) * glm_vals$se.fit
glm_vals$lower <- glm_vals$fit - conf_interval
glm_vals$upper <- glm_vals$fit + conf_interval

#prep for plotting
glm_predict <- data.frame(temp = test_predictors[[1]], prop = glm_vals$fit,lower=glm_vals$lower,upper=glm_vals$upper)
glm_predict <- glm_predict %>% mutate(stage="Larvae")

#plot
Fig_5d <- ggplot(data=glm_predict, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm),color=my_pal[1],linewidth=1)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Gulf Stream Flounder Larvae")

#modeling eggs
spp_cdf<-cdfs %>% filter(!data == "EcoMon" & spp == "Gulfstream Flounder")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm2<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
glm_vals2 <- predict(test_glm2, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
lower <- glm_vals2$fit - qnorm(0.975) * glm_vals2$se.fit
upper <- glm_vals2$fit + qnorm(0.975) * glm_vals2$se.fit

#prep for plotting
glm_predict2 <- data.frame(temp = test_predictors[[1]], prop = glm_vals2$fit,lower=lower,upper=upper)
glm_predict2 <- glm_predict2 %>% mutate(stage="Egg")


Fig_5c <- ggplot(data=glm_predict2, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm,color=data),linewidth=1)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Gulf Stream Flounder Eggs")+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))

#compare both models
glm_out <- rbind(glm_predict,glm_predict2)

#plot
Fig_5e <- ggplot(data=glm_out, aes(x=temp,y=prop,color=stage))+
  geom_line(linewidth=1.5)+
  geom_ribbon(aes(ymin = lower, ymax = upper,fill=stage), alpha = 0.4, linetype = 0,show.legend = F) +
  theme_Publication()+
  scale_color_manual(values = c(my_pal[8],my_pal[12]))+
  scale_fill_manual(values = c(my_pal[8],my_pal[12]))+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Gulf Stream Flounder: Eggs & Larvae")

Fig_5_top <- (Fig_5a + Fig_5b)
Fig_5_bot <- (Fig_5c / Fig_5d)| Fig_5e
Fig_5 <- Fig_5_top / Fig_5_bot  + plot_annotation(tag_levels = "A")

ggsave("Figure_5.png", Fig_5, width=11, height=8.5,units=c("in"))

# Figure 6 - Fourspot--------------------------------------------------
Fig_6a <- boot_summary %>%
  filter(!data == "EcoMon" & spp == "Fourspot Flounder")%>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI),linewidth=1.5,alpha=0.5)+
  geom_point(data = cdfs %>% filter(!data == "EcoMon" & spp == "Fourspot Flounder"),aes(x=temp_bin,y=q,color=data),size=2.25)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Single Parameter Quotient (SPQ)",
       title="Fourspot Flounder Eggs")+
  theme_Publication()+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))

Fig_6b <- boot_summary %>%
  filter(data == "EcoMon" & spp == "Fourspot Flounder")%>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI),linewidth=1.5,alpha=0.5,show.legend=F)+
  geom_point(data = cdfs %>% filter(data == "EcoMon" & spp == "Fourspot Flounder"),aes(x=temp_bin,y=q,color=data),size=2.25,
             show.legend=F)+
  scale_color_manual(values = c(my_pal[1]))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Single Parameter Quotient (SPQ)",title="Fourspot Flounder Larvae")+
  theme_Publication()

Fig_6_top <- Fig_6a+Fig_6b

#modeling larvae
test_predictors<-list(temp_bin=seq(-0.5,30,0.5))
spp_cdf<-cdfs %>% filter(data == "EcoMon" & spp == "Fourspot Flounder")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
glm_vals <- predict(test_glm, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
conf_interval <- qnorm(0.975) * glm_vals$se.fit
glm_vals$lower <- glm_vals$fit - conf_interval
glm_vals$upper <- glm_vals$fit + conf_interval

#prep for plotting
glm_predict <- data.frame(temp = test_predictors[[1]], prop = glm_vals$fit,lower=glm_vals$lower,upper=glm_vals$upper)
glm_predict <- glm_predict %>% mutate(stage="Larvae")

#plot
Fig_6d <- ggplot(data=glm_predict, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm),color=my_pal[1],linewidth=1)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Fourspot Flounder Larvae")

#modeling eggs
spp_cdf<-cdfs %>% filter(!data == "EcoMon" & spp == "Fourspot Flounder")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm2<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
glm_vals2 <- predict(test_glm2, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
lower <- glm_vals2$fit - qnorm(0.975) * glm_vals2$se.fit
upper <- glm_vals2$fit + qnorm(0.975) * glm_vals2$se.fit

#prep for plotting
glm_predict2 <- data.frame(temp = test_predictors[[1]], prop = glm_vals2$fit,lower=lower,upper=upper)
glm_predict2 <- glm_predict2 %>% mutate(stage="Egg")


Fig_6c <- ggplot(data=glm_predict2, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm,color=data),linewidth=1)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Fourspot Flounder Eggs")+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))

#compare both models
glm_out <- rbind(glm_predict,glm_predict2)

#plot
Fig_6e <- ggplot(data=glm_out, aes(x=temp,y=prop,color=stage))+
  geom_line(linewidth=1.5)+
  geom_ribbon(aes(ymin = lower, ymax = upper,fill=stage), alpha = 0.4, linetype = 0,show.legend = F) +
  theme_Publication()+
  scale_color_manual(values = c(my_pal[8],my_pal[12]))+
  scale_fill_manual(values = c(my_pal[8],my_pal[12]))+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Fourspot Flounder: Eggs & Larvae")

Fig_6_top <- (Fig_6a + Fig_6b)
Fig_6_bot <- (Fig_6c / Fig_6d)| Fig_6e
Fig_6 <- Fig_6_top / Fig_6_bot  + plot_annotation(tag_levels = "A")

ggsave("Figure_6.png", Fig_6, width=11, height=8.5,units=c("in"))

# Figure 7 - Butterfish--------------------------------------------------
Fig_7a <- boot_summary %>%
  filter(!data == "EcoMon" & spp == "Butterfish")%>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI),linewidth=1.5,alpha=0.5)+
  geom_point(data = cdfs %>% filter(!data == "EcoMon" & spp == "Butterfish"),aes(x=temp_bin,y=q,color=data),size=2.25)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Single Parameter Quotient (SPQ)",
       title="Butterfish Eggs")+
  theme_Publication()+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))

Fig_7b <- boot_summary %>%
  filter(data == "EcoMon" & spp == "Butterfish")%>%
  ggplot(aes(x=temp_bin,colour = data))+
  geom_hline(yintercept = 1,linetype="dashed")+
  geom_errorbar(aes(ymin = lower_CI,ymax =upper_CI),linewidth=1.5,alpha=0.5,show.legend=F)+
  geom_point(data = cdfs %>% filter(data == "EcoMon" & spp == "Butterfish"),aes(x=temp_bin,y=q,color=data),size=2.25,
             show.legend=F)+
  scale_color_manual(values = c(my_pal[1]))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="Single Parameter Quotient (SPQ)",title="Butterfish Larvae")+
  theme_Publication()

Fig_7_top <- Fig_7a+Fig_7b

#modeling larvae
test_predictors<-list(temp_bin=seq(-0.5,30,0.5))
spp_cdf<-cdfs %>% filter(data == "EcoMon" & spp == "Butterfish")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
glm_vals <- predict(test_glm, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
conf_interval <- qnorm(0.975) * glm_vals$se.fit
glm_vals$lower <- glm_vals$fit - conf_interval
glm_vals$upper <- glm_vals$fit + conf_interval

#prep for plotting
glm_predict <- data.frame(temp = test_predictors[[1]], prop = glm_vals$fit,lower=glm_vals$lower,upper=glm_vals$upper)
glm_predict <- glm_predict %>% mutate(stage="Larvae")

#plot
Fig_7d <- ggplot(data=glm_predict, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm),color=my_pal[1],linewidth=1)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Butterfish Larvae")

#modeling eggs
spp_cdf<-cdfs %>% filter(!data == "EcoMon" & spp == "Butterfish")
spp_cdf <- spp_cdf %>% mutate(round_cdf = round(cdf_q_norm,digits = 4))
test_glm2<-glm(data=spp_cdf,formula=round_cdf~temp_bin, family = "quasibinomial")
glm_vals2 <- predict(test_glm2, newdata=test_predictors, type = "response", se.fit=T)

# Calculate the upper and lower bounds for the confidence intervals
lower <- glm_vals2$fit - qnorm(0.975) * glm_vals2$se.fit
upper <- glm_vals2$fit + qnorm(0.975) * glm_vals2$se.fit

#prep for plotting
glm_predict2 <- data.frame(temp = test_predictors[[1]], prop = glm_vals2$fit,lower=lower,upper=upper)
glm_predict2 <- glm_predict2 %>% mutate(stage="Egg")


Fig_7c <- ggplot(data=glm_predict2, aes(x=temp,y=prop))+
  geom_line(linewidth=1)+
  geom_line(data=spp_cdf,aes(x=temp_bin,y=cdf_q_norm,color=data),linewidth=1)+
  scale_color_manual(values = c(my_pal[10],my_pal[4]),labels=c("EcoMon Eggs","NYOS Eggs"))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4) +
  theme_Publication()+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Butterfish Eggs")+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))

#compare both models
glm_out <- rbind(glm_predict,glm_predict2)

#plot
Fig_7e <- ggplot(data=glm_out, aes(x=temp,y=prop,color=stage))+
  geom_line(linewidth=1.5)+
  geom_ribbon(aes(ymin = lower, ymax = upper,fill=stage), alpha = 0.4, linetype = 0,show.legend = F) +
  theme_Publication()+
  scale_color_manual(values = c(my_pal[8],my_pal[12]))+
  scale_fill_manual(values = c(my_pal[8],my_pal[12]))+
  theme(legend.position = "inside",legend.position.inside = c(0.15,0.8))+
  labs(x = expression(paste("Sea Surface Temperature (",degree,"C)")),y="SPQ Cumulative Fraction",title="Butterfish: Eggs & Larvae")

Fig_7_top <- (Fig_7a + Fig_7b)
Fig_7_bot <- (Fig_7c / Fig_7d)| Fig_7e
Fig_7 <- Fig_7_top / Fig_7_bot  + plot_annotation(tag_levels = "A")

ggsave("Figure_7.png", Fig_7, width=11, height=8.5,units=c("in"))

