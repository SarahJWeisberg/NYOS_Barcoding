# Wed Aug  7 13:52:00 2024 ------------------------------

#read in data
fish<-read.csv(here("data/Barcoding_Results_Full.csv"))
stations<-read.csv(here("data/Station_Info.csv"))

#convert Cruise_ID, Station_IN to factor
stations$Cruise_ID<-as.factor(stations$Cruise_ID)
fish$Cruise.ID<-as.factor(fish$Cruise.ID)
stations$Station_ID<-as.factor(stations$Station_ID)
fish$Station.ID<-as.factor(fish$Station.ID)

#use CTD SST data
stations <- stations %>% rename(SST = SST_CTD) %>% mutate(SST = as.numeric(SST)) %>%
  drop_na(SST)
  
#filter just SST data
stations_filter <- stations %>% dplyr::select(Cruise_ID,Station_ID,SST) %>% 
  mutate(cruise_station = paste(Cruise_ID,Station_ID,sep = "_")) %>%
  select(cruise_station,SST)

#remove all columns but genus, species, common name, station, type
fish_filter<- fish %>% 
  dplyr::select(c("Genus.species","Common_Name","Cruise.ID","Station.ID","Stage")) %>% 
  rename("Cruise_ID" = "Cruise.ID") %>% rename("Station_ID"="Station.ID")
  
#filter out 2103, juveniles
fish_filter<-fish_filter %>% filter(Cruise_ID != "2103") %>% filter(Stage!="Juvenile")


#plot temp data - all
stations %>%
  ggplot(aes(x=SST)) +
  geom_histogram(binwidth = 0.5, alpha=0.5, color="#A58AFF",fill="#A58AFF")

#tally all
fish_tally <- fish_filter %>% group_by(Common_Name,Stage) %>% tally()

#single species example: silver hake
# #silver hake eggs
# silver<-fish_filter %>% filter(Common_Name == "Silver Hake" & Stage == "Egg") 
# silver_tally<-silver %>% group_by(Cruise_ID,Station_ID) %>% tally()
# 
# silver_tally <- left_join(stations_filter,silver_tally,by=c("Cruise_ID","Station_ID")) %>%
#   rename("SilverHake" = "n") %>% mutate(SilverHake = ifelse(is.na(SilverHake),0,SilverHake))
# 
# ggplot() +
#   geom_histogram(data=stations,aes(x=SST, color=Cruise_ID, fill=Cruise_ID),binwidth = 0.5, alpha=0.5)+
#   geom_histogram(data=silver,aes(x=SST),binwidth = 0.5, alpha=0.5,color="black")

#SST distribution
temp_dist_NYOS <- stations_filter %>% mutate(temp_bin = round(SST*2)/2)

bins<- temp_dist_NYOS %>% dplyr::select(temp_bin) %>% filter(!is.na(temp_bin)) %>% unique()

#plotstations
# temp_dist_NYOS%>%
# ggplot(aes(x=SST, color=Cruise_ID, fill=Cruise_ID)) +
#   geom_histogram(binwidth = 0.5, alpha=0.5)

#tally, calculate proportion
temp_dist_NYOS <- temp_dist_NYOS %>% group_by(temp_bin) %>% tally %>% filter(is.na(temp_bin) == F) %>% 
  mutate(total=sum(n),prop_samples=n/total) %>% dplyr::select(temp_bin,prop_samples,n) %>% rename(n_temp = n)

#generate clean dataframe
fish_filter<- fish_filter %>% mutate(cruise_station = paste(Cruise_ID,Station_ID,sep = "_")) %>%
  filter(Stage == "Egg", !Common_Name =="") %>% select(Common_Name,cruise_station) %>% group_by(cruise_station, Common_Name) %>% tally()

#select for spp of interest
fish_filter <- fish_filter %>% filter(Common_Name %in% c("Silver Hake","Gulfstream Flounder","Fourspot Flounder","Butterfish"))
fish_filter <- left_join(fish_filter,stations_filter,by="cruise_station")

rm(fish,fish_tally,stations,bins)


# modeling ----------------------------------------------------------------


# #double logistic fitting
# library(minpack.lm)
# double_logistic <- function(x, A1, k1, x01, A2, k2, x02) {
#   A1 / (1 + exp(-k1 * (x - x01))) + A2 / (1 + exp(-k2 * (x - x02)))
# }
# 
# A1 <- 0.1
# k1 <- 60
# x01 <- 14.5
# A2 <- 0.1
# k2 <- (-60) 
# x02<- 18
# 
# plot(x=x,y=double_logistic(x, A1, k1, x01, A2, k2, x02))
# plot(x=Lewis_silver$temp_bin,y=Lewis_silver$q_norm)
# 
# x<-seq(5,27,by=0.5)
# 
# # Example data
# df <- data.frame(x = seq(1, 100, length.out = 100), y = double_logistic(seq(1, 100, length.out = 100), 1, 0.1, 20, 0.5, 0.2, 60) + rnorm(100, sd = 0.05))
# 
# # Fit the model
# fit <- nls(y ~ double_logistic(x, A1, k1, x01, A2, k2, x02), data = df,
#            start = list(A1 = 1, k1 = 0.1, x01 = 20, A2 = 0.5, k2 = 0.2, x02 = 60))
# 
# 
# fit<-nls(data=Lewis_silver,q_norm~double_logistic(temp_bin, A1, k1, x01, A2, k2, x02),
#        start = list(A1 = 0.1, k1 = 20, x01 = 12.5, A2 = 0.1, k2 = (-20), x02 = 18))
# 
# # Summary of the fit
# summary(fit)
# 
# # Plot
# plot(df$x, df$y, main = "Double Logistic Fit", xlab = "x", ylab = "y")
# lines(df$x, predict(fit, newdata = df), col = "red")


# bootstrapping q=1 line --------------------------------------------------

# #modified from Hare, Richardson code
# silver_q<-cdfs %>% filter(spp == "Silver Hake")
# bins<-temp_dist_NYOS$temp_bin
# 
# #make vectors to match matlab code
# silver_tally <- silver_tally %>% filter(!is.na(SST))
# abund_sub<-silver_tally$SilverHake
# stations_filter <- stations_filter %>% filter(!is.na(SST))
# env_sub<-stations_filter$SST
# 
# # Parameters# Parameterstemp_bin
# boot <- 10000  # Number of bootstrap iterations (replace with your value)
# ci_l <- 0.05 # Lower confidence interval level (replace with your value)
# ci_u <- 0.95 # Upper confidence interval level (replace with your value)
# 
# # Initialize matrices to store results
# q_rand <- matrix(0, nrow = boot, ncol = length(bins))
# sta_bins_rand <- numeric(length(bins))
# sta_bins_pro_rand <- numeric(length(bins))
# abund_sub_bins_rand <- numeric(length(bins))
# abund_sub_bins_pro_rand <- numeric(length(bins))
# 
# set.seed(19)
# 
# for (jj in seq_len(boot)) {
#   # Generate random bootstrap samples
#   bootstrap <- runif(length(abund_sub))
#   
#   # Sort bootstrap values and reorder abund_sub accordingly
#   IX <- order(bootstrap)
#   abund_sub_boot <- abund_sub[IX]
#   
#   for (yy in seq_len(length(bins))) {
#     # Find indices within the bin range
#     a <- which(env_sub > bins[yy] & env_sub <= bins[yy + 1])
#     
#     # Bin statistics
#     sta_bins_rand[yy] <- length(a)
#     sta_bins_pro_rand[yy] <- length(a) / length(env_sub)
#     
#     # Calculate bootstrap statistics based on 'type'
#       abund_sub_bins_rand[yy] <- sum(abund_sub_boot[a])
#       abund_sub_bins_pro_rand[yy] <- sum(abund_sub_boot[a]) / sum(abund_sub_boot)
#   }
#   
#   # Calculate quotient
#   q_rand[jj, ] <- abund_sub_bins_pro_rand / sta_bins_pro_rand
# }
# 
# # Compute confidence intervals for the quotients
# q_rand_ci_lower <- numeric(length(bins))
# q_rand_ci_upper <- numeric(length(bins))
# 
# for (yy in seq_len(length(bins))) {
#   if (any(is.nan(q_rand[, yy]))) {
#     q_rand_ci_lower[yy] <- NA
#     q_rand_ci_upper[yy] <- NA
#   } else {
#     # Calculate empirical cumulative distribution function (ECDF)
#     ecdf_result <- ecdf(q_rand[, yy])
#     
#     # Interpolation for confidence intervals
#     q_rand_ci_lower[yy] <- quantile(q_rand[, yy], ci_l)
#     q_rand_ci_upper[yy] <- quantile(q_rand[, yy], ci_u)
#   }
# }
# 
# # Output results
# results <- list(
#   q_rand = q_rand,
#   q_rand_ci_lower = q_rand_ci_lower,
#   q_rand_ci_upper = q_rand_ci_upper
# )
# 
# plot(x=temp_dist_NYOS$temp_bin,y=results[["q_rand_ci_lower"]])
# 
# print(results)
# 
# silver_test<-as.data.frame(cbind(results[["q_rand_ci_upper"]],silver_q$temp_bin,silver_q$q))
# colnames(silver_test) <- c("upper","temp","q")
# ggplot(silver_test)+
#   geom_point(aes(x=temp,y=q))+
#   geom_point(aes(x=temp,y=upper),col="blue")
