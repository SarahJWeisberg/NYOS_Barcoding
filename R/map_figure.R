
#code to create Figure 1 for Chapter 4

#need to run 
source(here("R/Lewis_data.R")) #egg data from Lewis_et_al_2016

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
map_pal <- paletteer_d("nbapalettes::warriors_city2")


#run NYOS mapping code
#data_paper.Rmd
Fig_1c  <- 

stations_lewis <- db_trim %>% select(cruise_station,SST) %>% mutate(data="Lewis")
#load in countries for plotting 
world <- ne_countries(scale = "medium", returnclass = "sf")

ecomon <- read.csv("data/EcoMon_Plankton_Data.csv")

#convert data and create year, month, day columns: 

ecomon$date <- as.Date(ecomon$date, format =  "%d-%b-%y") #b here is the code for month when its words 
ecomon$month <- month(ecomon$date)
ecomon$year <- year(ecomon$date)
ecomon$day <- day(ecomon$date)

ecomon$ID_points <- as.character(rownames(ecomon))


#make ecomon spatial with the same defined coordinate system
ecomon_sp <- st_as_sf(x = ecomon, 
                      coords = c("lon", "lat"),
                      crs = 4326)

ecomon_sp <- ecomon_sp %>% mutate(season = ifelse(month %in% c(12,1,2),"winter",
                                                ifelse(month %in% c(3,4,5), "spring",
                                                       ifelse(month %in% c(6,7,8),"summer",
                                                              "fall"))))

ecomon <- ecomon %>% mutate(cruise_station = paste(cruise_name,station,sep = "_"))

lewis_sp <- left_join(stations_lewis,ecomon,by="cruise_station") %>% filter(!is.na(lat))

lewis_sp <- st_as_sf(x = lewis_sp, 
                      coords = c("lon", "lat"),
                      crs = 4326)
lewis_sp <- lewis_sp %>% mutate(season = ifelse(month %in% c(12,1,2),"winter",
                                                ifelse(month %in% c(3,4,5), "spring",
                                                       ifelse(month %in% c(6,7,8),"summer",
                                                              "fall"))))
# #Call in bathymetric data
# #This part of script converts the bathymetric contours and accompanying data into a cropped dataframe. 
bathy <- getNOAA.bathy(-80, -65, 33.5, 46, resolution=1); bathydf <- as.xyz(bathy) 
bathydf$col <- NA
bathydf$col <- ifelse(bathydf$V3 >= -100, ifelse(bathydf$V3 < -15, "b","a"),"c")
bathydf$col <- factor(bathydf$col, levels=c("a",""))
# ####################
# #Separate land from ocean.
# #This part of script helps to properly overlay the coastline data on the background bathymetric layer since `marmap` can only provide very
# #coarse coastline resolution. 
bathy_sea <- bathydf; bathy_sea$V3[bathy_sea$V3 > 1] <- NA

ecomon_egg_map <- ggplot() + geom_sf(data = world) + 
  geom_raster(data=bathy_sea, aes(x=V1, y=V2),fill="#f1f5f8") +
  #   #This creates bathymetric contour lines at: 50, 100, 500, 1000, 2000, and 3000 m depth thresholds
  # geom_contour(data=bathy_sea, aes(x= V1, y=V2, z=V3), 
  #              breaks=c(0,-50,-100,-500,-1000,-2000,-3000), colour="black",linewidth = 0.1) +
  geom_sf(data = lewis_sp, size = 0.5,aes(color=season),show.legend = F) + 
  scale_color_manual(values = c(map_pal[3],map_pal[4],map_pal[2],map_pal[1]))+
  geom_sf(data = world,fill="lightgrey") +
  coord_sf(xlim=c(-78, -65), ylim=c(35,45), expand = F)+
  # geom_text(aes(x=x, y=y, label=text),
  #           data=data.frame(x=c(rep(-65.25,4)),
  #                           y=c(42.45, 42,41.5, 41),
  #                           text=c("100 m","1000 m","2000 m","3000 m")),
            #color=c(rep("black",4)),size=2.5)+
  labs(x=NULL,y=NULL)+
  #labs(title = "EcoMon Eggs")+
  theme_Publication()+
  theme(axis.text.x=element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  annotate("text",x=-70,y=36,label = "Stations with Egg Data",size=4)
ecomon_egg_map

Fig_1a <- ggplot() + geom_sf(data = world) + 
  geom_raster(data=bathy_sea, aes(x=V1, y=V2),fill="#f1f5f8") +
  #   #This creates bathymetric contour lines at: 50, 100, 500, 1000, 2000, and 3000 m depth thresholds
  geom_contour(data=bathy_sea, aes(x= V1, y=V2, z=V3), 
               breaks=c(0,-50,-100,-500,-1000,-2000,-3000), colour="black",linewidth = 0.1) +
  geom_sf(data = ecomon_sp, size = 0.5,aes(color=season))+
  geom_sf(data = world,fill="lightgrey") +
  coord_sf(xlim=c(-78, -65), ylim=c(34,45), expand = F)+
  labs(x=expression(paste("Longitude (",degree,"W)")),
       y=expression(paste("Latitude (",degree,"N)")))+
  scale_color_manual(values = c(map_pal[2],map_pal[3],map_pal[4],map_pal[1]))+
  # geom_text(aes(x=x, y=y, label=text),
  #           data=data.frame(x=c(rep(-65.25,3)),
  #                           y=c(42,41.5, 41),
  #                           text=c("1000 m","2000 m","3000 m")),
            #color=c(rep("black",3)),size=2.5)+
  theme_Publication()+
  theme(legend.position = "bottom")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  annotation_custom(ggplotGrob(ecomon_egg_map),xmin=-72, xmax=-65, ymin=34, ymax=39)+
  guides(color = guide_legend(override.aes = list(size=4)))+
  ggtitle("EcoMon Data")
Fig_1a+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

# Fig_1 <- Fig_1a / Fig_1b &theme(legend.position = "bottom")
# Fig_1 <- Fig_1 + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')
# Fig_1
ggsave(filename = "Figure1a.png", plot=Fig_1a, width = 8.5, height=11, units=c("in"))
ggsave(filename = "Figure1b.png", plot=Fig_1b, width = 8.5, height=11, units=c("in"))

