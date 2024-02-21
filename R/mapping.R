#Title: Barcoding Results Mapping

# Purpose: This script creates a map of the DNA barcoding results
#           for ichthyoplankton samples collected on NYOS cruises

# DataFiles: ''

# Author: S. Weisberg
# Contact details: sarah.j.weisberg@stonybrook.edu

# Tue Jun  7 10:32:55 2022 ------------------------------
library(RColorBrewer)
library(Di)

#map eggs
basemap+
  geom_point(data=eggs_map,mapping=aes(x=Longitude, y=Latitude, color = Species, size = n),position = "jitter")+
  scale_colour_brewer(palette = "Paired")

basemap+
  geom_point(data=eggs_map,mapping=aes(x=Longitude, y=Latitude, color = Common_Name, size = n), position = "jitter")+
  scale_color_manual(values = getPalette(colourCount)) 

#map larvae
basemap+
  geom_point(data=larvae_map,mapping=aes(x=Longitude, y=Latitude, color = Common_Name, size = n),position = "jitter")+
  scale_colour_brewer(palette = "Paired")

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
pie(rep(1, 25), col = c25)


# install.packages("Polychrome")
library(Polychrome)

# build-in color palette
unique<-length(unique(tallied_filter$Common_Name))
Glasbey = glasbey.colors(unique)
swatch(Glasbey)
