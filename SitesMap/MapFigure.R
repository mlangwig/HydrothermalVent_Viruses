## Modified from original author, Dr. Mirna Vazquez Rosas Landa
## email: mirnavrl@austin.utexas.edu

#read input
BR <- readr::read_tsv("input.tsv")
input <- readr::read_tsv("input_VentPlume.txt")

#raster_df <- readr::read_tsv("salinity.txt")

#libraries
library(ggmap)
library(RColorBrewer)
library(pals)
library(sf)
library(dplyr)

# read in shape file
#sf_map <- st_read("data/water_bodies_carto.shp") #read in shape file. Pay attention to CRS for coordinate system
#sf_map <- sf_map %>% st_transform(crs = 4326) #make coordinate ref system lat / long

#convert lat and long & salinity to numeric for plotting
input$Latitude<-as.numeric(input$Latitude)
input$Longitude<-as.numeric(input$Longitude)
#BR$AvgSalinity<-as.numeric(BR$AvgSalinity)

#color scheme
#mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
mycolors <- pals::kelly(n=7)

# map = shape file...extension .shp
# raster = salinity...extension .geotiff
# brick for average

dev.off()
#plotting
mapworld <- borders("world",colour = "gray80",fill="white") 
mp<-ggplot() + mapworld + ylim(-60,90) 
overview <- mp + 
  geom_point(data = input, aes(x = Longitude, y = Latitude, color=factor(Site2),
                            shape=Site)) +
  scale_color_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1",
                              "#8c510a", "#000000")) + #choose colors of sites
  #scale_size(limits = c(2, 29), breaks = c(2, 10, 20, 27, 28)) +
  scale_shape_manual(values = c(15,16,17)) +
  theme_gray() +
  labs(x= "Longitude", y = "Latitude") + #Lat and Long labels
  guides(color = guide_legend(title = "Site")) + #change legend title
  #coord_cartesian(xlim=c(-123,-121.5),ylim=c(37,38.5)) + #set limits of map
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12),
        legend.key=element_rect(fill="white"),
        legend.title = element_blank(),
        panel.border = element_rect(color = "black", 
                                    fill = NA, 
                                    linewidth = .5))
overview

ggsave(overview, filename = "Overview_Map_PlumeVent.png", width = 12, height = 6)

# Zoom in Lau

dat_Lau <- BR[BR$Site2=="Eastern Lau Spreading Center",]

Lau <- mp + 
  geom_point(data = dat_Lau, aes(x = Longitude, y = Latitude, color=factor(Site))) +
  #scale_color_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1")) + #choose colors of sites
  #scale_size(limits = c(2, 29), breaks = c(2, 10, 20, 27, 28)) +
  theme_gray() +
  labs(x= "Longitude", y = "Latitude") + #Lat and Long labels
  ggtitle("Eastern Lau Spreading Center") +
  guides(color = guide_legend(title = "Site")) + #change legend title
  coord_cartesian(xlim=c(-176.61,-176.18),ylim=c(-20.7,-22.3)) + #set limits of map
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 6))
Lau

# Zoom in Brothers

dat_BV <- BR[BR$Site2=="Brothers volcano",]

BV <- mp + 
  geom_point(data = dat_BV, aes(x = Longitude, y = Latitude, color=factor(Site))) +
  #scale_color_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1")) + #choose colors of sites
  #scale_size(limits = c(2, 29), breaks = c(2, 10, 20, 27, 28)) +
  theme_gray() +
  labs(x= "Longitude", y = "Latitude") + #Lat and Long labels
  ggtitle("Brothers volcano") +
  guides(color = guide_legend(title = "Site")) + #change legend title
  coord_cartesian(xlim=c(179.075,179.05),ylim=c(-34.89,-34.85)) + #set limits of map
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 6))
BV

# Zoom in Guaymas

dat_Guay <- BR[BR$Site2=="Guaymas Basin",]

Guay <- mp + 
  geom_point(data = dat_Guay, aes(x = Longitude, y = Latitude, color=factor(Site))) +
  #scale_color_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1")) + #choose colors of sites
  #scale_size(limits = c(2, 29), breaks = c(2, 10, 20, 27, 28)) +
  theme_gray() +
  labs(x= "Longitude", y = "Latitude") + #Lat and Long labels
  ggtitle("Guaymas Basin") +
  guides(color = guide_legend(title = "Site")) + #change legend title
  coord_cartesian(xlim=c(-111.41,-111.405),ylim=c(27.0, 27.015)) + #set limits of map
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 6))
Guay

# Zoom in EPR

dat_EPR <- BR[BR$Site2=="East Pacific Rise",]

EPR <- mp + 
  geom_point(data = dat_EPR, aes(x = Longitude, y = Latitude, color=factor(Site))) +
  #scale_color_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1")) + #choose colors of sites
  #scale_size(limits = c(2, 29), breaks = c(2, 10, 20, 27, 28)) +
  theme_gray() +
  labs(x= "Longitude", y = "Latitude") + #Lat and Long labels
  ggtitle("East Pacific Rise") +
  guides(color = guide_legend(title = "Site")) + #change legend title
  coord_cartesian(xlim=c(-104.295,-104.275),ylim=c(9.85, 9.7)) + #set limits of map
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 6))
EPR

# Zoom in MAR

dat_MAR <- BR[BR$Site2=="Mid-Atlantic Ridge",]

MAR <- mp + 
  geom_point(data = dat_MAR, aes(x = Longitude, y = Latitude, color=factor(Site))) +
  #scale_color_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1")) + #choose colors of sites
  #scale_size(limits = c(2, 29), breaks = c(2, 10, 20, 27, 28)) +
  theme_gray() +
  labs(x= "Longitude", y = "Latitude") + #Lat and Long labels
  ggtitle("Mid-Atlantic Ridge") +
  guides(color = guide_legend(title = "Site")) + #change legend title
  coord_cartesian(xlim=c(-33.91,-32.2),ylim=c(37.3, 36.1)) + #set limits of map
  theme(axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 6))
MAR

library(patchwork)

patch_plot <- overview/(Guay + EPR + MAR)/(Lau + BV)

patch_plot +
  plot_annotation(tag_levels = 'A')

ggsave(patch_plot, filename = "Map.png", width = 16, height = 15)

### map with raster
# ggplot() + 
#   geom_sf(data = sf_map) +
#   geom_point(data = BR, aes(x = Longitude, y = Latitude, color=factor(Site), 
#                             size  = AvgSalinity)) +
#   #scale_color_manual(values=c("#4F508C","#B56478","#CE9A28","#28827A", "#3F78C1")) + #choose colors of sites
#   #scale_size(limits = c(2, 29), breaks = c(2, 10, 20, 27, 28)) +
#   theme_gray() +
#   labs(x= "Longitude", y = "Latitude") + #Lat and Long labels
#   guides(color = guide_legend(title = "Site"), 
#          size = guide_legend(title = "Average Salinity")) + #change legend title
#   coord_cartesian(xlim=c(-123,-121.5),ylim=c(37,38.5)) + #set limits of map
#   theme(axis.text = element_text(size = 10), 
#         axis.title = element_text(size = 12))