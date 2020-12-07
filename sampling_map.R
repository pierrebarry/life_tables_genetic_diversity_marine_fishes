#------------------------------------------------------------------------#
#                                                                        #
#                          SAMPLING MAP                                  #
#                                                                        #
#------------------------------------------------------------------------#

Salut !
# Load packages ----
library(oceanmap)
library(rgdal)
library(png)
library(openxlsx)
library(readxl)
library(RColorBrewer)
library(ggplot2)
library("rnaturalearth")
library("rnaturalearthdata")
library(sf)
library(ggsn)

# Draw map ----
color_med_atl=data.frame(Location=c("Gulf of Lion","Costa Calida","Algarve","Bay of Biscay"),
                         Col=brewer.pal(n = 4, name = "RdBu"))
load(file="Data/sampling.Rdata")
sampling$lon=as.numeric(sampling$lon)
sampling$lat=as.numeric(sampling$lat)
for (i in 1:nrow(sampling)){
  if (sampling$LOCATION[i]=="Mar Menor"){
    sampling$LOCATION[i]="Costa Calida"
  }
}
sampling$LOCATION=factor(sampling$LOCATION)
sites <- st_as_sf(data.frame(longitude = sampling[is.na(sampling$lon)==FALSE,]$lon, 
                             latitude = sampling[is.na(sampling$lat)==FALSE,]$lat, 
                             loc=sampling[is.na(sampling$lon)==FALSE,]$LOCATION,
                             detail_loc=sampling[is.na(sampling$lon)==FALSE,]$DETAILED_LOCATION,
                             sp=sampling[is.na(sampling$lon)==FALSE,]$SPECIES_CODE),
                  coords = c("longitude", "latitude"),
                  crs = 4326,
                  agr = "constant"
)

col=c()
for (i in 1:length(sites$loc)){
  col[i]=as.character(color_med_atl$Col[which(sites$loc[i]==color_med_atl$Location)])
}
sites$col=as.character(col)
world <- ne_countries(scale='medium',returnclass = 'sf')
lon <- c (-10, 10)
lat <- c (35, 47.5)

europe <- ggplot(data = world) +
  geom_sf(fill="grey",lwd=0) +
  coord_sf(xlim = c(-100,100), ylim = c(-25,75), expand = FALSE)+
    geom_rect(xmin = -10, xmax = 10, ymin = 35, ymax = 47.5, 
              fill = NA, colour = "black", size = 0.5) +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

all<-ggplot(data = world) +
  geom_sf(fill="grey") +
  geom_sf(data = sites, size = 2, shape = 1, col=col,aes(label=detail_loc)) +
  coord_sf(xlim = c(-10, 10), ylim = c(35,47.5), expand = FALSE)+
  annotate(geom = "text", x = 5, y = 38, label = "Mediterranean \n Sea", 
           fontface = "italic", color = "grey22", size = 5.5) +
  annotate(geom = "text", x = -7.5, y = 45, label = "Atlantic \n Ocean", 
           fontface = "italic", color = "grey22", size = 5.5) +
  annotate(geom = "text", x = 5, y = 42.35, label = "Gulf of \n Lion", 
            size = 4) +
  annotate(geom = "text", x = 0.75, y = 37.25, label = "Costa \n Calida", 
           size = 4) +
  annotate(geom = "text", x = -8, y = 36.5, label = "Algarve", 
           size = 4) +
  annotate(geom = "text", x = -2.85, y = 44.5, label = "Bay of \n Biscay", 
           size = 4) +
  theme(plot.tag.position = 'topleft')+
  xlab("")+
  ylab("")+
  scalebar(dist = 250, 
           dist_unit = "km", 
           model = 'WGS84',
           st.size = 2.5,
           transform=TRUE,
           x.min=8,
           x.max=9,
           y.min=36,
           y.max=45)+
  north(x.min=-9,x.max=-5,
        y.min=46,y.max=47)+
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA),
        legend.position = "none")

# Export map ----
pdf(paste(wd,"/figures/sampling_map.pdf",sep=""),width=5.5,height=5.5)
map<-ggplot() +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  annotation_custom(ggplotGrob(all), xmin = 0, xmax = 100, ymin = 0, 
                    ymax = 100) +
  annotation_custom(ggplotGrob(europe), xmin = 2.5, xmax = 37.5, ymin = 80, 
                    ymax = 100) + 
  theme_void()
print(map)
dev.off()

# Species for notebook -----
p<-vector('list',length(levels(sites$sp)))

for (j in 1:length(levels(sites$sp))){
  
  sites <- st_as_sf(data.frame(longitude = sampling[is.na(sampling$lon)==FALSE,]$lon, 
                               latitude = sampling[is.na(sampling$lat)==FALSE,]$lat, 
                               loc=sampling[is.na(sampling$lon)==FALSE,]$LOCATION,
                               detail_loc=sampling[is.na(sampling$lon)==FALSE,]$DETAILED_LOCATION,
                               sp=sampling[is.na(sampling$lon)==FALSE,]$SPECIES_CODE),
                    coords = c("longitude", "latitude"),
                    crs = 4326,
                    agr = "constant"
  )
  col=c()
  for (i in 1:length(sites$loc)){
    col[i]=as.character(color_med_atl$Col[which(sites$loc[i]==color_med_atl$Location)])
  }
  sites$col=as.character(col)
  
  shape=c()
  
  for (i in 1:length(sites$loc)){
    if (sites$sp[i]==levels(sites$sp)[j]){
      shape[i]=19
    } else {
      shape[i]=1
    }
  }
  sites$shape=as.numeric(shape)
  
  size=c()
  
  for (i in 1:length(sites$loc)){
    if (sites$sp[i]==levels(sites$sp)[j]){
      size[i]=4
    } else {
      size[i]=2
    }
  }
  sites$size=as.numeric(size)
  
  img<-readPNG(paste("figures/",levels(sites$sp)[j],".png",sep=""))
  g <- rasterGrob(img, interpolate=TRUE)
  
  p[[j]] <- local({
    all<-ggplot(data = world) +
      geom_sf(fill="grey") +
      geom_sf(data = sites, size = sites$size, shape = sites$shape, col=col,aes(label=detail_loc)) +
      coord_sf(xlim = c(-10, 10), ylim = c(35,47.5), expand = FALSE)+
      annotate(geom = "text", x = 5, y = 38, label = "Mediterranean \n Sea", 
               fontface = "italic", color = "grey22", size = 5.5) +
      annotate(geom = "text", x = -7.5, y = 45, label = "Atlantic \n Ocean", 
               fontface = "italic", color = "grey22", size = 5.5) +
      annotate(geom = "text", x = 5, y = 42.35, label = "Gulf of \n Lion", 
               size = 4) +
      annotate(geom = "text", x = 0.75, y = 37.25, label = "Costa \n Calida", 
               size = 4) +
      annotate(geom = "text", x = -8, y = 36.5, label = "Algarve", 
               size = 4) +
      annotate(geom = "text", x = -2.85, y = 44.5, label = "Bay of \n Biscay", 
               size = 4) +
      theme(plot.tag.position = 'topleft')+
      xlab("")+
      ylab("")+
      scalebar(dist = 250, 
               dist_unit = "km", 
               model = 'WGS84',
               st.size = 3,
               transform=TRUE,
               x.min=8,
               x.max=9,
               y.min=36,
               y.max=45)+
      north(x.min=-9,x.max=-5,
            y.min=46,y.max=47)+
      theme(panel.grid.major = element_blank(),
            panel.background = element_rect(fill = "transparent"),
            panel.border = element_rect(fill = NA),
            legend.position = "none")+
      annotation_custom(g, 
                        ymin=46, 
                        ymax=47.5,
                        xmin=-10, 
                        xmax=-5) 
    
  })
  
  
}

names(p)=levels(sites$sp)
species_map<-p

# Map each location ----
lon <- c (3.25, 4.15)
lat <- c (43.15, 43.65)
gulf_of_lion<-ggplot(data = world) +
  geom_sf(fill="grey") +
  geom_sf(data = sites, size = 1, shape = 19, aes(col = loc)) +
  coord_sf(xlim = lon, ylim = lat, expand = FALSE)+
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

lon <- c (-1.5, -0.15)
lat <- c (37.5,38.25)
mar_menor<-ggplot(data = world) +
  geom_sf(fill="grey") +
  geom_sf(data = sites, size = 1, shape = 19, aes(col = loc)) +
  coord_sf(xlim = lon, ylim = lat, expand = FALSE)+
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

lon <- c (-9, -7.75)
lat <- c (36.75,37.25)
algarve<-ggplot(data = world) +
  geom_sf(fill="grey") +
  geom_sf(data = sites, size = 1, shape = 19, aes(col = loc)) +
  coord_sf(xlim = lon, ylim = lat, expand = FALSE)+
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )

lon <- c (-2.25, -1)
lat <- c (43.25,44.75)
bay_of_biscay<-ggplot(data = world) +
  geom_sf(fill="grey") +
  geom_sf(data = sites, size = 2, shape = 19, aes(col = loc)) +
  coord_sf(xlim = lon, ylim = lat, expand = FALSE)+
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
  )



pdf(paste(wd,"/figures/sampling_map_loc.pdf",sep=""),width=5,height=5)
ggplot() +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  annotation_custom(ggplotGrob(all), xmin = 25, xmax = 75, ymin = 25, 
                    ymax = 75) +
  annotation_custom(ggplotGrob(gulf_of_lion), xmin = 75, xmax = 100, ymin = 75, 
                    ymax = 100) +
  annotation_custom(ggplotGrob(mar_menor), xmin = 75, xmax = 100, ymin = 0, 
                    ymax = 25) +
  annotation_custom(ggplotGrob(algarve), xmin = 0, xmax = 25, ymin = 0, 
                    ymax = 25) +
  annotation_custom(ggplotGrob(bay_of_biscay), xmin = 0, xmax = 25, ymin = 75, 
                    ymax = 100) +
  theme_void()
dev.off()

#------
