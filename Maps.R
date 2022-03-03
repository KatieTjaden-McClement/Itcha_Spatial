# Map making
# Katie Tjaden-McClement
# Feb 72 2022

list.of.packages <- c("raster",           # If you are working with Rasters this is essential
                      "osmdata",          # Access the OpenStreet map interface -> it is constantly updated
                      "sf",               # The go to package for spatial operations
                      "bcmaps",           # The BC specific resources
                      "spex",             # ?
                      "dplyr",            # Tidyverse package for data manipulation
                      "ggplot2",          # Plotting tools
                      "rgeos",            # Spatial operation backend package
                      "rgdal",
                      "spatialEco",       # Handy tools for spatial ecology
                      "MODISTools",       # Tool for extracting NDVI
                      "leaflet",          # Tool for interactive maps!
                      "tidyverse",
                      "tmap",
                      "cowplot",
                      "ggspatial",
                      "ggmap")         

lapply(list.of.packages, require, character.only = TRUE)

##### Plotting camera locations #####
# Set up
# Get camera locations
sta <- as_tibble(read.csv("stations.csv", header=T))
head(sta)

# Plot them
plot(sta$longitude, sta$latitude, asp=T, las=1,
     xlab="Longitude", ylab="Latitude")


# The following code determines the best UTM for your data

mlong <- mean(sta$longitude); mlat <- mean(sta$latitude)

# UTM finder function
lonlat2UTM <-  function(lonlat) {
  utm <-  (floor((lonlat[1] + 180) / 6) %% 60) + 1
  if(lonlat[2] > 0) {
    utm + 32600
  } else{
    utm + 32700
  }
}

crs.utm <- lonlat2UTM(c(mlong,mlat))

## Convert your camera points into a shapefile using simple features (sf)
sta.wgs <- st_as_sf(sta, coords = c("longitude", "latitude"), crs=4326) # Note the CRS tells R that the projection is WGS1984

## Convert the lat long to UTM - st_transform
sta.utm <- st_transform(sta.wgs,crs=crs.utm) 

### Buffer these points by 18km to create an area of interest(AOI) 
tmp <- st_buffer(sta.utm, 18000)

# Get the bounding box coordinates
aoi.utm <- st_as_sfc(st_bbox(tmp))
aoi.wgs <- st_bbox(st_transform(aoi.utm, 4326))
#remove(tmp)

# Check your aoi is sensible
plot(st_geometry(aoi.utm))  # Plot the box
plot(st_geometry(sta.utm), add=T, col="red", pch=19) # Add your camera stations

#st_write(aoi.utm, paste0("Exported_data/", "aoi_utm.shp"), append=F)

# Nicer station map
#Map background of BC
mad_map_BC <- get_map(getbb("British Columbia"), maptype = "satellite")

#Map background of the location
myLocation_cams <- c(-125.66908, 52.03613, -123.89290, 53.00564)

#add background
mad_map_cams <- get_map(location=myLocation_cams, maptype = "roadmap")

g1 <- ggmap(mad_map_cams) +
  geom_sf(data=sta.wgs, inherit.aes = F) +
  ggtitle("Itcha-Ilgachuz Camera Trapping Project") +
  labs(x = "Longitude", y = "Latitude")
g1


##### Plotting I-I seasonal ranges #####
II_ranges <- st_read(dsn = "Spatial_Layers/Itcha_seasonal_ranges.shp")
plot(st_geometry(II_ranges), col = "red")

ranges.utm <- st_transform(II_ranges, crs = crs(aoi.utm)) #UTM zone 10
ranges.wgs <- st_transform(II_ranges, crs = crs(sta.wgs))

tmp <- st_buffer(ranges.utm, 30000)

aoi_ranges.utm <- st_as_sfc(st_bbox(tmp))
aoi_ranges.wgs <- st_bbox(st_transform(aoi_ranges.utm, 4326))
aoi_ranges.wgs

#Map background of the location
myLocation_ranges <- c(-126.27771, 51.32935, -122.80277, 53.63412)

#add background
mad_map_ranges <- get_map(location=myLocation_ranges, maptype = "roadmap")

g2 <- ggmap(mad_map_ranges) +
  geom_sf(data=ranges.wgs, inherit.aes = F,
          aes(fill = BCHab_code), alpha = 0.5) +
  ggtitle("Itcha-Ilgachuz caribou range") +
  labs(x = "Longitude", y = "Latitude", 
       fill = "Seasonal range") +
  annotation_scale(location = "bl", width_hint = 0.3) +
  theme(legend.position = c(1.2, 0.8))
g2

# add inset map of BC with red outline of aoi

# get BC outline
mapBC <- bc_bound()
#mapBC <- st_transform(mapBC, crs = crs(sta.wgs))
plot(st_geometry(mapBC))

inset <- ggplot(mapBC) +
  geom_sf(fill = "white", colour = "black") +
  geom_sf(data = aoi_ranges.utm,
          color = "red", fill = "white") +
  theme_void()

# plot together, need to play around with plot window size, x,y coords
g2_inset <- ggdraw() +
  draw_plot(g2, x = -0.15, y = 0) +
  draw_plot(inset, x = 0.63, y = 0.1, width = 0.4, height = 0.4)
g2_inset

ggsave("II_ranges.png", bg = "white")


