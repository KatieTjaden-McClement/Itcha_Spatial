### Extracting spatial covariates for Itcha
### Katie Tjaden-McClement
### 22 Feb 2022

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
                      "ggmap")         

lapply(list.of.packages, require, character.only = TRUE)

##### Cam locations and set AOI #####

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
myLocation_BH <- c(-125.66908, 52.03613, -123.89290, 53.00564) #Bighorn location BB
#add background
mad_map_BH <- get_map(location=myLocation_BH, maptype = "roadmap")

g1 <- ggmap(mad_map_BH) +
  geom_sf(data=sta.wgs, inherit.aes = F) +
  ggtitle("Itcha-Ilgachuz Camera Trapping Project")
g1


##### Extract elevation data #####

# Download Digital Elevation Model for AOI
DEM_raster <- cded_raster(aoi = aoi.utm)
plot(DEM_raster)

# Check the stations are there!
plot(st_geometry(sta.wgs), add=T)

# Extract value from raster at each point and add it to your dataframe
sta$Elevation <- raster::extract(DEM_raster, sta.wgs)


##### Seasonal range polygons #####
ogrListLayers("Spatial_Layers/II_Seasonal_Ranges.gdb")

II_ranges <- st_read(dsn = "Spatial_Layers/II_Seasonal_Ranges.gdb", 
                     layer = "BCHabitat_Herd_Boundaries_20190730")
summary(II_ranges)

II_ranges <- filter(II_ranges, Herd_Name == "Itcha-Ilgachuz")

plot(st_geometry(II_ranges))

tm_shape(II_ranges) +
  tm_polygons()
# not working now... try to export layer from Arc as shp


##### Cutblock and Fire polygons and road layer #####

# Cutblock layer from gov data warehouse
cut <- st_read(dsn = "Spatial_Layers/CUT_BLOCKS/CNS_CUT_BL_polygon.shp")
cut <- st_transform(cut, crs = crs(aoi.utm)) #UTM zone 10
cut <- st_intersection(cut, aoi.utm) # clip to AOI

# Fire layer
fire <- st_read(dsn = "Spatial_Layers/FIRE_POLYS/H_FIRE_PLY_polygon.shp")
fire <- st_transform(fire, crs = crs(aoi.utm)) #UTM zone 10
fire <- st_intersection(fire, aoi.utm) # clip to AOI

# Roads
roads <- st_read(dsn = "Spatial_Layers/Cariboo_Consolidated_Roads_with_data_corrections_20211124.gdb/",
                 layer = "Cariboo_Consolidated_Roads")

roads <- st_transform(roads, crs = crs(aoi.utm))
roads <- st_intersection(roads, aoi.utm)

# Plot all layers with cams
par(mfrow=c(1,1))
plot(st_geometry(aoi.utm))  # Plot the box
plot(st_geometry(roads), add = T, col = "grey") #add roads
plot(st_geometry(fire), add = T, col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5),
     border = "yellow") #add fire
plot(st_geometry(cut), add = T, col = rgb(red = 0, green = 1, blue = 0, alpha = 0.5), 
     border = "green") #add cutblocks
plot(st_geometry(sta.utm), add=T, col="red", pch=19) #add cameras
#everything looks good!!

fire_40 <- filter(fire, FIRE_YEAR > (2022-40)) # down to 125 polys from 275
cut_40 <- filter(cut, HARVESTYR > (2022-40)) # down to 6016 from 6072


##### Buffer cams and extract covariates ####

# Create 500 m buffers around camera
cams_buff <- st_buffer(sta.utm, dist = 500)
plot(st_geometry(cams_buff))

# clip cutblock polygons to camera buffers
buff_cut <- st_intersection(cams_buff, cut_40)
plot(st_geometry(buff_cut))
head(buff_cut)

length(unique(buff_cut$station_id)) #32 camera sites
nrow(buff_cut) # 46 fire polygons overlapping with camera buffers

# Need to remove overlapping polygons:
buff_cut <- buff_cut %>% 
  group_by(station_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
plot(st_geometry(buff_cut)) #good

buff_cut <- as_Spatial(buff_cut) %>% 
  st_as_sf()

buff_cut <- st_intersection(cams_buff, buff_cut)
plot(st_geometry(buff_cut))

# Calculate area cut within each cam buffer
buff_cut$area_cut <- st_area(buff_cut$geometry)

cut_summary <- buff_cut %>% 
  group_by(station_id) %>% 
  summarise(area_cut = sum(area_cut)) %>% 
  mutate(prop_cut = as.numeric(area_cut/(pi*500^2)))
hist(cut_summary$prop_cut)
# can merge prop_cut column with station csv later, adding zeros for missing stations (no cutblocks)

### Repeat with burnt areas ###
# clip fire polygons to camera buffers
buff_fire <- st_intersection(cams_buff, fire_40)
plot(st_geometry(buff_fire))
head(buff_fire)

length(unique(buff_fire$station_id)) #32 camera sites
nrow(buff_fire) # 46 fire polygons overlapping with camera buffers

# Need to remove overlapping polygons:
buff_fire <- buff_fire %>% 
  group_by(station_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
plot(st_geometry(buff_fire)) #good

buff_fire <- as_Spatial(buff_fire) %>% 
  st_as_sf()

buff_fire <- st_intersection(cams_buff, buff_fire)
plot(st_geometry(buff_fire))

# Calculate area cut within each cam buffer
buff_fire$area_fire <- st_area(buff_fire$geometry)

fire_summary <- buff_fire %>% 
  group_by(station_id) %>% 
  summarise(area_fire = sum(area_fire)) %>% 
  mutate(prop_fire = as.numeric(area_fire/(pi*500^2)))
hist(fire_summary$prop_fire)

