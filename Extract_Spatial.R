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
                      "ggmap",
                      "lubridate")         

lapply(list.of.packages, require, character.only = TRUE)

theme_set(theme_classic())

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

### Puntzi lat, long, and utm coords for Carolyn:
# sta.utm_puntzi <- filter(sta.utm, grepl("ITCHA4", station_id)) %>% 
#   arrange(station_id)
# puntzi_utm <- as_tibble(st_coordinates(sta.utm_puntzi)) %>% 
#   rename(UTM.E = X, UTM.N = Y)
# 
# sta.wgs_puntzi <- filter(sta.wgs, grepl("ITCHA4", station_id)) %>% 
#   arrange(station_id)
# puntzi_lat_long <- as_tibble(st_coordinates(sta.wgs_puntzi)) %>% 
#   rename(Longitude = X, Latitude = Y)
# 
# puntzi_utm_latlong <- tibble(station_id = sta.utm_puntzi$station_id, 
#                              puntzi_lat_long,
#                              puntzi_utm)
# write.csv(puntzi_utm_latlong, file = "Outputs/puntzi_utm_latlon.csv")

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

### Distances between cameras ###
#install.packages("spatstat")
# library(spatstat)
# 
# sta.utm_sp <- filter(sta.utm, !grepl("ITCHA4", station_id)) %>% 
#   as_Spatial()
# aoi_sp <- as_Spatial(aoi)
# plot(aoi)
# 
# bbox <- st_bbox(aoi)
# 
# sta_owin <- as.owin(bbox)
# 
# pts <- coordinates(sta.utm_sp)
# p <- ppp(pts[,1], pts[,2], window = sta_owin)
# 
# test <- nndist(p, k = 1:5)/1000
# test <- as.data.frame(test)
# 
# test$station_id <- sta.utm_sp$station_id
# write_csv(test, file = "station_distances.csv")
# went through and manually removed extra distance measures for side/corner stations or 
# those that have been re-set in excel and did summary stats

##### Extract elevation data #####

# Download Digital Elevation Model for AOI
DEM_raster <- cded_raster(aoi = aoi.utm)
plot(DEM_raster)

# Check the stations are there!
plot(st_geometry(sta.wgs), add=T)

# Extract value from raster at each point and add it to your dataframe
sta$Elevation <- raster::extract(DEM_raster, sta.wgs)

##### Extract ruggedness data #####

vrm3 <- vrm(DEM_raster, s = 3) # double check exactly what s = 3 means
plot(vrm3)

sta$vrm <- raster::extract(vrm3, sta.wgs)

##### Extract NDVI data #####
bands <- mt_bands(product = "MOD13Q1") # MODIS/Terra Vegetation Indices (NDVI/EVI) 16-Day L3 Global 250m SIN Grid
head(bands)

other_bands <- mt_bands(product = "MCD12Q1")
head(other_bands)

tmp <- sta %>% 
  dplyr::select("station_id", "longitude", "latitude")
colnames(tmp) <- c("site_name", "lon", "lat")
tmp <- tmp[,c(1,3,2)]
str(tmp)

# list available dates for a product at a center of aoi
dates <- mt_dates(product = "MOD13Q1", lat = mlat, lon = mlong)
head(dates)
tail(dates)

# get NDVI
itcha_ndvi <- mt_batch_subset(product = "MOD13Q1",
                              df=tmp,
                              band = "250m_16_days_NDVI",
                              start = "2020-09-01",
                              end = "2021-10-31",
                              km_lr = 0,
                              km_ab = 0,
                              internal = TRUE)
hist(itcha_ndvi$value)

# Get Quality Assurance data for NDVI values:
itcha_vi_qa <- mt_batch_subset(product = "MOD13Q1",
                               df = tmp,
                               band = "250m_16_days_VI_Quality",
                               start = "2020-09-01",
                               end = "2021-10-31",
                               km_lr = 0,
                               km_ab = 0,
                               internal = TRUE)
hist(itcha_vi_qa$value)
# need to convert from decimal to binary to correspond to QA metrics outlined in user guide

itcha_vi_qa <- rename(itcha_vi_qa, "qa_value" = "value")

# add to itcha_ndvi data
itcha_ndvi$qa_value <- itcha_vi_qa$qa_value

# Remove stuff that is below zero, which are indicative of water
itcha_ndvi_filtered_old <- filter(itcha_ndvi, value > 0)
# removes 4784-4297 = 487 values

# Valid range of values is -2000 to 10,000 and fill value is -3000, so just remove those
itcha_ndvi <- filter(itcha_ndvi, value >= -2000) # removes 26 values

itcha_ndvi$calendar_date <- strptime(itcha_ndvi$calendar_date, "%Y-%m-%d")

# Plot NDVI through time
itcha_ndvi <- itcha_ndvi %>% 
  mutate(grid = word(site, 1, sep = "-"),
         calendar_date = as.POSIXct(calendar_date))

ggplot(itcha_ndvi, aes(x = calendar_date, y = value, colour = grid)) +
  geom_point(alpha = 0.1) + 
  geom_line(stat = "smooth", method = "loess", 
            alpha = 0.8, size = 1)

# summarize to monthly station NDVI values
sta_ndvi <- itcha_ndvi %>%
  mutate(month = format_ISO8601(calendar_date, precision = "ym")) %>% 
  group_by(grid, site, month) %>% 
  summarise(ndvi = mean(value),
            ndvi_scaled = ndvi*0.0001)
# no missing values after new filtering criteria (keeps negative values - likely snow and ice at high elevation)

ggplot(sta_ndvi, aes(x = ym(month), 
                     y = ndvi_scaled, colour = grid)) +
  geom_point(alpha = 0.1) + 
  geom_line(stat = "smooth", method = "loess", 
            alpha = 0.8, size = 1) +
  labs(x = "Date", y = "NDVI", colour = "Camera Grid")

#write.csv(sta_ndvi, "station_ndvi_Sept2020-Oct2021.csv")

##### Cutblock and Fire polygons and road layer #####

# Cutblock layer from gov data warehouse
cut <- st_read(dsn = "Spatial_Layers/CUT_BLOCKS_ranges.gdb")
cut <- st_transform(cut, crs = crs(aoi.utm)) #UTM zone 10

# Fire layer
fire <- st_read(dsn = "Spatial_Layers/FIRE_POLYS_ranges.gdb")
fire <- st_transform(fire, crs = crs(aoi.utm)) #UTM zone 10

# Roads
roads <- st_read(dsn = "Spatial_Layers/Cariboo_Consolidated_Roads_with_data_corrections_20211124.gdb/",
                 layer = "Cariboo_Consolidated_Roads")

roads <- st_transform(roads, crs = crs(aoi.utm))

# need to remove roads that don't actually exist:
roads <- filter(roads, TRANSPORT_LINE_TYPE_CODE != "X")


##### Distance to water #####

# downloaded data from https://maps.canada.ca/czs/index-en.html
lakes <- st_read(dsn = "Spatial_Layers/Watercourses/waterbody_2.shp") %>% 
  select("feature_id", "geometry")
lakes <- st_transform(lakes, crs = crs(aoi.utm))
plot(st_geometry((lakes)))

rivers <- st_read(dsn = "Spatial_Layers/Watercourses/watercourse_1.shp") %>% 
  select("feature_id", "geometry")
rivers <- st_transform(rivers, crs = crs(aoi.utm))
plot(st_geometry(rivers), add = T)

plot(st_geometry(sta.utm), add = T)

# bind together
lakes_polyline <- st_cast(lakes, to = "LINESTRING")
plot(st_geometry(lakes_polyline))

water <- rbind(lakes_polyline, rivers)
plot(st_geometry(water))

# find nearest water feature to each cam
closest_water <- st_nearest_feature(sta.utm, water)

closest_water <- water[closest_water, ]
plot(st_geometry(sta.utm))
plot(st_geometry(closest_water), add = T) #looks good

# calculate distance to closest water feature
dist_water <- st_distance(sta.utm, closest_water, by_element = T)
dist_water_km <- as.numeric(dist_water/1000)
hist(dist_water_km)

dist_water <- tibble(station_id = sta.utm$station_id, dist_water_km)

# how accurate are these water sources? could be lots of tiny streams and stuff...

## compare to waterbody and watercourse layers from bc data catalougue (freshwater atlas):
# https://catalogue.data.gov.bc.ca/dataset/freshwater-atlas-stream-network
# https://catalogue.data.gov.bc.ca/dataset/freshwater-atlas-lakes

fwa_lakes <- st_read("Spatial_Layers/Watercourses/GeoBC/FWA_LAKES_POLY.gdb") %>% 
  st_transform(crs = crs(aoi.utm)) %>% 
  select(WATERBODY_POLY_ID) %>% 
  rename(ID = WATERBODY_POLY_ID) %>% 
  st_cast("MULTILINESTRING")
str(fwa_lakes)

fwa_streams <- st_read("Spatial_Layers/Watercourses/GeoBC/FWA_STREAM_NETWORKS_SP.gdb") %>% 
  st_transform(crs = crs(aoi.utm)) %>% 
  select(LINEAR_FEATURE_ID) %>% 
  rename(ID = LINEAR_FEATURE_ID) %>% 
  st_zm() # remove unneeded z dimension
str(fwa_streams)

fwa_water <- rbind(fwa_lakes, fwa_streams)
plot(st_geometry((fwa_water))) # even denser coverage....

# could remove streams below a certian magnitude
# e.g. order 1 streams (the majority of streams in the dataset) only run during wet periods

# find nearest water feature to each cam
closest_fwa_water <- st_nearest_feature(sta.utm, fwa_water)

closest_fwa_water <- fwa_water[closest_fwa_water, ]
plot(st_geometry(sta.utm))
plot(st_geometry(closest_fwa_water), add = T) #looks good

# calculate distance to closest water feature
dist_fwa_water <- st_distance(sta.utm, closest_fwa_water, by_element = T)
dist_fwa_water_km <- as.numeric(dist_fwa_water/1000)
hist(dist_fwa_water_km)

### Wetlands ###
# lots of issues with data quality (see emails with Robin Steenweg)
# FWA seems to be best resource, let's see how it looks
fwa_wetlands <- st_read("Spatial_Layers/Watercourses/GeoBC/FWA_WETLANDS_POLY") %>% 
  st_transform(crs = crs(sta.utm))
head(fwa_wetlands)

plot(st_geometry(fwa_wetlands))
plot(st_geometry(sta.utm), col = "red", add = T)
# looks not bad (but you need to expand map to see clearly otherwise it looks a bit crazy)

### BEC Zones ###
bec <- st_read("Spatial_Layers/BEC/BEC_BIOGEOCLIMATIC_POLY.gdb") %>% 
  st_transform(crs(aoi.utm))

plot(st_geometry(bec))

bec_plot <- st_intersection(bec, aoi.utm) %>% 
  st_transform(crs(sta.wgs))
# plot(bec_plot["ZONE_NAME"])
# plot(st_geometry(sta.utm), add = T)

ggmap(mad_map_cams) +
  geom_sf(data = bec_plot, aes(fill = ZONE_NAME), 
          inherit.aes = F) +
  geom_sf(data = sta.wgs, inherit.aes = F)
# most grids are almost entirely within one or 2 bec zones...

# extract zone and subzones at cams
bec_cams <- st_intersection(sta.utm, bec) %>% 
  select(station_id, ZONE_NAME, SUBZONE_NAME) %>%
  as.data.frame() %>% 
  select(-geometry)
# all subzones either very dry very cold, very dry cold, or undifferentiated

# Plot all layers with cams
# par(mfrow=c(1,1))
# plot(st_geometry(aoi.utm))  # Plot the box
# plot(st_geometry(roads), add = T, col = "grey") #add roads
# plot(st_geometry(fire), add = T, col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5),
#      border = "yellow") #add fire
# plot(st_geometry(cut), add = T, col = rgb(red = 0, green = 1, blue = 0, alpha = 0.5), 
#      border = "green") #add cutblocks
# plot(st_geometry(sta.utm), add=T, col="red", pch=19) #add cameras
# #everything looks good!!

fire_40 <- filter(fire, FIRE_YEAR > (2022-40)) 
cut_40 <- filter(cut, HARVEST_YEAR > (2022-40)) 

cut_20 <- filter(cut, HARVEST_YEAR > (2022-20))
cut_20_40 <- filter(cut, between(HARVEST_YEAR, (2022-40), (2022-20)))

##### Buffer cams and extract covariates ####

# Create 500 m buffers around camera
cams_buff <- st_buffer(sta.utm, dist = 500)
plot(st_geometry(cams_buff))

### Cutblocks < 40 years old
# clip cutblock polygons to camera buffers
buff_cut_40 <- st_intersection(cams_buff, cut_40)
plot(st_geometry(buff_cut_40))
head(buff_cut_40)

length(unique(buff_cut_40$station_id)) #71 camera sites
nrow(buff_cut_40) # 157 cutblock polygons overlapping with camera buffers

# Need to remove overlapping polygons:
buff_cut_40 <- buff_cut_40 %>% 
  group_by(station_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
#plot(st_geometry(buff_cut)) #good

buff_cut_40 <- as_Spatial(buff_cut_40) %>% 
  st_as_sf()

buff_cut_40 <- st_intersection(cams_buff, buff_cut_40)
plot(st_geometry(buff_cut_40))

# Calculate area cut within each cam buffer
buff_cut_40$area_cut <- st_area(buff_cut_40$geometry)

cut_40_summary <- buff_cut_40 %>% 
  group_by(station_id) %>% 
  summarise(area_cut = sum(area_cut)) %>% 
  mutate(prop_cut_40 = as.numeric(area_cut/(pi*500^2)))
hist(cut_40_summary$prop_cut_40)
# can merge prop_cut column with station csv later, adding zeros for missing stations (no cutblocks)

### Cutblocks < 20 years old
# clip cutblock polygons to camera buffers
buff_cut_20 <- st_intersection(cams_buff, cut_20)
plot(st_geometry(buff_cut_20))
head(buff_cut_20)

length(unique(buff_cut_20$station_id)) #25 camera sites
nrow(buff_cut_20) # 64 cutblock polygons overlapping with camera buffers

# Need to remove overlapping polygons:
buff_cut_20 <- buff_cut_20 %>% 
  group_by(station_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
#plot(st_geometry(buff_cut)) #good

buff_cut_20 <- as_Spatial(buff_cut_20) %>% 
  st_as_sf()

buff_cut_20 <- st_intersection(cams_buff, buff_cut_20)
plot(st_geometry(buff_cut_20))

# Calculate area cut within each cam buffer
buff_cut_20$area_cut <- st_area(buff_cut_20$geometry)

cut_20_summary <- buff_cut_20 %>% 
  group_by(station_id) %>% 
  summarise(area_cut = sum(area_cut)) %>% 
  mutate(prop_cut_20 = as.numeric(area_cut/(pi*500^2)))
hist(cut_20_summary$prop_cut_20, breaks = 20)

### Cutblocks 20-40 years old
# clip cutblock polygons to camera buffers
buff_cut_20_40 <- st_intersection(cams_buff, cut_20_40)
plot(st_geometry(buff_cut_20_40))
head(buff_cut_20_40)

length(unique(buff_cut_20_40$station_id)) #59 camera sites
nrow(buff_cut_20_40) # 93 cutblock polygons overlapping with camera buffers

# Need to remove overlapping polygons:
buff_cut_20_40 <- buff_cut_20_40 %>% 
  group_by(station_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
#plot(st_geometry(buff_cut)) #good

buff_cut_20_40 <- as_Spatial(buff_cut_20_40) %>% 
  st_as_sf()

buff_cut_20_40 <- st_intersection(cams_buff, buff_cut_20_40)
plot(st_geometry(buff_cut_20_40))

# Calculate area cut within each cam buffer
buff_cut_20_40$area_cut <- st_area(buff_cut_20_40$geometry)

cut_20_40_summary <- buff_cut_20_40 %>% 
  group_by(station_id) %>% 
  summarise(area_cut = sum(area_cut)) %>% 
  mutate(prop_cut_20_40 = as.numeric(area_cut/(pi*500^2)))
hist(cut_20_40_summary$prop_cut_20_40, breaks = 20)

### Repeat with burnt areas ###
# clip fire polygons to camera buffers
buff_fire <- st_intersection(cams_buff, fire_40)
#plot(st_geometry(buff_fire))
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
#plot(st_geometry(buff_fire)) #good

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

### Roads
buff_roads <- st_intersection(cams_buff, roads)
#plot(st_geometry(buff_roads))

# Calculate road density within each cam buffer
buff_roads$road_length <- st_length(buff_roads$geometry)

road_summary <- buff_roads %>% 
  group_by(station_id) %>% 
  summarise(road_length = sum(road_length)) %>% 
  mutate(road_dens_m = as.numeric(road_length*1e6/(pi*500^2)),
         road_dens_km = road_dens_m/1000)
hist(road_summary$road_dens_km)

### distance from cam to nearest cutblock edge
# find nearest cutblock
closest_cut <- st_nearest_feature(sta.utm, cut_40)

closest_cut <- cut_40[closest_cut, ]
plot(st_geometry(sta.utm))
plot(st_geometry(closest_cut), add = T) #looks good

# calculate distance to closest cut
dist_cut <- st_distance(sta.utm, closest_cut, by_element = T)
dist_cut_km <- as.numeric(dist_cut/1000)


### Make tibble with extracted covariates
station_sp_cov <- sta %>% 
  select(station_id, Elevation, vrm) %>% 
  left_join(bec_cams, by = "station_id") %>% 
  left_join(fire_summary, by = "station_id") %>% 
  left_join(cut_40_summary, by = "station_id") %>% 
  left_join(cut_20_summary, by = "station_id") %>%
  left_join(cut_20_40_summary, by = "station_id") %>%
  left_join(road_summary, by = "station_id") %>% 
  select(station_id, prop_cut_40, prop_cut_20, prop_cut_20_40, prop_fire, 
         road_dens_km, Elevation, vrm, ZONE_NAME, SUBZONE_NAME) %>% 
  mutate(dist_cut_km = dist_cut_km,
         dist_water_km = dist_water_km)

station_sp_cov[is.na(station_sp_cov)] <- 0

# convert to long format to plot this faceted over covariates:
station_sp_cov_long <- station_sp_cov %>% 
  select(-ZONE_NAME, -SUBZONE_NAME, -prop_cut_20_40,
         -prop_cut_20, -vrm) %>% 
  pivot_longer(cols = (-station_id),
               names_to = "covariate",
               values_to = "value") %>% 
  mutate(grid_id = case_when(grepl("ITCHA1", station_id) ~ "Ilgachuz",
                             grepl("ITCHA2", station_id) ~ "Saddle",
                             grepl("ITCHA3", station_id) ~ "Itcha",
                             grepl("ITCHA4", station_id) ~ "Puntzi",
                             grepl("ITCHA5", station_id) ~ "Chezacut",
                             grepl("ITCHA6", station_id) ~ "Satah"))

ggplot(station_sp_cov_long, aes(x = value)) + 
  geom_histogram(fill = "darkblue") +
  facet_wrap(~covariate, scales = "free_x")
ggsave("Outputs/covariate_hists.png", width = 8, height = 6)


ggplot(station_sp_cov_long, aes(x = value, fill = grid_id)) + 
  geom_density(alpha = 0.4) +
  facet_wrap(~covariate, scales = "free")

write.csv(station_sp_cov, file = "Outputs/station_spatial_covariates.csv")


##### Seasonal range disturbance summary #####
II_ranges <- st_read(dsn = "Spatial_Layers/Itcha_seasonal_ranges.shp")
II_ranges <- st_transform(II_ranges, crs = crs(aoi.utm)) #UTM zone 10

plot(st_geometry(II_ranges), col = "red")

# get total area of each range
range_area_summary <- II_ranges %>% 
  mutate(area = st_area(geometry)) %>% 
  group_by(BCHab_code) %>% 
  summarise(total_area_m2 = sum(area),
            total_area_km2 = as.numeric(total_area_m2)/1e6)

# Get BB for seasonal ranges with 1 km buffer
tmp <- st_buffer(II_ranges, 1000)
bbox_ranges <- st_as_sfc(st_bbox(tmp))

#st_write(bbox_ranges, dsn = "bbox_ranges.shp")

# need to download cutblock and fire layers with wider aoi to capture full range (raw roads layer is good)

### Cutblocks
# Wider cutblock layer
cut_ranges <- st_read("Spatial_Layers/CUT_BLOCKS_ranges.gdb")
cut_ranges <- st_transform(cut_ranges, crs = crs(aoi.utm))
#plot(st_geometry(cut_ranges), add = T) #looks good

min(cut_ranges$HARVEST_YEAR) #1960
hist(cut_ranges$HARVEST_YEAR)

# clip cutblock polygons seasonal ranges
ranges_cut <- st_intersection(II_ranges, cut_ranges)
#plot(st_geometry(ranges_cut), col = "green")
head(ranges_cut)

# Need to remove overlapping polygons:
ranges_cut <- ranges_cut %>% 
  group_by(BCHab_code) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
#plot(st_geometry(ranges_cut))

ranges_cut <- as_Spatial(ranges_cut) %>% 
  st_as_sf()

ranges_cut <- st_intersection(II_ranges, ranges_cut)
#plot(st_geometry(ranges_cut))

# Calculate area and proportion cut in each range
ranges_cut$area_cut <- st_area(ranges_cut$geometry)
head(ranges_cut)

cut_summary_ranges <- ranges_cut %>% 
  group_by(BCHab_code) %>% 
  summarise(area_cut_m2 = sum(area_cut),
            area_cut_km2 = as.numeric(area_cut_m2)/1e6) 
cut_summary_ranges

### Fire
# Wider fire layer
fire_ranges <- st_read("Spatial_Layers/FIRE_POLYS_ranges.gdb")
fire_ranges <- st_transform(fire_ranges, crs = crs(aoi.utm))
#plot(st_geometry(fire_ranges), col = "red")

min(fire_ranges$FIRE_YEAR) #1921
hist(fire_ranges$FIRE_YEAR)

# clip fire polygons seasonal ranges
ranges_fire <- st_intersection(II_ranges, fire_ranges)
#plot(st_geometry(ranges_fire), col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5))
head(ranges_fire)

# Need to remove overlapping polygons:
ranges_fire <- ranges_fire %>% 
  group_by(BCHab_code) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
#plot(st_geometry(ranges_fire), col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5))

ranges_fire <- as_Spatial(ranges_fire) %>% 
  st_as_sf()

ranges_fire <- st_intersection(II_ranges, ranges_fire)
#plot(st_geometry(ranges_fire), col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5))

# Calculate area and proportion cut in each range
ranges_fire$area_fire <- st_area(ranges_fire$geometry)
head(ranges_cut)

fire_summary_ranges <- ranges_fire %>% 
  group_by(BCHab_code) %>% 
  summarise(area_fire_m2 = sum(area_fire),
            area_fire_km2 = as.numeric(area_fire_m2)/1e6) 
fire_summary_ranges

### Roads
roads <- st_read(dsn = "Spatial_Layers/Cariboo_Consolidated_Roads_with_data_corrections_20211124.gdb/",
                 layer = "Cariboo_Consolidated_Roads")
roads <- st_transform(roads, crs = crs(aoi.utm))
# need to remove roads that don't actually exist:
roads <- filter(roads, TRANSPORT_LINE_TYPE_CODE != "X")

range_roads <- st_intersection(II_ranges, roads)
plot(st_geometry(range_roads))

# Calculate road density within each cam buffer
range_roads$road_length <- st_length(range_roads$geometry)

road_summary_ranges <- range_roads %>% 
  group_by(BCHab_code) %>% 
  summarise(road_length_m = sum(road_length),
            road_length_km = road_length_m/1000)


### Total disturbance

# buffer roads by 50m
roads_buffered <- st_buffer(roads, 50)
roads_buffered_range <- st_buffer(range_roads, 50)
#plot(st_geometry(roads_buffered), col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5))

# merge all disturbance layers
roads_tomerge <- select(roads_buffered_range, BCHab_code)
cut_tomerge <- select(ranges_cut, BCHab_code)
fire_tomerge <- select(ranges_fire, BCHab_code)

all_dist_ranges <- rbind(roads_tomerge, cut_tomerge, fire_tomerge)

# remove overlap
all_dist_ranges <- all_dist_ranges %>% 
  group_by(BCHab_code) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
#plot(st_geometry(all_dist_ranges), col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5))

all_dist_ranges <- as_Spatial(all_dist_ranges) %>% 
  st_as_sf()

all_dist_ranges <- st_intersection(II_ranges, all_dist_ranges)

total_dist_summary_ranges <- all_dist_ranges %>% 
  mutate(all_dist_area = st_area(geometry)) %>% 
  group_by(BCHab_code) %>% 
  summarise(all_dist_area_m2 = sum(all_dist_area),
            all_dist_area_km2 = as.numeric(all_dist_area_m2/1e6))

# put it all together
disturbance_summary_ranges <- cbind(range_area_summary, cut_summary_ranges,
                                    fire_summary_ranges, road_summary_ranges,
                                    total_dist_summary_ranges) %>% 
  as_tibble() %>% 
  select(c("BCHab_code", "total_area_km2", "area_cut_km2", 
           "area_fire_km2", "road_length_km", "all_dist_area_km2")) %>% 
  mutate(prop_cut = area_cut_km2/total_area_km2,
         prop_fire = area_fire_km2/total_area_km2,
         road_dens = road_length_km/total_area_km2,
         prop_total_dist = all_dist_area_km2/total_area_km2)
disturbance_summary_ranges

write.csv(disturbance_summary_ranges, "disturbance_summary_ranges.csv")

##### Disturbance within grids #####
# quantify levels fo disturbance within each camera grid

# buffer camera grids by 1.3km to get 200km^2 grid areas
grids <- sta.utm %>% 
  mutate(grid_id = word(station_id, 1, sep = "-")) %>% 
  group_by(grid_id) %>% 
  summarise() %>% 
  st_cast("POLYGON") %>% 
  st_convex_hull() %>% 
  st_buffer(1300)
plot(st_geometry(grids))
head(grids)

grids_area <- grids %>% 
  mutate(area = st_area(grids$geometry)/1e6) %>% 
  as_tibble() %>% 
  select(-geometry)
grids_area # all very close to 200, good

# Road density for each grid
roads_grids <- st_intersection(grids, roads)

roads_grids$road_length <- st_length(roads_grids$geometry)

road_summary_grid <- roads_grids %>% 
  group_by(grid_id) %>% 
  summarise(road_length_m = sum(road_length),
            road_length_km = as.numeric(road_length_m/1000))

# Fire in each grid
fire_grids <- st_intersection(grids, fire_ranges) %>% 
  select(grid_id, FIRE_LABEL) %>% 
  rename(ID = FIRE_LABEL)

fire_grids <- fire_grids %>% 
  group_by(grid_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
plot(st_geometry(fire_grids), col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5))

fire_grids <- as_Spatial(fire_grids) %>% 
  st_as_sf()

fire_grids <- st_intersection(grids, fire_grids)

fire_summary_grid <- fire_grids %>% 
  mutate(area_fire = st_area(geometry)) %>% 
  group_by(grid_id) %>% 
  summarise(area_fire_m2 = sum(area_fire),
            area_fire_km2 = as.numeric(area_fire_m2)/1e6) 
fire_summary_grid

# Cutblocks in each grid
cut_grids <- st_intersection(grids, cut_ranges) %>% 
  select(grid_id, VEG_CONSOLIDATED_CUT_BLOCK_ID) %>% 
  rename(ID = VEG_CONSOLIDATED_CUT_BLOCK_ID)

cut_grids <- cut_grids %>% 
  group_by(grid_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
plot(st_geometry(cut_grids), col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5))

cut_grids <- as_Spatial(cut_grids) %>% 
  st_as_sf()

cut_grids <- st_intersection(grids, cut_grids)

cut_summary_grid <- cut_grids %>% 
  mutate(area_cut = st_area(geometry)) %>% 
  group_by(grid_id) %>% 
  summarise(area_cut_m2 = sum(area_cut),
            area_cut_km2 = as.numeric(area_cut_m2)/1e6) 
cut_summary_grid

# clip road buff layer to grid areas:
roads_buff_grids <- st_intersection(grids, roads_buffered) %>% 
  select(grid_id, CCR_ID) %>% 
  rename(ID = CCR_ID)

# merge layers
head(roads_buff_grids)
head(fire_grids)
head(cut_grids)

roads_buff_grids <- select(roads_buff_grids, grid_id)

all_dist_grids <- rbind(roads_buff_grids, fire_grids, cut_grids)
head(all_dist_grids)

# remove overlap:
all_dist_grids <- all_dist_grids %>% 
  group_by(grid_id) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>%
  ungroup() %>%
  # union polygons
  st_union()
plot(st_geometry(all_dist_grids), col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5))

all_dist_grids <- as_Spatial(all_dist_grids) %>% 
  st_as_sf()

all_dist_grids <- st_intersection(grids, all_dist_grids)
plot(st_geometry(all_dist_grids), col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5))

dist_summary_grids <- all_dist_grids %>% 
  mutate(dist_area = st_area(geometry)) %>% 
  group_by(grid_id) %>% 
  summarise(dist_area_m2 = sum(dist_area),
            dist_area_km2 = as.numeric(dist_area_m2/1e6)) %>% 
  left_join(grids_area, by = "grid_id") %>% 
  mutate(prop_dist = dist_area_km2/area)
as_tibble(dist_summary_grids)

# put it all together... but some dataframes don't have all grids so can't just cbind
disturbance_summary_grids <- cbind(dist_summary_grids, fire_summary_grid) %>% 
  left_join(as_tibble(cut_summary_grid), by = "grid_id") %>% 
  left_join(as_tibble(road_summary_grid), by = "grid_id") %>% 
  select(grid_id, area, dist_area_km2, area_fire_km2, area_cut_km2, road_length_km, prop_dist) %>% 
  mutate(prop_fire = as.numeric(area_fire_km2/area),
         prop_cut = as.numeric(area_cut_km2/area),
         dens_road = as.numeric(road_length_km/area),
         prop_dist = as.numeric(prop_dist)) %>% 
  as_tibble() %>%
  select(-geometry.x)
disturbance_summary_grids[is.na(disturbance_summary_grids)] <- 0
disturbance_summary_grids

write_csv(disturbance_summary_grids, "disturbance_summary_grids.csv")
