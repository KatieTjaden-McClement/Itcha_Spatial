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

##### Extract NDVI data #####
bands <- mt_bands(product = "MOD13Q1") # MODIS/Terra Vegetation Indices (NDVI/EVI) 16-Day L3 Global 250m SIN Grid
head(bands)

tmp <- sta %>% 
  select("station_id", "longitude", "latitude")
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

# Remove stuff that is below zero, which are indicative of water
itcha_ndvi <- filter(itcha_ndvi, value > 0)

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
#plot(st_geometry(buff_cut)) #good

buff_cut <- as_Spatial(buff_cut) %>% 
  st_as_sf()

buff_cut <- st_intersection(cams_buff, buff_cut)
#plot(st_geometry(buff_cut))

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
#plot(st_geometry(buff_fire))

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

### Make tibble with extracted covariates
station_sp_cov <- select(sta, station_id)
station_sp_cov <- left_join(station_sp_cov, fire_summary, by = "station_id")
station_sp_cov <- left_join(station_sp_cov, cut_summary, by = "station_id")
station_sp_cov <- left_join(station_sp_cov, road_summary, by = "station_id")
station_sp_cov <- left_join(station_sp_cov, sta, by = "station_id")

station_sp_cov <- select(station_sp_cov, station_id, prop_cut, prop_fire, road_dens_km, Elevation)
station_sp_cov[is.na(station_sp_cov)] <- 0

hist(station_sp_cov$prop_cut)
hist(station_sp_cov$prop_fire)
hist(station_sp_cov$road_dens_km)
hist(station_sp_cov$Elevation)

write.csv(station_sp_cov, file = "Station_spatial_covariates.csv")


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
