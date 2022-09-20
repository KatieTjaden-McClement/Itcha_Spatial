## Exploring Primary Productivity within Itcha and across regions
# Does the vegetation response following disturbance support DMAC?

library(tidyverse)
library(sf)
library(leaflet)
library(raster)
library(MODISTools)
library(lubridate)

# load station data for reference 
sta <- as_tibble(read.csv("stations.csv", header=T))
sta.wgs <- st_as_sf(sta, coords = c("longitude", "latitude"), crs=4326)

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

sta.utm <- st_transform(sta.wgs,crs=crs.utm) #crs.utm is from Extract_Spatial script

# area of interest the included entire Itcha Ilgachuz range
aoi <- st_read("bbox_ranges/")
plot(st_geometry(aoi))
plot(st_geometry(sta.utm), add = T)

aoi.wgs <- st_bbox(st_transform(aoi, 4326))

# trying out different packages for accessing daily-scale NDVI/EVI from Landsat
# MODIS 16-day coverage is too coarse and is missing entries

### rLandsat package
#devtools::install_github("atlanhq/rLandsat")
library(rLandsat)

# see what Landsat 8 products are available
landsat_search(min_date = "2017-09-01", max_date = "2019-09-01",
               country = "Canada")

### getSpatialData package
# https://rpubs.com/ials2un/getlandsat

#devtools::install_github("16EAGLE/getSpatialData")
library(getSpatialData)

set_aoi(st_geometry(aoi)) # sets aoi for session
view_aoi()

login_USGS(username = "katietm")

get_products()
evi_records <- getMODIS_records(time_range = c("2000-09-01", "2022-08-22"),
                                products = "modis_mod13q1_v6")
head(evi_records)
#view_records(evi_records)

plot(st_geometry(evi_records[466,]))
plot(st_geometry(sta.wgs), add = T)

itcha_evi_records_july <- evi_records %>% 
  mutate(month = month(start_time, label = T, abbr = T)) %>% 
  filter(month == "Jul")

#order_data(records = itcha_evi_records_july, wait_to_complete = T)
check_availability(itcha_evi_records_july)

set_archive(dir_archive = "MODIS_Data")

login_earthdata(username = "katietm")
services()

evi_files <- getMODIS_data(records = itcha_evi_records_july)

evi_files_test <- raster("MODIS_Data_datasets/modis_mod13q1_v6")
# files are in .hdr format, need to convert to geotiff

#devtools:::install_github("gearslaboratory/gdalUtils")
library(gdalUtils)

gdalinfo("MODIS_Data_datasets/modis_mod13q1_v6/MOD13Q1.A2001193.h10v03.006.2015143124411.hdf")

install.packages("terra")
library(terra)

hds_test <- terra::rast("MODIS_Data_datasets/modis_mod13q1_v6/MOD13Q1.A2001193.h10v03.006.2015143124411.hdf")

#unzip downloaded file
untar("Landsat_Data_datasets/landsat_8_c1/LC08_L1TP_049023_20210425_20210501_01_T1_LEVEL_sr_evi.tar.gz")

#read in evi data for 2021_04_25
evi_2021_04_25 <- raster("Landsat_Data/LC08_L1TP_049023_20210425_20210501_01_T1_sr_evi.tif")
evi_2021_04_25
plot(evi_2021_04_25)

# crop to aoi
crs(aoi)

evi_2021_04_25_cropped <- crop(evi_2021_04_25, aoi)
plot(evi_2021_04_25_cropped)
plot(st_geometry(sta.utm), add = T)

# crop to satah grid to plot - entire raster is too large
satah <- filter(sta.utm, grepl("ITCHA6", station_id))
satah_aoi <- st_buffer(satah, 2000)

evi_2021_04_25_satah <- crop(evi_2021_04_25, satah_aoi)
plot(evi_2021_04_25_satah)
plot(st_geometry(satah), add = T)

normalize <- function(x) {
  min <- raster::minValue(x)
  max <- raster::maxValue(x)
  return(255* (x - min) / (max - min))
}
evi_2021_04_25_satah <- normalize(evi_2021_04_25_satah)

pal <- colorNumeric(c("red", "orange", "yellow", "green", "darkgreen"), 
                    values(evi_2021_04_25_satah),
                    na.color = "transparent")

leaflet() %>% addTiles() %>% 
  addRasterImage(evi_2021_04_25_satah, opacity = 0.7,
                 colors = pal) %>% 
  addLegend(pal = pal, 
            values = values(evi_2021_04_25_satah))

# landsat has much finer spatial resolution (30m), 
# but coarser temporal resolution compared to MODIS - MODIS is better for what I'm trying to do

# look at evi over time for a specific cutblock
cut <- st_read(dsn = "Spatial_Layers/CUT_BLOCKS_ranges.gdb")
cut <- st_transform(cut, crs = crs(aoi.utm)) #UTM zone 10
head(cut)
plot(st_geometry(cut))

cut_2003 <- filter(cut, HARVEST_YEAR == 2003)
plot(st_geometry(sta.utm))
plot(st_geometry(cut_2003), add = T)

cut_2003_wgs <- st_transform(cut_2003, crs = crs(sta.wgs))

leaflet(cut_2003_wgs) %>% 
  addTiles() %>% 
  addPolygons(label = cut_2003_wgs$VEG_CONSOLIDATED_CUT_BLOCK_ID) %>% 
  addCircleMarkers(lng = sta$longitude, lat = sta$latitude,
                   color = "red")

cutblock1_2003 <- filter(cut, VEG_CONSOLIDATED_CUT_BLOCK_ID == "249034")
# fairly large cutblock in top left corner of Satah grid

cutblock1_centre <- st_geometry(st_centroid(cutblock1_2003)) %>% 
  st_transform(crs = crs(sta.wgs)) %>% 
  unlist()

dates <- mt_dates(product = "MOD13Q1", lat = cutblock1_centre[2],
                  lon = cutblock1_centre[1])
head(dates$calendar_date, 1)
tail(dates$calendar_date, 1)

cutblock1_evi <- mt_subset(product = "MOD13Q1",
                           lat = cutblock1_centre[2],
                           lon = cutblock1_centre[1],
                           band = "250m_16_days_EVI",
                           start = "2000-02-18",
                           end = "2022-07-12",
                           km_lr = 5,
                           km_ab = 5,
                           internal = TRUE)
head(cutblock1_evi)

cutblock1_evi <- filter(cutblock1_evi, value >= -2000)

cutblock1_evi$calendar_date <- as.POSIXct(cutblock1_evi$calendar_date)

cutblock1_evi <- cutblock1_evi %>% 
  mutate(month = month(calendar_date, label = T, abbr = T),
         year = year(calendar_date))
  
ggplot(cutblock1_evi, aes(x = calendar_date, y = value, 
                          group = month, colour = month)) +
  geom_vline(xintercept = as.POSIXct("2003-03-01"),
             colour = "red", linetype = "dashed") +
  geom_point() +
  geom_line() +
  scale_color_discrete(type = "viridis")

cutblock1_evi_july <- cutblock1_evi %>% 
  filter(month == "Jul") %>% 
  group_by(year, pixel) %>% 
  summarise(jul_evi = mean(value))
head(cutblock1_evi_july)

cutblock1_evi_july <- right_join(cutblock1_evi, cutblock1_evi_july,
                                 by = c("year", "pixel"))

cutblock1_evi_july <- cutblock1_evi_july %>% 
  filter(month == "Jul") %>%
  mutate(year_pixel = paste(year, pixel, sep = "_")) %>% 
  distinct(year_pixel, .keep_all = T)


cutblock1_evi_july_raster <- mt_to_raster(cutblock1_evi_july)
#some years have multiple rasters, with the second one having obvious errors

test <- filter(cutblock1_evi_july, calendar_date == "2000-07-27")
# not all pixels have valid evi estimates for the first observation of the month,
# so mt_to_raster is creating a raster from the later date for those few entries

# try converting to rasters before averaging july evi values:
cut1_evi_all_raster <- mt_to_raster(cutblock1_evi)

filter_to_july <- as_tibble(names(cut1_evi_all_raster)) %>% 
  mutate(date = ymd(gsub("X", "", value)),
         year = year(date),
         month = month(date, label = T, abbr = T)) %>% 
  filter(month == "Jul")

cut1_evi_jul_rasters <- raster::subset(cut1_evi_all_raster,
                                       filter_to_july$value)
cut1_evi_jul_rasters <- projectRaster(cut1_evi_jul_rasters, crs = crs(sta.utm))

cut1_evi_jul2000_mask <- mask(cut1_evi_jul_rasters$X2000.07.11, cutblock1_2003)
plot(cut1_evi_jul2000_mask)
plot(st_geometry(cutblock1_2003), add = T)

plot(cut1_evi_jul_rasters$X2000.07.11)
plot(cut1_evi_jul_rasters$X2021.07.12)
plot(st_geometry(cutblock1_2003), add = T)

cutblock1_wgs <- st_transform(cutblock1_2003, crs = crs(sta.wgs))
cutblock1_2002_evi_raster_wgs <- projectRaster(cutblock1_evi_july_raster$X2002.07.12,
                                     crs = crs(sta.wgs))
cutblock1_2004_evi_raster_wgs <- projectRaster(cutblock1_evi_july_raster$X2004.07.11,
                                               crs = crs(sta.wgs))
cutblock1_2022_evi_raster_wgs <- projectRaster(cutblock1_evi_july_raster$X2022.07.12,
                                               crs = crs(sta.wgs))

pal <- colorNumeric(c("red", "orange", "yellow", "green", "darkgreen"), 
                    values(cutblock1_evi_raster_wgs),
                    na.color = "transparent")

leaflet(cutblock1_wgs) %>% 
  addTiles() %>% 
  addPolygons() %>% 
  addRasterImage(cutblock1_evi_raster_wgs, opacity = 0.7,
                 colors = pal) %>% 
  addLegend(pal = pal, 
            values = values(cutblock1_2022_evi_raster_wgs))

# average across july evi valeus for each year
filter_to_july <- filter_to_july %>% 
  group_by(year, month) %>% 
  mutate(group = cur_group_id())

cutblock1_july_evi_means <- stackApply(cut1_evi_jul_rasters, 
                                       indices = filter_to_july$group,
                                       fun = mean) # returns rasterBrick object
plot(cutblock1_july_evi_means$index_8)

# check that it averaged values like expected:
cutblock1_july_evi_means$index_1@data@values[2000:2020]
cut1_evi_jul_rasters$X2000.07.11@data@values[2000:2020]
cut1_evi_jul_rasters$X2000.07.27@data@values[2000:2020]
# looks good!

# extract average evi values across the cutblock for each year
cutblock1_jul_evi <- raster::extract(cutblock1_july_evi_means,
                                     cutblock1_2003, method = "bilinear",
                                     fun = mean)
cutblock1_jul_evi <- cutblock1_jul_evi[,]
  
cutblock1_jul_evi <- cutblock1_jul_evi %>% 
  as_tibble() %>% 
  rename(july_evi = value) %>% 
  mutate(cutblock_id = cutblock1_2003$VEG_CONSOLIDATED_CUT_BLOCK_ID,
         year = c(2000:2022),
         year_since_cut = year - 2003)

ggplot(cutblock1_jul_evi, aes (x = year, y = july_evi)) +
  geom_vline(xintercept = 2003,
             color = "red", linetype = "dashed") +
  geom_point()

# convert above workflow into a function
# then can use *apply to run for each cutblock

extract_jul_evi_dist <- function(dist_shp){
  
  dist_centre <- st_geometry(st_centroid(dist_shp)) %>% 
    st_transform(crs = "+proj=longlat +datum=WGS84") %>% 
    unlist()
  
  dist_bbox <- st_bbox(dist_shp)
  ab <- ceiling((dist_bbox[4] - dist_bbox[2])/1000*1.2)
  lr <- ceiling((dist_bbox[3] - dist_bbox[1])/1000*1.2)
  
  dates <- mt_dates(product = "MOD13Q1", lat = dist_centre[2],
                    lon = dist_centre[1])
  
  dist_evi_all <- mt_subset(product = "MOD13Q1",
                            lat = dist_centre[2],
                            lon = dist_centre[1],
                            band = "250m_16_days_EVI",
                            start = head(dates$calendar_date, 1),
                            end = tail(dates$calendar_date, 1),
                            km_lr = lr,
                            km_ab = ab,
                            internal = TRUE)
  
  dist_evi_all <- filter(dist_evi_all, value >= -2000)
  
  dist_evi_all_raster <- mt_to_raster(dist_evi_all)
  
  filter_to_july <- as_tibble(names(dist_evi_all_raster)) %>% 
    mutate(date = ymd(gsub("X", "", value)),
           year = year(date),
           month = month(date, label = T, abbr = T)) %>% 
    filter(month == "Jul") %>% 
    group_by(year, month) %>% 
    mutate(group = cur_group_id())
  
  dist_evi_jul_rasters <- raster::subset(dist_evi_all_raster,
                                         filter_to_july$value)
  dist_evi_jul_rasters <- projectRaster(dist_evi_jul_rasters, 
                                        crs = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")
  
  dist_july_evi_means <- stackApply(dist_evi_jul_rasters, 
                                    indices = filter_to_july$group,
                                    fun = mean)
  
  dist_shp <- st_zm(dist_shp)
  dist_jul_evi <- raster::extract(dist_july_evi_means,
                                  dist_shp, method = "bilinear",
                                  fun = mean)
  
  dist_jul_evi <- dist_jul_evi %>% 
    as_tibble() 
  
  dist_jul_evi <-  pivot_longer(data = dist_jul_evi,
                                cols = names(dist_jul_evi), 
                                names_to = "group",
                                names_prefix = "index_", 
                                values_to = "july_evi") %>% 
    mutate(id = dist_shp$id,
           year = c(year(head(dates$calendar_date, 1)):
                      year(tail(dates$calendar_date, 1))))
  
}

cutblock1_2003 <- cutblock1_2003 %>% 
  rename(id = VEG_CONSOLIDATED_CUT_BLOCK_ID)

fun_test <- extract_jul_evi(dist_shp = cutblock1_2003)
# can add harvest year/other attributes back in based on id column after

# keep only cutblock actually in the I-I caribou range:
II_cut <- st_intersection(cut, II_ranges) # removes 7321 cutblocks
plot(st_geometry(II_cut))
summary(II_cut$AREA_HA)

II_cut <- II_cut %>% 
  rename(id = VEG_CONSOLIDATED_CUT_BLOCK_ID)

cut_2003 <- filter(II_cut, HARVEST_YEAR == 2003) # 264 cutblocks
cut_2003_list <- split(cut_2003, seq(nrow(cut_2003)))

cut_2003_jul_evi <- lapply(cut_2003_list, FUN = extract_jul_evi)
# Started around 6:30pm - think it got interuppted sometime ~ 1am bc laptop went to sleep,
# got through ~ 50 cutblocks
# Error in curl::curl_fetch_memory(url, handle = handle) : 
# Operation was aborted by an application callback

cut_2001 <- filter(II_cut, HARVEST_YEAR == 2001) # 305 cutblocks
cut_2002 <- filter(II_cut, HARVEST_YEAR == 2002) # 403 cutblocks
cut_2004 <- filter(II_cut, HARVEST_YEAR == 2004) # 232 cutblocks
cut_2005 <- filter(II_cut, HARVEST_YEAR == 2005) # 294 cutblocks

# this is going to take forever and probably downloading overlapping areas over and over again..
# try to download entire time series for entire range then extract cutblock values....

# tried download from USGS EarthExplorer



# try a larger spatial scale for a single temporal period

# look at area around chezacut
chez <- filter(sta.utm, grepl("ITCHA5", station_id))
chez_aoi <- st_convex_hull(st_union(st_buffer(chez, 2000)))
plot(st_geometry(chez_aoi))

chez_centre <- st_geometry(st_centroid(chez_aoi)) %>% 
  st_transform(crs = crs(sta.wgs)) %>% 
  unlist()

chez_evi_dates <- mt_dates(product = "MOD13Q1",
                            lat = chez_centre[2],
                            lon = chez_centre[1])

chez_evi_Jul2021 <- mt_subset(product = "MOD13Q1",
                               lat = chez_centre[2],
                               lon = chez_centre[1],
                               band = "250m_16_days_EVI",
                               start = "2021-07-01",
                               end = "2021-07-15",
                               km_lr = 20,
                               km_ab = 12,
                               internal = TRUE)

chez_evi_Jul2021_raster <- mt_to_raster(chez_evi_Jul2021, reproject = T)
chez_evi_Jul2021_raster <- projectRaster(chez_evi_Jul2021_raster$X2021.07.12,
                                          crs = crs(chez_aoi))

plot(chez_evi_Jul2021_raster)
plot(st_geometry(chez_aoi), add = T)

# crop to cutblocks in the buffered Satah grid
chez_cut <- st_crop(cut, chez_aoi)
plot(st_geometry(chez_cut))
plot(st_geometry(chez_aoi), add = T) # looks good

# create mask of evi values to cutblock areas
# satah_evi_cut <- mask(satah_evi_Jul2021_raster, satah_cut)
# plot(satah_evi_cut)
# plot(st_geometry(satah_aoi), add = T) #looks good

mean_cut_evi <- raster::extract(chez_evi_Jul2021_raster,
                                chez_cut, method = "simple",
                                fun = mean)
chez_cut$mean_evi <- mean_cut_evi

ggplot(chez_cut, aes(x = HARVEST_YEAR, y = mean_evi)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = mean(chez_intact@data@values, na.rm = T),
             show.legend = "intact forest",
             colour = "darkgreen")

chez_intact <- mask(chez_evi_Jul2021_raster$X2021.07.12, 
                    chez_cut, inverse = T)
plot(chez_intact)

mean(chez_intact@data@values, na.rm = T) #0.2790723

# do above for entire I-I range
st_area(aoi)/1e6 # 33,325 km^2 total area
# download evi data for area of 200km by 250km

aoi_centre <- st_geometry(st_centroid(aoi)) %>% 
  st_transform(crs = crs(sta.wgs)) %>% 
  unlist()

itcha_evi_Jul2021 <- mt_subset(product = "MOD13Q1",
                              lat = aoi_centre[2],
                              lon = aoi_centre[1],
                              band = "250m_16_days_EVI",
                              start = "2021-07-01",
                              end = "2021-07-15",
                              km_lr = 100, #can't dowload any larger area
                              km_ab = 100,
                              internal = TRUE)
itcha_evi_Jul2021 <- filter(itcha_evi_Jul2021, value >= -2000)

ii_evi_Jul2021_raster <- mt_to_raster(itcha_evi_Jul2021, reproject = T)
#jul_mean <- calc(ii_evi_Jul2021_raster, fun = mean)

ii_evi_Jul2021_raster <- projectRaster(ii_evi_Jul2021_raster$X2021.07.12,
                                         crs = crs(aoi.utm))

plot(ii_evi_Jul2021_raster)
plot(st_geometry(aoi), add = T)
plot(st_geometry(sta.utm), add = T)
# covers all camera grids with a buffer, but missing the top corners, check in relation to seasonal ranges

II_ranges <- st_read(dsn = "Spatial_Layers/Itcha_seasonal_ranges.shp")
II_ranges <- st_transform(II_ranges, crs = crs(aoi.utm))

plot(ii_evi_Jul2021_raster)
plot(st_geometry(II_ranges), add = T)
# not quite making it to the edges of the range, need to download strips on each side and merge then to fill out the whole range

# start with just one time point to test the coverage needed and merging method:
st_bbox(aoi.wgs)[3]

# try to cover right and then left sides of range - just need to  download and merge two rasters
itcha_evi_Jul2021_right <- mt_subset(product = "MOD13Q1",
                               lat = aoi_centre[2] + 0.2,
                               lon = st_bbox(aoi.wgs)[3] - 0.5, # xmax
                               band = "250m_16_days_EVI",
                               start = "2021-07-01",
                               end = "2021-07-15",
                               km_lr = 100, 
                               km_ab = 100,
                               internal = TRUE)
# started 11:58am, finished 12:28 pm

# remove fill values:
itcha_evi_Jul2021_right <- filter(itcha_evi_Jul2021_right, value >= -2000)

ii_evi_Jul2021_right_raster <- mt_to_raster(itcha_evi_Jul2021_right, 
                                            reproject = T)
ii_evi_Jul2021_right_raster <- projectRaster(ii_evi_Jul2021_right_raster$X2021.07.12,
                                       crs = crs(aoi.utm))
plot(ii_evi_Jul2021_right_raster)
plot(st_geometry(II_ranges), add = T)
#plot(ii_evi_Jul2021_raster, add = T)

# ii_evi_Jul2021_right_raster <- projectRaster(ii_evi_Jul2021_right_raster, ii_evi_Jul2021_raster)
# 
# ii_evi_Jul2021_right_merge <- raster::merge(ii_evi_Jul2021_raster, ii_evi_Jul2021_right_raster)
# plot(ii_evi_Jul2021_right_merge)
# 
# ii_evi_Jul2021_right_merge_plot <- as.data.frame(ii_evi_Jul2021_right_merge,
#                                                  xy = TRUE)
# 
# ggplot() +
#   geom_raster(ii_evi_Jul2021_right_merge_plot,
#             mapping = aes(x = x, y = y, fill = layer)) +
#   geom_sf(data = II_ranges, alpha = 0.5) +
#   coord_sf() 

itcha_evi_Jul2021_left <- mt_subset(product = "MOD13Q1",
                                     lat = aoi_centre[2] - 0.2,
                                     lon = st_bbox(aoi.wgs)[1] + 0.5, # xmin
                                     band = "250m_16_days_EVI",
                                     start = "2021-07-01",
                                     end = "2021-07-15",
                                     km_lr = 100, 
                                     km_ab = 100,
                                     internal = TRUE)
itcha_evi_Jul2021_left <- filter(itcha_evi_Jul2021_left, value >= -2000)

ii_evi_Jul2021_left_raster <- mt_to_raster(itcha_evi_Jul2021_left, 
                                            reproject = T)

ii_evi_Jul2021_left_raster <- projectRaster(ii_evi_Jul2021_left_raster$X2021.07.12,
                                             crs = crs(aoi.utm))

plot(ii_evi_Jul2021_left_raster)
plot(ii_evi_Jul2021_right_raster, add = T)
plot(st_geometry(II_ranges), add = T)

template <- projectRaster(from = ii_evi_Jul2021_left_raster,
                          to = ii_evi_Jul2021_right_raster,
                          alignOnly = TRUE)

ii_evi_Jul2021_left_raster_aligned <- projectRaster(from = ii_evi_Jul2021_left_raster, 
                                                    to = template)
plot(ii_evi_Jul2021_left_raster_aligned)
plot(ii_evi_Jul2021_right_raster, add = T)

ii_evi_Jul2021_merged <- raster::merge(ii_evi_Jul2021_right_raster, 
                                       ii_evi_Jul2021_left_raster_aligned)

plot(ii_evi_Jul2021_merged)
plot(st_geometry(II_ranges), add = T)

ii_evi_mask <- mask(ii_evi_Jul2021_merged, st_zm(II_ranges))
plot(ii_evi_mask)

# make above into function to feed in just July dates 
jul_dates <- mt_dates(product = "MOD13Q1", lat = aoi_centre[2],
                  lon = aoi_centre[1]) %>% 
  mutate(month = month(ymd(calendar_date))) %>% 
  filter(month == 7) %>% 
  pull(calendar_date)
head(jul_dates) #45 July dates

head(jul_dates$calendar_date, 1)
tail(jul_dates$calendar_date, 1)

get_july_evi_rast_ii_range <- function(dates){
  itcha_evi_left <- mt_subset(product = "MOD13Q1",
                              lat = aoi_centre[2] - 0.2,
                              lon = st_bbox(aoi.wgs)[1] + 0.5, # xmin
                              band = "250m_16_days_EVI",
                              start = dates,
                              end = dates,
                              km_lr = 100, 
                              km_ab = 100,
                              internal = TRUE)
  itcha_evi_left <- filter(itcha_evi_left, value >= -2000)
  
  ii_evi_left_raster <- mt_to_raster(itcha_evi_left, 
                                     reproject = T)
  ii_evi_left_raster <- projectRaster(ii_evi_left_raster, 
                                      crs = crs(aoi.utm))
  
  itcha_evi_right <- mt_subset(product = "MOD13Q1",
                               lat = aoi_centre[2] + 0.2,
                               lon = st_bbox(aoi.wgs)[3] - 0.5, # xmax
                               band = "250m_16_days_EVI",
                               start = dates,
                               end = dates,
                               km_lr = 100, 
                               km_ab = 100,
                               internal = TRUE)
  itcha_evi_right <- filter(itcha_evi_right, value >= -2000)
  
  ii_evi_right_raster <- mt_to_raster(itcha_evi_right, 
                                      reproject = T)
  ii_evi_right_raster <- projectRaster(ii_evi_right_raster, 
                                       crs = crs(aoi.utm))
  
  template <- projectRaster(from = ii_evi_left_raster,
                            to = ii_evi_right_raster,
                            alignOnly = TRUE)
  
  ii_evi_left_raster_aligned <- projectRaster(from = ii_evi_left_raster, 
                                              to = template)
  ii_evi_merged <- raster::merge(ii_evi_right_raster, 
                                 ii_evi_left_raster_aligned)
  ii_evi_mask <- mask(ii_evi_merged, st_zm(II_ranges))
}

test <- get_july_evi_rast_ii_range(dates = jul_dates[1])
plot(test)

ii_range_jul_evi_all <- lapply(jul_dates, FUN = get_july_evi_rast_ii_range)
names(ii_range_jul_evi_all) <- jul_dates

test_rasterstack <- do.call(stack, test2)
