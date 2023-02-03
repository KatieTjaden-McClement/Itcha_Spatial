## Exploring Primary Productivity within Itcha and across regions
# Does the vegetation response following disturbance support DMAC?

#mt_to_raster was replaced with mt_to_terra in version 1.1.4
#need to download older verison

devtools::install_version("MODISTools", "1.1.2")

library(tidyverse)
library(sf)
library(leaflet)
library(raster)
library(MODISTools)
library(lubridate)

theme_set(theme_classic())

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

# keep only cutblocks actually in the I-I caribou range:
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

get_evi_rast_ii_range <- function(dates){
  itcha_evi_left <- mt_subset(product = "MOD13Q1",
                              lat = aoi_centre[2] - 0.2,
                              lon = st_bbox(aoi.wgs)[1] + 0.5, # xmin
                              band = "250m_16_days_EVI",
                              start = dates,
                              end = dates,
                              km_lr = 100, 
                              km_ab = 100,
                              internal = TRUE)
  ii_evi_left_raster <- mt_to_raster(itcha_evi_left, 
                                     reproject = F)
  
  itcha_evi_right <- mt_subset(product = "MOD13Q1",
                               lat = aoi_centre[2] + 0.2,
                               lon = st_bbox(aoi.wgs)[3] - 0.5, # xmax
                               band = "250m_16_days_EVI",
                               start = dates,
                               end = dates,
                               km_lr = 100, 
                               km_ab = 100,
                               internal = TRUE)
  ii_evi_right_raster <- mt_to_raster(itcha_evi_right, 
                                      reproject = F)
  
  ii_evi_merged <- raster::merge(ii_evi_right_raster, 
                                 ii_evi_left_raster)
  
  ii_evi_utm <- projectRaster(ii_evi_merged, 
                              crs = crs(II_ranges))
  
  ii_evi_mask <- mask(ii_evi_utm, st_zm(II_ranges))
}

test <- get_evi_rast_ii_range(dates = jul_dates[1]) #start 3pm, finish 3:45
plot(test)

ii_range_jul_evi_all <- lapply(jul_dates, FUN = get_evi_rast_ii_range)
load("ii_range_jul_evi_all.RData")

plot(ii_range_jul_evi_all[[12]]) #looks good!
crs(ii_range_jul_evi_all[[1]]) #all good!

names(ii_range_jul_evi_all) <- jul_dates

#save(II_ranges, aoi_centre, aoi.wgs, file = "ii_ranges_evi_files.RData")

##### Extracting July EVI across II cutblocks #####

# unprocessed list output of get_july_evi_rast_ii_range function lapply loop
load("II_EVI_downloads/ii_range_jul_evi_all.RData")
head(ii_range_jul_evi_all)

plot(ii_range_jul_evi_all[[1]])

#convert list to rasterStack
ii_range_jul_evi_all_stack <- stack(ii_range_jul_evi_all)
head(ii_range_jul_evi_all_stack)
plot(ii_range_jul_evi_all_stack)

# need to average july evi rasters for each year 
# use a group id for each year from jul_dates then average within groups
merge_within_year_groups <- as_tibble(jul_dates) %>% 
  mutate(year = year(ymd(value)),
         month = month(ymd(value))) %>% 
  group_by(year, month) %>% 
  mutate(group = cur_group_id())

ii_range_jul_evi <- stackApply(ii_range_jul_evi_all_stack, 
                               indices = merge_within_year_groups$group,
                               fun = mean) # returns rasterBrick object
plot(ii_range_jul_evi$index_1)

plot(st_geometry(II_ranges))
plot(ii_range_jul_evi$index_1, add = T)

# plot 2000 and 2022 July EVI in the II range for the TWS poster
# ii_park <- filter(II_ranges, BCHab_code == "HEWSR")
# 
# png("evi_raster_2000.png")
# plot(ii_range_jul_evi$index_1, ext = extent(II_ranges),
#      axes = F, box = F)
# plot(st_geometry(sta.utm), add = T, pch = 19, cex = 0.15)
# dev.off()
# 
# 
# png("evi_raster_2021.png")
# plot(ii_range_jul_evi$index_22, ext = extent(II_ranges),
#      axes = F, box = F)
# plot(st_geometry(sta.utm), add = T, pch = 19, cex = 0.15)
# dev.off()
# 
# plot(ii_range_jul_evi$index_23, ext = extent(II_ranges),
#      axes = F, box = F)

# extract average yearly evi across cutblocks of different ages
cut <- st_read(dsn = "Spatial_Layers/CUT_BLOCKS_ranges.gdb")
cut <- st_transform(cut, crs = crs(aoi.utm))

# keep only cutblocks actually in the I-I caribou range:
ii_cut <- st_intersection(cut, II_ranges)

# extracting evi values across years for one cutblock:
ii_cut1 <- slice(ii_cut, 1)

plot(ii_range_jul_evi$index_1)
plot(st_geometry(ii_cut1), add = T)

cut1_yearly_evi <- raster::extract(ii_range_jul_evi,
                             st_zm(ii_cut1), method = "bilinear",
                             fun = mean)
cut1_yearly_evi <- cut1_yearly_evi[,]

cut1_yearly_evi <- cut1_yearly_evi %>% 
  as_tibble() %>% 
  rename(july_evi = value) %>% 
  mutate(cutblock_id = ii_cut1$VEG_CONSOLIDATED_CUT_BLOCK_ID,
         year = c(2000:2022),
         year_since_cut = year - unique(ii_cut1$HARVEST_YEAR))

ggplot(cut1_yearly_evi, aes(x = year_since_cut, y = july_evi)) +
  geom_point() # +
  # geom_smooth()

# make above into function:
extract_yearly_evi_cutblocks <- function(cutblocks, #list of cutblock polygons
                                         evi_raster){
  cut_yearly_evi <- raster::extract(evi_raster,
                                    st_zm(cutblocks), method = "bilinear",
                                    fun = mean)
  cut_yearly_evi <- cut_yearly_evi[,]
  
  cut_yearly_evi <- cut_yearly_evi %>% 
    as_tibble() %>% 
    rename(july_evi = value) %>% 
    mutate(cutblock_id = cutblocks$VEG_CONSOLIDATED_CUT_BLOCK_ID,
           year = c(2000:2022),
           year_since_cut = year - unique(cutblocks$HARVEST_YEAR))
  
  }

# put into lapply to output listm can then rbind into full long data frame

# test with list of first 10 cutblocks
cut_1_10_list <- slice_head(ii_cut, n = 10) 
cut_1_10_list <- split(cut_1_10_list, seq(nrow(cut_1_10_list)))

extract_yearly_evi_cutblocks_test <- lapply(cut_1_10_list, ii_range_jul_evi,
                                            FUN = extract_yearly_evi_cutblocks)
extract_yearly_evi_cutblocks_test <- do.call("rbind", extract_yearly_evi_cutblocks_test)

ggplot(extract_yearly_evi_cutblocks_test, aes(x = year_since_cut,
                                              y = july_evi)) +
  geom_point(aes(colour = as.character(cutblock_id))) +
  geom_smooth(aes(colour = as.character(cutblock_id)), se = F)

# lapply across all cutblocks
ii_cut_list <- split(ii_cut, seq(nrow(ii_cut)))

ii_cut_yearly_evi <- lapply(ii_cut_list, ii_range_jul_evi,
                            FUN = extract_yearly_evi_cutblocks)
ii_cut_yearly_evi_df <- do.call("rbind", ii_cut_yearly_evi)

ggplot(ii_cut_yearly_evi_df, aes(x = year_since_cut,
                                 y = july_evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest harvest",
       y = "Mean July EVI")


##### Extracting July EVI across II burnt areas #####

fire_all <- st_read(dsn = "Spatial_Layers/FIRE_POLYS_ranges.gdb")
fire_all <- st_transform(fire_all, crs = crs(aoi.utm))

ii_fire <- st_intersection(fire_all, II_ranges)
head(ii_fire)

extract_yearly_evi_fires <- function(fires #list of fire polygons
){
  fire_yearly_evi <- raster::extract(ii_range_jul_evi,
                                    st_zm(fires), method = "bilinear",
                                    fun = mean)
  fire_yearly_evi <- fire_yearly_evi[,]
  
  fire_yearly_evi <- fire_yearly_evi %>% 
    as_tibble() %>% 
    rename(july_evi = value) %>% 
    mutate(fire_id = fires$FIRE_LABEL,
           year = c(2000:2022),
           year_since_fire = year - unique(fires$FIRE_YEAR))
}

ii_fire_list <- split(ii_fire, seq(nrow(ii_fire)))

ii_fire_yearly_evi <- lapply(ii_fire_list, 
                             FUN = extract_yearly_evi_fires)
ii_fire_yearly_evi_df <- do.call("rbind", ii_fire_yearly_evi)

ggplot(ii_fire_yearly_evi_df, aes(x = year_since_fire,
                                 y = july_evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since fire",
       y = "Mean July EVI")

### See how intact forest compares to burnt and cut areas in EVI over time

# masking out cut and burnt areas from evi rasterstack
not_burnt <- mask(ii_range_jul_evi, st_zm(ii_fire), inverse = T)
plot(not_burnt$index_1) #looks good!

ii_intact <- mask(not_burnt, st_zm(ii_cut), inverse = T)
plot(ii_intact$index_1) #looks good!

intact_mean_yearly_evi <- data.frame(mean_jul_evi = cellStats(ii_intact, "mean")) 
intact_mean_yearly_evi$year <- c(2000:2022)

ggplot(intact_mean_yearly_evi, aes(x = year, y = mean_jul_evi)) +
  geom_point()

# something going wrong with methods and recentering everything to year since cut?
# plot all 2001 cutblocks in from 2000-2022

ii_cut_2002_ids <- ii_cut %>% 
  filter(HARVEST_YEAR == 2002) %>% 
  pull(VEG_CONSOLIDATED_CUT_BLOCK_ID)

ii_cut_2002 <- ii_cut %>% 
  filter(HARVEST_YEAR == 2002)

plot(st_geometry(ii_cut_2002))

ii_cut_2002_evi <- filter(ii_cut_yearly_evi_df, 
                          cutblock_id %in% ii_cut_2002_ids)

ggplot(ii_cut_2002_evi, aes(x = year_since_cut, y = july_evi)) +
  geom_point() +
  geom_smooth()

# snapshot of EVI in different aged cutblocks in a 2020 rather than longitudinal

ii_cutblock_evi_in_2022 <- filter(ii_cut_yearly_evi_df, year == 2022)

ggplot(ii_cutblock_evi_in_2022, aes(x = year_since_cut, y = july_evi)) +
  geom_point(size = 0.2, alpha = 0.2) +
  geom_smooth() +
  labs(y = "July 2022 EVI", x = "Cutblock age (years)")


## try masking out areas burnt in the last 40 years...


## try taking only the July 27/28 EVI raster for each year like Gagne 2016
jul_dates

# need to average july evi rasters for each year 
# use a group id for each year from jul_dates then average within groups
keep_second_year_groups <- as_tibble(jul_dates) %>% 
  mutate(year = year(ymd(value)),
         month = month(ymd(value))) %>% 
  group_by(year, month) %>% 
  mutate(group = cur_group_id()) %>% 
  slice_tail(n = 1) %>% 
  pull(value)

ii_endJul_evi <- ii_range_jul_evi_all[keep_second_year_groups]
ii_endJul_evi <- stack(ii_endJul_evi)

plot(ii_endJul_evi$X2018.07.28)

ii_cut_yearly_endJul_evi <- lapply(ii_cut_list, ii_endJul_evi,
                                   FUN = extract_yearly_evi_cutblocks)
ii_cut_yearly_endJul_evi_df <- do.call("rbind", ii_cut_yearly_endJul_evi)

ggplot(ii_cut_yearly_endJul_evi_df, aes(x = year_since_cut,
                                        y = july_evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  geom_smooth(data = ii_cut_yearly_evi_df, colour = "red") +
  ylim(c(0, 0.7)) +
  labs(x = "Years since forest harvest",
       y = "July 27/28 EVI")
# using the last July evi vs. average july evi makes basically 
# no difference to the average evi in cuts over the years



##### Delta EVI #####

# extract sept-beginning oct evi to get delta evi
sep_oct_dates <- mt_dates(product = "MOD13Q1", lat = aoi_centre[2],
                          lon = aoi_centre[1]) %>% 
  mutate(month = month(ymd(calendar_date)),
         year = year(ymd(calendar_date)),
         month_day = format(as.Date(calendar_date), "%m-%d")) %>% 
  filter(month %in% c(9,10),
         month_day < "10-20") %>% 
  pull(calendar_date)
head(sep_oct_dates) #69 dates

ii_range_fall_evi_all <- lapply(sept_oct_dates, FUN = get_evi_rast_ii_range)
load("II_EVI_downloads/ii_range_fall_evi_all.RData")

names(ii_range_fall_evi_all) <- sep_oct_dates

ii_range_fall_evi_stack <- stack(ii_range_fall_evi_all)
head(ii_range_fall_evi_stack)
mean(ii_range_fall_evi_stack)

plot(ii_range_fall_evi_stack$X2000.09.13)
plot(ii_range_fall_evi_stack$X2000.09.29)
plot(ii_range_fall_evi_stack$X2000.10.15)

dens(ii_range_fall_evi_stack$X2000.09.13)
hist(ii_range_fall_evi_stack$X2000.09.29, 
     col = "blue", opacity = 0.5, add = T)
hist(ii_range_fall_evi_stack$X2000.10.15, 
     col = "red", opacity = 0.5, add = T)

# Use average sept. EVI as in Serrouya et al. 2021 paper
sep_dates <- mt_dates(product = "MOD13Q1", lat = aoi_centre[2],
                       lon = aoi_centre[1]) %>% 
  mutate(month = month(ymd(calendar_date))) %>% 
  filter(month == 9)

# filter out october rasters
ii_range_sep_evi_all <- raster::subset(ii_range_fall_evi_stack,
                                   grep('\\.09\\.', 
                                        names(ii_range_fall_evi_stack), 
                                        value = T))
names(ii_range_sep_evi_all)

sep_year_groups <- as_tibble(sep_dates) %>% 
  mutate(year = year(ymd(calendar_date)),
         month = month(ymd(calendar_date))) %>% 
  group_by(year, month) %>% 
  mutate(group = cur_group_id())

ii_range_sep_evi <- stackApply(ii_range_sep_evi_all, 
                               indices = sep_year_groups$group,
                               fun = mean)
plot(ii_range_sep_evi$index_18)

# subtract average september evi from average july evi
ii_range_delta_evi <- ii_range_jul_evi - ii_range_sep_evi
plot(ii_range_delta_evi$layer.4)

# look at delta evi across cutblocks
ii_cut_yearly_delta_evi <- lapply(ii_cut_list, ii_range_delta_evi,
                                  FUN = extract_yearly_evi_cutblocks)
ii_cut_yearly_delta_evi_df <- do.call("rbind", ii_cut_yearly_delta_evi)

ii_cut_yearly_delta_evi_df <- ii_cut_yearly_delta_evi_df %>% 
  rename(delta_evi = july_evi)

ggplot(ii_cut_yearly_delta_evi_df, aes(x = year_since_cut,
                                       y = delta_evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest harvest",
       y = "Delta EVI")


##### Other caribou ranges #####
bc_caribou_ranges <- st_read("Spatial_Layers/BC_Caribou_Ranges/GCPB_CARIBOU_POPULATION_SP.gdb/")
head(bc_caribou_ranges)

plot(st_geometry(bc_caribou_ranges))

#filter to southern group of southern mountain caribou
sm_caribou_ranges <- filter(bc_caribou_ranges,
                            HERD_NAME %in% c("Purcells South", "Purcell Central",
                                             "South Selkirks", "Central Selkirks",
                                             "Monashee", "Columbia South",
                                             "Frisby Boulder", "Columbia North",
                                             "Groundhog", "Central Rockies", "Wells Gray South",
                                             "Wells Gray North", "Barkerville", "North Cariboo",
                                             "Narrow Lake", "George Mountain", "Hart Ranges",
                                             "Redrock-Prairie Creek"))
plot(st_geometry(sm_caribou_ranges))

# Start with purcells range - looks like most cutblocks and studied in kinley and Apps 2001
purcells_range <- filter(bc_caribou_ranges,
                         HERD_NAME %in% c("Purcells South", 
                                          "Purcell Central")) %>% 
  st_union() %>%  # merge south and central
  st_transform(crs = 32611) #UTM Zone 11N

plot(st_geometry(purcells_range))

purcells_bbox <- purcells_range %>% 
  st_buffer(1000)
  
purcells_bbox <- st_as_sfc(st_bbox(purcells_bbox))
crs(purcells_range)

plot(purcells_bbox)
plot(st_geometry(purcells_range), add = T)
#st_write(purcells_bbox, dsn = "purcells_bbox.shp")

purc_aoi_centre <- st_geometry(st_centroid(purcells_range)) %>% 
  st_transform(crs = crs(sta.wgs)) %>% 
  unlist()

purc_dates <- mt_dates(product = "MOD13Q1", lat = purc_aoi_centre[2],
                                    lon = purc_aoi_centre[1]) %>% 
  mutate(month = month(ymd(calendar_date))) %>% 
  filter(month %in% c(7, 9)) %>% 
  pull(calendar_date)

#save(purc_aoi_centre, purcells_range, file = "purcells_aoi_and_range.RData")
plot(purcells_range)
purcells_range_sf <- st_as_sf(purcells_range)

# write function for extracting evi
get_evi_purcells_range <- function(dates){
  evi_top <- mt_subset(product = "MOD13Q1",
                       lat = purc_aoi_centre[2] + 0.4,
                       lon = purc_aoi_centre[1] + 0.1, # xmin
                       band = "250m_16_days_EVI",
                       start = dates,
                       end = dates,
                       km_lr = 100, 
                       km_ab = 100,
                       internal = TRUE)
  evi_top_raster <- mt_to_raster(evi_top, 
                                 reproject = F)
  
  evi_bottom <- mt_subset(product = "MOD13Q1",
                          lat = purc_aoi_centre[2] - 0.4,
                          lon = purc_aoi_centre[1] - 0.15, # xmin
                          band = "250m_16_days_EVI",
                          start = dates,
                          end = dates,
                          km_lr = 70, 
                          km_ab = 70,
                          internal = TRUE)
  evi_bottom_raster <- mt_to_raster(evi_bottom, 
                                    reproject = F)
  
  evi_merged <- raster::merge(evi_top_raster, 
                              evi_bottom_raster)
  
  evi_utm <- projectRaster(evi_merged, 
                           crs = crs(purcells_range))
  
  ii_evi_mask <- raster::mask(evi_utm, purcells_range_sf)
}

purcells_range_evi_all <- lapply(purc_dates, FUN = get_evi_purcells_range)
#started 5pm Jan 30

#save(purcells_range_evi_all, file = "purcells_range_evi_all.RData")

load("Spatial_Layers/Productivity_comparison/Purcells/purcells_range_evi_all.RData")
plot(purcells_range_evi_all[[1]])

names(purcells_range_evi_all) <- purc_dates

purcells_evi_all_stack <- stack(purcells_range_evi_all)

# filter to july dates, average within
purc_jul_evi_all <- raster::subset(purcells_evi_all_stack,
                                   grep('\\.07\\.', 
                                        names(purcells_evi_all_stack), 
                                        value = T))

purc_jul_year_groups <- as_tibble(purc_dates) %>% 
  mutate(year = year(ymd(value)),
         month = month(ymd(value))) %>% 
  filter(month == 7) %>% 
  group_by(year) %>% 
  mutate(group = cur_group_id())

purc_jul_evi <- stackApply(purc_jul_evi_all, 
                           indices = purc_jul_year_groups$group,
                           fun = mean)
plot(purc_jul_evi$index_1)

# filter to Sept dates, average within
purc_sep_evi_all <- raster::subset(purcells_evi_all_stack,
                                   grep('\\.09\\.', 
                                        names(purcells_evi_all_stack), 
                                        value = T))

purc_sep_year_groups <- as_tibble(purc_dates) %>% 
  mutate(year = year(ymd(value)),
         month = month(ymd(value))) %>% 
  filter(month == 9) %>% 
  group_by(year) %>% 
  mutate(group = cur_group_id())

purc_sep_evi <- stackApply(purc_sep_evi_all, 
                           indices = purc_sep_year_groups$group,
                           fun = mean)
plot(purc_sep_evi$index_1)

# Delta EVI
purc_delta_evi <- purc_jul_evi - purc_sep_evi
plot(purc_delta_evi$layer.1)

### Extract Jul evi and delta evi across purcells cutblocks and fires
purc_cut_all <- st_read("Spatial_Layers/Productivity_comparison/Purcells/Purcells_Cons_Cutblocks/VEG_CONSOLIDATED_CUT_BLOCKS_SP.gdb") %>% 
  st_transform(crs = crs(purcells_range))
plot(st_geometry(purc_cut_all)) # looks good
plot(st_geometry(purcells_range), add = T)

purc_cut <- st_intersection(purc_cut_all, purcells_range)
# 7441 cutblocks
plot(st_geometry(purc_cut))

purc_cut_list <- split(purc_cut, seq(nrow(purc_cut)))

purc_cut_yearly_jul_evi <- lapply(purc_cut_list, purc_jul_evi,
                                  FUN = extract_yearly_evi_cutblocks)
purc_cut_yearly_jul_evi_df <- do.call("rbind", purc_cut_yearly_jul_evi)

ggplot(purc_cut_yearly_jul_evi_df, aes(x = year_since_cut,
                                       y = july_evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest harvest",
       y = "Mean July EVI")

# compare trends in Itcha vs. Purcells
ii_cut_yearly_evi_df$site <- "Itcha-Ilgachuz"
purc_cut_yearly_jul_evi_df$site <- "Purcells"

jul_evi_range_comp <- rbind(purc_cut_yearly_jul_evi_df,
                            ii_cut_yearly_evi_df)

ggplot(jul_evi_range_comp, aes(x = year_since_cut, y = july_evi,
                               colour = site, group = site)) +
  geom_jitter(size = 0.1, alpha = 0.05, width = 0.2) +
  geom_smooth()


