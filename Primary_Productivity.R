## Exploring Primary Productivity within Itcha and across regions
# Does the vegetation response following disturbance support DMAC?

#mt_to_raster was replaced with mt_to_terra in version 1.1.4
#need to download older verison

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(leaflet)
library(raster)
library(MODISTools)
library(rgdal)
library(bcmaps)
library(lubridate)
library(NatParksPalettes)
library(cowplot)

library(lme4)
library(merTools)

library(beepr)

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

### Download EVI data for Itcha-Ilgachuz caribou range

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

# remove fill values: - update - don't do this, it messes things up
#itcha_evi_Jul2021_right <- filter(itcha_evi_Jul2021_right, value >= -2000)

ii_evi_Jul2021_right_raster <- mt_to_raster(itcha_evi_Jul2021_right, 
                                            reproject = T)
ii_evi_Jul2021_right_raster <- projectRaster(ii_evi_Jul2021_right_raster$X2021.07.12,
                                       crs = crs(aoi.utm))
plot(ii_evi_Jul2021_right_raster)
plot(st_geometry(II_ranges), add = T)
#plot(ii_evi_Jul2021_raster, add = T)

itcha_evi_Jul2021_left <- mt_subset(product = "MOD13Q1",
                                     lat = aoi_centre[2] - 0.2,
                                     lon = st_bbox(aoi.wgs)[1] + 0.5, # xmin
                                     band = "250m_16_days_EVI",
                                     start = "2021-07-01",
                                     end = "2021-07-15",
                                     km_lr = 100, 
                                     km_ab = 100,
                                     internal = TRUE)

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

<<<<<<< HEAD
ii_range_jul_evi_all <- lapply(jul_dates, FUN = get_july_evi_rast_ii_range) # started 9:30 Sept. 19
save(ii_range_jul_evi_all, file = "ii_range_evi_all.RData")
=======
ii_range_jul_evi_all <- lapply(jul_dates, FUN = get_evi_rast_ii_range)
#load("II_EVI_downloads/ii_range_jul_evi_all.RData")

plot(ii_range_jul_evi_all[[12]]) #looks good!
crs(ii_range_jul_evi_all[[1]]) #all good!
>>>>>>> dbda7120feb87df66523e1acca57c11b7453167f

names(ii_range_jul_evi_all) <- jul_dates

#save(II_ranges, aoi_centre, aoi.wgs, file = "ii_ranges_evi_files.RData")

##### Extracting July EVI across II cutblocks #####

# unprocessed list output of get_july_evi_rast_ii_range function lapply loop
load("Spatial_Layers/Productivity_comparison/Itcha_Ilgachuz/ii_range_jul_evi_all.RData")
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
                                         evi_raster, #rasterbrick of yearly evi within range
                                         type #time of evi values/absolute of delta evi?
){
  cut_yearly_evi <- raster::extract(evi_raster,
                                    st_zm(cutblocks), method = "bilinear",
                                    fun = mean)
  cut_yearly_evi <- cut_yearly_evi[,]
  
  cut_yearly_evi <- cut_yearly_evi %>% 
    as_tibble() %>% 
    rename(evi = value) %>% 
    mutate(cutblock_id = cutblocks$VEG_CONSOLIDATED_CUT_BLOCK_ID,
           year = c(2000:2022),
           year_since_cut = year - unique(cutblocks$HARVEST_YEAR),
           evi_type = type)
  
}

# put into lapply to output list can then rbind into full long data frame

# test with list of first 10 cutblocks
cut_1_10_list <- slice_head(ii_cut, n = 10) 
cut_1_10_list <- split(cut_1_10_list, seq(nrow(cut_1_10_list)))

extract_yearly_evi_cutblocks_test <- lapply(cut_1_10_list, ii_range_jul_evi,
                                            "july_evi",
                                            FUN = extract_yearly_evi_cutblocks)
extract_yearly_evi_cutblocks_test <- do.call("rbind", extract_yearly_evi_cutblocks_test)

ggplot(extract_yearly_evi_cutblocks_test, aes(x = year_since_cut,
                                              y = evi)) +
  geom_point(aes(colour = as.character(cutblock_id))) +
  geom_smooth(aes(colour = as.character(cutblock_id)), se = F)

# lapply across all cutblocks
ii_cut_list <- split(ii_cut, seq(nrow(ii_cut)))

ii_cut_yearly_evi <- lapply(ii_cut_list, ii_range_jul_evi, "july_evi",
                            FUN = extract_yearly_evi_cutblocks)
ii_cut_yearly_evi_df <- do.call("rbind", ii_cut_yearly_evi)

ggplot(ii_cut_yearly_evi_df, aes(x = year_since_cut,
                                 y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest harvest",
       y = "Mean July EVI")


##### Extracting July EVI across II burnt areas #####

fire_all <- st_read(dsn = "Spatial_Layers/FIRE_POLYS_ranges.gdb")
fire_all <- st_transform(fire_all, crs = crs(aoi.utm))

ii_fire <- st_intersection(fire_all, II_ranges)
head(ii_fire)

extract_yearly_evi_fires <- function(fires, #list of fire polygons
                                     evi_raster,
                                     type
){
  fire_yearly_evi <- raster::extract(evi_raster,
                                     st_zm(fires), method = "bilinear",
                                     fun = mean)
  fire_yearly_evi <- fire_yearly_evi[,]
  
  fire_yearly_evi <- fire_yearly_evi %>% 
    as_tibble() %>% 
    rename(evi = value) %>% 
    mutate(fire_id = fires$FIRE_LABEL,
           year = c(2000:2022),
           year_since_fire = year - unique(fires$FIRE_YEAR),
           evi_type = type)
}

ii_fire_list <- split(ii_fire, seq(nrow(ii_fire)))

ii_fire_yearly_evi <- lapply(ii_fire_list, ii_range_jul_evi,
                             "july_evi",
                             FUN = extract_yearly_evi_fires) 
ii_fire_yearly_evi_df <- do.call("rbind", ii_fire_yearly_evi)

ggplot(ii_fire_yearly_evi_df, aes(x = year_since_fire,
                                 y = evi)) +
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

# plot all 2001 cutblocks in from 2000-2022
ii_cut_2002_ids <- ii_cut %>% 
  filter(HARVEST_YEAR == 2002) %>% 
  pull(VEG_CONSOLIDATED_CUT_BLOCK_ID)

ii_cut_2002 <- ii_cut %>% 
  filter(HARVEST_YEAR == 2002)
plot(st_geometry(ii_cut_2002))

ii_cut_2002_evi <- filter(ii_cut_yearly_evi_df, 
                          cutblock_id %in% ii_cut_2002_ids)

ggplot(ii_cut_2002_evi, aes(x = year_since_cut, y = evi)) +
  geom_point() +
  geom_smooth()

# snapshot of EVI in different aged cutblocks in a 2020 rather than longitudinal
ii_cutblock_evi_in_2022 <- filter(ii_cut_yearly_evi_df, year == 2022)

ggplot(ii_cutblock_evi_in_2022, aes(x = year_since_cut, y = evi)) +
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
                                  "delta_evi",
                                  FUN = extract_yearly_evi_cutblocks)
ii_cut_yearly_delta_evi_df <- do.call("rbind", ii_cut_yearly_delta_evi)

ggplot(ii_cut_yearly_delta_evi_df, aes(x = year_since_cut,
                                       y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest harvest",
       y = "Delta EVI")

# delta evi across burns
ii_fire_yearly_delta_evi <- lapply(ii_fire_list, ii_range_delta_evi,
                                  "delta_evi",
                                  FUN = extract_yearly_evi_fires)
ii_fire_yearly_delta_evi_df <- do.call("rbind", ii_fire_yearly_delta_evi)

ggplot(ii_fire_yearly_delta_evi_df, aes(x = year_since_fire,
                                        y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest fire",
       y = "Delta EVI")

##### Other caribou ranges #####
bc_caribou_ranges <- st_read("Spatial_Layers/Productivity_comparison/BC_Caribou_Ranges/GCPB_CARIBOU_POPULATION_SP.gdb/")
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

##### Purcells #####
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
#started 5pm Jan 30, finished by 9:30am Feb 3

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

### Extract Jul evi across purcells cutblocks and fires
purc_cut_all <- st_read("Spatial_Layers/Productivity_comparison/Purcells/Purcells_Cons_Cutblocks/VEG_CONSOLIDATED_CUT_BLOCKS_SP.gdb") %>% 
  st_transform(crs = crs(purcells_range))
# plot(st_geometry(purc_cut_all)) # looks good
# plot(st_geometry(purcells_range), add = T)

purc_cut <- st_intersection(purc_cut_all, purcells_range)
# 7441 cutblocks
plot(st_geometry(purc_cut))

purc_cut_list <- split(purc_cut, seq(nrow(purc_cut)))

purc_cut_yearly_jul_evi <- lapply(purc_cut_list, purc_jul_evi, "july_evi",
                                  FUN = extract_yearly_evi_cutblocks)
purc_cut_yearly_jul_evi_df <- do.call("rbind", purc_cut_yearly_jul_evi)

ggplot(purc_cut_yearly_jul_evi_df, aes(x = year_since_cut,
                                       y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest harvest",
       y = "Mean July EVI")

## Fires
purc_fire_all <- st_read("Spatial_Layers/Productivity_comparison/Purcells/Purcells_Fires/PROT_HISTORICAL_FIRE_POLYS_SP.gdb/") %>% 
  st_transform(crs = crs(purcells_range))
plot(st_geometry(purc_fire_all)) # looks good
plot(st_geometry(purcells_range), add = T)

purc_fire <- st_intersection(purc_fire_all, purcells_range)
# 600 fires
plot(st_geometry(purc_fire))

purc_fire_list <- split(purc_fire, seq(nrow(purc_fire)))

purc_fire_yearly_jul_evi <- lapply(purc_fire_list, purc_jul_evi, "july_evi",
                                   FUN = extract_yearly_evi_fires)
purc_fire_yearly_jul_evi_df <- do.call("rbind", purc_fire_yearly_jul_evi)

ggplot(purc_fire_yearly_jul_evi_df, aes(x = year_since_fire,
                                        y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since fire",
       y = "Mean July EVI")

### Extract delta evi across purcells cutblocks and fires
# cutblocks
purc_cut_yearly_delta_evi <- lapply(purc_cut_list, purc_delta_evi, "delta_evi",
                                  FUN = extract_yearly_evi_cutblocks)
purc_cut_yearly_delta_evi_df <- do.call("rbind", purc_cut_yearly_delta_evi)

ggplot(purc_cut_yearly_delta_evi_df, aes(x = year_since_cut,
                                       y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest harvest",
       y = "Delta EVI")

ii_cut_yearly_delta_evi_df$site <- "Itcha-Ilgachuz"
purc_cut_yearly_delta_evi_df$site <- "Purcells"

cut_delta_evi_range_comp <- rbind(purc_cut_yearly_delta_evi_df,
                                  ii_cut_yearly_delta_evi_df)

ggplot(cut_delta_evi_range_comp, aes(x = year_since_cut, y = evi,
                                   colour = site, group = site)) +
  geom_jitter(size = 0.1, alpha = 0.05, width = 0.2) +
  geom_smooth() +
  scale_color_manual(values = natparks.pals("Saguaro", 2)) +
  labs(x = "Years since forest harvest",
       y = "Delta EVI",
       colour = "")

# fires
purc_fire_yearly_delta_evi <- lapply(purc_fire_list, purc_delta_evi, "delta_evi",
                                    FUN = extract_yearly_evi_fires)
purc_fire_yearly_delta_evi_df <- do.call("rbind", purc_fire_yearly_delta_evi)

ggplot(purc_fire_yearly_delta_evi_df, aes(x = year_since_fire,
                                         y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest fire",
       y = "Delta EVI")


##### Barkerville #####

bv_range <- filter(bc_caribou_ranges,
                   HERD_NAME == "Barkerville") %>%
  st_transform(crs = 32610) #UTM zone 10N
plot(st_geometry(bv_range))

bv_bbox <- bv_range %>% 
  st_buffer(1000)

bv_bbox <- st_as_sfc(st_bbox(bv_bbox))
crs(bv_range)

plot(bv_bbox)
plot(st_geometry(bv_range), add = T)
#st_write(bv_bbox, dsn = "Spatial_Layers/Productivity_comparison/Barkerville/bv_bbox.shp")

bv_aoi_centre <- st_geometry(st_centroid(bv_range)) %>% 
  st_transform(crs = crs(sta.wgs)) %>% 
  unlist()

bv_dates <- mt_dates(product = "MOD13Q1", lat = bv_aoi_centre[2],
                       lon = bv_aoi_centre[1]) %>% 
  mutate(month = month(ymd(calendar_date))) %>% 
  filter(month %in% c(7, 9)) %>% 
  pull(calendar_date)

# save(bv_aoi_centre, bv_range, file = "Spatial_Layers/Productivity_comparison/Barkerville/bv_aoi_and_range.RData")
# plot(bv_range)
# bv_range_sf <- st_as_sf(bv_range)

# write function for extracting evi
get_evi_bv_range <- function(dates){
  evi <- mt_subset(product = "MOD13Q1",
                       lat = bv_aoi_centre[2], 
                       lon = bv_aoi_centre[1], # xmin
                       band = "250m_16_days_EVI",
                       start = dates,
                       end = dates,
                       km_lr = 90, 
                       km_ab = 50,
                       internal = TRUE)
  evi_raster <- mt_to_raster(evi, 
                             reproject = F)
  
  evi_utm <- projectRaster(evi_raster, 
                           crs = crs(bv_range))
  
  bv_evi_mask <- raster::mask(evi_utm, bv_range)
}

bv_range_evi_all <- lapply(bv_dates, FUN = get_evi_bv_range)
#save(bv_range_evi_all, file = "bv_range_evi_all.RData")

load("Spatial_Layers/Productivity_comparison/Barkerville/bv_range_evi_all.RData")
plot(bv_range_evi_all[[1]])

names(bv_range_evi_all) <- bv_dates

bv_evi_all_stack <- stack(bv_range_evi_all)

# filter to july dates, average within
bv_jul_evi_all <- raster::subset(bv_evi_all_stack,
                                   grep('\\.07\\.', 
                                        names(bv_evi_all_stack), 
                                        value = T))

bv_jul_year_groups <- as_tibble(bv_dates) %>% 
  mutate(year = year(ymd(value)),
         month = month(ymd(value))) %>% 
  filter(month == 7) %>% 
  group_by(year) %>% 
  mutate(group = cur_group_id())

bv_jul_evi <- stackApply(bv_jul_evi_all, 
                           indices = bv_jul_year_groups$group,
                           fun = mean)
plot(bv_jul_evi$index_1)
plot(st_geometry(bv_cut), add = T)

# filter to Sept dates, average within
bv_sep_evi_all <- raster::subset(bv_evi_all_stack,
                                   grep('\\.09\\.', 
                                        names(bv_evi_all_stack), 
                                        value = T))

bv_sep_year_groups <- as_tibble(bv_dates) %>% 
  mutate(year = year(ymd(value)),
         month = month(ymd(value))) %>% 
  filter(month == 9) %>% 
  group_by(year) %>% 
  mutate(group = cur_group_id())

bv_sep_evi <- stackApply(bv_sep_evi_all, 
                           indices = bv_sep_year_groups$group,
                           fun = mean)
plot(bv_sep_evi$index_1)

# Delta EVI
bv_delta_evi <- bv_jul_evi - bv_sep_evi
plot(bv_delta_evi$layer.1)

### Extract Jul evi across bv cutblocks and fires
bv_cut_all <- st_read("Spatial_Layers/Productivity_comparison/Barkerville/BV_Cons_Cutblocks/VEG_CONSOLIDATED_CUT_BLOCKS_SP.gdb/") %>% 
  st_transform(crs = crs(bv_range))
# plot(st_geometry(bv_cut_all)) # looks good
# plot(st_geometry(bv_range), add = T)

bv_cut <- st_intersection(bv_cut_all, bv_range)
# 3946 cutblocks
#plot(st_geometry(bv_cut))

bv_cut_list <- split(bv_cut, seq(nrow(bv_cut)))

bv_cut_yearly_jul_evi <- lapply(bv_cut_list, bv_jul_evi, "july_evi",
                                  FUN = extract_yearly_evi_cutblocks)
bv_cut_yearly_jul_evi_df <- do.call("rbind", bv_cut_yearly_jul_evi)

ggplot(bv_cut_yearly_jul_evi_df, aes(x = year_since_cut,
                                       y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest harvest",
       y = "Mean July EVI")

# fires
bv_fire_all <- st_read("Spatial_Layers/Productivity_comparison/Barkerville/BV_Fires/PROT_HISTORICAL_FIRE_POLYS_SP/") %>% 
  st_transform(crs = crs(bv_range)) %>% 
  rename(FIRE_LABEL = FIRELABEL)
plot(st_geometry(bv_fire_all)) # looks good
plot(st_geometry(bv_range), add = T)

bv_fire <- st_intersection(bv_fire_all, bv_range)
# 77 fires
plot(st_geometry(bv_fire))

bv_fire_list <- split(bv_fire, seq(nrow(bv_fire)))

bv_fire_yearly_jul_evi <- lapply(bv_fire_list, bv_jul_evi, "july_evi",
                                   FUN = extract_yearly_evi_fires)
bv_fire_yearly_jul_evi_df <- do.call("rbind", bv_fire_yearly_jul_evi)

ggplot(bv_fire_yearly_jul_evi_df, aes(x = year_since_fire,
                                        y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since fire",
       y = "Mean July EVI")

### Extract delta evi across bv cutblocks and fires
# cutblocks
bv_cut_yearly_delta_evi <- lapply(bv_cut_list, bv_delta_evi, "delta_evi",
                                    FUN = extract_yearly_evi_cutblocks)
bv_cut_yearly_delta_evi_df <- do.call("rbind", bv_cut_yearly_delta_evi)

ggplot(bv_cut_yearly_delta_evi_df, aes(x = year_since_cut,
                                         y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest harvest",
       y = "Delta EVI")

# fires
bv_fire_yearly_delta_evi <- lapply(bv_fire_list, bv_delta_evi, "delta_evi",
                                     FUN = extract_yearly_evi_fires)
bv_fire_yearly_delta_evi_df <- do.call("rbind", bv_fire_yearly_delta_evi)

ggplot(bv_fire_yearly_delta_evi_df, aes(x = year_since_fire,
                                          y = evi)) +
  geom_point(size = 0.25, alpha = 0.15) +
  geom_smooth() +
  labs(x = "Years since forest fire",
       y = "Delta EVI")

##### Compare EVI across ranges ####

## add inset map with range locations below legend
focal_ranges <- bc_caribou_ranges %>% 
  filter(HERD_NAME %in% c("Itcha-Ilgachuz", "Barkerville",
                          "Purcells South", "Purcell Central")) %>% 
  mutate(caribou_range = if_else(HERD_NAME %in% c("Purcells South", "Purcell Central"),
                                 "Purcells",
                                 HERD_NAME))
bc <- bcmaps::bc_bound()

inset_map <- ggplot(bc) +
  geom_sf(size = 0.3) +
  geom_sf(data = focal_ranges, 
          aes(fill = caribou_range),
          size = 0.3) +
  scale_fill_viridis_d() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

### July EVI

# compare trends in Itcha vs. Purcells
ii_cut_yearly_evi_df$site <- "Itcha-Ilgachuz"
purc_cut_yearly_jul_evi_df$site <- "Purcells"
bv_cut_yearly_jul_evi_df$site <- "Barkerville"

cut_jul_evi_range_comp <- rbind(purc_cut_yearly_jul_evi_df,
                                ii_cut_yearly_evi_df,
                                bv_cut_yearly_jul_evi_df)
# save(cut_jul_evi_range_comp, 
#      file = "Spatial_Layers/Productivity_comparison/cut_jul_evi_range_comp.RData")

ggplot(cut_jul_evi_range_comp, aes(x = year_since_cut, y = evi,
                                   colour = site, group = site)) +
  geom_jitter(size = 0.05, alpha = 0.05, width = 0.2) +
  geom_smooth() +
  geom_ribbon(stat = "smooth", method = "gam", size = 0.1,
              fill = NA, color = "black") +
  scale_color_viridis_d() +
  labs(x = "Years since forest harvest",
       y = "Mean July EVI",
       colour = "")

# fire
purc_fire_yearly_jul_evi_df$site <- "Purcells"
ii_fire_yearly_evi_df$site <- "Itcha-Ilgachuz"
bv_fire_yearly_jul_evi_df$site <- "Barkerville"

fire_jul_evi_range_comp <- rbind(purc_fire_yearly_jul_evi_df,
                                 ii_fire_yearly_evi_df,
                                 bv_fire_yearly_jul_evi_df)
# save(fire_jul_evi_range_comp, 
#      file = "Spatial_Layers/Productivity_comparison/fire_jul_evi_range_comp.RData")

ggplot(fire_jul_evi_range_comp, aes(x = year_since_fire, y = evi,
                                    colour = site, group = site)) +
  geom_jitter(size = 0.1, alpha = 0.05, width = 0.2) +
  geom_smooth() +
  scale_color_viridis_d() +
  labs(x = "Years since forest fire",
       y = "Mean July EVI",
       colour = "")

cut_jul_evi_range_comp_merge <- cut_jul_evi_range_comp %>% 
  rename(dist_id = cutblock_id,
         year_since_dist = year_since_cut) %>% 
  mutate(dist_type = "Cutblock")

fire_jul_evi_range_comp_merge <- fire_jul_evi_range_comp %>% 
  rename(dist_id = fire_id,
         year_since_dist = year_since_fire) %>% 
  mutate(dist_type = "Fire")

all_jul_evi_comp <- rbind(cut_jul_evi_range_comp_merge,
                          fire_jul_evi_range_comp_merge)

jul_evi_dist_comp_plot <- ggplot(all_jul_evi_comp, aes(x = year_since_dist, y = evi,
                                                       colour = site, group = site)) +
  geom_jitter(size = 0.1, alpha = 0.05, width = 0.2) +
  geom_smooth() +
  geom_ribbon(stat = "smooth", method = "gam", size = 0.1,
              fill = NA, color = "white") +
  scale_colour_viridis_d() +
  labs(x = "Years since disturbance",
       y = "Mean July EVI",
       colour = "Caribou range") +
  theme(legend.position = c(1.2, 0.8),
        plot.margin = margin(0.5, 5, 0.5, 0.5, "cm")) +
  facet_wrap(~dist_type, scales = "free_x")
jul_evi_dist_comp_plot

jul_evi_dist_comp_plot_inset <- ggdraw() +
  draw_plot(jul_evi_dist_comp_plot, x = 0, y = 0) +
  draw_plot(inset_map, x = 0.66, y = 0.15, width = 0.45, height = 0.45)
jul_evi_dist_comp_plot_inset

ggsave("Outputs/plots/jul_evi_dist_comp_plot2.jpg", 
       width = 8, height = 4)

### Delta EVI
# cutblocks
ii_cut_yearly_delta_evi_df$site <- "Itcha-Ilgachuz"
purc_cut_yearly_delta_evi_df$site <- "Purcells"
bv_cut_yearly_delta_evi_df$site <- "Barkerville"

cut_delta_evi_range_comp <- rbind(purc_cut_yearly_delta_evi_df,
                                  ii_cut_yearly_delta_evi_df,
                                  bv_cut_yearly_delta_evi_df)

# fire
ii_fire_yearly_delta_evi_df$site <- "Itcha-Ilgachuz"
purc_fire_yearly_delta_evi_df$site <- "Purcells"
bv_fire_yearly_delta_evi_df$site <- "Barkerville"

fire_delta_evi_range_comp <- rbind(purc_fire_yearly_delta_evi_df,
                                   ii_fire_yearly_delta_evi_df,
                                   bv_fire_yearly_delta_evi_df)

# combine and plot
cut_delta_evi_range_comp_merge <- cut_delta_evi_range_comp %>% 
  rename(dist_id = cutblock_id,
         year_since_dist = year_since_cut) %>% 
  mutate(dist_type = "Cutblock")

fire_delta_evi_range_comp_merge <- fire_delta_evi_range_comp %>% 
  rename(dist_id = fire_id,
         year_since_dist = year_since_fire) %>% 
  mutate(dist_type = "Fire")

all_delta_evi_comp <- rbind(cut_delta_evi_range_comp_merge,
                            fire_delta_evi_range_comp_merge)

delta_evi_dist_comp_plot <- ggplot(all_delta_evi_comp, aes(x = year_since_dist, y = evi,
                                                           colour = site, group = site)) +
  geom_jitter(size = 0.1, alpha = 0.05, width = 0.2) +
  geom_smooth() +
  geom_ribbon(stat = "smooth", method = "gam", size = 0.1,
              fill = NA, color = "black") +
  scale_color_viridis_d() +
  labs(x = "Years since disturbance",
       y = "Delta EVI",
       colour = "Caribou range") +
  facet_wrap(~dist_type, scales = "free_x")
delta_evi_dist_comp_plot

ggsave("Outputs/plots/delta_evi_dist_comp_plot.jpg", 
       width = 8, height = 4)


#### Athabasca herd ranges ####
ab_caribou_ranges <- st_read("Spatial_Layers/Productivity_comparison/AB_CaribouRange/Caribou_Range.shp")
head(ab_caribou_ranges)
plot(st_geometry(ab_caribou_ranges))
unique(ab_caribou_ranges$LOCALRANGE)

atha_range <- filter(ab_caribou_ranges, 
                     LOCALRANGE %in% c("East Side Athabasca",
                                       "West Side Athabasca"))
plot(st_geometry(atha_range))

atha_centre <- atha_range %>% 
  st_union() %>% 
  st_transform(crs(sta.wgs)) %>% 
  st_centroid(atha_range) %>% 
  as_vector()

atha_utm_crs <- lonlat2UTM(atha_centre)

atha_range <- st_transform(atha_range, crs = atha_utm_crs)
plot(st_geometry(atha_range))

atha_bbox <- atha_range %>% 
  st_buffer(1000)

atha_bbox <- st_as_sfc(st_bbox(atha_bbox))
crs(atha_range)

plot(atha_bbox)
plot(st_geometry(atha_range), add = T)
st_write(atha_bbox, dsn = "Spatial_Layers/Productivity_comparison/Athabasca/atha_bbox.shp")

atha_aoi_centre <- st_union(atha_range) 
atha_aoi_centre <- st_geometry(st_centroid(atha_aoi_centre)) %>% 
  st_transform(crs = crs(sta.wgs)) %>% 
  unlist()

atha_dates <- mt_dates(product = "MOD13Q1", lat = atha_aoi_centre[2],
                     lon = atha_aoi_centre[1]) %>% 
  mutate(month = month(ymd(calendar_date))) %>% 
  filter(month %in% c(7, 9)) %>% 
  pull(calendar_date)

#save(atha_aoi_centre, atha_range, file = "Spatial_Layers/Productivity_comparison/Athabasca/atha_aoi_and_range.RData")
# plot(atha_range)
# atha_range_sf <- st_as_sf(atha_range)

# write function for extracting evi
get_evi_atha_range <- function(dates){
  evi_right <- mt_subset(product = "MOD13Q1",
                             lat = atha_aoi_centre[2] + 0.3, 
                             lon = atha_aoi_centre[1] + 1.1, # xmin
                             band = "250m_16_days_EVI",
                             start = dates,
                             end = dates,
                             km_lr = 100, 
                             km_ab = 100,
                             internal = TRUE)
  evi_right_raster <- mt_to_raster(evi_right, 
                                  reproject = F)
  
  evi_left <- mt_subset(product = "MOD13Q1",
                        lat = atha_aoi_centre[2] -0.7, 
                        lon = atha_aoi_centre[1] - 0.6,
                        band = "250m_16_days_EVI",
                        start = dates,
                        end = dates,
                        km_lr = 100, 
                        km_ab = 100,
                        internal = TRUE)
  evi_left_raster <- mt_to_raster(evi_left, 
                                   reproject = F)
  
  evi_middle <-  mt_subset(product = "MOD13Q1",
                       lat = atha_aoi_centre[2], 
                       lon = atha_aoi_centre[1] - 0.1, # xmin
                       band = "250m_16_days_EVI",
                       start = dates,
                       end = dates,
                       km_lr = 70, 
                       km_ab = 100,
                       internal = TRUE)
  evi_middle_raster <- mt_to_raster(evi_middle, 
                                    reproject = F)
  
  evi_merged <- raster::merge(evi_right_raster, 
                              evi_left_raster,
                              evi_middle_raster)
  
  evi_utm <- projectRaster(evi_merged, 
                           crs = crs(atha_range))
  
  atha_evi_mask <- raster::mask(evi_utm, atha_range)
}

atha_range_evi_all <- lapply(atha_dates, FUN = get_evi_atha_range)
plot(st_geometry(atha_range))
plot(atha_range_evi_all[[1]], add = T)

<<<<<<< HEAD
load("Spatial_Layers/Productivity_comparison/Athabasca/atha_range_evi_all.RData")
plot(atha_range_evi_all[[1]])

names(atha_range_evi_all) <- atha_dates[1:91]

atha_evi_all_stack <- stack(atha_range_evi_all)

# filter to july dates, average within
atha_jul_evi_all <- raster::subset(atha_evi_all_stack,
                                 grep('\\.07\\.', 
                                      names(atha_evi_all_stack), 
                                      value = T))

atha_jul_year_groups <- as_tibble(atha_dates) %>% 
  mutate(year = year(ymd(value)),
         month = month(ymd(value))) %>% 
  filter(month == 7) %>% 
  group_by(year) %>% 
  mutate(group = cur_group_id())

atha_jul_evi <- stackApply(atha_jul_evi_all, 
                         indices = atha_jul_year_groups$group,
                         fun = mean)
plot(atha_jul_evi$index_1)
plot(st_geometry(atha_range), add = T)

# filter to Sept dates, average within
atha_sep_evi_all <- raster::subset(atha_evi_all_stack,
                                 grep('\\.09\\.', 
                                      names(atha_evi_all_stack), 
                                      value = T))

atha_sep_year_groups <- as_tibble(atha_dates[1:91]) %>% 
  mutate(year = year(ymd(value)),
         month = month(ymd(value))) %>% 
  filter(month == 9) %>% 
  group_by(year) %>% 
  mutate(group = cur_group_id())

atha_sep_evi <- stackApply(atha_sep_evi_all, 
                         indices = atha_sep_year_groups$group,
                         fun = mean)
plot(atha_sep_evi$index_1)

# Delta EVI
atha_delta_evi <- atha_jul_evi - atha_sep_evi
plot(atha_delta_evi$layer.1)

### NTEMS Cutblock Data
ntems_cut <- raster("Spatial_Layers/NTEMS_Forest_Harvest_1985-2020/CA_Forest_Harvest_1985-2020.tif")
plot(ntems_cut)
crs(ntems_cut)

atha_bbox_nad83 <- st_transform(atha_bbox, crs = crs(ntems_cut)) %>% 
  st_as_sf()
atha_cut <- crop(ntems_cut, extent(atha_bbox_nad83))
plot(atha_cut)

atha_cut <- projectRaster(atha_cut, crs = crs(atha_range)); beep(8)
plot(atha_cut)
plot(st_geometry(atha_range), add = T)

atha_cut_mask <- mask(atha_cut, atha_range)
plot(atha_cut_mask)

# could vectorize these cutblock polygons and then use old function or try to remake it using rasters for evi and cutblocks...
atha_cut_vec <- rasterToPolygons(atha_cut_mask, na.rm = T, dissolve = F); beep(8)
# started 11:23 am... didn't finish by end of day

### Try with terra package instead, seems to be quicker
ntems_cut <- rast("Spatial_Layers/NTEMS_Forest_Harvest_1985-2020/CA_Forest_Harvest_1985-2020.tif")
plot(ntems_cut)
crs(ntems_cut)
# test <- clamp(ntems_cut, lower = 1, values = F)
# hist(test)

atha_bbox_nad83 <- vect(atha_bbox_nad83)
atha_range_sv <- vect(atha_range)

atha_cut_crop <- crop(ntems_cut, atha_bbox_nad83)
plot(atha_cut_crop)

atha_cut_crop <- terra::project(atha_cut_crop, crs(atha_range),
                                method = "near") # year data is categorigal, not continuous

atha_cut_mask <- mask(atha_cut_crop, atha_range_sv)
plot(atha_cut_mask)

hist(atha_cut_mask$`CA_Forest_Harvest_1985-2020`)

# remove 0s
atha_cut_mask <- clamp(atha_cut_mask, lower = 1, upper = 3000,
                       values = F)
plot(atha_cut_mask)
plot(atha_range_sv, add = T)

hist(atha_cut_mask$`CA_Forest_Harvest_1985-2020`)

# turn into polygons:
atha_cut_poly <- as.polygons(atha_cut_mask, values = T,
                             dissolve = T, na.rm = T)
atha_cut_poly <- atha_cut_poly %>% 
  rename(year_cut = `CA_Forest_Harvest_1985-2020`)
plot(atha_cut_poly)
plot(atha_range_sv, add = T) # looks good
# one polygon feature for each year... disaggregate

atha_cut_poly <- disagg(atha_cut_poly)
atha_cut_poly #20051 cutblock polygons... 
plot(atha_cut_poly[3000])/1e6

# exclude "cutblocks" that are only one grid cell - likely roads/not actually cutblocks?
atha_cut_poly$area <- expanse(atha_cut_poly)
mean(atha_cut_poly$area) # ~3 ha

atha_cut_poly <- atha_cut_poly %>% 
  dplyr::filter(area > 1000) %>% 
  mutate(cutblock_id = row_number())
atha_cut_poly # 11143 cutblocks
plot(atha_cut_poly)
sum(atha_cut_poly$area)/1e6 # 601 km2

# atha_cut_list <- split(atha_cut_poly, "area")
# atha_cut_list

# convert to spatraster
atha_jul_evi_sr <- rast(atha_jul_evi)
atha_jul_evi_sr

# extract evi across cutblock polygons:
atha_cut_yearly_evi <- terra::extract(atha_jul_evi_sr, atha_cut_poly,
                                 fun = mean, method = "bilinear")

# rename evi "index x" to years, covert to long format, add year since cut info
colnames(atha_cut_yearly_evi)[2:24] <- c(2000:2022)

atha_cut_yearly_evi <- atha_cut_yearly_evi %>% 
  rename(cutblock_id = ID) %>% 
  pivot_longer(cols = -cutblock_id, names_to = "year_evi",
               values_to = "evi") %>% 
  left_join(as_tibble(atha_cut_poly), by = "cutblock_id") %>% 
  mutate(year_since_cut = as.numeric(year_evi) - year_cut,
         type = "jul_evi")

hist(atha_cut_yearly_evi$year_since_cut)

ggplot(atha_cut_yearly_evi, aes(x = year_since_cut,
                                y = jul_evi)) +
  #geom_point(size = 0.1, alpha = 0.05) +
  geom_smooth()

# should weight cutblock influence by size in the model somehow...

save(atha_range_evi_all, file = "atha_range_evi_all.RData")


##### Snake-Sahtaneh #####
snake <- vect("GCPB_CARIBOU_POPULATION_SP.gdb/") %>% 
  filter(HERD_NAME == "Snake-Sahtaneh")
plot(snake)

snake_centroid <- snake %>% 
  terra::project(y = crs(sta.wgs)) %>% 
  terra::centroids() %>%
  crds()
snake_centroid[[1]]

snake_utm_crs <- lonlat2UTM(c(snake_centroid[[1]], snake_centroid[[2]]))
# UTM zone 10 (same as itcha)

snake_range <- terra::project(x = snake, y = "epsg:32610")
plot(snake_range)

snake_dates <- mt_dates(product = "MOD13Q1", lat = snake_centroid[[2]],
                       lon = snake_centroid[[1]]) %>% 
  mutate(month = month(ymd(calendar_date))) %>% 
  filter(month %in% c(7, 9)) %>% 
  pull(calendar_date)

save(snake_centroid, snake_range, file = "Spatial_Layers/Productivity_comparison/Snake_Sahtaneh/snake_aoi_and_range.RData")
# plot(atha_range)
# atha_range_sf <- st_as_sf(atha_range)

# write function for extracting evi
get_evi_snake_range <- function(dates){
  evi_top <- mt_subset(product = "MOD13Q1",
                         lat = snake_centroid[[2]]+ 0.3, 
                         lon = snake_centroid[[1]], # xmin
                         band = "250m_16_days_EVI",
                         start = dates,
                         end = dates,
                         km_lr = 100, 
                         km_ab = 50,
                         internal = TRUE)
  evi_top_raster <- mt_to_terra(evi_top, 
                                   reproject = F)
  # evi_top_utm <- terra::project(evi_top_raster,
  #                        "epsg:32610")
  # 
  # terra::plot(evi_top_utm)
  # terra::plot(snake_range, add = T)
  
  evi_bottom <- mt_subset(product = "MOD13Q1",
                        lat = snake_centroid[[2]]-0.4, # can probably go to like -0.6
                        lon = snake_centroid[[1]], # xmin
                        band = "250m_16_days_EVI",
                        start = dates,
                        end = dates,
                        km_lr = 100,
                        km_ab = 60,
                        internal = TRUE)
  evi_bottom_raster <- mt_to_terra(evi_bottom,
                                  reproject = F)
  # evi_bottom_utm <- terra::project(evi_bottom_raster,
  #                                  "epsg:32610")
  # 
  # plot(snake_range)
  # plot(evi_top_utm, add = T)
  # plot(evi_bottom_utm, add = T)
  # plot(snake_range, add = T)
 
  evi_merged <- terra::merge(evi_top_raster,
                             evi_bottom_raster)

  evi_utm <- terra::project(evi_merged,
                     "epsg:32610")

  snake_evi_mask <- terra::mask(evi_utm, snake_range)
  # plot(snake_evi_mask)
  # plot(snake_range, add = T)
}

snake_range_evi_all <- lapply(snake_dates, FUN = get_evi_snake_range)
# started 12:30 pm July 24, finished by 9:30am July 28

save(snake_range_evi_all, file = "snake_range_evi_all.RData")



##### Modelling #####
load("Spatial_Layers/Productivity_comparison/cut_jul_evi_range_comp.RData")
load("Spatial_Layers/Productivity_comparison/fire_jul_evi_range_comp.RData")

#cutblocks

# ii_yearly_cut_model <- lmer(evi ~ as.character(year_since_cut) + (1 | cutblock_id),
#                             data = ii_cut_yearly_evi_df)

# try all ranges together with interactions:
#filter to years that all ranges have data for:
cut_evi_model_data <- cut_jul_evi_range_comp %>% 
  filter(!is.na(evi), year_since_cut <= 50,
         year_since_cut >= -20)

# run all ranges together adding fixed effect for range?
yearly_cut_model <- lmer(evi ~ as.character(year_since_cut)*site + (1 | cutblock_id),
                         data = cut_evi_model_data)
summary(yearly_cut_model)

plot(yearly_cut_model) # looks cloud-y

cut_evi_pred_int <- predictInterval(yearly_cut_model, 
                                    newdata = cut_evi_model_data,
                                    which = "fixed") # all data is too big

cut_evi_model_data %>% 
  mutate(pred_evi = predict(yearly_cut_model)) %>% 
  group_by(site, year_since_cut) %>% 
  summarise(pred_evi = mean(pred_evi)) %>% 
  ggplot(aes(x = year_since_cut, y = pred_evi, 
             group = site, colour = site)) +
  geom_point() +
  geom_line()

# predict yearly evi values

# from Clayton Lamb's github:
boot.dat <- data.frame()
for(i in 1:500){
  dat <- cut_evi_model_data %>%
    group_by(cutblock_id) %>%
    sample_frac(1, replace = TRUE) #pulls the same number, but with replacement so thats the bootstrapping part
  
  # I don't think you need to run the model to bootstrap, just need to resample...
  # mod <- dat %>%
  #   lmer(evi ~ as.character(year_since_cut)*site + (1 | cutblock_id),
  #        data = .)
  # 
  # dat$pred <-predict(mod)
  
  boot.dat <- dat %>%
    group_by(site, year_since_cut) %>%
    summarise(mean_evi = mean(evi, na.rm=TRUE)) %>% #why is this evi, not pred??
    mutate(iter=i) %>%
    rbind(boot.dat)
  
  print(i)
}

cut_boot_ests <- boot.dat %>%
  group_by(site, year_since_cut) %>%
  summarise(mean=mean(mean_evi, na.rm=TRUE),
            lower=quantile(mean_evi, 0.05),
            upper=quantile(mean_evi,0.95)) %>% 
  mutate(dist_type = "Cutblock") %>% 
  rename(year_since_dist = year_since_cut)

cut_boot_ests %>%
  ggplot(aes(x=year_since_dist, y=mean, 
             group = site, colour = site)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width=0.01, alpha=0.2) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d() +
  labs(x = "Time since cut (years)", y = "Mean EVI",
       colour = "")

## fire
range(fire_jul_evi_range_comp$year_since_fire)
fire_evi_model_data <- fire_jul_evi_range_comp %>% 
  filter(!is.na(evi), year_since_fire <= 90,
         year_since_fire >= -18)

yearly_fire_model <- lmer(evi ~ as.character(year_since_fire)*site 
                          + (1 | fire_id),
                         data = fire_evi_model_data)
summary(yearly_fire_model)

plot(yearly_fire_model) # looks cloud-y

# fire_evi_pred_int <- predictInterval(yearly_fire_model, 
#                                     newdata = fire_evi_model_data,
#                                     which = "fixed") # all data is too big

fire_evi_model_data %>% 
  mutate(pred_evi = predict(yearly_fire_model)) %>% 
  group_by(site, year_since_fire) %>% 
  summarise(pred_evi = mean(pred_evi)) %>% 
  ggplot(aes(x = year_since_fire, y = pred_evi, 
             group = site, colour = site)) +
  geom_point() +
  geom_line()

# predict yearly evi values

# from Clayton Lamb's github:
boot.dat_fire <- data.frame()
for(i in 1:500){
  dat <- fire_evi_model_data %>%
    group_by(fire_id) %>%
    sample_frac(1, replace = TRUE) #pulls the same number, but with replacement so thats the bootstrapping part
  
  # I don't think you need to run the model to bootstrap, just need to resample...
  # mod <- dat %>%
  #   lmer(evi ~ as.character(year_since_fire)*site + (1 | fireblock_id),
  #        data = .)
  # 
  # dat$pred <-predict(mod)
  
  boot.dat_fire <- dat %>%
    group_by(site, year_since_fire) %>%
    summarise(mean_evi = mean(evi, na.rm=TRUE)) %>% #why is this evi, not pred??
    mutate(iter=i) %>%
    rbind(boot.dat_fire)
  
  print(i)
}

fire_boot_ests <- boot.dat_fire %>%
  group_by(site, year_since_fire) %>%
  summarise(mean=mean(mean_evi, na.rm=TRUE),
            lower=quantile(mean_evi, 0.05),
            upper=quantile(mean_evi,0.95)) %>% 
  mutate(dist_type = "Fire") %>% 
  rename(year_since_dist = year_since_fire)

fire_boot_ests %>%
  ggplot(aes(x=year_since_dist, y=mean, 
             group = site, colour = site)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width=0.01, alpha=0.4) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d() +
  labs(x = "Time since fire (years)", y = "Mean EVI",
       colour = "")

# plot together
dist_boot_ests <- rbind(cut_boot_ests, fire_boot_ests)

jul_evi_dist_boot_plot <- ggplot(data = dist_boot_ests,
       aes(x=year_since_dist, y=mean, 
           group = site, colour = site)) +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                width=0.01, alpha=0.4) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d() +
  labs(x = "Time since disturbance (years)", y = "Mean EVI",
       colour = "") +
  theme(legend.position = c(1.15, 0.8),
        plot.margin = margin(0.5, 5, 0.5, 0.5, "cm")) +
  facet_wrap(~dist_type, nrow = 1, scales = "free_x")

jul_evi_dist_boot_comp_plot_inset <- ggdraw() +
  draw_plot(jul_evi_dist_boot_plot, x = 0, y = 0) +
  draw_plot(inset_map, x = 0.68, y = 0.1, width = 0.42, height = 0.42)
jul_evi_dist_boot_comp_plot_inset

ggsave("Outputs/plots/jul_evi_dist_boot_comp_plot.jpg", 
       width = 9, height = 4.5)
