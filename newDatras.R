# devtools::install_github("ices-tools-prod/icesDatras")
rm(list =ls())
library(icesDatras)
library(ggplot2)
library(viridis)
library(sf)
library(dplyr)

# Calculate distance in kilometers between two points
earth_distance <- function (long1, lat1, long2, lat2) {
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6371 # Mean radius (km) of the earth
  d <- R * c
  return(d)
}

# getDatrasDataOverview(surveys = "BITS")

# 
# hh_bits <- getDATRAS(record = "HH", "BITS", years = 1991:2017, quarters = c(1,4))
# saveRDS(hh_bits, file = "~/git/ices-dk/LFI/hh_bits.rds")
# hl_bits <- getDATRAS(record = "HL", "BITS", years = 1991:2017, quarters = c(1,4))
# saveRDS(hl_bits, file = "~/git/ices-dk/LFI/hl_bits.rds")
# ca_bits <- getDATRAS(record = "CA", "BITS", years = 1991:2017, quarters = c(1,4))
# saveRDS(ca_bits, file = "~/git/ices-dk/LFI/ca_bits.rds")

clean_hh <- readRDS("~/git/ices-dk/LFI/hh_bits.rds") %>%
  filter(Quarter == 1, 
         Year >= 2001,
         HaulVal == "V", 
         DayNight == "D") %>% 
  mutate_each(funs(replace(., . < 0, NA))) %>% 
  mutate(newDistance = NA, # dummy variable for later
         distanceType = NA,
    GroundSpeed = GroundSpeed * 1852 / 60, # Convert from nautical miles/hour to meters/minute
    LatLongDistance = earth_distance(long1 = ShootLong,
                                     lat1 = ShootLat,
                                     long2 = HaulLong,
                                     lat2 = HaulLat) * 1000,
    LatLongDistance = ifelse(LatLongDistance > 1 &
                             LatLongDistance/HaulDur > 50 &
                             LatLongDistance/HaulDur < 500,
                             LatLongDistance,
                             NA),
    rawDistance = Distance,
    SpeedHaulDur = GroundSpeed * HaulDur)

clean_hh_dist %>%
  select(newDistance,
         HaulDur,
         WingSpread,
         DoorSpread,
         Depth) %>% 
  summarise_each(funs(sum(is.na(.))))

# Add the missing values
clean_hh_dist <- clean_hh
clean_hh_dist$distanceType[is.na(clean_hh_dist$newDistance)] <- "rawDistance"
clean_hh_dist$newDistance[is.na(clean_hh_dist$newDistance)] <- clean_hh_dist$rawDistance[is.na(clean_hh_dist$newDistance)]
clean_hh_dist$distanceType[is.na(clean_hh_dist$newDistance)] <- "LatLongDistance"
clean_hh_dist$newDistance[is.na(clean_hh_dist$newDistance)] <- clean_hh_dist$LatLongDistance[is.na(clean_hh_dist$newDistance)]
clean_hh_dist$distanceType[is.na(clean_hh_dist$newDistance)] <- "SpeedHaulDur"
clean_hh_dist$newDistance[is.na(clean_hh_dist$newDistance)] <- clean_hh_dist$SpeedHaulDur[is.na(clean_hh_dist$newDistance)]

needSpeed <- clean_hh_dist[is.na(clean_hh_dist$newDistance),]

## Use linear models to estimate distance
needSpeedShip <- unique(needSpeed$Ship)
withSpeedShip <- unique(clean_hh_dist$Ship[!is.na(clean_hh_dist$GroundSpeed)])
#
# ID ships that have some missing GroundSpeed records
canSpeedShip <- needSpeedShip[needSpeedShip %in% withSpeedShip]
canSpeed <- clean_hh_dist[clean_hh_dist$Ship %in% canSpeedShip, ]
#
# Split the data into the two types: only a few missing and completely missing
canSpeedNA <- canSpeed[is.na(canSpeed$GroundSpeed) & is.na(canSpeed$newDistance),]
canSpeedOK <- canSpeed[!is.na(canSpeed$GroundSpeed),]

# Linear model to check speed
if(length(canSpeedShip)>1){
  cs1 <- lm(GroundSpeed ~ Ship + Year + Depth, canSpeedOK, na.action = "na.fail")
}else{
  cs1 <- lm(GroundSpeed ~ Year + Depth, canSpeedOK, na.action = "na.fail")  
}

cs2 <- lm(GroundSpeed ~ Year + Depth, canSpeedOK, na.action = "na.fail")
# 

clean_hh_dist <- bind_rows(
  clean_hh_dist %>% 
    filter(Ship %in% canSpeedShip) %>% 
    mutate(lmShipSpeedHaulDur = predict(cs1, .) * HaulDur
    ),
  clean_hh_dist %>% 
    filter(!Ship %in% canSpeedShip)
) %>%
  mutate( lmYear = predict(cs2, .) * HaulDur)

clean_hh_dist$distanceType[is.na(clean_hh_dist$newDistance)] <- "lmShipHaulSpeed"
clean_hh_dist$newDistance[is.na(clean_hh_dist$newDistance) &
                          clean_hh_dist$Ship %in% canSpeedShip] <- clean_hh_dist$lmShipSpeedHaulDur[is.na(clean_hh_dist$newDistance) &
                                                                                               clean_hh_dist$Ship %in% canSpeedShip]

clean_hh_dist$distanceType[is.na(clean_hh_dist$newDistance)] <- "lmYear"
clean_hh_dist$newDistance[is.na(clean_hh_dist$newDistance)] <- clean_hh_dist$lmYear[is.na(clean_hh_dist$newDistance)]


################
## Wingspread ##
################
# clean_hh_ws <- clean_hh
# clean_hh_ws$newWingSpread <- NA
# Linear model: WingSpread ~ a * log(Depth)) + b
lmWingForm <- lm(WingSpread ~ I(log(Depth)), data = clean_hh_dist)

table(clean_hh_dist$Country, is.na(clean_hh_dist$Distance))

ggplot(clean_hh_dist, aes(DoorSpread, WingSpread, color = Country)) +
  geom_point() +
  facet_grid(~Country)

# lmWingDat <- clean_hh_ws[is.na(clean_hh_ws$WingSpread) &
#                            !is.na(clean_hh_ws$Depth),]
# lmWingDatLen <- nrow(lmWingDat)
# lmWingBoot <- matrix(bootCase(lmWingForm, function(x) predict(x, lmWingDat), B = nb), nb, lmWingDatLen)

clean_hh_dist <- bind_rows(
  clean_hh_dist %>% 
    filter(is.na(WingSpread)) %>% 
    mutate(newWingSpread = predict(lmWingForm, .),
           distanceType = paste0(distanceType, "_lmDepth")),
  clean_hh_dist %>% 
    filter(!is.na(WingSpread)) %>%
    mutate(newWingSpread = WingSpread,
           distanceType = paste0(distanceType, "_raw"))
) %>% 
  mutate(sweptarea = (newWingSpread * newDistance)/1000) # meters to square km 

# clean_hh_dist %>%
#   select(newDistance,
#          newWingSpread,
#          sweptarea,
#          HaulDur,
#          WingSpread,
#          DoorSpread,
#          Depth) %>% 
#   summarise_each(funs(sum(is.na(.))))

ggplot(clean_hh_dist, aes(x = newWingSpread,
                          y = newDistance)) +
  geom_point(aes(color = distanceType), alpha = .2)

ggplot(clean_hh_dist, aes(x = GroundSpeed,
                          y = newDistance)) +
  geom_point(aes(color = Year), alpha = 0.4)  +
  facet_wrap(~distanceType)


ggplot(clean_hh_dist, aes(x = Year, y = sweptarea)) +
  geom_jitter(alpha=I(1/4), aes(fill = distanceType,
                                color = distanceType), shape = 21, alpha = 0.4) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2")


# 
# source("https://raw.githubusercontent.com/briatte/ggcorr/master/ggcorr.R")
# 
# hh_cor <- clean_hh_dist %>% 
#   dplyr::select(LatLongDistance, 
#          rawDistance, 
#          lmShipSpeedHaulDur, 
#          lmYear, 
#          SpeedHaulDur)
#   
# ggcorr(hh_cor, label = TRUE)


devtools::install_github("tidyverse/ggplot2")
library(ggplot2)
library(sf)

ggplot(nc) + geom_sf()

ggplot(clean_hh_dist, aes(x = ShootLat, y = ShootLong)) +
  geom_sf()
  geom_histogram() +
  facet_grid(Gear ~ Year)


  ggplot(clean_hh_dist, aes(x = Distance, y = WingSpread)) +
  geom_point(stat = "identity", aes(fill = Year,
                                    color = Year), 
             shape = 21,
             alpha = 0.5) +
  scale_fill_viridis(discrete = FALSE) +
  scale_color_viridis(discrete = FALSE) +
  labs(title = "Initial swept area",
       subtitle = paste0("HaulVal = 'V'"),
       x = "Distance", y = "Wing spread", fill = "Year") +
  theme_minimal() +
  theme(legend.key = element_blank())







### Combine biological and survey code to make biomass per haul per species

### LFI per trawl
  
### Plot biomass per species per length class over time
