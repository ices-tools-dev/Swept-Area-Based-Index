rm(list = ls())
###############
# Script Info #
###############
# PURPOSE: Swept Area calculations from DATRAS
# AUTHOR: Scott Large 2015
# REVIEWED/EXTENDED BY: Axel G. Rossberg (Cefas) 2015

############
# PACKAGES #
############
#
library(devtools)
devtools::install_github("ices-dk/rICES")
library(rICES)
# library(ggplot2)
library(data.table)
library(reshape2)
library(arm)
library(car)
library(DMwR) #for lofactor(..)
#
##################
# Load functions #
##################
#
# Calculate distance in kilometers between two points
earth.dist <- function (long1, lat1, long2, lat2) {
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
#
na.false <- function(x) {return(replace(x, which(is.na(x)), FALSE))}
#
#############
# Load data #
#############
#
# If FALSE, mostly suppress CI computation
need.CI <- FALSE
#
# Number of bootstrap replicates
if(need.CI){
  nb <- 1000
}else{
  nb <- 3 ## just to make code run through
}
#
# Confidence interval range
CV <- .95
#
# Set seed for reproducable results
setSeed <- set.seed(627)
#
HH <- getDATRAS(record = "HH", 
                survey = "NS-IBTS",
                startyear = 1965,
                endyear = 2014,
                quarters = c(1:4),
                parallel = T,
                cores = 4)
#
###################################
# Clean up raw data from getHHfun #
###################################
#
# Remove extra columns
if(any(colnames(HH) == "V1")) HH[, ":=" (V1 = NULL)]
#
# Remove extra spaces from character strings
cnames <- colnames(HH)
for(cname in cnames) {
  set(HH, j = cname, value = gsub("[[:space:]]", "", HH[[cname]]))
}
#
# Change appropriate to numeric
numCols <- c("Quarter", "SweepLngt", "HaulNo", "Year", "month",
             "Day", "TimeShot", "Stratum", "HaulDur", "ShootLat",
             "ShootLong", "HaulLat", "HaulLong", "Depth", "StdSpecRecCode",
             "Netopening", "Rigging", "Tickler", "Distance", "Warplngt",
             "Warpdia", "DoorSurface", "DoorWgt", "DoorSpread", "WingSpread",
             "Buoyancy", "KiteDim", "WgtGroundRope", "TowDir", "GroundSpeed",
             "SpeedWater", "SurCurDir", "SurCurSpeed", "BotCurDir",  "BotCurSpeed",
             "WindDir", "WindSpeed",  "SwellDir",  "SwellHeight",  "SurTemp",
             "BotTemp", "SurSal", "BotSal")
HH[, (numCols) := lapply(.SD, as.numeric), .SDcols = numCols]
#
# Change -9 values to NA
for(cname in cnames){
  set(HH, i = which(HH[[cname]] == -9), j = cname, value = NA)
}
#
# Select only valid hauls during the day and select important columns
hauls <- HH[HH$HaulVal=='V' & DayNight == "D", ]
keepers <- c("StNo","Quarter", "StatRec", "Depth", "HaulDur", "Distance",
             "GroundSpeed", "Ship", "Year", "WingSpread", "HaulLat",
             "HaulLong", "ShootLat", "ShootLong")
others <- colnames(hauls)[!colnames(hauls) %in% keepers]
hauls[, c(others):= NULL]
#
#Convert from nautical miles/hour to meters/minute
hauls[, GroundSpeed := GroundSpeed * 1852 / 60]
#
###################
# DISTANCE HAULED #
###################
#
## Calculate haversine distance with shoot and haul coordinates ##
hauls[, LatLongDistance := earth.dist(long1 = ShootLong,
                                      lat1 = ShootLat,
                                      long2 = HaulLong,
                                      lat2 = HaulLat) * 1000]
#
## RAW DISTANCE ##
hauls[!is.na(Distance), c("newDist", "qualityDistance") :=
        list(Distance, "rawDistance")]
#
## HAVERSINE DISTANCE ##
# if haul and shoot coordinates are the same (i.e., equal to zero) then leave as NA
# Also ingore hauls with funny implied speed.
hauls[is.na(Distance) & !is.na(LatLongDistance) & LatLongDistance > 1 &
        LatLongDistance/HaulDur>50 & LatLongDistance/HaulDur<500 ,
      c("newDist", "qualityDistance") :=
        list(LatLongDistance, "LatLongDistance")]
#
## DURATION X SPEED ##
# HaulDur x GroundSpeed
hauls[!is.na(GroundSpeed) & is.na(newDist),
      c("newDist", "qualityDistance") :=
        list(GroundSpeed * HaulDur, "SpeedHaulDur")]
#
# Hauls that don't have raw distance, shoot/haul coordinates, or GroundSpeed can be estimated
# using linear models to predict GroundSpeed. Two different types of missing data are found:
# 1) ships that have no GroundSpeed records and we need to estimate from other ships, and
# 2) ships that have only a few missing GroundSpeed records and we can use ship as a factor in the lm
#
needSpeed <- hauls[is.na(newDist) & is.na(GroundSpeed) & !is.na(HaulDur),]

if(nrow(needSpeed)>0){
  
  needSpeedShip <- unique(needSpeed$Ship)
  withSpeedShip <- unique(hauls$Ship[!is.na(hauls$GroundSpeed)])
  #
  #
  # ID ships that have some missing GroundSpeed records
  canSpeedShip <- needSpeedShip[needSpeedShip %in% withSpeedShip]
  canSpeed <- hauls[Ship %in% canSpeedShip, ]
  #
  # Split the data into the two types: only a few missing and completely missing
  canSpeedNA <- canSpeed[is.na(GroundSpeed) & is.na(newDist),]
  canSpeedOK <- canSpeed[!is.na(GroundSpeed),]
  #
  if(length(canSpeedShip)>1){
    cs1 <- lm(GroundSpeed ~ Ship + Quarter + Year + Depth, canSpeedOK)
  }else{
    cs1 <- lm(GroundSpeed ~ Quarter + Year + Depth, canSpeedOK)  
  }
  cs1 <- step(cs1)  # simplify model based on AIC
  lmShipSpeedHaulDurLen <- nrow(canSpeedNA)
  lmShipSpeedHaulDurBoot <- matrix(bootCase(cs1, function(x) predict(x, canSpeedNA) * canSpeedNA$HaulDur,
                                            B = nb),
                                   nb,
                                   lmShipSpeedHaulDurLen)
  #
  hauls[Ship %in% canSpeedShip & is.na(GroundSpeed) & is.na(newDist),
        c("newDist", "newDistCIlower", "newDistCIupper", "qualityDistance") :=
          list(apply(lmShipSpeedHaulDurBoot, 2, median),
               apply(lmShipSpeedHaulDurBoot, 2, quantile, (1-CV)/2, na.rm = T),
               apply(lmShipSpeedHaulDurBoot, 2, quantile, (1+CV)/2, na.rm = T),
               "lmShipSpeedHaulDur")]
  ##
  stillNeedSpeed <- hauls[is.na(newDist),]
  cs2 <- lm(GroundSpeed ~ Quarter + Year, canSpeedOK)
  cs2 <- step(cs2)  # simplify model based on AIC
  lmQYearSpeedHaulDurLen <- nrow(stillNeedSpeed)
  lmQYearSpeedHaulDurBoot <- matrix(bootCase(cs2, function(x) predict(x, stillNeedSpeed) * stillNeedSpeed$HaulDur,
                                             B = nb),
                                    nb,
                                    lmQYearSpeedHaulDurLen)
  #
  hauls[is.na(newDist),
        c("newDist", "newDistCIlower", "newDistCIupper", "qualityDistance") :=
          list(apply(lmQYearSpeedHaulDurBoot, 2, median),
               apply(lmQYearSpeedHaulDurBoot, 2, quantile, (1-CV)/2, na.rm = T),
               apply(lmQYearSpeedHaulDurBoot, 2, quantile, (1+CV)/2, na.rm = T),
               "lmQYearSpeedHaulDur")]
  #
  hauls[,
        c("newDist", "newDistCIlower", "newDistCIupper") :=
          list(newDist / 1000,
               newDistCIlower / 1000,
               newDistCIupper / 1000)]
  
} # if(nrow(needSpeed)>0)

#
######################
# WING & DOOR SPREAD #
######################
#
if(any(!is.na(hauls$WingSpread))){
  # Linear model: WingSpread ~ a * log(Depth)) + b
  lmWingForm <- lm(WingSpread ~ I(log(Depth)), hauls)
  lmWingDat <- hauls[is.na(WingSpread) & !is.na(Depth),]
  lmWingDatLen <- nrow(lmWingDat)
  lmWingBoot <- matrix(bootCase(lmWingForm, function(x) predict(x, lmWingDat), B = nb), nb, lmWingDatLen)
  #
  ## RAW WINGSPREAD ##
  whichHasWingSpread <- which(!is.na(hauls$WingSpread))
  whichIsOutlier <- whichHasWingSpread[na.false(lofactor(log(hauls$WingSpread[whichHasWingSpread]),5)>4)]
  hauls[!is.na(WingSpread), c("newWingSpread", "qualityWingSpread") :=
          list(WingSpread, "rawWingSpread"),]
  hauls[whichIsOutlier,c("newWingSpread", "qualityWingSpread") := list(NA, NA)]
  hauls[qualityWingSpread=="rawWingSpread" ]
  #
  ## LINEAR MODEL WINGSPREAD ##
  hauls[is.na(WingSpread) & !is.na(Depth),
        c("newWingSpread","newWingSpreadCIlower", "newWingSpreadCIupper", "qualityWingSpread") :=
          list(apply(lmWingBoot, 2, median),
               apply(lmWingBoot, 2, quantile, (1-CV)/2, na.rm = T),
               apply(lmWingBoot, 2, quantile, (1+CV)/2, na.rm = T),
               "lmWingSpread")]
  #
  ## MEDIAN WINGSPREAD ##
  #
  rawWingMat <- hauls$WingSpread[!is.na(hauls$WingSpread)]
  #?#? The CI seems to be the CI of the median (it shrinks with sample size), but not CI for how well we know WS. 
  #?#? Same with bootstrap estimats above! 
  hauls[is.na(newWingSpread), c("newWingSpread","newWingSpreadCIlower", "newWingSpreadCIupper", "qualityWingSpread") :=
          list(median(rawWingMat),
               quantile(replicate(nb, median(sample(rawWingMat, rep = TRUE))),  (1-CV)/2),
               quantile(replicate(nb, median(sample(rawWingMat, rep = TRUE))),  (1+CV)/2),
               "medianWingSpread")]
  
}else{ ## case where no WingSpread is available:
  hauls[, c("newWingSpread","newWingSpreadCIlower", "newWingSpreadCIupper", "qualityWingSpread") :=
          list(1,
               1,
               1,
               "NoWingSpread")]
}
#?#? What units are we converting to here? 
hauls[,
      c("newWingSpread", "newWingSpreadCIlower", "newWingSpreadCIupper") :=
        list(newWingSpread / 1000,
             newWingSpreadCIupper / 1000,
             newWingSpreadCIlower / 1000)]

##############
# SWEPT AREA #
##############
#
if(need.CI){
  #?#? CI are not simply multiplicative like this...
  hauls[, c("newSweptArea", "CIupper", "CIlower", "qualityCode") := list(newWingSpread * newDist,
                                                                         newWingSpreadCIupper * newDistCIupper,
                                                                         newWingSpreadCIlower * newDistCIlower,
                                                                         paste0(qualityDistance, qualityWingSpread)),]
}else{
  hauls[, c("newSweptArea", "qualityCode") := list(newWingSpread * newDist,
                                                   paste0(qualityDistance, qualityWingSpread)),]  
}
# save(hauls, file = "~/Data Products/cleanHH_v001.rdat")
# Explore data graphically
ggplot(hauls, aes(x = Year, y = newSweptArea)) +
  geom_jitter(alpha=I(1/4), aes(color = qualityCode)) +
  theme_bw()
