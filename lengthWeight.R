rm(list = ls())
###############
# Script Info #
###############
# PURPOSE: Weight-length relationships from DATRAS
# AUTHOR: Scott Large 2015
# REVIEWED/EXTENDED BY: Axel G. Rossberg (Cefas) 2015
#
############
# PACKAGES #
############
#
library(devtools)
install_github("ices-dk/rICES")
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
na.false <- function(x) {return(replace(x, which(is.na(x)), FALSE))}
#
# HL <- getDATRAS(record = "HL", 
#                 survey = "NS-IBTS",
#                 startyear = 1965,
#                 endyear = 2014,
#                 quarters = c(1:4),
#                 parallel = T,
#                 cores = 4)

# save(HL, file = "02072015_NS-IBTS_HL.rdata")
dir()
load("02072015_NS-IBTS_HL.rdata")

###################################
# Clean up raw data from getHLfun #
###################################
#
# Remove extra columns
if(any(colnames(HL) == "V1")) HL[, ":=" (V1 = NULL)]
#
# Remove extra spaces from character strings
cnames <- colnames(HL)
for(cname in cnames) {
  set(HL, j = cname, value = gsub("[[:space:]]", "", HL[[cname]]))
}
#
# Change appropriate to numeric
numCols <- c("Quarter", "SweepLngt", "HaulNo", "Year", "SpecCode", "SpecVal",
             "TotalNo", "CatIdentifier", "NoMeas", "SubFactor", "SubWgt",
             "CatCatchWgt", "LngtCode", "LngtClass", "HLNoAtLngt",
             "DateofCalculation","Valid_Aphia")    
HL[, (numCols) := lapply(.SD, as.numeric), .SDcols = numCols]
#
# Change -9 values to NA
for(cname in cnames){
  set(HL, i = which(HL[[cname]] == -9), j = cname, value = NA)
}
#
# # Select only important columns
# keepers <- c("StNo","Quarter", "StatRec", "Depth", "HaulDur", "Distance",
#              "GroundSpeed", "Ship", "Year", "WingSpread", "HaulLat",
#              "HaulLong", "ShootLat", "ShootLong")
# others <- colnames(hauls)[!colnames(hauls) %in% keepers]
# hauls[, c(others):= NULL]
#
#Convert from nautical miles/hour to meters/minute


