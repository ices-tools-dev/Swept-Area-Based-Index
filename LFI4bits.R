rm(list = ls())
####################
# LFI for BITS #
####################
# Purpose: To download data from DATRAS, extract oultiers, clean data, 
# calculate length-weight parameters,
# and best LFI time series for BITS survey 
# Authors: Scott Large, Colin Millar and Adriana Villamor
# Date: Februray 2018

#~~~~~~~~~~~~~~~#
# Load packages #
#~~~~~~~~~~~~~~~#

library(icesDatras)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(icesVocab)
library(reshape2)
library(lubridate)
library(gamm4)

#~~~~~~~~~~~~~~~#
# Download data #
#~~~~~~~~~~~~~~~#


hh_bits <- getDATRAS(record = "HH", "BITS", years = 1991:2017, quarters = 1)

hl_bits <- getDATRAS(record = "HL", "BITS", years = 1991:2017, quarters = 1)

ca_bits<-getDATRAS(record="CA",survey =  "BITS", years = 1991:2017, quarters = 1)

speclist <- getCodeList("SpecWoRMS")


hl_bits <- left_join(hl_bits, 
                     speclist %>% 
                       select(Key, 
                              Description),
                     by = c("Valid_Aphia" = "Key")) %>%
  rename(Species = Description)



ca_bits <- left_join(ca_bits, 
                     speclist %>% 
                       select(Key, 
                              Description),
                     by = c("Valid_Aphia" = "Key")) %>%
  rename(Species = Description)

##Define species list includded in the LFI calculation
LFIspecies<- c("Gadus morhua", 
               "Platichthys flesus",
               "Pleuronectes platessa",
               "Scophthalmus maximus",
               "Merlangius merlangus")

ca_bits <- ca_bits%>%
  filter(Species %in% LFIspecies)
hl_bits <- hl_bits%>%
  filter(Species %in% LFIspecies)

# Transform LngtClass with LngtCode "." and "0" to cm!

ca_bits<-rbind(ca_bits%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               ca_bits%>%filter(!LngtCode%in%c(".", "0")))

hl_bits<-rbind(hl_bits%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               hl_bits%>%filter(!LngtCode%in%c(".", "0")))


#need to remove common columns, "RecordType" and "DateofCalculation" before merging

ca_bits <- ca_bits %>% select(-(RecordType)) %>% select(-(DateofCalculation))  
hh_bits <- hh_bits %>% select(-(RecordType)) %>% select(-(DateofCalculation))

bits <- left_join(hl_bits, ca_bits)

#Only valid and Day hauls
bits <- left_join(bits, hh_bits)%>%
  filter( DayNight =="D",HaulVal =="V")


bits[bits == -9] <- NA
bits[bits == 0] <- NA

#ca_bits[ca_bits == NaN] <- NA  
bits[bits == Inf] <- NA

#save bitsq4 for speed
save(bits, file= "bitsq1.RData")
load("bitsq1.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plotting data to check and extract outliers #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#adding an id column to track back outliers in the regresions

bits$id <- 1:nrow(bits) 

#select data over the top 85% of the range of values for check,

#Outliers1: Check extreme Individual Weights
bits%>% ggplot(aes(IndWgt,LngtClass) )+geom_point(aes(colour= Country))+
  facet_wrap(~Species,scales = "free")+theme(text = element_text(size=8))

ggsave("IndWgtLngtClassQ1.tiff", units="in", width=5, height=5, dpi=300)

a <- bits %>% filter(Species == "Gadus morhua")%>%
  filter(!is.na(IndWgt))
b <- bits %>%filter(Species == "Merlangius merlangus")%>%
  filter(!is.na(IndWgt))
c <- bits %>%filter(Species == "Platichthys flesus")%>%
  filter(!is.na(IndWgt)) 
d <- bits %>% filter(Species == "Pleuronectes platessa")%>%
  filter(!is.na(IndWgt)) 
e <- bits %>% filter(Species == "Scophthalmus maximus")%>%
  filter(!is.na(IndWgt))


Outl1a <- a %>% filter( IndWgt > max(a$IndWgt, na.rm=TRUE)*0.85)
Outl1b <- b %>% filter( IndWgt > max(b$IndWgt, na.rm=TRUE)*0.85)
Outl1c <- c %>% filter( IndWgt > max(c$IndWgt, na.rm=TRUE)*0.85)
Outl1d <- d %>% filter( IndWgt > max(d$IndWgt, na.rm=TRUE)*0.85)
Outl1e <- e %>% filter( IndWgt > max(e$IndWgt, na.rm=TRUE)*0.85)
Outl1 <- rbind(Outl1a, Outl1b, Outl1c, Outl1d, Outl1e)
Outl1$outlier <- "IndWgt"

#Check points out of the regression

  # fit the model
  gam1 <- gam(log(IndWgt) ~ log(LngtClass),
              data = a,
              drop.unused.levels = FALSE)
  
    
  unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 2.5))
  gam.check(gam1)
  
  a$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 2.5
  
  a$threshold <- ifelse(a$LngtClass < 20, 5, 1) * 2.5
  a$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > a$threshold
  

  a$id[a$outlier == TRUE]  
  
  lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = a)
  ggsave("CodlogFitQ1.tiff", units="in", width=5, height=5, dpi=300)
  lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = a)
  ggsave("CodFitQ1.tiff", units="in", width=5, height=5, dpi=300)
  
 
  # fit the model
  gam1 <- gam(log(IndWgt) ~ log(LngtClass),
              data = b,
              drop.unused.levels = FALSE)
  
  
  unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 2.5))
  gam.check(gam1)

  b$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 2.5
  
  b$threshold <- ifelse(b$LngtClass < 25, 5, 1) * 2.5
  b$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > b$threshold
  
  
  b$id[b$outlier == TRUE]  
  
  lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = b)
  ggsave("HakelogFitQ1.tiff", units="in", width=5, height=5, dpi=300)
  lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = b)
  ggsave("HakeFitQ1.tiff", units="in", width=5, height=5, dpi=300)
  
  # fit the model
  gam1 <- gam(log(IndWgt) ~ log(LngtClass),
              data = c,
              drop.unused.levels = FALSE)
  
  
  unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 5))
  gam.check(gam1)
  

  c$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 5
  
  #c$threshold <- ifelse(c$LngtClass < 20, 4, 1) * 4.5
  #c$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > c$threshold
  
  
  c$id[c$outlier == TRUE]  
  
  lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = c)
  ggsave("FlologFitQ1.tiff", units="in", width=5, height=5, dpi=300)
  lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = c)
  ggsave("FloFitQ1.tiff", units="in", width=5, height=5, dpi=300)
  
  
  # fit the model
  gam1 <- gam(log(IndWgt) ~ log(LngtClass),
              data = d,
              drop.unused.levels = FALSE)
  
  
  unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 5))
  gam.check(gam1)
  
  
  d$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 5
  
  #d$threshold <- ifelse(d$LngtClass < 30, 5, 1) * 4.5
  #d$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > d$threshold
  
  
  d$id[d$outlier == TRUE]  
  
  lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = d)
  ggsave("PlalogFitQ1.tiff", units="in", width=5, height=5, dpi=300)
  lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = d)
  ggsave("PlaFitQ1.tiff", units="in", width=5, height=5, dpi=300)
  
  # fit the model
  gam1 <- gam(log(IndWgt) ~ log(LngtClass),
              data = e,
              drop.unused.levels = FALSE)
  
  
  unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 3.5))
  gam.check(gam1)
  
  
  e$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 3.5
  
  e$threshold <- ifelse(e$LngtClass < 20, 4, 1) * 3.5
  e$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > e$threshold
  
  
  e$id[e$outlier == TRUE]  
  
  lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = e)
  ggsave("TurlogFitQ1.tiff", units="in", width=5, height=5, dpi=300)
  lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = e)
  ggsave("TurFitQ1.tiff", units="in", width=5, height=5, dpi=300)
  
#finish exporting outliers for IndWgt
  Outl2a <- a %>% filter(outlier ==TRUE)%>%select(-c(outlier, threshold))
  Outl2b <- b %>% filter(outlier ==TRUE)%>%select(-c(outlier, threshold))
  Outl2c <- c %>% filter(outlier ==TRUE)%>%select(-c(outlier, threshold))
  Outl2d <- d %>% filter(outlier ==TRUE)%>%select(-c(outlier, threshold))
  Outl2e <- e %>% filter(outlier ==TRUE)%>%select(-c(outlier, threshold))
  Outl2 <- rbind(Outl2a, Outl2b, Outl2c, Outl2d, Outl2e)
  Outl2$outlier <- "IndWgt"
  
a <- a %>% select(-c(outlier, threshold))
b <- b %>% select(-c(outlier, threshold))
c <- c %>% select(-c(outlier))
d <- d %>% select(-c(outlier))
e <- e %>% select(-c(outlier, threshold))

#Extreme length classes by year
  
bits%>% ggplot(aes(LngtClass, Year) )+geom_point(aes(colour= Country))+
  facet_wrap(~Species,scales = "free") + theme(text = element_text(size=8))

ggsave("LngtClassYearQ1.tiff", units="in", width=5, height=5, dpi=300)


Outl3a <- a %>% filter( LngtClass > max(a$LngtClass, na.rm=TRUE)*0.95)
Outl3b <- b %>% filter( LngtClass > max(b$LngtClass, na.rm=TRUE)*0.95)
Outl3c <- c %>% filter( LngtClass > max(c$LngtClass, na.rm=TRUE)*0.95)
Outl3d <- d %>% filter( LngtClass > max(d$LngtClass, na.rm=TRUE)*0.95)
Outl3e <- e %>% filter( LngtClass > max(e$LngtClass, na.rm=TRUE)*0.95)
Outl3 <- rbind(Outl3a, Outl3b, Outl3c, Outl3d, Outl3e)
Outl3$outlier <- "LngtClass"

#Extreme number of individuals at length
bits%>% ggplot(aes(LngtClass, HLNoAtLngt) )+geom_point(aes(colour= Country))+ 
  facet_wrap(~Species,scales = "free") + theme(text = element_text(size=8))

ggsave("LngtClassHLNoATLngtQ1.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

Outl4a <- a %>% filter( HLNoAtLngt > max(a$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl4b <- b %>% filter( HLNoAtLngt > max(b$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl4c <- c %>% filter( HLNoAtLngt > max(c$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl4d <- d %>% filter( HLNoAtLngt > max(d$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl4e <- e %>% filter( HLNoAtLngt > max(e$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl4 <- rbind(Outl4a, Outl4b, Outl4c, Outl4d, Outl4e)
Outl4$outlier <- "HLNoAtLngt"

#DoorSpread and WingSpread consistency and extreme values 

bits%>% ggplot(aes(DoorSpread,HaulNo) )+geom_point(aes(colour= Country))+
  facet_wrap(~Country,scales = "free")

ggsave("DoorSpreadHaulQ1.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

Outl5a <- bits %>% filter(DoorSpread > max(bits$DoorSpread, na.rm=TRUE)*0.85)
Outl5a <- Outl5a[!duplicated(Outl5a[c("Country", "Year", "Gear", "HaulNo")]),]

#weird series of doorspread = 50 in Pol and also in RUS
#How could I make this filter more general?
Outl5b <- bits %>% filter(Country ==c("POL", "RUS", "EST"), DoorSpread == c(50,68,70))
Outl5b <- Outl5b[!duplicated(Outl5b[c("Country", "Year", "Gear", "HaulNo")]),]
Outl5 <- rbind(Outl5a, Outl5b)
Outl5$outlier <- "DoorSpread"

bits%>% ggplot(aes(WingSpread,HaulNo) )+geom_point(aes(colour= Country))+
  facet_wrap(~Country,scales = "free")
ggsave("WingSpreadHaulQ1.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

#Something here for wingspread, poland and estonia, wrong
Outl6 <- bits %>% filter(Country ==c("POL", "RUS", "EST"), WingSpread > 0)
Outl6 <- Outl6[!duplicated(Outl6[c("Country", "Year", "Gear", "HaulNo")]),]
Outl6$outlier <- "WingSpread"


#Distance
#there are still some zeros that have not become NAs here, no idea whytf!!
bits%>% ggplot(aes(Distance,HaulNo) )+geom_point(aes(colour= Country))+
  facet_wrap(~Country,scales = "free")

ggsave("DistanceQ1.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

##I need to do df for each country, an extract 0.85 upper, and those less than 50 meters

Outl7a <- bits %>% filter(Country == "DEN")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7a <- Outl7a[!duplicated(Outl7a[c("Year", "Gear", "HaulNo")]),]

Outl7b <- bits %>% filter(Country == "GFR")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7b <- Outl7b[!duplicated(Outl7b[c("Year", "Gear", "HaulNo")]),]

Outl7c <- bits %>% filter(Country == "LAT")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7c <- Outl7c[!duplicated(Outl7c[c("Year", "Gear", "HaulNo")]),]

Outl7d <- bits %>% filter(Country == "LTU")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7d <- Outl7d[!duplicated(Outl7d[c("Year", "Gear", "HaulNo")]),]

Outl7e <- bits %>% filter(Country == "POL")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7e <- Outl7e[!duplicated(Outl7e[c("Year", "Gear", "HaulNo")]),]

Outl7f <- bits %>% filter(Country == "RUS")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7f <- Outl7f[!duplicated(Outl7f[c("Year", "Gear", "HaulNo")]),]

Outl7g <- bits %>% filter(Country == "SWE")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7g <- Outl7g[!duplicated(Outl7g[c("Year", "Gear", "HaulNo")]),]

Outl7h <- bits %>% filter(Country == "EST")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7h <- Outl7h[!duplicated(Outl7h[c("Year", "Gear", "HaulNo")]),]

Outl7 <- rbind(Outl7a, Outl7b, Outl7c, Outl7d, Outl7e, Outl7f, Outl7g, Outl7h)
Outl7$outlier <- "Distance"

# When Distance is NA, we can calculate Distance2 with the coordinates given

bits$HaulLat[bits$HaulLat == -9] <- NA
bits$HaulLong[bits$HaulLong == -9] <- NA
bits$ShootLat[bits$ShootLat == -9] <- NA
bits$ShootLong[bits$ShootLong == -9] <- NA

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

#Distance2 is the distance as calculated with the submitted coordinates
bits$Distance2 <- 1000*earth_distance(bits$HaulLong, bits$HaulLat,bits$ShootLong, bits$ShootLat)

bits%>% ggplot(aes(Distance2,HaulNo) )+geom_point(aes(colour= Country))+
                                facet_wrap(~Country, scales = "free")
ggsave("Distance2Q1.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')


#Check coordinates resulting in hauls longer than 10 km Distance2 
#and less than 50m

Outl8a <- bits %>% filter(Distance2 >10000)
Outl8a <- Outl8a[!duplicated(Outl8a[c("Country","Year", "Gear", "HaulNo")]),]
Outl8b <- bits %>% filter(Distance2 < 50)
Outl8b <- Outl8b[!duplicated(Outl8b[c("Country","Year", "Gear", "HaulNo")]),]
Outl8 <- rbind(Outl8a, Outl8b)%>%select(-(Distance2))
Outl8$outlier <- "Coordinates"

#relation between reported and calculated distance, less data

bits%>% ggplot(aes(Distance,Distance2) )+geom_point(aes(colour= Country))+
  facet_wrap(~Country, scales = "free")
ggsave("DistancevsCoordinatesQ1.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

#here, I should do some extraction as the one of IndWgt, later.

#GroundSpeed
bits%>% ggplot(aes(GroundSpeed, HaulNo))+geom_point(aes(colour= Country))+
  facet_wrap(~Country, scales = "free")
ggsave("GroundSpeedHaulQ1.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

#Extract GroundSpeed less than 2.0 and more than 7
Outl9a <- bits %>% filter(GroundSpeed >7)
Outl9a <- Outl8a[!duplicated(Outl8a[c("Country","Year", "Gear", "HaulNo")]),]
Outl9b <- bits %>% filter(GroundSpeed < 2)
Outl9b <- Outl8b[!duplicated(Outl8b[c("Country","Year", "Gear", "HaulNo")]),]
Outl9 <- rbind(Outl9a, Outl9b)%>% select(-(Distance2))
Outl9$outlier <- "GroundSpeed"

bits <- bits %>% select(-Distance2)

bits%>% ggplot(aes(Distance,(HaulDur*GroundSpeed)))+
  geom_point(aes(colour= Country))+facet_wrap(~Country)

#Haul Duration
bits%>% ggplot(aes(HaulDur, HaulNo))+
  geom_point(aes(color=Country))+facet_wrap(~Country)

ggsave("HaulDurQ1.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

Outl10 <- bits %>% filter(HaulDur > max(HaulDur, na.rm=TRUE)*0.75)
Outl10 <- Outl10[!duplicated(Outl9[c("Country","Year", "Gear", "HaulNo")]),]
Outl10$outlier <- "HaulDur" 

#merge all outliers and export excel

OutlQ1 <- rbind(Outl1,Outl2,Outl3,Outl4,Outl5,Outl6,Outl7,Outl8,Outl9,Outl10 )

###################
#Repeat everything with quarter 4 changing name of saved plots
###################

hh_bits <- getDATRAS(record = "HH", "BITS", years = 1991:2017, quarters = 4)

hl_bits <- getDATRAS(record = "HL", "BITS", years = 1991:2017, quarters = 4)

ca_bits<-getDATRAS(record="CA",survey =  "BITS", years = 1991:2017, quarters = 4)

speclist <- getCodeList("SpecWoRMS")


hl_bits <- left_join(hl_bits, 
                     speclist %>% 
                       select(Key, 
                              Description),
                     by = c("Valid_Aphia" = "Key")) %>%
  rename(Species = Description)



ca_bits <- left_join(ca_bits, 
                     speclist %>% 
                       select(Key, 
                              Description),
                     by = c("Valid_Aphia" = "Key")) %>%
  rename(Species = Description)

##Define species list includded in the LFI calculation
LFIspecies<- c("Gadus morhua", 
               "Platichthys flesus",
               "Pleuronectes platessa",
               "Scophthalmus maximus",
               "Merlangius merlangus")

ca_bits <- ca_bits%>%
  filter(Species %in% LFIspecies)
hl_bits <- hl_bits%>%
  filter(Species %in% LFIspecies)

# Transform LngtClass with LngtCode "." and "0" to cm!

ca_bits<-rbind(ca_bits%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               ca_bits%>%filter(!LngtCode%in%c(".", "0")))

hl_bits<-rbind(hl_bits%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               hl_bits%>%filter(!LngtCode%in%c(".", "0")))


#need to remove common columns, "RecordType" and "DateofCalculation" before merging

ca_bits <- ca_bits %>% select(-(RecordType)) %>% select(-(DateofCalculation))  
hh_bits <- hh_bits %>% select(-(RecordType)) %>% select(-(DateofCalculation))

bits <- left_join(hl_bits, ca_bits)

#Only valid and Day hauls
bits <- left_join(bits, hh_bits)%>%
  filter( DayNight =="D",HaulVal =="V")


bits[bits == -9] <- NA
bits[bits == 0] <- NA

#ca_bits[ca_bits == NaN] <- NA  
bits[bits == Inf] <- NA

#save bitsq4 for speed
#save(bits, file= "bitsq4.RData")
#load("bitsq4.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plotting data to check and extract outliers #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#adding an id column to track back outliers in the regresions

bits$id <- 1:nrow(bits) 

#select data over the top 85% of the range of values for check,

#Outliers1: Check extreme Individual Weights
bits%>% ggplot(aes(IndWgt,LngtClass) )+geom_point(aes(colour= Country))+
  facet_wrap(~Species,scales = "free")+theme(text = element_text(size=8))

ggsave("IndWgtLngtClassQ4.tiff", units="in", width=5, height=5, dpi=300)

a <- bits %>% filter(Species == "Gadus morhua")%>%
  filter(!is.na(IndWgt))
b <- bits %>%filter(Species == "Merlangius merlangus")%>%
  filter(!is.na(IndWgt))
c <- bits %>%filter(Species == "Platichthys flesus")%>%
  filter(!is.na(IndWgt)) 
d <- bits %>% filter(Species == "Pleuronectes platessa")%>%
  filter(!is.na(IndWgt)) 
e <- bits %>% filter(Species == "Scophthalmus maximus")%>%
  filter(!is.na(IndWgt))


Outl1a <- a %>% filter( IndWgt > max(a$IndWgt, na.rm=TRUE)*0.85)
Outl1b <- b %>% filter( IndWgt > max(b$IndWgt, na.rm=TRUE)*0.85)
Outl1c <- c %>% filter( IndWgt > max(c$IndWgt, na.rm=TRUE)*0.85)
Outl1d <- d %>% filter( IndWgt > max(d$IndWgt, na.rm=TRUE)*0.85)
Outl1e <- e %>% filter( IndWgt > max(e$IndWgt, na.rm=TRUE)*0.85)
Outl1 <- rbind(Outl1a, Outl1b, Outl1c, Outl1d, Outl1e)
Outl1$outlier <- "IndWgt"

#Check points out of the regression

# fit the model
gam1 <- gam(log(IndWgt) ~ log(LngtClass),
            data = a,
            drop.unused.levels = FALSE)


unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 2.5))
gam.check(gam1)

a$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 2.5

a$threshold <- ifelse(a$LngtClass < 20, 5, 1) * 2.5
a$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > a$threshold


a$id[a$outlier == TRUE]  

lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = a)
ggsave("CodlogFitQ4.tiff", units="in", width=5, height=5, dpi=300)
lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = a)
ggsave("CodFitQ4.tiff", units="in", width=5, height=5, dpi=300)


# fit the model
gam1 <- gam(log(IndWgt) ~ log(LngtClass),
            data = b,
            drop.unused.levels = FALSE)


unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 2.5))
gam.check(gam1)

b$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 2.5

b$threshold <- ifelse(b$LngtClass < 25, 5, 1) * 2.5
b$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > b$threshold


b$id[b$outlier == TRUE]  

lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = b)
ggsave("HakelogFitQ4.tiff", units="in", width=5, height=5, dpi=300)
lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = b)
ggsave("HakeFitQ4.tiff", units="in", width=5, height=5, dpi=300)

# fit the model
gam1 <- gam(log(IndWgt) ~ log(LngtClass),
            data = c,
            drop.unused.levels = FALSE)


unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 5))
gam.check(gam1)


c$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 5

#c$threshold <- ifelse(c$LngtClass < 20, 4, 1) * 4.5
#c$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > c$threshold


c$id[c$outlier == TRUE]  

lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = c)
ggsave("FlologFitQ4.tiff", units="in", width=5, height=5, dpi=300)
lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = c)
ggsave("FloFitQ4.tiff", units="in", width=5, height=5, dpi=300)


# fit the model
gam1 <- gam(log(IndWgt) ~ log(LngtClass),
            data = d,
            drop.unused.levels = FALSE)


unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 5))
gam.check(gam1)


d$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 5

d$threshold <- ifelse(d$LngtClass < 30, 6, 1) *5
d$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > d$threshold


d$id[d$outlier == TRUE]  

lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = d)
ggsave("PlalogFitQ4.tiff", units="in", width=5, height=5, dpi=300)
lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = d)
ggsave("PlaFitQ4.tiff", units="in", width=5, height=5, dpi=300)

# fit the model
gam1 <- gam(log(IndWgt) ~ log(LngtClass),
            data = e,
            drop.unused.levels = FALSE)


unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 3.5))
gam.check(gam1)


e$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 3.5

e$threshold <- ifelse(e$LngtClass < 20, 4, 1) * 3.5
e$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > e$threshold


e$id[e$outlier == TRUE]  

lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = e)
ggsave("TurlogFitQ4.tiff", units="in", width=5, height=5, dpi=300)
lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = e)
ggsave("TurFitQ4.tiff", units="in", width=5, height=5, dpi=300)

#finish exporting outliers for IndWgt
Outl2a <- a %>% filter(outlier ==TRUE)%>%select(-c(outlier, threshold))
Outl2b <- b %>% filter(outlier ==TRUE)%>%select(-c(outlier, threshold))
Outl2c <- c %>% filter(outlier ==TRUE)%>%select(-c(outlier))
Outl2d <- d %>% filter(outlier ==TRUE)%>%select(-c(outlier, threshold))
Outl2e <- e %>% filter(outlier ==TRUE)%>%select(-c(outlier, threshold))
Outl2 <- rbind(Outl2a, Outl2b, Outl2c, Outl2d, Outl2e)
Outl2$outlier <- "IndWgt"

a <- a %>% select(-c(outlier, threshold))
b <- b %>% select(-c(outlier, threshold))
c <- c %>% select(-c(outlier))
d <- d %>% select(-c(outlier, threshold))
e <- e %>% select(-c(outlier, threshold))

#Extreme length classes by year

bits%>% ggplot(aes(LngtClass, Year) )+geom_point(aes(colour= Country))+
  facet_wrap(~Species,scales = "free") + theme(text = element_text(size=8))

ggsave("LngtClassYearQ4.tiff", units="in", width=5, height=5, dpi=300)


Outl3a <- a %>% filter( LngtClass > max(a$LngtClass, na.rm=TRUE)*0.95)
Outl3b <- b %>% filter( LngtClass > max(b$LngtClass, na.rm=TRUE)*0.95)
Outl3c <- c %>% filter( LngtClass > max(c$LngtClass, na.rm=TRUE)*0.95)
Outl3d <- d %>% filter( LngtClass > max(d$LngtClass, na.rm=TRUE)*0.95)
Outl3e <- e %>% filter( LngtClass > max(e$LngtClass, na.rm=TRUE)*0.95)
Outl3 <- rbind(Outl3a, Outl3b, Outl3c, Outl3d, Outl3e)
Outl3$outlier <- "LngtClass"

#Extreme number of individuals at length
bits%>% ggplot(aes(LngtClass, HLNoAtLngt) )+geom_point(aes(colour= Country))+ 
  facet_wrap(~Species,scales = "free") + theme(text = element_text(size=8))

ggsave("LngtClassHLNoATLngtQ4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

Outl4a <- a %>% filter( HLNoAtLngt > max(a$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl4b <- b %>% filter( HLNoAtLngt > max(b$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl4c <- c %>% filter( HLNoAtLngt > max(c$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl4d <- d %>% filter( HLNoAtLngt > max(d$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl4e <- e %>% filter( HLNoAtLngt > max(e$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl4 <- rbind(Outl4a, Outl4b, Outl4c, Outl4d, Outl4e)
Outl4$outlier <- "HLNoAtLngt"

#DoorSpread and WingSpread consistency and extreme values 

bits%>% ggplot(aes(DoorSpread,HaulNo) )+geom_point(aes(colour= Country))+
  facet_wrap(~Country,scales = "free")

ggsave("DoorSpreadHaulQ4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

Outl5a <- bits %>% filter(DoorSpread > max(bits$DoorSpread, na.rm=TRUE)*0.85)
Outl5a <- Outl5a[!duplicated(Outl5a[c("Country", "Year", "Gear", "HaulNo")]),]

#weird series of doorspread = 50 in Pol and also in RUS
#How could I make this filter more general?
Outl5b <- bits %>% filter(Country ==c("POL", "RUS", "EST"), DoorSpread == c(50,68,70))
Outl5b <- Outl5b[!duplicated(Outl5b[c("Country", "Year", "Gear", "HaulNo")]),]
Outl5 <- rbind(Outl5a, Outl5b)
Outl5$outlier <- "DoorSpread"

bits%>% ggplot(aes(WingSpread,HaulNo) )+geom_point(aes(colour= Country))+
  facet_wrap(~Country,scales = "free")
ggsave("WingSpreadHaulQ4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

#Something here for wingspread, poland and estonia, wrong
Outl6 <- bits %>% filter(Country ==c("POL", "RUS", "EST"), WingSpread > 0)
Outl6 <- Outl6[!duplicated(Outl6[c("Country", "Year", "Gear", "HaulNo")]),]
Outl6$outlier <- "WingSpread"


#Distance
#there are still some zeros that have not become NAs here, no idea whytf!!
bits%>% ggplot(aes(Distance,HaulNo) )+geom_point(aes(colour= Country))+
  facet_wrap(~Country,scales = "free")

ggsave("DistanceQ4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

##I need to do df for each country, an extract 0.85 upper, and those less than 50 meters

Outl7a <- bits %>% filter(Country == "DEN")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7a <- Outl7a[!duplicated(Outl7a[c("Year", "Gear", "HaulNo")]),]

Outl7b <- bits %>% filter(Country == "GFR")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7b <- Outl7b[!duplicated(Outl7b[c("Year", "Gear", "HaulNo")]),]

Outl7c <- bits %>% filter(Country == "LAT")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7c <- Outl7c[!duplicated(Outl7c[c("Year", "Gear", "HaulNo")]),]

Outl7d <- bits %>% filter(Country == "LTU")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7d <- Outl7d[!duplicated(Outl7d[c("Year", "Gear", "HaulNo")]),]

Outl7e <- bits %>% filter(Country == "POL")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7e <- Outl7e[!duplicated(Outl7e[c("Year", "Gear", "HaulNo")]),]

Outl7f <- bits %>% filter(Country == "RUS")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7f <- Outl7f[!duplicated(Outl7f[c("Year", "Gear", "HaulNo")]),]

Outl7g <- bits %>% filter(Country == "SWE")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7g <- Outl7g[!duplicated(Outl7g[c("Year", "Gear", "HaulNo")]),]

Outl7h <- bits %>% filter(Country == "EST")%>%
  filter(Distance > max(Distance, na.rm=TRUE)*0.95)
Outl7h <- Outl7h[!duplicated(Outl7h[c("Year", "Gear", "HaulNo")]),]

Outl7 <- rbind(Outl7a, Outl7b, Outl7c, Outl7d, Outl7e, Outl7f, Outl7g, Outl7h)
Outl7$outlier <- "Distance"

# When Distance is NA, we can calculate Distance2 with the coordinates given

bits$HaulLat[bits$HaulLat == -9] <- NA
bits$HaulLong[bits$HaulLong == -9] <- NA
bits$ShootLat[bits$ShootLat == -9] <- NA
bits$ShootLong[bits$ShootLong == -9] <- NA

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

#Distance2 is the distance as calculated with the submitted coordinates
bits$Distance2 <- 1000*earth_distance(bits$HaulLong, bits$HaulLat,bits$ShootLong, bits$ShootLat)

bits%>% ggplot(aes(Distance2,HaulNo) )+geom_point(aes(colour= Country))+
  facet_wrap(~Country, scales = "free")
ggsave("Distance2Q4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')


#Check coordinates resulting in hauls longer than 10 km Distance2 
#and less than 50m

Outl8a <- bits %>% filter(Distance2 >10000)
Outl8a <- Outl8a[!duplicated(Outl8a[c("Country","Year", "Gear", "HaulNo")]),]
Outl8b <- bits %>% filter(Distance2 < 50)
Outl8b <- Outl8b[!duplicated(Outl8b[c("Country","Year", "Gear", "HaulNo")]),]
Outl8 <- rbind(Outl8a, Outl8b)%>%select(-(Distance2))
Outl8$outlier <- "Coordinates"

#relation between reported and calculated distance, less data

bits%>% ggplot(aes(Distance,Distance2) )+geom_point(aes(colour= Country))+
  facet_wrap(~Country, scales = "free")
ggsave("DistancevsCoordinatesQ4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

#here, I should do some extraction as the one of IndWgt, later.

#GroundSpeed
bits%>% ggplot(aes(GroundSpeed, HaulNo))+geom_point(aes(colour= Country))+
  facet_wrap(~Country, scales = "free")
ggsave("GroundSpeedHaulQ4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

#Extract GroundSpeed less than 2.0 and more than 7
Outl9a <- bits %>% filter(GroundSpeed >7)
Outl9a <- Outl8a[!duplicated(Outl8a[c("Country","Year", "Gear", "HaulNo")]),]
Outl9b <- bits %>% filter(GroundSpeed < 2)
Outl9b <- Outl8b[!duplicated(Outl8b[c("Country","Year", "Gear", "HaulNo")]),]
Outl9 <- rbind(Outl9a, Outl9b)%>% select(-(Distance2))
Outl9$outlier <- "GroundSpeed"

bits <- bits %>% select(-Distance2)

bits%>% ggplot(aes(Distance,(HaulDur*GroundSpeed)))+
  geom_point(aes(colour= Country))+facet_wrap(~Country)

#Haul Duration
bits%>% ggplot(aes(HaulDur, HaulNo))+
  geom_point(aes(color=Country))+facet_wrap(~Country)

ggsave("HaulDurQ4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

#Outl10 <- bits %>% filter(HaulDur > max(HaulDur, na.rm=TRUE)*0.75)
#Outl10 <- Outl10[!duplicated(Outl9[c("Country","Year", "Gear", "HaulNo")]),]
#Outl10$outlier <- "HaulDur" 

#merge all outliers and export excel

OutlQ4 <- rbind(Outl1,Outl2,Outl3,Outl4,Outl5,Outl6,Outl7,Outl8,Outl9 )


#merge outliers from both quarters
BITSoutliers <- rbind(OutlQ1,OutlQ4)
write.csv(BITSoutliers, file = "BITSoutliers.csv")
BITSoutliers_Denmark <- BITSoutliers %>% filter(Country == "DEN")
write.csv(BITSoutliers_Denmark, file = "BITSoutliers_Denmark.csv")
BITSoutliers_Estonia <- BITSoutliers %>% filter(Country == "EST")
write.csv(BITSoutliers_Estonia, file = "BITSoutliers_Estonia.csv")
BITSoutliers_Germany <- BITSoutliers %>% filter(Country == "GFR")
write.csv(BITSoutliers_Germany, file = "BITSoutliers_Germany.csv")
BITSoutliers_Latvia <- BITSoutliers %>% filter(Country == "LAT")
write.csv(BITSoutliers_Latvia, file = "BITSoutliers_Latvia.csv")
BITSoutliers_Lithuania <- BITSoutliers %>% filter(Country == "LTU")
write.csv(BITSoutliers_Lithuania, file = "BITSoutliers_Lithuania.csv")
BITSoutliers_Poland <- BITSoutliers %>% filter(Country == "POL")
write.csv(BITSoutliers_Poland, file = "BITSoutliers_Poland.csv")
BITSoutliers_Russia <- BITSoutliers %>% filter(Country == "RUS")
write.csv(BITSoutliers_Russia, file = "BITSoutliers_Russia.csv")
BITSoutliers_Sweden <- BITSoutliers %>% filter(Country == "SWE")
write.csv(BITSoutliers_Sweden, file = "BITSoutliers_Sweden.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fitting model for length-weight relationship#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#If I use Year, Species, Country and Gear, I will only fit ca 50000 missing values
#If I fit by Species and year, I fit 120000 missing values
#but still only cod goes further down than 2007.
#Fitting only by species gives back the whole time series for all 5 species, but
#Following Colin's idea:
#fitting a model that allows the slope and power (Wt = aL^b) to change with year, 
#in a smooth way with yearly randomness about the smooth trend, by species
#old lm 
#bitslm <- bits %>% 
#  group_by(Species)%>%
#  filter(!is.na(IndWgt), IndWgt > 0) %>%
#  do(lm = lm(log(IndWgt) ~ log(LngtClass), data = ., 
#             singular.ok = T, 
#             na.action = na.exclude))


##Get the coefficients
#bitslm <- bitslm %>% 
#  tidy(lm)
#bitslm <- bitslm %>% select(term,
#                                  estimate) %>% 
#  spread(term, estimate)%>%
#  rename(intercept = `(Intercept)`, 
#         slope = `log(LngtClass)`)


plot(bits$IndWgt, bits$LngtClass)
plot(log(bits$IndWgt), log(bits$LngtClass))

sum(is.na(bits$IndWgt)) #320466

cod <-bits%>% 
  filter(Species=="Gadus morhua")%>%
  select(Year, Species, IndWgt, LngtClass)
sum(is.na(cod$IndWgt))

# cod is the data.frame with Individual weight, length and year.
cod <- within(cod, {
  fYear_a <- fYear_b <- factor(Year, levels = sort(unique(cod$Year)))
  Year_a <- Year_b <- Year
  #date = lubridate::ymd(paste(Year, Month, Day))
  #yday = yday(date)
})

# fit the model
gam1<- gam(log(IndWgt) ~ 1 + s(fYear_a, bs = "re") + s(Year_a, k = 13) +
            log(LngtClass) + s(fYear_b, by = log(LngtClass), bs = "re") +
            s(Year_b, k = 13, by = log(LngtClass)),
          data = cod, select = TRUE, family = gaussian(),
          drop.unused.levels = FALSE)

#year_k is the maximum degrees of freedom for the year smoothers,
#I usually set this to half of the number of years that you have data.
#You need to make a dataframe that has lots of copies of the year column 
#because you are fitting lots of smoothers on year

#Plotting the Model
plot(gam1) 
#se stands for standard error Bands

summary(gam1)

cod2 <- cod
cod2$IndWgt <-NA

pred <- predict(gam1, cod2)

library(memisc)

cod2 <- to.data.frame(pred, as.vars = 0,name="IndWgt2")
cod <- merge(cod, cod2, all = TRUE)


# Calculate distance with the coordinates
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

bits$Distance2 <- 1000*(earth_distance(bits$ShootLong, bits$ShootLat, bits$HaulLong, bits$HaulLat))

#subsitute NAs in Distance for this calculation

bits <- transform(bits, Distance = ifelse(!is.na(Distance), Distance, Distance2))

sum(is.na(bits$Distance))
sum(is.na(bits$ShootLat))
sum(is.na(bits$Distance2))   


p <- ggplot(bits, aes(Distance, Year))
q<- p + geom_point() 
q

bits$Distance[bits$Distance < 50] <- NA

p <- ggplot(bits, aes(Distance, Year))
q<- p + geom_point() 
q

bits$Distance[bits$Distance < 10] <- NA

p <- ggplot(bits, aes(DoorSpread, Year))
q<- p + geom_point() 
q

bits$DoorSpread[bits$DoorSpread < 10] <- NA


bits%>% ggplot(aes(DoorSpread,Ship) )+geom_point()
#is everything ok here or the lower values are too low?
bits%>% ggplot(aes(DoorSpread,Distance) )+geom_point()

#The amount of DoorSpread = 50 is kind of weird, to what is it due?
#Too big distances are suspicious?

bits%>% ggplot(aes(Distance,HaulNo) )+geom_point()
bits%>% ggplot(aes(Distance,CatCatchWgt) )+geom_point()


sum(is.na(bits$IndWgt))
bits <- left_join(bits, bitslm)
bits$logIndWgt <- bits$intercept+bits$slope*log(bits$LngtClass)
bits$IndWgtFit <- 10^(bits$logIndWgt)
bits <- transform(bits, IndWgt = ifelse(!is.na(IndWgt), IndWgt, IndWgtFit))

sum(is.na(bits$IndWgt))
bits%>% ggplot(aes(LngtClass,IndWgt) )+geom_point()

bits%>% ggplot(aes(Year,HLNoAtLngt) )+geom_point()+ facet_wrap(~Species,scales = "free")

# there are a few crazy numbers, but I leave them all as possible, will check with submitter PENDING
bits$WgAtLngt <- bits$IndWgt*bits$HLNoAtLngt

#Calculations of Biomass at length per km2

bits$sweptarea <- bits$DoorSpread*bits$Distance

bits$biomdens <- bits$WgAtLngt/bits$sweptarea

bits%>% ggplot(aes(Year,WgAtLngt) )+geom_point()+facet_wrap(~Species)
# I will assume that is a crazy value:
bits$WgAtLngt[bits$WgAtLngt > 40] <- NA
bits%>% ggplot(aes(Year,biomdens) )+geom_point()+facet_wrap(~Species)


#Check Duration and distance of hauls:

bits%>% ggplot(aes(Distance,HaulDur) )+geom_point()+facet_wrap(~Year)
#some weird points here and there, but I will keep them by now

bits%>% ggplot(aes(Year,biomdens) )+geom_point()+facet_wrap(~Species)

#create different LFI time series for different L<sub>LF from 20 to 50 cm

  bits <- bits %>%
    mutate(length_bin = case_when(
      LngtClass > 0 &
        LngtClass <= 20 ~ 20,
      LngtClass > 20 &
        LngtClass <= 30 ~ 30,
      LngtClass > 30 &
        LngtClass <= 40 ~ 40,
      LngtClass > 40 &
        LngtClass <= 50 ~ 50,
      LngtClass >50 ~ 60,
      TRUE ~ NA_real_))

bits <-  bits[!is.na(bits$biomdens),]  

#all this calculates the different time series of LFI= (B  for L > L<sub> lf)/Btotal
# for each value of L<sub>lf from 20 to 50
# probably I could make it nicer with some loops and functions,
#but it seems to work ;-)

lfi_calc <-bits %>%
  group_by(Year, length_bin) %>%
  summarise(lfi = sum(biomdens))

 lfi_20<- lfi_calc %>%
  group_by(Year)%>%
  filter(length_bin > 20)%>%
  summarise(lfi20 = sum(lfi)) 

lfi_30<- lfi_calc %>%
  group_by(Year)%>%
  filter(length_bin > 30)%>%
  summarise(lfi30 = sum(lfi)) 
  
lfi_40<- lfi_calc %>%
  group_by(Year)%>%
  filter(length_bin > 40)%>%
  summarise(lfi40 = sum(lfi)) 
 
lfi_50<- lfi_calc %>%
  group_by(Year)%>%
  filter(length_bin > 50)%>%
  summarise(lfi50 = sum(lfi)) 
   
lfiperyear <- lfi_calc %>% group_by(Year) %>% summarise(lfitot= sum(lfi))  

lfi_20 <- left_join(lfi_20,lfiperyear)     
lfi_30 <- left_join(lfi_30,lfiperyear) 
lfi_40 <- left_join(lfi_40,lfiperyear) 
lfi_50 <- left_join(lfi_50,lfiperyear) 

lfi_20$lfi20 <- lfi_20$lfi20/lfi_20$lfitot
lfi_30$lfi30 <- lfi_30$lfi30/lfi_30$lfitot
lfi_40$lfi40 <- lfi_40$lfi40/lfi_40$lfitot
lfi_50$lfi50 <- lfi_50$lfi50/lfi_50$lfitot

lfi_timeseries <- left_join(lfi_20,lfi_30)
lfi_timeseries <- left_join(lfi_timeseries,lfi_40)
lfi_timeseries <- left_join(lfi_timeseries,lfi_50)

lfi_timeseries <- subset(lfi_timeseries, select = -lfitot) 
lfi_timeseries <- melt(lfi_timeseries, id="Year")

ggplot(data=lfi_timeseries,
       aes(x=Year, y=value, colour=variable)) +
  geom_line()

#Fitting a 5th degree polynomial finction to each LFI time series 

model20 <- lm(lfi_20$lfi20 ~ poly(lfi_20$Year,5))

summary(model20)

model30 <- lm(lfi_30$lfi30 ~ poly(lfi_30$Year,5))

summary(model30)

model40 <- lm(lfi_40$lfi40 ~ poly(lfi_40$Year,5))

summary(model40)

model50 <- lm(lfi_50$lfi50 ~ poly(lfi_50$Year,5))

summary(model50)

#model40 and model50 give same results, check!

ggplot(data=lfi_20,
       aes(x=Year, y=lfi20)) +
  geom_point()+geom_line()
ggplot(data=lfi_30,
       aes(x=Year, y=lfi30)) +
  geom_point()+geom_line()
ggplot(data=lfi_40,
       aes(x=Year, y=lfi40)) +
  geom_point()+geom_line()
ggplot(data=lfi_50,
       aes(x=Year, y=lfi50)) +
  geom_point()+geom_line()


