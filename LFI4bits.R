rm(list = ls())
####################
# LFI for BITS #
####################
# Purpose: To download data from DATRAS, clean data, calculate length-weight parameters,
# and best LFI time series for BITS survey 
# Authors: Scott Large, Szymon Smoliński and Adriana Villamor
# Date: January 2018

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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plotting data to check and extract outliers #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

bits$id <- 1:nrow(bits) 

#select data over the top 85% of the range of values for check,

#Outliers1: Check extreme Individual Weights
bits%>% ggplot(aes(IndWgt,LngtClass) )+geom_point(aes(colour= Country))+
  facet_wrap(~Species,scales = "free")

ggsave("IndWgtQ1.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

a <- bits %>% filter(Species == "Gadus morhua")%>%
  select(Year, Species, IndWgt, LngtClass, id) %>%
  filter(!is.na(IndWgt)) %>%filter(IndWgt > 0)
b <- bits %>%filter(Species == "Merlangius merlangus")%>%
  select(Year, Species, IndWgt, LngtClass, id)%>%
  filter(!is.na(IndWgt)) %>%filter(IndWgt > 0)
c <- bits %>%filter(Species == "Platichthys flesus")%>%
  select(Year, Species, IndWgt, LngtClass, id)%>%
  filter(!is.na(IndWgt)) %>%filter(IndWgt > 0)
d <- bits %>% filter(Species == "Pleuronectes platessa")%>%
  select(Year, Species, IndWgt, LngtClass, id)%>%
  filter(!is.na(IndWgt)) %>%filter(IndWgt > 0)
e <- bits %>% filter(Species == "Scophthalmus maximus")%>% 
  select(Year, Species, IndWgt, LngtClass, id)%>%
  filter(!is.na(IndWgt)) %>%filter(IndWgt > 0)


Outl1a <- a %>% filter( IndWgt > max(a$IndWgt, na.rm=TRUE)*0.85)
Outl1b <- b %>% filter( IndWgt > max(b$IndWgt, na.rm=TRUE)*0.85)
Outl1c <- c %>% filter( IndWgt > max(c$IndWgt, na.rm=TRUE)*0.85)
Outl1d <- d %>% filter( IndWgt > max(d$IndWgt, na.rm=TRUE)*0.85)
Outl1e <- e %>% filter( IndWgt > max(e$IndWgt, na.rm=TRUE)*0.85)

#save bitsq4 for speed
save(bits, file= "bitsq4.RData")
load("bitsq4.RData")

#Check points out of the regression

  a <- within(a, {
    fYear_a <- fYear_b <- factor(Year, levels = sort(unique(a$Year)))
    Year_a <- Year_b <- Year
  })
  
  # fit the model
  gam1 <- gam(log(IndWgt) ~ log(LngtClass),
              data = a,
              drop.unused.levels = FALSE)
  
    
  unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 2.5))
  gam.check(gam1)
  
  a$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 3
  
  a$threshold <- ifelse(a$LngtClass < 20, 4, 1) * 2.5
  a$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > a$threshold
  

  a$id[a$outlier == TRUE]  
  
  lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = a)
  #ggsave
  lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = a)
  #ggsave
  
  b <- within(b, {
    fYear_a <- fYear_b <- factor(Year, levels = sort(unique(b$Year)))
    Year_a <- Year_b <- Year
   })
  
  # fit the model
  gam1 <- gam(log(IndWgt) ~ log(LngtClass),
              data = b,
              drop.unused.levels = FALSE)
  
  
  unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 2.5))
  gam.check(gam1)

  b$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 3
  
  b$threshold <- ifelse(b$LngtClass < 20, 4, 1) * 2.5
  b$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > b$threshold
  
  
  b$id[a$outlier == TRUE]  
  
  lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = b)
  #ggsave
  lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = b)
  #ggsave
  
  c <- within(c, {
    fYear_a <- fYear_b <- factor(Year, levels = sort(unique(c$Year)))
    Year_a <- Year_b <- Year
  })
  
  # fit the model
  gam1 <- gam(log(IndWgt) ~ log(LngtClass),
              data = c,
              drop.unused.levels = FALSE)
  
  
  unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 2.5))
  gam.check(gam1)
  

  c$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 3
  
  c$threshold <- ifelse(c$LngtClass < 20, 4, 1) * 2.5
  c$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > c$threshold
  
  
  c$id[a$outlier == TRUE]  
  
  lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = c)
  #ggsave
  lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = c)
  #ggsave
  
  
  d <- within(d, {
    fYear_a <- fYear_b <- factor(Year, levels = sort(unique(d$Year)))
    Year_a <- Year_b <- Year
  })
  
  # fit the model
  gam1 <- gam(log(IndWgt) ~ log(LngtClass),
              data = d,
              drop.unused.levels = FALSE)
  
  
  unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 2.5))
  gam.check(gam1)
  
  
  d$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 3
  
  d$threshold <- ifelse(d$LngtClass < 20, 4, 1) * 2.5
  d$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > d$threshold
  
  
  d$id[a$outlier == TRUE]  
  
  lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = d)
  #ggsave
  lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = d)
  #ggsave
  
  e <- within(e, {
    fYear_a <- fYear_b <- factor(Year, levels = sort(unique(e$Year)))
    Year_a <- Year_b <- Year
  })
  
  # fit the model
  gam1 <- gam(log(IndWgt) ~ log(LngtClass),
              data = e,
              drop.unused.levels = FALSE)
  
  
  unname(which(abs(residuals(gam1, type = "scaled.pearson")) > 2.5))
  gam.check(gam1)
  
  
  e$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > 3
  
  e$threshold <- ifelse(e$LngtClass < 20, 4, 1) * 2.5
  e$outlier <- abs(residuals(gam1, type = "scaled.pearson")) > e$threshold
  
  
  d$id[a$outlier == TRUE]  
  
  lattice::xyplot(log(LngtClass) ~ log(IndWgt), groups = outlier, data = e)
  #ggsave
  lattice::xyplot(LngtClass ~ IndWgt, groups = outlier, data = e)
  #ggsave
  
  
  #finish exporting outliers for IndWgt
  
  


#Extreme length classes by year
  
bits%>% ggplot(aes(LngtClass, Year) )+geom_point(aes(colour= Country))+ 
    facet_wrap(~Species,scales = "free")

#ggsave

a <- bits %>% filter(Species == "Gadus morhua")
b <- bits %>%filter(Species == "Merlangius merlangus")
c <- bits %>%filter(Species == "Platichthys flesus")
d <- bits %>% filter(Species == "Pleuronectes platessa")
e <- bits %>% filter(Species == "Scophthalmus maximus")



Outl2a <- a %>% filter( LngtClass > max(a$LngtClass, na.rm=TRUE)*0.95)
Outl2b <- b %>% filter( LngtClass > max(b$LngtClass, na.rm=TRUE)*0.95)
Outl2c <- c %>% filter( LngtClass > max(c$LngtClass, na.rm=TRUE)*0.95)
Outl2d <- d %>% filter( LngtClass > max(d$LngtClass, na.rm=TRUE)*0.95)
Outl2e <- e %>% filter( LngtClass > max(e$LngtClass, na.rm=TRUE)*0.95)


#Extreme number of individuals at length
bits%>% ggplot(aes(LngtClass, HLNoAtLngt) )+geom_point(aes(colour= Country))+ 
  facet_wrap(~Species,scales = "free")

ggsave("outl2q4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

Outl2a <- a %>% filter( HLNoAtLngt > max(a$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl2b <- b %>% filter( HLNoAtLngt > max(b$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl2c <- c %>% filter( HLNoAtLngt > max(c$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl2d <- d %>% filter( HLNoAtLngt > max(d$HLNoAtLngt, na.rm=TRUE)*0.85)
Outl2e <- e %>% filter( HLNoAtLngt > max(e$HLNoAtLngt, na.rm=TRUE)*0.85)


#DoorSpread and WingSpread consistency and extreme values 

bits%>% ggplot(aes(DoorSpread,HaulNo) )+geom_point(aes(colour= Country))
#ggsave

ggsave("outl9q4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

Outl4a <- bits %>% filter(DoorSpread > max(bits$DoorSpread, na.rm=TRUE)*0.85)
Outl4a <- Outl4a[!duplicated(Outl4a[c("Country", "Year", "Gear", "HaulNo")]),]

#weird series of doorspread = 50 in Pol and also in RUS
#How could I make this filter more general?
Outl4b <- bits %>% filter(Country =="POL", DoorSpread == 50)
Outl4b <- Outl4b[!duplicated(Outl4b[c("Country", "Year", "Gear", "HaulNo")]),]
Outl4c <- bits %>% filter(Country =="RUS", DoorSpread == c(68,70))
Outl4c <- Outl4c[!duplicated(Outl4c[c("Country", "Year", "Gear", "HaulNo")]),]

bits%>% ggplot(aes(WingSpread,HaulNo) )+geom_point(aes(colour= Country))
ggsave("outl6q4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

#Something here for wingspread, poland and estonia, wrong
Outl4d <- bits %>% filter(Country ==c("POL", "RUS", "EST"), WingSpread > 0)
Outl4d <- Outl4d[!duplicated(Outl4d[c("Country", "Year", "Gear", "HaulNo")]),]


#Distance
#there are still some zeros that have not become NAs here, no idea whytf!!
sum(is.na(bits$Distance)) 
bits%>% ggplot(aes(Distance,HaulNo) )+geom_point()+facet_wrap(~Country)
ggsave("outl5q4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

##I need to do df for each country, an extract 0.85 upper, and those less than 50 meters

Outl5a <- bits %>% group_by("Country")%>%
  filter(Distance > max(bits$Distance, na.rm=TRUE)*0.85)
Outl5a <- Outl5a[!duplicated(Outl5a[c("Country", "Year", "Gear", "HaulNo")]),]

# Calculate distance with the coordinates given

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

bits%>% ggplot(aes(Distance2,HaulNo) )+geom_point()+facet_wrap(~Country)
#ggsave

#outliers all hauls more than 10 km distance2 and those less than 50m

#relation between reported and calculated distance, less data

bits%>% ggplot(aes(Distance,Distance2) )+geom_point()+facet_wrap(~Country)


#GroundSpeed
bits%>% ggplot(aes(GroundSpeed, HaulNo))+geom_point()+facet_wrap(~Country)
#Extract GroundSpeed less than 2.0
ggsave("outl7q4.tiff", units="in", width=5, height=4, dpi=300, compression = 'lzw')

bits%>% ggplot(aes(Distance,(HaulDur*GroundSpeed)))+geom_point()+facet_wrap(~Country)

#Haul Duration
bits%>% ggplot(aes(HaulDur, HaulNo))+geom_point()+facet_wrap(~Country)
bits%>% ggplot(aes(HaulDur, Distance))+geom_point()+facet_wrap(~Country)
#I think this outliers in Denmark are the same as before


#merge all outliers and export excel

outliers1 <- rbind(Outl1a, Outl1b, Outl1c, Outl1d, Outl1e, Outl2a, Outl2b, Outl2c,
                  Outl2d, Outl2e, Outl3a, Outl3b, Outl3c, Outl3d, Outl3e, Outl4a,
                  Outl4b, Outl4c, Outl4d, Outl5a)

#Repeat everything with quarter 4 changing name of saved plots


#merge outliers from both quarters
outliers <- rbind(outliers1,outliers2)
write.csv(outliers, file = "BITSouliers.csv")

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


