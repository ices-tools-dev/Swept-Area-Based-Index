rm(list = ls())
####################
# LFI for BITS #
####################
# Purpose: To calculate swept area, length-weight parameters and best LFI time series
# for BITS survey 
# Authors: Scott Large, Szymon Smoliński and Adriana Villamor
# Date: December 2017

#~~~~~~~~~~~~~~~#
# Load packages #
#~~~~~~~~~~~~~~~#

library(icesDatras)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
#devtools::install_github("ices-tools-prod/icesVocab")
library(icesVocab)
library(reshape2)


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

# Change -9 values to NA

ca_bits[ca_bits == -9] <- NA

ca_bits[ca_bits == Inf] <- NA

hl_bits[hl_bits == -9] <- NA

hl_bits[hl_bits == Inf] <- NA

#need to remove common columns, "RecordType" and "DateofCalculation"

ca_bits <- ca_bits %>% select(-(RecordType)) %>% select(-(DateofCalculation))  
hh_bits <- hh_bits %>% select(-(RecordType)) %>% select(-(DateofCalculation))

bits <- left_join(hl_bits, ca_bits)


bits <- left_join(bits, hh_bits)%>%
  filter( DayNight =="D",HaulVal =="V")

#plotting data to check outliers
sum(is.na(bits$IndWgt))
bits%>% ggplot(aes(IndWgt,LngtClass) )+geom_point()
Outl1 <- bits %>% filter(IndWgt> 7500)

write.csv(Outl1, "DATRASoutliers.csv")

sum(is.na(bits$LngtClass))
bits%>% ggplot(aes(LngtClass, Year) )+geom_point()+ facet_wrap(~Species,scales = "free")
bits%>% ggplot(aes(LngtClass, HLNoAtLngt) )+geom_point()+ facet_wrap(~Species,scales = "free")

#select the last top quarter for check, 

#bits2 <- bits[!is.na(bits$HLNoAtLngt),]
#bits2 <- bits2%>%
#  group_by(Species)
#n <- 2
#out2 <- bits2[bits2$HLNoAtLngt > quantile(bits2$HLNoAtLngt,prob=1-n/100),]

#Any of these is doing what I want!!

#out2 <- bits2 %>% group_by(Species) %>%filter(HLNoAtLngt > quantile(bits2$HLNoAtLngt, 0.98))
#out2%>% ggplot(aes(LngtClass, HLNoAtLngt) )+geom_point()+ facet_wrap(~Species,scales = "free")         

Outl <- bits %>% filter(Species == "Gadus morhua", HLNoAtLngt> 4000)
Outl <- bits %>% filter(Species == "Merlangius merlangus", HLNoAtLngt> 1000)
Outl <- bits %>% filter(Species == "Platichthys flesus", HLNoAtLngt> 4000)
Outl <- bits %>% filter(Species == "Pleuronectes platessa", HLNoAtLngt> 200)

bits%>% ggplot(aes(CatCatchWgt, Year) )+geom_point()+ facet_wrap(~Species,scales = "free")
#same, select upper third for check, PENDING
Outl <- bits %>% filter(Species == "Gadus morhua", CatCatchWgt > 4e+06)
Outl <- Outl[!duplicated(Outl[c("Country", "Year", "Gear", "HaulNo")]),] 

Outl <- bits %>% filter(Species == "Merlangius merlangus", CatCatchWgt> 5e+05)
Outl <- Outl[!duplicated(Outl[c("Country", "Year", "Gear", "HaulNo")]),]

Outl <- bits %>% filter(Species == "Platichthys flesus", CatCatchWgt> 2500000)
Outl <- Outl[!duplicated(Outl[c("Country", "Year", "Gear", "HaulNo")]),]

Outl <- bits %>% filter(Species == "Pleuronectes platessa", CatCatchWgt> 3e+05)
Outl <- Outl[!duplicated(Outl[c("Country", "Year", "Gear", "HaulNo")]),]

Outl <- bits %>% filter(Species == "Scophthalmus maximus", CatCatchWgt > 20000)
Outl <- Outl[!duplicated(Outl[c("Country", "Year", "Gear", "HaulNo")]),]

sum(is.na(bits$DoorSpread))
bits%>% ggplot(aes(DoorSpread,HaulNo) )+geom_point()
Outl2 <- bits %>% filter(DoorSpread > 300)

bits%>% ggplot(aes(Ship,CatCatchWgt) )+geom_point()+ facet_wrap(~Year,scales = "free")
bits%>% ggplot(aes(HaulNo,CatCatchWgt) )+geom_point()+ facet_wrap(~Year,scales = "free")



##To check relation between length and weight, fit the linear models
#If I use Year, Species, Country and Gear, I will only fit ca 50000 missing values
#If I fit by Species and year, I fit 120000 missing values
#but still only cod goes further down than 2007
#Species, Country and Gear improves a bit the time series 
#fitting only by species gives back the whole time series for all 5 species
#Colin?


bitslm <- bits %>%
  group_by(Species)%>%
  filter(!is.na(IndWgt), IndWgt > 0, IndWgt < 7500) %>%
  do(lm = lm(log(IndWgt) ~ log(LngtClass), data = ., 
             singular.ok = T, 
             na.action = na.exclude))

##Get the coefficients
bitslm <- bitslm %>% 
  tidy(lm)
bitslm <- bitslm %>% select(term,
                                  estimate) %>% 
  spread(term, estimate)%>%
  rename(intercept = `(Intercept)`, 
         slope = `log(LngtClass)`)

##Plot slopes distribution

bitslm%>%ggplot(aes(slope))+geom_histogram()+facet_wrap(~Species,scales = "free")
#nothing to see when fitting only by Species

# better this: dat <-  dat %>% mutate(x = replace(x, x<0, NA))

bits[bits == -9] <- NA


sum(is.na(bits$Distance)) #96189 out of 349839
bits%>% ggplot(aes(Distance,DoorSpread) )+geom_point()
#for ouliers also distance > 70000
 
bits$Distance[bits$Distance > 7000] <- NA
bits$DoorSpread[bits$DoorSpread > 300] <- NA
bits$IndWgt[bits$IndWgt > 7500] <- NA

sum(is.na(bits$LngtClass))
sum(is.na(bits$DoorSpread))
sum(is.na(bits$Netopening))

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
sum(is.na(bits$Distance2))   


p <- ggplot(bits, aes(Distance, Year))
q<- p + geom_point() 
q

bits$Distance[bits$Distance > 10000] <- NA

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


