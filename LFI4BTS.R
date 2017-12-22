rm(list = ls())
####################
# LFI for BTS#
####################
# Purpose: To calculate best LFI time series for BTS, calculating sweptarea from
#standard beamtrawl lengths as reported in
#http://ices.dk/marine-data/Documents/DATRAS%20Manuals/WGBEAM_Manual.pdf
# Authors: Adriana Villamor
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

#~~~~~~~~~~~~~~~#
# Download data #
#~~~~~~~~~~~~~~~#

hh_bts <- getDATRAS(record = "HH", "BTS", years = 1987:2017, quarters = 3)

hl_bts <- getDATRAS(record = "HL", "BTS", years = 1987:2017, quarters = 3)

ca_bts<-getDATRAS(record="CA",survey =  "BTS", years = 1987:2017, quarters = 3)

speclist <- getCodeList("SpecWoRMS")

# speclist <- rbind(speclist,
#                   read.csv("Aphias_WoRMS_add.csv"))
hl_bts <- left_join(hl_bts, 
                     speclist %>% 
                       select(Key, 
                              Description),
                     by = c("Valid_Aphia" = "Key")) %>%
  rename(Species = Description)



ca_bts <- left_join(ca_bts, 
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

ca_bts <- ca_bts%>%
  filter(Species %in% LFIspecies)
hl_bts <- hl_bts%>%
  filter(Species %in% LFIspecies)

# Transform LngtCode . and 0 to LngtClass in cm!

ca_bts<-rbind(ca_bts%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               ca_bts%>%filter(!LngtCode%in%c(".", "0")))

hl_bts<-rbind(hl_bts%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               hl_bts%>%filter(!LngtCode%in%c(".", "0")))

# Change -9 values to NA

ca_bts[ca_bits == -9] <- NA

#ca_bits[ca_bits == NaN] <- NA  
ca_bts[ca_bits == Inf] <- NA


##Fit the linear models
ca_btslm <- ca_bts %>% 
  group_by(Country, Year, Species)%>%
  filter(!is.na(IndWgt), IndWgt > 0) %>%
  do(lm = lm(log(IndWgt) ~ log(LngtClass), data = ., 
             singular.ok = T, 
             na.action = na.exclude))


##Get the coefficients
ca_btslm <- ca_btslm %>% 
  tidy(lm)
ca_btslm <- ca_btslm %>% select(term,
                                  estimate) %>% 
  spread(term, estimate)%>%
  rename(intercept = `(Intercept)`, 
         slope = `log(LngtClass)`)

##Plot slopes distribution
ca_btslm%>%ggplot(aes(slope))+geom_histogram()+facet_wrap(~Species,scales = "free")

#Only Platichthys flesus seems has more spread slopes, but that´s ok

ca_btslm%>%group_by(Year, Species)%>%summarise(medianslope=(median(slope, na.rm = T)))%>%
  ggplot(aes(Year, medianslope,color=Species))+geom_line(size=1)+scale_color_brewer(palette = "Dark2")+theme_bw()


#regression seems quite fine, merging the three pieces together
#only valid and Day hauls

hl_bts <- hl_bts[ ,-1]  
hh_bts <- hh_bts[ ,-1]        

bts <- left_join(hl_bts, ca_bts)

bts <- left_join(bts, hh_bts)
  
bts <- bts %>%filter( DayNight == "D", IndWgt > 0, HaulVal =="V")

# better this: dat <-  dat %>% mutate(x = replace(x, x<0, NA))
bts[bts == -9] <- NA

sum(is.na(bts$Distance)) #2195 out of 14300

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

bts$Distance2 <- 1000*(earth_distance(bts$ShootLong, bts$ShootLat, bts$HaulLong, bts$HaulLat))

p <- ggplot(bts, aes(Distance, Distance2))
q<- p + geom_point() 
q

#subsitute NAs in Distance for this calculation

bts <- transform(bts, Distance = ifelse(!is.na(Distance), Distance, Distance2))

sum(is.na(bts$Distance))
sum(is.na(bts$Distance2))  

plot(bts$Distance, bts$Year)
#remove two outliers

bts <- bts[!bts$Distance > 15000, ]

sum(is.na(bts$Netopening)) 

#I will add values of Beam Trawl length by country according to Annex 1 of 
#http://ices.dk/marine-data/Documents/DATRAS%20Manuals/WGBEAM_Manual.pdf

#BEL 4, FR 4, GFR 7, NED 8, ENG 4

bts <- bts %>%
  mutate(Netopening = case_when(
    Country == "BEL" ~ 4,
    Country == "FR"~ 4,
    Country =="GFR"~ 7,
    Country =="ENG" ~ 4,
    TRUE ~ NA_real_))

sum(is.na(bts$Netopening))

bts$sweptarea <- bts$Netopening*bts$Distance

bts$WgtAtLngt <- bts$IndWgt*bts$HLNoAtLngt

bts$biomdens <- bts$WgtAtLngt/bts$sweptarea

#create different LFI for different L<sub>LF from 20 to 50

bts <- bts %>%
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

bts <-  bts[!is.na(bts$biomdens),]  

#all this calculates the different time series of LFI= (B  for L > L<sub> lf)/Btotal
# for each value of L<sub>lf from 20 to 50
# probably I could make it nicer with some loops and functions,
#but it seems to work ;-)

lfi_calc <-bts %>%
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

# too few data!, model does not run

model20 <- lm(lfi_20$lfi20 ~ poly(lfi_20$Year,5))

summary(model20)

model30 <- lm(lfi_30$lfi30 ~ poly(lfi_30$Year,5))

summary(model30)

model40 <- lm(lfi_40$lfi40 ~ poly(lfi_40$Year,5))

summary(model40)

model50 <- lm(lfi_50$lfi50 ~ poly(lfi_50$Year,5))

summary(model50)

#best fit is lfi when Llf = 30, kind of stable since 2006 

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
