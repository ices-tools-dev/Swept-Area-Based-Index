rm(list = ls())
####################
# LFI for the BITS #
####################
# Purpose: To calculate swept area and length-weight parameters for 
# Authors: Scott Large and Szymon Smoliński
# Date: 31 March 2017

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

#From 1993 to 2007, almost all data are cod, after that all 5 species are present, 
#should not be a problem right?

hh_bits <- getDATRAS(record = "HH", "BITS", years = 1991:2017, quarters = 1)

hl_bits <- getDATRAS(record = "HL", "BITS", years = 1991:2017, quarters = 1)

ca_bits<-getDATRAS(record="CA",survey =  "BITS", years = 1991:2017, quarters = 1)

speclist <- getCodeList("SpecWoRMS")

# speclist <- rbind(speclist,
#                   read.csv("Aphias_WoRMS_add.csv"))
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

# Transform LngtCode . and 0 to LngtClass in cm!

ca_bits<-rbind(ca_bits%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               ca_bits%>%filter(!LngtCode%in%c(".", "0")))

hl_bits<-rbind(hl_bits%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               hl_bits%>%filter(!LngtCode%in%c(".", "0")))

# Change -9 values to NA

ca_bits[ca_bits == -9] <- NA

#ca_bits[ca_bits == NaN] <- NA  
ca_bits[ca_bits == Inf] <- NA


##Fit the linear models
ca_bitslm <- ca_bits %>% 
  group_by(Country, Ship, Year, Gear,Species)%>%
  filter(!is.na(IndWgt), IndWgt > 0) %>%
  do(lm = lm(log(IndWgt) ~ log(LngtClass), data = ., 
             singular.ok = T, 
             na.action = na.exclude))


##Get the coefficients
ca_bitslm <- ca_bitslm %>% 
  tidy(lm)
ca_bitslm <- ca_bitslm %>% select(term,
                                  estimate) %>% 
  spread(term, estimate)%>%
  rename(intercept = `(Intercept)`, 
         slope = `log(LngtClass)`)

##Plot slopes distribution
ca_bitslm%>%ggplot(aes(slope))+geom_histogram()+facet_wrap(~Species,scales = "free")

#Only Platichthys flesus seems has more spread slopes, but that´s ok

ca_bitslm%>%group_by(Year, Species)%>%summarise(medianslope=(median(slope, na.rm = T)))%>%
  ggplot(aes(Year, medianslope,color=Species))+geom_line(size=1)+scale_color_brewer(palette = "Dark2")+theme_bw()


#regression seems quite fine, merging the three pieces together
#only valid and Day hauls

hl_bits <- hl_bits[ ,-1]  
hh_bits <- hh_bits[ ,-1]        

bits <- left_join(hl_bits, ca_bits)

bits <- left_join(bits, hh_bits) %>%
  filter( DayNight =="D", IndWgt > 0 ,HaulVal =="V")

# better this: dat <-  dat %>% mutate(x = replace(x, x<0, NA))
bits[bits == -9] <- NA

sum(is.na(bits$Distance)) #2195 out of 14300

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

bits$Distance2 <- 1000*(earth_distance(bits$ShootLong, bits$ShootLat, bits$HaulLong, bits$HaulLat))

p <- ggplot(bits, aes(Distance, Distance2))
q<- p + geom_point() 
q

#subsitute NAs in Distance for this calculation

bits <- transform(bits, Distance = ifelse(!is.na(Distance), Distance, Distance2))

sum(is.na(bits$Distance))
sum(is.na(bits$Distance2))  #572 out of 14300, is that ok? 




#remove outliers, up to where?

bits <- bits[!bits$Distance > 15000, ]


p <- ggplot(bits, aes(Distance, Year))
q<- p + geom_point() 
q

#do you think this is an outlier?

#FROM here, I will go to the end of the script and try two things:
#1:same LFI with DoorSpread instead of NetOpening
#2: back with NetOpening, substituting missing data 
#with the mean value per year


sum(is.na(bits$Netopening)) #5410, a bit too much

p <- ggplot(bits, aes(Netopening, Ship))
q<- p + geom_point() 
q
#2 outliers in DAN2 I will transform them to NA

bits$Netopening[bits$Netopening > 25] <- NA

sum(is.na(bits$DoorSpread)) #5725, even more

p <- ggplot(bits, aes(DoorSpread, Ship))
q<- p + geom_point() 
q

# one oulier here, out?

p <- ggplot(bits, aes(Netopening, Year))
q<- p + geom_point() 
q

p <- ggplot(bits, aes(Distance, Year))
q<- p + geom_point() 
q


#Still 26HF and HAF have very low values of Netopening, what to do with these?

#mean Netopening per ship and year shitttHEREEE

#temp <- aggregate(bits$Netopening,
#          list(Ship = bits$Ship, Year= bits$Year),
#          mean)

#there are too many ships with missing data, what to do?

p <- ggplot(bits, aes(DoorSpread, Ship)) #also here, same data I think
q<- p + geom_point() 
q
  
bits$sweptarea <- bits$Netopening*bits$Distance

bits$WgtAtLngt <- bits$IndWgt*bits$HLNoAtLngt

bits$biomdens <- bits$WgtAtLngt/bits$sweptarea

#trying to find more informative graphs

p <- ggplot(bits, aes(LngtClass, biomdens, color = Species))
q<- p + geom_point() 
r<- q + facet_grid(Species ~ Year)
r

#in order not to have to run all this everytime:

save(bits, file = "bits.RData")

load("bits.RData")
#create different LFI for different L<sub>LF from 10 to 50


# ahhhhh!!!

  bits <- bits %>%
    mutate(length_bin = case_when(
      LngtClass > 0 &
        LngtClass <= 10 ~ 10,
      LngtClass > 11 & 
            LngtClass <= 20 ~ 20,
      LngtClass > 21 &
        LngtClass <= 30 ~ 30,
      LngtClass > 31 &
        LngtClass <= 40 ~ 40,
      LngtClass > 41 &
        LngtClass <= 50 ~ 50,
      LngtClass >51 ~ 60,
      TRUE ~ NA_real_))


bits2 <-  bits[!is.na(bits$biomdens),]  
  
  
#lfi_all <- data.frame(matrix(ncol = 6, nrow = 26 ))
                        
#colnames(lfi_all)<- c("Year","lfi10","lfi20","lfi30","lfi40","lfi50")

#lfi_all$Year <- unique(bits$Year, na.rm =TRUE)
  
#lfi_all$lfi10 <- bits %>%
#    group_by(Year)%>%
#    filter(bits$length_bin >10) %>%
#    sum(biomdens)
  
#lfi_all$lfi10 <- bits %>% 
#  group_by(Year) %>%
#  sum(bits$biomdens, na.rm =TRUE)
  #

lfi_calca <-bits2 %>%
  group_by(Year, length_bin) %>%
  summarise(lfi = sum(biomdens))


lfi_10<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 10)%>%
  summarise(lfi10 = sum(lfi))
 

 lfi_20<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 20)%>%
  summarise(lfi20 = sum(lfi)) 

lfi_30<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 30)%>%
  summarise(lfi30 = sum(lfi)) 
  
lfi_40<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 40)%>%
  summarise(lfi40 = sum(lfi)) 
 
lfi_50<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 50)%>%
  summarise(lfi50 = sum(lfi)) 
   
lfiperyear <- lfi_calca %>% group_by(Year) %>% summarise(lfitot= sum(lfi))  

lfi_10 <- left_join(lfi_10,lfiperyear) 
lfi_20 <- left_join(lfi_20,lfiperyear)     
lfi_30 <- left_join(lfi_30,lfiperyear) 
lfi_40 <- left_join(lfi_40,lfiperyear) 
lfi_50 <- left_join(lfi_50,lfiperyear) 

lfi_10$lfi10 <- lfi_10$lfi10/lfi_10$lfitot
lfi_20$lfi20 <- lfi_20$lfi20/lfi_20$lfitot
lfi_30$lfi30 <- lfi_30$lfi30/lfi_30$lfitot
lfi_40$lfi40 <- lfi_40$lfi40/lfi_40$lfitot
lfi_50$lfi50 <- lfi_50$lfi50/lfi_50$lfitot

lfi_timeseries <- left_join(lfi_10,lfi_20)
lfi_timeseries <- left_join(lfi_timeseries,lfi_30)
lfi_timeseries <- left_join(lfi_timeseries,lfi_40)
lfi_timeseries <- left_join(lfi_timeseries,lfi_50)

lfi_timeseries <- subset(lfi_timeseries, select = -lfitot) 
lfi_timeseries <- melt(lfi_timeseries, id="Year")

ggplot(data=lfi_timeseries,
       aes(x=Year, y=value, colour=variable)) +
  geom_line()

model10 <- lm(lfi_10$lfi10 ~ poly(lfi_10$Year,5))

summary(model10)

model20 <- lm(lfi_20$lfi20 ~ poly(lfi_20$Year,5))

summary(model20)

model30 <- lm(lfi_30$lfi30 ~ poly(lfi_30$Year,5))

summary(model30)

model40 <- lm(lfi_40$lfi40 ~ poly(lfi_40$Year,5))

summary(model40)

model50 <- lm(lfi_50$lfi50 ~ poly(lfi_50$Year,5))

summary(model50)

#best is lfi when Llf = 50!, but maybe taking into account p-values, 
#residuals and dg might be better Llf = 40? 

ggplot(data=lfi_50,
       aes(x=Year, y=lfi50)) +
  geom_point()+geom_line()

ggplot(data=lfi_10,
       aes(x=Year, y=lfi10)) +
  geom_point()+geom_line()
ggplot(data=lfi_20,
       aes(x=Year, y=lfi20)) +
  geom_point()+geom_line()
ggplot(data=lfi_30,
       aes(x=Year, y=lfi30)) +
  geom_point()+geom_line()
ggplot(data=lfi_40,
       aes(x=Year, y=lfi40)) +
  geom_point()+geom_line()
# Test the function
my_summary(df1, 1)



write.csv(lfi, file = "LFI.csv")


#better

#########################################
#december 12: with DoorSpread instead of NetOpening:
#################################################
sum(is.na(bits$DoorSpread)) #5725, even more

p <- ggplot(bits, aes(DoorSpread, Ship))
q<- p + geom_point() 
q

# one oulier here, out?

bits$DoorSpread[bits$DoorSpread > 400] <- NA

p <- ggplot(bits, aes(Distance, Year))
q<- p + geom_point() 
q

#i will remove the values over 7500

bits$Distance[bits$Distance > 7500] <- NA

p <- ggplot(bits, aes(DoorSpread, Ship)) #also here, same data I think
q<- p + geom_point() 
q
#also small numbers should be substituted by NAs

bits$DoorSpread[bits$DoorSpread < 10] <- NA


bits$sweptarea <- bits$DoorSpread*bits$Distance

bits$WgtAtLngt <- bits$IndWgt*bits$HLNoAtLngt

bits$biomdens <- bits$WgtAtLngt/bits$sweptarea

#trying to find more informative graphs

p <- ggplot(bits, aes(LngtClass, biomdens, color = Species))
q<- p + geom_point() 
r<- q + facet_grid(Species ~ Year)
r


bits <- bits %>%
  mutate(length_bin = case_when(
    LngtClass > 0 &
      LngtClass <= 10 ~ 10,
    LngtClass > 11 & 
      LngtClass <= 20 ~ 20,
    LngtClass > 21 &
      LngtClass <= 30 ~ 30,
    LngtClass > 31 &
      LngtClass <= 40 ~ 40,
    LngtClass > 41 &
      LngtClass <= 50 ~ 50,
    LngtClass >51 ~ 60,
    TRUE ~ NA_real_))


sum(is.na(bits$biomdens))

bits2 <-  bits[!is.na(bits$biomdens),]  


#lfi_all <- data.frame(matrix(ncol = 6, nrow = 26 ))

#colnames(lfi_all)<- c("Year","lfi10","lfi20","lfi30","lfi40","lfi50")

#lfi_all$Year <- unique(bits$Year, na.rm =TRUE)

#lfi_all$lfi10 <- bits %>%
#    group_by(Year)%>%
#    filter(bits$length_bin >10) %>%
#    sum(biomdens)

#lfi_all$lfi10 <- bits %>% 
#  group_by(Year) %>%
#  sum(bits$biomdens, na.rm =TRUE)
#

lfi_calca <-bits2 %>%
  group_by(Year, length_bin) %>%
  summarise(lfi = sum(biomdens))


lfi_10<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 10)%>%
  summarise(lfi10 = sum(lfi))


lfi_20<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 20)%>%
  summarise(lfi20 = sum(lfi)) 

lfi_30<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 30)%>%
  summarise(lfi30 = sum(lfi)) 

lfi_40<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 40)%>%
  summarise(lfi40 = sum(lfi)) 

lfi_50<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 50)%>%
  summarise(lfi50 = sum(lfi)) 

lfiperyear <- lfi_calca %>% group_by(Year) %>% summarise(lfitot= sum(lfi))  

lfi_10 <- left_join(lfi_10,lfiperyear) 
lfi_20 <- left_join(lfi_20,lfiperyear)     
lfi_30 <- left_join(lfi_30,lfiperyear) 
lfi_40 <- left_join(lfi_40,lfiperyear) 
lfi_50 <- left_join(lfi_50,lfiperyear) 

lfi_10$lfi10 <- lfi_10$lfi10/lfi_10$lfitot
lfi_20$lfi20 <- lfi_20$lfi20/lfi_20$lfitot
lfi_30$lfi30 <- lfi_30$lfi30/lfi_30$lfitot
lfi_40$lfi40 <- lfi_40$lfi40/lfi_40$lfitot
lfi_50$lfi50 <- lfi_50$lfi50/lfi_50$lfitot

lfi_timeseries <- left_join(lfi_10,lfi_20)
lfi_timeseries <- left_join(lfi_timeseries,lfi_30)
lfi_timeseries <- left_join(lfi_timeseries,lfi_40)
lfi_timeseries <- left_join(lfi_timeseries,lfi_50)

lfi_timeseries <- subset(lfi_timeseries, select = -lfitot) 
lfi_timeseries <- melt(lfi_timeseries, id="Year")

ggplot(data=lfi_timeseries,
       aes(x=Year, y=value, colour=variable)) +
  geom_line()

model10 <- lm(lfi_10$lfi10 ~ poly(lfi_10$Year,5))

summary(model10)

model20 <- lm(lfi_20$lfi20 ~ poly(lfi_20$Year,5))

summary(model20)

model30 <- lm(lfi_30$lfi30 ~ poly(lfi_30$Year,5))

summary(model30)

model40 <- lm(lfi_40$lfi40 ~ poly(lfi_40$Year,5))

summary(model40)

model50 <- lm(lfi_50$lfi50 ~ poly(lfi_50$Year,5))

summary(model50)



ggplot(data=lfi_50,
       aes(x=Year, y=lfi50)) +
  geom_point()+geom_line()

ggplot(data=lfi_10,
       aes(x=Year, y=lfi10)) +
  geom_point()+geom_line()
ggplot(data=lfi_20,
       aes(x=Year, y=lfi20)) +
  geom_point()+geom_line()
ggplot(data=lfi_30,
       aes(x=Year, y=lfi30)) +
  geom_point()+geom_line()
ggplot(data=lfi_40,
       aes(x=Year, y=lfi40)) +
  geom_point()+geom_line()
# Test the function
my_summary(df1, 1)




#######################################
#substituting NAs of Netopening with mean of the year
############################################################


sum(is.na(bits$Netopening)) #5410, a bit too much

p <- ggplot(bits, aes(Netopening, Ship))
q<- p + geom_point() 
q
#2 outliers in DAN2 I will transform them to NA

bits$Netopening[bits$Netopening > 25] <- NA

sum(is.na(bits$Netopening)) #5410, a bit too much

p <- ggplot(bits, aes(Netopening, Ship))
q<- p + geom_point() 
q

p <- ggplot(bits, aes(Netopening, Year))
q<- p + geom_point() 
q

#remove also the small ones?


##this substitutes NAs with the mean of all years...


bits$Netopening[is.na(bits$Netopening)] <- round(mean(bits$Netopening, na.rm = TRUE))

sum(is.na(bits$Netopening))

p <- ggplot(bits, aes(Netopening, Year))
q<- p + geom_point() 
q

#graphs change far too much, 
#but still the selected time series is for Llf = 40

bits$sweptarea <- bits$Netopening*bits$Distance

bits$WgtAtLngt <- bits$IndWgt*bits$HLNoAtLngt

bits$biomdens <- bits$WgtAtLngt/bits$sweptarea


#in order not to have to run all this everytime:


# ahhhhh!!!

bits <- bits %>%
  mutate(length_bin = case_when(
    LngtClass > 0 &
      LngtClass <= 10 ~ 10,
    LngtClass > 11 & 
      LngtClass <= 20 ~ 20,
    LngtClass > 21 &
      LngtClass <= 30 ~ 30,
    LngtClass > 31 &
      LngtClass <= 40 ~ 40,
    LngtClass > 41 &
      LngtClass <= 50 ~ 50,
    LngtClass >51 ~ 60,
    TRUE ~ NA_real_))


bits2 <-  bits[!is.na(bits$biomdens),]  



lfi_calca <-bits2 %>%
  group_by(Year, length_bin) %>%
  summarise(lfi = sum(biomdens))


lfi_10<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 10)%>%
  summarise(lfi10 = sum(lfi))


lfi_20<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 20)%>%
  summarise(lfi20 = sum(lfi)) 

lfi_30<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 30)%>%
  summarise(lfi30 = sum(lfi)) 

lfi_40<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 40)%>%
  summarise(lfi40 = sum(lfi)) 

lfi_50<- lfi_calca %>%
  group_by(Year)%>%
  filter(length_bin > 50)%>%
  summarise(lfi50 = sum(lfi)) 

lfiperyear <- lfi_calca %>% group_by(Year) %>% summarise(lfitot= sum(lfi))  

lfi_10 <- left_join(lfi_10,lfiperyear) 
lfi_20 <- left_join(lfi_20,lfiperyear)     
lfi_30 <- left_join(lfi_30,lfiperyear) 
lfi_40 <- left_join(lfi_40,lfiperyear) 
lfi_50 <- left_join(lfi_50,lfiperyear) 

lfi_10$lfi10 <- lfi_10$lfi10/lfi_10$lfitot
lfi_20$lfi20 <- lfi_20$lfi20/lfi_20$lfitot
lfi_30$lfi30 <- lfi_30$lfi30/lfi_30$lfitot
lfi_40$lfi40 <- lfi_40$lfi40/lfi_40$lfitot
lfi_50$lfi50 <- lfi_50$lfi50/lfi_50$lfitot

lfi_timeseries <- left_join(lfi_10,lfi_20)
lfi_timeseries <- left_join(lfi_timeseries,lfi_30)
lfi_timeseries <- left_join(lfi_timeseries,lfi_40)
lfi_timeseries <- left_join(lfi_timeseries,lfi_50)

lfi_timeseries <- subset(lfi_timeseries, select = -lfitot) 
lfi_timeseries <- melt(lfi_timeseries, id="Year")

ggplot(data=lfi_timeseries,
       aes(x=Year, y=value, colour=variable)) +
  geom_line()

model10 <- lm(lfi_10$lfi10 ~ poly(lfi_10$Year,5))

summary(model10)

model20 <- lm(lfi_20$lfi20 ~ poly(lfi_20$Year,5))

summary(model20)

model30 <- lm(lfi_30$lfi30 ~ poly(lfi_30$Year,5))

summary(model30)

model40 <- lm(lfi_40$lfi40 ~ poly(lfi_40$Year,5))

summary(model40)

model50 <- lm(lfi_50$lfi50 ~ poly(lfi_50$Year,5))

summary(model50)

#best is lfi when Llf = 50!, but maybe taking into account p-values, 
#residuals and dg might be better Llf = 40? 

ggplot(data=lfi_50,
       aes(x=Year, y=lfi50)) +
  geom_point()+geom_line()

ggplot(data=lfi_10,
       aes(x=Year, y=lfi10)) +
  geom_point()+geom_line()
ggplot(data=lfi_20,
       aes(x=Year, y=lfi20)) +
  geom_point()+geom_line()
ggplot(data=lfi_30,
       aes(x=Year, y=lfi30)) +
  geom_point()+geom_line()
ggplot(data=lfi_40,
       aes(x=Year, y=lfi40)) +
  geom_point()+geom_line()
# Test the function
my_summary(df1, 1)




##############################################################
#November 27th, try to do the same thing with ROCKALL survey, with weigth at length 
# and with area as effort measures


hh_rockall <- getDATRAS(record = "HH", "ROCKALL", years = 2001:2017, quarters = 3)

hl_rockall <- getDATRAS(record = "HL", "ROCKALL", years = 2001:2017, quarters = 3)

ca_rockall<-getDATRAS(record="CA",survey =  "ROCKALL", years = 2001:2017, quarters = 3)

speclist <- getCodeList("SpecWoRMS")

# speclist <- rbind(speclist,
#                   read.csv("Aphias_WoRMS_add.csv"))
hl_rockall <- left_join(hl_rockall, 
                     speclist %>% 
                       select(Key, 
                              Description),
                     by = c("Valid_Aphia" = "Key")) %>%
  rename(Species = Description)



ca_rockall <- left_join(ca_rockall, 
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

ca_rockall <- ca_rockall%>%
  filter(Species %in% LFIspecies)
hl_rockall <- hl_rockall%>%
  filter(Species %in% LFIspecies)

# Transform LngtCode . and 0 to LngtClass in cm!

ca_rockall<-rbind(ca_rockall%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               ca_rockall%>%filter(!LngtCode%in%c(".", "0")))

hl_rockall<-rbind(hl_rockall%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               hl_rockall%>%filter(!LngtCode%in%c(".", "0")))


#this has already been done, the fit is good, with the exception of Platichthys flessus

##Fit the linear models
ca_rockalllm <- ca_rockall %>% 
  filter(!is.na(IndWgt)) %>%
  group_by(Country, Ship, Year, Gear,Species)%>%
  do(lm = lm(log(IndWgt) ~ log(LngtClass), data = ., 
             singular.ok = T, 
             na.action = na.exclude))

##Get the coefficients
ca_rockalllm <- ca_rockalllm %>% 
  tidy(lm)
ca_rockalllm <- ca_rockalllm %>% select(term,
                                  estimate) %>% 
  spread(term, estimate)%>%
  rename(intercept = `(Intercept)`, 
         slope = `log(LngtClass)`)

##Plot slopes distribution
ca_rockalllm%>%ggplot(aes(slope))+geom_histogram()+facet_wrap(~Species,scales = "free")
ca_bit2lm%>%ggplot(aes(slope))+geom_histogram()+facet_wrap(~Species,scales = "free")


ca_rockalllm%>%group_by(Year, Species)%>%summarise(medianslope=(median(slope, na.rm = T)))%>%
  ggplot(aes(Year, medianslope,color=Species))+geom_line(size=1)+scale_color_brewer(palette = "Dark2")+theme_bw()


#here


hl_rockall <- hl_rockall[ ,-1]  
hh_rockall <- hh_rockall[ ,-1]        

rockall <- left_join(hl_rockall, ca_rockall)


rockall <- left_join(rockall, hh_rockall) %>%
  filter( DayNight =="D", HaulVal =="V")


rockall$WgtAtLngt <- rockall$IndWgt*rockall$HLNoAtLngt
rockall$cpue <- rockall$HLNoAtLngt/(rockall$Netopening*rockall$Distance)

#trying to find more informative graphs

p <- ggplot(rockall, aes(LngtClass, cpue, color = Year))
q<- p + geom_point() 
r<- q + facet_grid(Species ~ Year)
r

#better
