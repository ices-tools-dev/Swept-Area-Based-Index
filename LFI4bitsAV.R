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

#remove two outliers

bits <- bits[!bits$Distance > 15000, ]
bits <- bits[!bits$Distance2 > 15000, ]

sum(is.na(bits$Netopening)) #5410, a bit too much

p <- ggplot(bits, aes(Netopening, Ship))
q<- p + geom_point() 
q
#2 outliers in DAN2 I will transform them to NA

bits$Netopening[bits$Netopening > 25] <- NA

p <- ggplot(bits, aes(Netopening, Ship))
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

#create different LFI for different L<sub>LF from 10 to 50


# ahhhhh!!!

  lfi10 <- bits%>%
  group_by(Year)%>%
  if ( LngtClass < 10) sum(bits$biomdens)





  
else if (condition == 20){
    group_by(Year)%>%
    do(lf= for(LngClass <= 20), sum(biomdens)/sum(biomdens)))
    %>% return(lf)
  }}

  df1 %>%
    group_by(!!group_var) %>%
    summarise(seats=n())
}

# Test the function
my_summary(df1, 1)



write.csv(lfi, file = "LFI.csv")


#better

## Length classes consistency

bits%>%ggplot(aes(LngtClass,WgtAtLngt,color=Species))+geom_point() + facet_grid(Year ~ Species)

bits%>%ggplot(aes(Year,WgtAtLngt,color=LngtClass))+geom_point() + facet_grid(bits$Species)


bits%>%ggplot(aes(LngtClass,IndWgt,color=Species))+geom_point()




#save(ca_bits,file="ca_bits.RData" )
load("ca_bits.RData")

dim(ca_bits%>%filter(IndWgt==-9))[1]/dim(ca_bits)[1]*100## percent of IndWght gaps

##Add latin names
speclist <- getCodeList("SpecWoRMS")
speclist <-rbind(speclist,read.csv("Aphias_WoRMS_add.csv"))
ca_bits <-left_join(ca_bits, speclist%>%select(Key, Description), by=c("Valid_Aphia"="Key"))%>%rename(Species=Description)
##Define species list includded in the LFI calculation
LFIspecies<- c("Gadus morhua", "Platichthys flesus", "Pleuronectes platessa",  "Scophthalmus maximus" , "Merlangius merlangus"  )

##Fit the linear models
ca_bitslm<-ca_bits%>%filter(Species%in%LFIspecies)%>%filter(IndWgt>0)%>%filter(!is.na(IndWgt))%>%group_by(Year,Gear, AreaCode, Species)%>%
  do(lm=lm(log(IndWgt) ~ log(LngtClass), data = ., singular.ok = T, na.action=na.exclude))

ca_bit2lm<-ca_bit2%>%filter(Species%in%LFIspecies)%>%filter(IndWgt>0)%>%filter(!is.na(IndWgt))%>%group_by(Year,Gear, AreaCode, Species)%>%
  do(lm=lm(log(IndWgt) ~ log(LngtClass), data = ., singular.ok = T, na.action=na.exclude))







#### To do: filter W-L models with more than (10, 30?) observations. Maybe additional consistency check +-1 sd of all slopes in the database? 
table(ca_bits$AreaCode, ca_bits$Species[ca_bits$Species %in% spec])

###Check records with no slope estimates 
View(ca_bitslm[which(is.na(ca_bitslm$slope)),])
View(ca_bitslm[which(is.na(ca_bitslm$slope)),]%>%group_by(Year)%>%summarise(n()))

ca_bits%>%filter(Year==1991&Gear=="P20"&AreaCode=="37G8"&Species=="Gadus morhua")
log(53)
ca_bits%>%filter(Year==1991&Gear=="P20"&AreaCode=="37G8"&Species=="Gadus morhua")%>%
  lm(log(LngtClass)~log(IndWgt), data=.)


#save.image("LFI.RData")
load("LFI.RData")





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
