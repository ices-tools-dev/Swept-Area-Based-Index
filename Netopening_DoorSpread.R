#Purpose: to get sensible measures of DoorSpread by Netopening by year,
#almost missing in cgfs survey
#Author: Adriana Villamor, January 2018



plot(bits$Netopening,bits$DoorSpread)
plot(bits$DoorSpread,bits$Netopening)


library(icesDatras)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
#devtools::install_github("ices-tools-prod/icesVocab")
library(icesVocab)
library(reshape2)

#I put togeher as many values of Netopening and DoorSpread as possible:

hh_bits <- getDATRAS(record = "HH", "BITS", years = 1991:2017, quarters = 1)
hh_nsibts  <- getDATRAS(record="HH", survey= "NS-IBTS",years = 1965:2017 , quarters = 1)
hh_cgfs <- getDATRAS(record = "HH", "FR-CGFS", years = 1988:2017, quarters = 4)
hh_ibts <- getDATRAS(record = "HH", "SWC-IBTS", years = 1985:2017, quarters = 1)

hh_bits <- hh_bits %>%
  select(Survey, Year, Gear, Country, Stratum, HaulVal, 
         DayNight, Netopening, DoorSpread)
hh_nsibts <- hh_nsibts %>%
  select(Survey, Year, Gear, Country, Stratum, HaulVal, 
         DayNight, Netopening, DoorSpread)
hh_cgfs <- hh_cgfs %>%
  select(Survey, Year, Gear, Country, Stratum, HaulVal, 
         DayNight, Netopening, DoorSpread)
hh_ibts <- hh_ibts %>%
  select(Survey, Year, Gear, Country, Stratum, HaulVal, 
         DayNight, Netopening, DoorSpread)
hh_all <- left_join(hh_nsibts, hh_bits)
hh_all <- left_join(hh_all, hh_cgfs)
hh_all <- left_join(hh_all, hh_ibts)

#only valid and Day hauls
hh_all <- hh_all %>%
  filter( DayNight =="D", HaulVal =="V")


hh_all[hh_all == -9] <- NA

hh_all[hh_all == Inf] <- NA
plot(hh_all$DoorSpread, hh_all$Netopening)

#log of DoorSpread vs Netropening are linear, negative correlation
plot(log(hh_all$DoorSpread), hh_all$Netopening)

#Fit the model:

hh_alllm <- hh_all %>% 
  group_by(Year)%>%
  filter(!is.na(DoorSpread), DoorSpread > 0) %>%
  do(lm = lm(log(DoorSpread) ~ Netopening, data = ., 
             singular.ok = T, 
             na.action = na.exclude))

##Get the coefficients
hh_alllm <- hh_alllm %>% 
  tidy(lm)
hh_alllm <- hh_alllm %>% select(term,
                                  estimate) %>% 
  spread(term, estimate)%>%
  rename(intercept = `(Intercept)`, 
         slope = `Netopening`)

##Plot slopes distribution
hh_alllm%>%ggplot(aes(slope))+geom_histogram()


#test this parameters in cgfs,
hh_cgfs <- getDATRAS(record = "HH", "FR-CGFS", years = 1988:2017, quarters = 4)

hl_cgfs <- getDATRAS(record = "HL", "FR-CGFS", years = 1988:2017, quarters = 4)

ca_cgfs<-getDATRAS(record="CA",survey =  "FR-CGFS", years = 1988:2017, quarters = 4)

speclist <- getCodeList("SpecWoRMS")


hl_cgfs <- left_join(hl_cgfs, 
                     speclist %>% 
                       select(Key, 
                              Description),
                     by = c("Valid_Aphia" = "Key")) %>%
  rename(Species = Description)



ca_cgfs <- left_join(ca_cgfs, 
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

ca_cgfs <- ca_cgfs%>%
  filter(Species %in% LFIspecies)
hl_cgfs <- hl_cgfs%>%
  filter(Species %in% LFIspecies)

# Transform LngtClass with LngtCode "." and "0" to cm!

ca_cgfs<-rbind(ca_cgfs%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               ca_cgfs%>%filter(!LngtCode%in%c(".", "0")))

hl_cgfs<-rbind(hl_cgfs%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               hl_cgfs%>%filter(!LngtCode%in%c(".", "0")))

# Change -9 values to NA

ca_cgfs[ca_cgfs == -9] <- NA

ca_cgfs[ca_cgfs == Inf] <- NA

hl_cgfs <- hl_cgfs[ ,-1]  
hh_cgfs <- hh_cgfs[ ,-1]        

cgfs <- left_join(hl_cgfs, ca_cgfs)

#only valid and Day hauls
cgfs <- left_join(cgfs, hh_cgfs) %>%
  filter( DayNight =="D",HaulVal =="V")

# better this: dat <-  dat %>% mutate(x = replace(x, x<0, NA))

cgfs[cgfs == -9] <- NA
ca_nsibts  <- getDATRAS(record="CA", survey= "NS-IBTS",years = 1965:2017 , quarters = 1)

speclist <- getCodeList("SpecWoRMS")

ca_nsibts <- ca_nsibts[ ,-1]  

all<- ca_nsibts

all <- left_join(all, 
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

all <- all%>%
  filter(Species %in% LFIspecies)

all<-rbind(all%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
           all%>%filter(!LngtCode%in%c(".", "0")))

all[all == -9] <- NA

all[all == Inf] <- NA
alllm <- all %>% 
  group_by(Species, Year)%>%
  filter(!is.na(IndWgt), IndWgt > 0) %>%
  do(lm = lm(log(IndWgt) ~ log(LngtClass), data = ., 
             singular.ok = T, 
             na.action = na.exclude))

##Get the coefficients
alllm <- alllm %>% 
  tidy(lm)
alllm <- alllm %>% select(term,
                          estimate) %>% 
  spread(term, estimate)%>%
  rename(intercept = `(Intercept)`, 
         slope = `log(LngtClass)`)

#Put the parameters together

cgfs <- left_join(cgfs, alllm)
sum(is.na(cgfs$IndWgt))
sum(is.na(cgfs$Netopening))
sum(is.na(cgfs$DoorSpread))

plot(cgfs$Netopening, cgfs$DoorSpread)
plot(cgfs$Netopening, log(cgfs$DoorSpread))


#calculate biomdens,

cgfs$logIndWgt <- cgfs$intercept+cgfs$slope*log(cgfs$LngtClass)
cgfs$IndWgt <- 10^(cgfs$logIndWgt)
cgfs$WgAtLngt <- cgfs$IndWgt*cgfs$HLNoAtLngt

cgfs <- cgfs[ ,-c(84,85)]
cgfs <- left_join(cgfs, hh_alllm)

cgfs$logDoorSpread <- cgfs$intercept+cgfs$slope*(cgfs$Netopening)
cgfs$DoorSpread <- 10^(cgfs$logDoorSpread)
cgfs$sweptarea <- cgfs$Distance*cgfs$DoorSpread
sum(is.na(cgfs$DoorSpread))

cgfs$biomdens <- cgfs$WgAtLngt/cgfs$sweptarea

#create different LFI time series for different L<sub>LF from 20 to 50 cm

cgfs <- cgfs %>%
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

cgfs <-  cgfs[!is.na(cgfs$biomdens),]  

#all this calculates the different time series of LFI= (B  for L > L<sub> lf)/Btotal
# for each value of L<sub>lf from 20 to 50
# probably I could make it nicer with some loops and functions,
#but it seems to work ;-)

lfi_calc <-cgfs %>%
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

#best fit is lfi when Llf = 30, still the same graph
#as with Netopening, but improved R squared and p 

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







