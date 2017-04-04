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

# hh_bits <- getDATRAS(record = "HH", "BITS", years = 2001:2017, quarters = 1)
# saveRDS(hh_bits, file = "~/git/ices-dk/LFI/hh_bits.rds")


# ca_bits<-getDATRAS(record="CA",survey =  "BITS", years = 2001:2017, quarters = 1)

ca_bits <- readRDS("~/git/ices-dk/LFI/ca_bits.rds")

speclist <- getCodeList("SpecWoRMS")

# speclist <- rbind(speclist,
#                   read.csv("Aphias_WoRMS_add.csv"))
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

##Fit the linear models
ca_bitslm <- ca_bits %>% 
  filter(Species %in% LFIspecies) %>% 
  filter(IndWgt > 0) %>% 
  filter(!is.na(IndWgt)) %>%
  group_by(Year, Gear, AreaCode, Species)%>%
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
# ca_bitslm %>% group_by()%>% mutate()


# ca_bitslm
# ## Script check (first row of df)
# ca_bits%>%filter(Year==1991&Gear=="FOT"&AreaCode=="38G3"&Species=="Gadus morhua" )%>%
#   lm(log(IndWgt)~log(LngtClass), data = .) 







## Length classes consistency
ca_bits%>%ggplot(aes(LngtClass,IndWgt,color=LngtCode))+geom_point()
ca_bits<-rbind(ca_bits%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               ca_bits%>%filter(!LngtCode%in%c(".", "0")))

ca_bits%>%filter(IndWgt>0)%>%ggplot(aes(LngtClass,IndWgt,color=LngtCode))+geom_jitter(alpha=0.2)

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

##Get the coefficients
ca_bitslm<-ca_bitslm%>%tidy(lm)
ca_bitslm<-ca_bitslm%>%select(term,estimate)%>%spread(term,estimate)%>%rename(intercept=`(Intercept)`, slope=`log(LngtClass)`)
ca_bitslm%>%group_by()%>%mutate()


ca_bitslm
## Script check (first row of df)
ca_bits%>%filter(Year==1991&Gear=="FOT"&AreaCode=="38G3"&Species=="Gadus morhua" )%>%
  lm(log(IndWgt)~log(LngtClass), data = .) 

##Plot slopes distribution
ca_bitslm%>%ggplot(aes(slope))+geom_histogram()+facet_wrap(~Species,scales = "free")

ca_bitslm%>%group_by(Year, Species)%>%summarise(medianslope=(median(slope, na.rm = T)))%>%
  ggplot(aes(Year, medianslope,color=Species))+geom_line(size=1)+scale_color_brewer(palette = "Dark2")+theme_bw()

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
