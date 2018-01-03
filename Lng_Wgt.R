#Purpose: to get sensible measures of IndWgt by LengthClass for each species,
#as in FR-CGFS and SWC-IBTS there are no IndWgt values at all!!
#Author: Adriana Villamor, December 2017


library(icesDatras)
library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
#devtools::install_github("ices-tools-prod/icesVocab")
library(icesVocab)
library(reshape2)


#will use data from North Sea, the closest survey, I think
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

plot(all$LngtClass, all$IndWgt)
plot(log(all$LngtClass), log(all$IndWgt))

#althought the distribution of the slopes looks better with, species, 
#country and year, I lost so many year data that the >LFI time series remains short. 
#so I decided to fit only species and year.

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

##Plot slopes distribution
alllm%>%ggplot(aes(slope))+geom_histogram()+facet_wrap(~Species,scales = "free")


#now, test this parameters in the FR-CGFS survey

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

#Put the parameters together

cgfs <- left_join(cgfs, alllm)
sum(is.na(cgfs$IndWgt))
sum(is.na(cgfs$Netopening))
sum(is.na(cgfs$DoorSpread))

plot(cgfs$Netopening, cgfs$DoorSpread)
plot(cgfs$Netopening, log(cgfs$DoorSpread))

#regression to get DoorSpread from Netopening:
cgfslm <- cgfs %>% 
  filter(!is.na(DoorSpread), DoorSpread > 0) %>%
  filter(!is.na(Netopening), Netopening > 0)
  do(lm = lm(log(Netopening) ~ log(DoorSpread), data = ., 
             singular.ok = T, 
             na.action = na.exclude))

##Get the coefficients
cgfslm <- cgfslm %>% 
  tidy(lm)
cgfslm <- cgfslm %>% select(term,
                          estimate) %>% 
  spread(term, estimate)%>%
  rename(intercept = `(Intercept)`, 
         slope = `log(Netopening)`)

##Plot slopes distribution
cgfslm%>%ggplot(aes(slope))+geom_histogram()+facet_wrap(~Year,scales = "free")

#This is not working, don´t know why...continue down


#calculate biomdens, should remove logs!! 

cgfs$logIndWgt <- cgfs$intercept+cgfs$slope*log(cgfs$LngtClass)
cgfs$IndWgt <- 10^(cgfs$logIndWgt)
cgfs$WgAtLngt <- cgfs$IndWgt*cgfs$HLNoAtLngt

plot(cgfs$DoorSpread, cgfs$Netopening)

#should do something with Netopening and DoorSpread...


#I use Netopening because DoorSpread is almost NAs
#but this is probably an issue

cgfs$sweptarea <- cgfs$Distance*cgfs$Netopening

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

#best fit is lfi when Llf = 30 

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

#best fit is LFI when Llf is 40, time series form 1997

#check with SWC-IBTS

#now, test this parameters in the SWC-IBTS survey
#I tried to use the parameters by country, but at the end the regressions 
#are not significant, probably because variability is so small (?).

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




hh_ibts <- getDATRAS(record = "HH", "SWC-IBTS", years = 1985:2017, quarters = 1)

hl_ibts <- getDATRAS(record = "HL", "SWC-IBTS", years = 1985:2017, quarters = 1)

ca_ibts<-getDATRAS(record="CA",survey =  "SWC-IBTS", years = 1985:2017, quarters = 1)

speclist <- getCodeList("SpecWoRMS")


hl_ibts <- left_join(hl_ibts, 
                     speclist %>% 
                       select(Key, 
                              Description),
                     by = c("Valid_Aphia" = "Key")) %>%
  rename(Species = Description)



ca_ibts <- left_join(ca_ibts, 
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

ca_ibts <- ca_ibts%>%
  filter(Species %in% LFIspecies)
hl_ibts <- hl_ibts%>%
  filter(Species %in% LFIspecies)

# Transform LngtClass with LngtCode "." and "0" to cm!

ca_ibts<-rbind(ca_ibts%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               ca_ibts%>%filter(!LngtCode%in%c(".", "0")))

hl_ibts<-rbind(hl_ibts%>%filter(LngtCode%in%c(".", "0"))%>%mutate(LngtClass=LngtClass/10),
               hl_ibts%>%filter(!LngtCode%in%c(".", "0")))

# Change -9 values to NA

ca_ibts[ca_ibts == -9] <- NA

ca_ibts[ca_ibts == Inf] <- NA

hl_ibts <- hl_ibts[ ,-1]  
hh_ibts <- hh_ibts[ ,-1]        

ibts <- left_join(hl_ibts, ca_ibts)

#only valid and Day hauls
ibts <- left_join(ibts, hh_ibts) %>%
  filter( DayNight =="D",HaulVal =="V")

# better this: dat <-  dat %>% mutate(x = replace(x, x<0, NA))

ibts[ibts == -9] <- NA

#Put the parameters together

ibts <- left_join(ibts, alllm)
sum(is.na(ibts$IndWgt))
sum(is.na(ibts$Netopening))
sum(is.na(ibts$DoorSpread))

plot(ibts$DoorSpread, ibts$Netopening)
plot(log(ibts$DoorSpread), log(ibts$Netopening))

#calculate biomdens 
ibts$logIndWgt <- ibts$intercept+ibts$slope*log(ibts$LngtClass)
ibts$IndWgt <- 10^(ibts$logIndWgt)
ibts$WgAtLngt <- ibts$IndWgt*ibts$HLNoAtLngt

#In this survey there are enough data on DoosSpread, so i use it

ibts$sweptarea <- ibts$Distance*ibts$DoorSpread

ibts$biomdens <- ibts$WgAtLngt/ibts$sweptarea

#create different LFI time series for different L<sub>LF from 20 to 50 cm

ibts <- ibts %>%
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

ibts <-  ibts[!is.na(ibts$biomdens),]  

#all this calculates the different time series of LFI= (B  for L > L<sub> lf)/Btotal
# for each value of L<sub>lf from 20 to 50
# probably I could make it nicer with some loops and functions,
#but it seems to work ;-)

lfi_calc <-ibts %>%
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

#lfi when llf is 20

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




#other approaches:

library(plyr)
# Break up d by species, then fit the specified model to each piece and
# return a list
models <- dlply(all, "Species", function(df) 
  lm(IndWgt ~ LngtClass, data = all))

# Apply coef to each model and return a data frame
ldply(models, coef)

# Print the summary of each model, (all summaries are the same, weird)
l_ply(models, summary, .print = TRUE)
#don't like it

#other approach:

library(lme4)
library(lattice)
xyplot(IndWgt ~ LngtClass, groups=Species, data=all, type='l')

fits <- lmList(IndWgt ~ LngtClass | Species, data=all)
fits

## boo, don't know what is better!

##gamm
library(RODBC); library(ggplot2)
library(gamm4)

ggplot(droplevels(all[is.na(all$IndWgt)==FALSE, ]), aes(x=LngtClass, y=IndWgt))+
  geom_point()+facet_wrap(~Species , drop = TRUE)+
  stat_smooth()

# work out random part first

g1 <- gamm4(IndWgt~s(LngtClass, by=Species),
            data=all, na.action=na.omit)

##booo

#other, smooth.splines

all.sp <- all %>% group_by(Species)

plot(all$IndWgt~all$LngtClass, pch = 19, xlab = "x", ylab = "y",
     main = "Smoothing Splines with df = 8")

spl8 <- smooth.spline(x = all$LngtClass, y = all$IndWgt, df = 8)
lines(spl8, col = 4)


ggplot(droplevels(all[is.na(all$IndWgt)==FALSE, ]), aes(x=LngtClass, y=IndWgt))+
  geom_point()+stat_smooth()+facet_wrap(~Species , drop = TRUE)

model.1 <- nls(IndWgt ~ SSlogis(LngtClass, ASym, xmid, scal), data= all)
coef.sig<-coef(summary(model.1))[,1]
est.p<-coef.sig[1]/(1+exp((coef.sig[2]-all$LngtClass)/coef.sig[3]))
points(all$LngtClass,est.p,col=2)

coef.sig<-coef(summary(model.1))[,1]
est.p<-coef.sig[1]/(1+exp((coef.sig[2]-all$LngtClass)/coef.sig[3]))
points(all$LngtClass,est.p,col=2)

###this looks neat:

library(mgcv)
model<-gam(LngtClass~IndWgt , data=all)

summary_model <- summary(model)
summary_model

p_table <- data.frame(summary_model$p.table)
p_table <- within(p_table, {lci <- Estimate - qnorm(0.975) * Std..Error
uci <- Estimate + qnorm(0.975) * Std..Error})
p_table

#using plot.gam() with shift=intercept and trans=exp (for a poisson model) does the job. I can then 
#add the original data using points() 

plot(model,pages=1,residuals=TRUE,all.terms=TRUE,shade=TRUE,shade.col=2)


ggplot(droplevels(all[is.na(all$IndWgt)==FALSE, ]), aes(x=LngtClass, y=IndWgt))+
  geom_point()+facet_wrap(~Species , drop = TRUE)

