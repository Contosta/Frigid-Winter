##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script calculates assesses changes in the timing (start and end), duration, and condition of winter 
#as defined by temperatures below freezing and snow-covered days for the manuscript, "Defining frigid winter illuminates its loss across seasonally snow-covered 
#areas of eastern North America" published in Environmental Research Letters (doi:10.1088/1748-9326/ab54f3)

#Code was developed by N. Casson, A. Contosta, and S. Nelson

#the dataset required to run this script ("metfin.csv") is located in this repository (with associated metadata)

####################################################################################
####################################################################################
####################################################################################

#Step 1: Initial set up (call libraries, import and compile climate data)

#call libraries

library(data.table)
library(dplyr)
library(Kendall)
library(matrixStats)
library(stringr)
library(tidyverse)
library(trend)
library(varhandle)
library(zoo)

#set working directory

setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\NSRC_Winter_Climate\\Data")

#read table
metfin <- read.table("metfin.csv", head = TRUE, sep = ",", na.strings=c("NA", "NAN"))

####################################################################################
#Flag temperature below freezing (TMAX and TMIN) and snow-covered days (SCD) that 
#define frigid winter
####################################################################################

###########
#frost days
###########

metfin$frostday = ifelse(metfin$TMINfin < 0, 1, 0)

############
#freeze days
############

metfin$freezeday = ifelse(metfin$TMAXfin < 0, 1, 0)

##################
#snow covered days
##################

metfin$modSCD = ifelse(metfin$modSWE > 0, 1, 0)

##################
#bare ground days
##################

metfin$modBGD = ifelse(metfin$modSWE == 0, 1, 0)

####################################################################################
#Determine start, end, and duration of frigid winter based on temperature
#defined as #of days when there is at least four freeze days within a seven day period
####################################################################################

#calculate cumulative freeze days and frost days over a seven day period
setDT(metfin)[, freezeday7 := Reduce(`+`, shift(freezeday, 0:6))]
setDT(metfin)[, frostday7 := Reduce(`+`, shift(frostday, 0:6))]

#select where freezeday7 is equal to 4 and the day of year is at least the second week in November (makes sure that at least one
#week is in cumulative freezeday7 calculation
startfz4 = metfin[which(metfin$freezeday7 == 4 & metfin$doy2 > 8),]

#subset data.frame to only include variables needed
DTfz = startfz4[ , c("SiteYr", "doy2", "freezeday7")]
names(DTfz) = c("SiteYr", "doyfz", "freezeday7")

#determine first day of year when there were four freeze days over a seven day period
FREEZbeg = DTfz %>%
  group_by(SiteYr) %>%
  filter(doyfz == min(doyfz)) %>%
  arrange(SiteYr, doyfz, freezeday7)
  names(FREEZbeg) = c("SiteYr", "FREEZbeg", "freezeday7")
  #remove duplicate row for 305801 in 1924 (not sure why there are two, but they have the same response
  FREEZbeg = FREEZbeg[!duplicated(FREEZbeg[ , 1]), ]

#determine last day of year when there were four freeze days over a seven day period
FREEZend = DTfz %>%
  group_by(SiteYr) %>%
  filter(doyfz == max(doyfz)) %>%
  arrange(SiteYr, doyfz, freezeday7)
  names(FREEZend) = c("SiteYr", "FREEZend", "freezeday7")

#join first and last days when there were four freeze days over a seven day period to metfin data.frame

metfin_mod = full_join(metfin, FREEZbeg, by=c('SiteYr'))

metfin_mod2 = full_join(metfin_mod, FREEZend, by=c('SiteYr'))

#flag instances where the first / last dates occur after / prior to the first day of November (data QA step)
metfin_mod2$startdate = ifelse(metfin_mod2$FREEZbeg > metfin_mod2$doy2, 0, 1)
metfin_mod2$enddate = ifelse(metfin_mod2$FREEZend < metfin_mod2$doy2, 1, 0)

#create subset where there are four frost days over a seven day period, and startdate and enddate are flagged as "1"
endft4 = metfin_mod2[which(metfin_mod2$frostday7 == 4 & metfin_mod2$startdate == 1 &metfin_mod2$enddate == 1),]

#create data.frame with freezeday7 and frostday7 to determine last day of cold period

DTft = endft4[ c("SiteYr", "doy2", "frostday7", "freezeday7")]
names(DTft)<-c("SiteYr", "doyft", "frostday7", "freezeday7")

#determine last day of year when there were four frost days over a seven day period

FROSend <- DTft %>%
  group_by(SiteYr) %>%
  filter(doyft == max(doyft)) %>%
  arrange(SiteYr, doyft, frostday7)
  names(FROSend) = c("SiteYr", "FROSend", "frostday7", "freezeday7")

#join last days when there were four frost days (end of cold period) over a seven day period to metfin data.frame

metfin_mod3 = full_join(metfin_mod2, FROSend, by=c('SiteYr'))

#rename as metfin.1

metfin.1 = metfin_mod3[order(metfin_mod3$Site_ID, metfin_mod3$WYear, metfin_mod3$doy2), ]

#remove columns no longer needed
metfin.1 = metfin.1[ , !names(metfin.1) %in% c("freezeday7.x", "frostday7.x", "freezeday7.y", "doyfz.y", "freezeday7.x.x", "startdate",
           "enddate", "frostday7.y", "freezeday7.y.y", "FREEZend")]

#calculate duration of cold period
metfin.1$COLDdur = metfin.1$FROSend - metfin.1$FREEZbeg

#indentify sustained cold period
metfin.1$COLD = ifelse(metfin.1$doy2 >= metfin.1$FREEZbeg & metfin.1$doy2 <= metfin.1$FROSend, 1, 0)

####################################################################################
#Determine start, end, and duration of frigid winter based on snow cover
#defined as #of days when there is at least four snow-covered days within a seven day period
####################################################################################

#calculate cumulative snow covered (SCD) and bare ground (BGD) over a seven day period
setDT(metfin.1)[, snowday7 := Reduce(`+`, shift(modSCD, 0:6))]
setDT(metfin.1)[, bareday7 := Reduce(`+`, shift(modBGD, 0:6))]

#select where snowday7 and bareday7 are equal to 4 and the day of year is at least the second week in November (makes sure that at least one
#week is in cumulative freezeday7 calculation
startsnow4 <-metfin.1[which(metfin.1$snowday7 == 4 & metfin.1$doy2 > 8), ]
endsnow4 <-metfin.1[which(metfin.1$bareday7 == 3 & metfin.1$doy2 > 8), ]

#subset to only include data needed for calculation
DTstartsnow <- startsnow4[ , c("SiteYr", "doy2")]
names(DTstartsnow)<-c("SiteYr", "doysnow")

#determine first day of snow covered season
mod.snowbeg <- DTstartsnow %>%
  group_by(SiteYr) %>%
  filter(doysnow == min(doysnow)) %>%
  arrange(SiteYr, doysnow)
  names(mod.snowbeg) = c("SiteYr", "modSWEbeg")

#subset to only include data needed for calculation
DTendsnow <- endsnow4[ , c("SiteYr", "doy2")]
names(DTendsnow)<-c("SiteYr", "doysnow")

#determine last day of snow covered season

mod.snowend <- DTendsnow %>%
  group_by(SiteYr) %>%
  filter(doysnow == max(doysnow)) %>%
  arrange(SiteYr, doysnow)
  names(mod.snowend) = c("SiteYr", "modSWEend")

#join SWEbeg and SWEend with metfin.1
metfin_mod6 <- full_join(metfin.1, mod.snowbeg, by=c("SiteYr"))
metfin_mod7 <- full_join(metfin_mod6, mod.snowend, by=c("SiteYr"))

metfin.2 = metfin_mod7[order(metfin_mod7$Site_ID, metfin_mod7$WYear, metfin_mod7$doy2), ]

#omit columns no longer needed
metfin.2 = metfin.2[ , !names(metfin.2) %in% c("snowday7", "bareday7")]

#calculate duration of snow-covered period
metfin.2$modSWEdur = ifelse(metfin.2$modSWEend - metfin.2$modSWEbeg > 0, metfin.2$modSWEend - metfin.2$modSWEbeg, 0)

#identify sustained snow-covered period
metfin.2$modSWEper = ifelse(metfin.2$doy2 >= metfin.2$modSWEbeg & metfin.2$doy2 <= metfin.2$modSWEend, 1, 0)

####################################################################################
#summarize number of frost days, snow-covered days, and average winter temperatures
#for four different definitions of winter: hiberal (dormant season), frigid, 
#astronomical, and meteorological. Determine maximum number of continuous frost days 
#and snow-covered days within hibernal winter
####################################################################################

########
#set up
########

#Create blank data.frame to which summary stats will be added
WYear = unique(metfin$WYear)
Site_ID = unique(metfin$Site_ID)

#paste LAT and LON together because there are duplicates of each
LAT_LON = unique(paste(metfin$LAT, metfin$LON, sep = ","))

#split the string to obtain correct number of sites (n = 37) for LAT and LON
spl_LAT_LON <- data.frame(do.call(rbind, str_split(LAT_LON, ",")))
names(spl_LAT_LON) <- c("LAT", "LON")

#unfactor the LAT and LON variables (convert categorical / factor variable back into numeric format)
LAT = unfactor(spl_LAT_LON$LAT)
LON = unfactor(spl_LAT_LON$LON)

#create empty data.frame into which summary statistics will be added
calc.ind = data.frame(matrix(data = NA, nrow =  length(WYear) * length(Site_ID), ncol = 4))
names(calc.ind) = c("WYear", "Site_ID", "LAT", "LON")

#populate data.frame with values for WYear (winter year)
calc.ind$WYear = rep(WYear, nrow(calc.ind) / length(WYear))
#sort data by WYear
calc.ind = calc.ind[order(calc.ind$WYear), ]
#populate data.frame with site ID, LAT, and LON
calc.ind$Site_ID = rep(Site_ID, nrow(calc.ind) / length(Site_ID))
calc.ind$LAT = rep(LAT, nrow(calc.ind) / length(LAT))
calc.ind$LON = rep(LON, nrow(calc.ind) / length(LON))

#create column for unique site x winter year combination
calc.ind$SiteYr = paste(calc.ind$Site_ID, calc.ind$WYear, sep = " ")

#add column to specify subregions within study area as defined by longitude
calc.ind$zone = ifelse(calc.ind$LON < -87, "west", 
                ifelse(calc.ind$LON >= -87 & calc.ind$LON <= -78, "central",
                ifelse(calc.ind$LON > -78, "east", NA)))

#order by winter year and site 
calc.ind = calc.ind[order(calc.ind$WYear, calc.ind$Site_ID), ]

#add columns to metfin.2 to identify winter using different seasonal definitions
#DOR = dormant period, CAL = calendar or meteorological winter, AST = astronomical winter

#Dormant period
metfin.2$DOR = ifelse(metfin.2$Month > 4, 0, 1)
#Calendar period
metfin.2$CAL = ifelse(metfin.2$Month == 12 | metfin.2$Month == 1 | metfin.2$Month == 2, 1, 0)
#Astronomical period
metfin.2$AST = ifelse(metfin.2$Month == 1 | metfin.2$Month == 2 | metfin.2$Month == 3, 1, 0)

##############################
#Cold Days (e.g., frost days)#
##############################

frost_dor = metfin.2[ , c("SiteYr", "frostday")]
frost_dor <- na.omit(frost_dor)

#count over entire "dormant" season
count.frost_dor = frost_dor %>%
  group_by(SiteYr) %>%
  summarize(count.frost_dor = sum(frostday))

#count over "frigid winter"

frost_win = metfin.2[metfin.2$COLD == 1, c("SiteYr", "frostday", "COLDdur")]
frost_win <- na.omit(frost_win)

count.frost_win = frost_win %>%
  group_by(SiteYr) %>%
  summarize(count.frost_win = sum(frostday))

#count over "astronomical winter"

frost_ast = metfin.2[metfin.2$AST == 1, c("SiteYr", "frostday")]
frost_ast <- na.omit(frost_ast)

count.frost_ast = frost_ast %>%
  group_by(SiteYr) %>%
  summarize(count.frost_ast = sum(frostday))

#count over "calendar winter"

frost_cal = metfin.2[metfin.2$CAL == 1, c("SiteYr", "frostday")]
frost_cal <- na.omit(frost_cal)

count.frost_cal = frost_cal %>%
  group_by(SiteYr) %>%
  summarize(count.frost_cal = sum(frostday))

#duration over entire "dormant" season

frost.rle.mean_dor = frost_dor %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$frostday==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Mean=if(length(tmp)==0) 0
             else mean(tmp)) })

frost.rle.max_dor = frost_dor %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$frostday==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Max=if(length(tmp)==0) 0
             else max(tmp)) })

frost.dur.mean_dor = frost.rle.mean_dor %>%
  group_by(SiteYr) %>%
  summarize(frost.dur.mean_dor = mean(Mean))

frost.dur.max_dor = frost.rle.max_dor %>%
  group_by(SiteYr) %>%
  summarize(frost.dur.max_dor = mean(Max))

#duration over "cold" period

frost.rle.mean_win = frost_win %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$frostday==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Mean=if(length(tmp)==0) 0
             else mean(tmp)) })

frost.rle.max_win = frost_win %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$frostday==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Max=if(length(tmp)==0) 0
             else max(tmp)) })

frost.dur.mean_win = frost.rle.mean_win %>%
  group_by(SiteYr) %>%
  summarize(frost.dur.mean_win = mean(Mean))

frost.dur.max_win = frost.rle.max_win %>%
  group_by(SiteYr) %>%
  summarize(frost.dur.max_win = mean(Max))

#duration over "astronomical" winter

frost.rle.mean_ast = frost_ast %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$frostday==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Mean=if(length(tmp)==0) 0
             else mean(tmp)) })

frost.rle.max_ast = frost_ast %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$frostday==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Max=if(length(tmp)==0) 0
             else max(tmp)) })

frost.dur.mean_ast = frost.rle.mean_ast %>%
  group_by(SiteYr) %>%
  summarize(frost.dur.mean_ast = mean(Mean))

frost.dur.max_ast = frost.rle.max_ast %>%
  group_by(SiteYr) %>%
  summarize(frost.dur.max_ast = mean(Max))

#duration over "calendar" winter

frost.rle.mean_cal = frost_cal %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$frostday==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Mean=if(length(tmp)==0) 0
             else mean(tmp)) })

frost.rle.max_cal = frost_cal %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$frostday==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Max=if(length(tmp)==0) 0
             else max(tmp)) })

frost.dur.mean_cal = frost.rle.mean_cal %>%
  group_by(SiteYr) %>%
  summarize(frost.dur.mean_cal = mean(Mean))

frost.dur.max_cal = frost.rle.max_cal %>%
  group_by(SiteYr) %>%
  summarize(frost.dur.max_cal = mean(Max))

#join frost stats togther and then add to calc.ind

frost.1 = full_join(count.frost_dor, count.frost_win, by = "SiteYr")
frost.2 = full_join(frost.1, count.frost_ast, by = "SiteYr")
frost.3 = full_join(frost.2, count.frost_cal, by = "SiteYr")
frost.4 = full_join(frost.3, frost.dur.mean_dor, by = "SiteYr")
frost.5 = full_join(frost.4, frost.dur.max_dor, by = "SiteYr")
frost.6 = full_join(frost.5, frost.dur.mean_win, by = "SiteYr")
frost.7 = full_join(frost.6, frost.dur.max_win, by = "SiteYr")
frost.8 = full_join(frost.7, frost.dur.mean_ast, by = "SiteYr")
frost.9 = full_join(frost.8, frost.dur.max_ast, by = "SiteYr")
frost.10 = full_join(frost.9, frost.dur.mean_cal, by = "SiteYr")
frost.11 = full_join(frost.10, frost.dur.max_cal, by = "SiteYr")

calc.ind.1 = full_join(calc.ind, frost.11, by = "SiteYr")

###################
#Snow-Covered Days#
###################

modSCD_dor = metfin.2[ , c("SiteYr", "modSCD")]
modSCD_dor <- na.omit(modSCD_dor)

#count over entire "dormant" season
count.modSCD_dor = modSCD_dor %>%
  group_by(SiteYr) %>%
  summarize(count.modSCD_dor = sum(modSCD))

#count over "frigid winter"

modSCD_win = metfin.2[metfin.2$COLD == 1, c("SiteYr", "modSCD", "COLDdur")]
modSCD_win <- na.omit(modSCD_win)

count.modSCD_win = modSCD_win %>%
  group_by(SiteYr) %>%
  summarize(count.modSCD_win = sum(modSCD))

#count over "astronomical winter"

modSCD_ast = metfin.2[metfin.2$AST == 1, c("SiteYr", "modSCD")]
modSCD_ast <- na.omit(modSCD_ast)

count.modSCD_ast = modSCD_ast %>%
  group_by(SiteYr) %>%
  summarize(count.modSCD_ast = sum(modSCD))

#count over "calendar winter"

modSCD_cal = metfin.2[metfin.2$CAL == 1, c("SiteYr", "modSCD")]
modSCD_cal <- na.omit(modSCD_cal)

count.modSCD_cal = modSCD_cal %>%
  group_by(SiteYr) %>%
  summarize(count.modSCD_cal = sum(modSCD))

#duration over entire "dormant" season

modSCD.rle.mean_dor = modSCD_dor %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$modSCD==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Mean=if(length(tmp)==0) 0
             else mean(tmp)) })

modSCD.rle.max_dor = modSCD_dor %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$modSCD==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Max=if(length(tmp)==0) 0
             else max(tmp)) })

modSCD.dur.mean_dor = modSCD.rle.mean_dor %>%
  group_by(SiteYr) %>%
  summarize(modSCD.dur.mean_dor = mean(Mean))

modSCD.dur.max_dor = modSCD.rle.max_dor %>%
  group_by(SiteYr) %>%
  summarize(modSCD.dur.max_dor = mean(Max))

#duration over "cold" period

modSCD.rle.mean_win = modSCD_win %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$modSCD==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Mean=if(length(tmp)==0) 0
             else mean(tmp)) })

modSCD.rle.max_win = modSCD_win %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$modSCD==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Max=if(length(tmp)==0) 0
             else max(tmp)) })

modSCD.dur.mean_win = modSCD.rle.mean_win %>%
  group_by(SiteYr) %>%
  summarize(modSCD.dur.mean_win = mean(Mean))

modSCD.dur.max_win = modSCD.rle.max_win %>%
  group_by(SiteYr) %>%
  summarize(modSCD.dur.max_win = mean(Max))

#duration over "astronomical" winter

modSCD.rle.mean_ast = modSCD_ast %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$modSCD==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Mean=if(length(tmp)==0) 0
             else mean(tmp)) })

modSCD.rle.max_ast = modSCD_ast %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$modSCD==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Max=if(length(tmp)==0) 0
             else max(tmp)) })

modSCD.dur.mean_ast = modSCD.rle.mean_ast %>%
  group_by(SiteYr) %>%
  summarize(modSCD.dur.mean_ast = mean(Mean))

modSCD.dur.max_ast = modSCD.rle.max_ast %>%
  group_by(SiteYr) %>%
  summarize(modSCD.dur.max_ast = mean(Max))

#duration over "calendar" winter

modSCD.rle.mean_cal = modSCD_cal %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$modSCD==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Mean=if(length(tmp)==0) 0
             else mean(tmp)) })

modSCD.rle.max_cal = modSCD_cal %>%
  group_by(SiteYr) %>%
  do({tmp <- with(rle(.$modSCD==1), lengths[values])
  data.frame(SiteYr=.$SiteYr, Max=if(length(tmp)==0) 0
             else max(tmp)) })

modSCD.dur.mean_cal = modSCD.rle.mean_cal %>%
  group_by(SiteYr) %>%
  summarize(modSCD.dur.mean_cal = mean(Mean))

modSCD.dur.max_cal = modSCD.rle.max_cal %>%
  group_by(SiteYr) %>%
  summarize(modSCD.dur.max_cal = mean(Max))

#join modSCD stats togther and then add to calc.ind

modSCD.1 = full_join(count.modSCD_dor, count.modSCD_win, by = "SiteYr")
modSCD.2 = full_join(modSCD.1, count.modSCD_ast, by = "SiteYr")
modSCD.3 = full_join(modSCD.2, count.modSCD_cal, by = "SiteYr")
modSCD.4 = full_join(modSCD.3, modSCD.dur.mean_dor, by = "SiteYr")
modSCD.5 = full_join(modSCD.4, modSCD.dur.max_dor, by = "SiteYr")
modSCD.6 = full_join(modSCD.5, modSCD.dur.mean_win, by = "SiteYr")
modSCD.7 = full_join(modSCD.6, modSCD.dur.max_win, by = "SiteYr")
modSCD.8 = full_join(modSCD.7, modSCD.dur.mean_ast, by = "SiteYr")
modSCD.9 = full_join(modSCD.8, modSCD.dur.max_ast, by = "SiteYr")
modSCD.10 = full_join(modSCD.9, modSCD.dur.mean_cal, by = "SiteYr")
modSCD.11 = full_join(modSCD.10, modSCD.dur.max_cal, by = "SiteYr")

calc.ind.2 = full_join(calc.ind.1, modSCD.11, by = "SiteYr")

#####################
#AVERAGE TEMPERATURES
#####################

TMINfin_dor = metfin.2[ , c("SiteYr", "TMINfin")]
TMINfin_dor <- na.omit(TMINfin_dor)

#avg over entire "dormant" season
TMIN_dor = TMINfin_dor %>%
  group_by(SiteYr) %>%
  summarize(TMIN_dor = mean(TMINfin))

#avg over "frigid winter"

TMINfin_win = metfin.2[metfin.2$COLD == 1, c("SiteYr", "TMINfin", "COLDdur")]
TMINfin_win <- na.omit(TMINfin_win)

TMIN_win = TMINfin_win %>%
  group_by(SiteYr) %>%
  summarize(TMIN_win = mean(TMINfin))

#avg over "astronomical winter"

TMINfin_ast = metfin.2[metfin.2$AST == 1, c("SiteYr", "TMINfin")]
TMINfin_ast <- na.omit(TMINfin_ast)

TMIN_ast = TMINfin_ast %>%
  group_by(SiteYr) %>%
  summarize(TMIN_ast = mean(TMINfin))

#avg over "calendar winter"

TMINfin_cal = metfin.2[metfin.2$CAL == 1, c("SiteYr", "TMINfin")]
TMINfin_cal <- na.omit(TMINfin_cal)

TMIN_cal = TMINfin_cal %>%
  group_by(SiteYr) %>%
  summarize(TMIN_cal = mean(TMINfin))

#join TMIN stats togther and then add to calc.ind

TMIN.1 = full_join(TMIN_dor, TMIN_win, by = "SiteYr")
TMIN.2 = full_join(TMIN.1, TMIN_ast, by = "SiteYr")
TMIN.3 = full_join(TMIN.2, TMIN_cal, by = "SiteYr")

calc.ind.3 = full_join(calc.ind.2, TMIN.3, by = "SiteYr")

#add statistics for start, and, and duration of frigid winter as defined by coldness or snow cover
FREEZbeg = FREEZbeg[ , c("SiteYr", "FREEZbeg")]
FROSend = FROSend[ , c("SiteYr", "FROSend")]
COLD = full_join(FREEZbeg, FROSend, by = "SiteYr")
COLD$COLDdur = COLD$FROSend - COLD$FREEZbeg

mod.snowbeg = mod.snowbeg[ , c("SiteYr", "modSWEbeg")]
mod.snowend = mod.snowend[ , c("SiteYr", "modSWEend")]
modSWE = full_join(mod.snowbeg, mod.snowend, by = "SiteYr")
modSWE$modSWEdur = modSWE$modSWEend - modSWE$modSWEbeg

WIN = full_join(COLD, modSWE, by = "SiteYr")

calc.ind.4 = full_join(calc.ind.3, WIN, by = "SiteYr")

####################################################################################
#Calculate change over time (Mann Kendall tau, p, sen slope) in the start, end, and 
#duration of frigid winter
####################################################################################

#rename calc.ind.3 data.frame to calc.fin
calc.fin =  calc.ind.4

#write function for sen slope

# Sen slope
sen = function(x, y){
  xx = outer(x, x, "-")
  yy = outer(y, y, "-")
  zz = yy/xx
  
  s1 = zz[lower.tri(zz)]
  s2 = subset(s1, s1!="NA" & s1!=Inf & s1!=-Inf)
  slope = median(s2)
  
  i1 = y - slope*x
  i2 = subset(i1, i1!="NA")
  intercept = median(i2)
  gi = data.frame(slope, intercept)
  return(gi)
}


#order data by Site_ID and WYear

calc.fin = calc.fin[order(calc.fin$Site_ID), ]

#omit duplicate SiteYrs from the beginning and the end of the dataframe to determine when to start / end model subsets
start.rows <- !duplicated(calc.fin$Site_ID)
end.rows <- !duplicated(calc.fin$Site_ID, fromLast = TRUE)

#extract original row indexes associated with the start and end of an Site_ID
#add one row onto the start and subtract one row from the end to remove potential artifacts with the first and last year of the record
sr.ind <- seq_along(calc.fin$Site_ID)[start.rows] + 1
er.ind <- seq_along(calc.fin$Site_ID)[end.rows] - 1

#check that the number of events is the same for the start and end
print(length(sr.ind))
print(length(er.ind))

#extract response, predictor, and co-variables (plus site-level info on co-variables)
inds = calc.fin[ , -(1:6)]

pred = calc.fin[ , "WYear"]

vars = (calc.fin[ , c(2:4, 6)])
vars.1 = vars[!duplicated(vars[1]), ]

#create containers to hold results of analysis
p.value = data.frame(matrix(nrow = length(sr.ind), ncol = ncol(inds)))
slope = data.frame(matrix(nrow = length(sr.ind), ncol = ncol(inds)))

#begin loop
for (h in seq(1,length(sr.ind))) {
  
  for (i in 1:ncol(inds)) {
    
    nn = data.frame(inds[sr.ind[h]:er.ind[h],i], pred[sr.ind[h]:er.ind[h]])
    nn. = nn[complete.cases(nn),]
    if(nrow(nn.)>0 ){
      p.value[h,i] = mk.test(nn.[,1])$p.value
      slope[h,i] = sen(nn.[,2], nn.[,1])$slope
    }
  }
}

#add names to each dataframe
names(p.value) = paste("pval", names(inds), sep = "_")
names(slope) = paste("slope", names(inds), sep = "_")

#compile tau, p, and sen slope into one master summary table
sum.tab.all = cbind(vars.1, p.value, slope)

#write calc.fin and sum.tab.all tables to run in other scripts
write.table(calc.fin, "calcfin_tdc.csv", col.names = T, row.names = F, sep = ",")
write.table(sum.tab.all, "winter_tdc.csv", col.names = T, row.names = F, sep = ",")
