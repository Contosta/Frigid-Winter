##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#This script produces the figures for the manuscript "Defining frigid winter illuminates its loss across seasonally snow-covered
#areas of eastern North America" published in Environmental Research Letters (doi:10.1088/1748-9326/ab54f3)

#Code was developed by N. Casson, A. Contosta, and S. Nelson

#the datasets required to run this script ("calcfin_tdc.csv", "metfin.csv" and "sumtab_tdc.csv") is located in this repository (with associated metadata)
#and was produced by running the "Winter_Timing_Duration_Condition.r" script along with the file
#"metfin.csv"

#Two spatial data files are also required and are located in this repository. They are:
# "gpr_000b11a_e.shp"; Source:Boundary Files, 2011 Census. Statistics Canada Catalogue no. 92-160-X.
# "cb_2016_us_state_500k.shp"; Source: US Census Bureau, Geography Division

#One image file ("legend.png") are required to create the map legends and is located in this repository

####################################################################################
#Initial set up
####################################################################################

##set wd
setwd("C:\\Users\\alix\\Box Sync\\UNH\\Projects\\NSRC_Winter_Climate\\Data")

##load libraries
library(tidyverse)
library(trend) 
library(reshape)
library(reshape2)
library(tidyr)
library(ggmap)
#library(ggplot2); may not be required with tidyverse 
library(maps)
library(mapproj)
library(mapdata)
library(rgeos)
library(maptools)
library(sp)
library(raster)
library(rgdal)
library(cowplot)
#library(dplyr); may not be required with tidyverse


#read in data
calcfin <- read.csv("calcfin_tdc.csv")
sumtab <- read.csv("winter_tdc.csv")
metfin <- read.csv("metfin.csv")

#classify sites into three zones
calcfin$ZONE <- ifelse(calcfin$LON< -87, "West", 
                     ifelse(calcfin$LON> -78, "East", "Central"))
calcfin$ZONE <- factor(calcfin$ZONE, levels=c("West", "Central", "East"))
metfin$ZONE <- ifelse(metfin$LON< -87, "West", 
                     ifelse(metfin$LON> -78, "East", "Central"))
metfin$ZONE <- factor(metfin$ZONE, levels=c("West", "Central", "East"))


#order zones from west to east
col_idx  <-  grep("ZONE", names(calcfin))
calcfin  <-  calcfin[, c(col_idx, (1:ncol(calcfin))[-col_idx])]

#classify dates into decades
calcfin$decade <- ifelse(calcfin$WYear>1909&calcfin$WYear<1920, 1910, 
                        ifelse(calcfin$WYear>1919&calcfin$WYear<1930, 1920,
                        ifelse(calcfin$WYear>1929&calcfin$WYear<1940, 1930,
                        ifelse(calcfin$WYear>1939&calcfin$WYear<1950, 1940,
                        ifelse(calcfin$WYear>1949&calcfin$WYear<1960, 1950,
                        ifelse(calcfin$WYear>1959&calcfin$WYear<1970, 1960,
                        ifelse(calcfin$WYear>1969&calcfin$WYear<1980, 1970,
                        ifelse(calcfin$WYear>1979&calcfin$WYear<1990, 1980,
                        ifelse(calcfin$WYear>1989&calcfin$WYear<2000, 1990,
                        ifelse(calcfin$WYear>1999&calcfin$WYear<2010, 2000,
                        ifelse(calcfin$WYear>2009&calcfin$WYear<2020, 2010,NA)))))))))))

metfin$decade <- ifelse(metfin$WYear>1909&metfin$WYear<1920, 1910, 
                        ifelse(metfin$WYear>1919&metfin$WYear<1930, 1920,
                        ifelse(metfin$WYear>1929&metfin$WYear<1940, 1930,
                        ifelse(metfin$WYear>1939&metfin$WYear<1950, 1940,
                        ifelse(metfin$WYear>1949&metfin$WYear<1960, 1950,
                        ifelse(metfin$WYear>1959&metfin$WYear<1970, 1960,
                        ifelse(metfin$WYear>1969&metfin$WYear<1980, 1970,
                        ifelse(metfin$WYear>1979&metfin$WYear<1990, 1980,
                        ifelse(metfin$WYear>1989&metfin$WYear<2000, 1990,
                        ifelse(metfin$WYear>1999&metfin$WYear<2010, 2000,
                        ifelse(metfin$WYear>2009&metfin$WYear<2020, 2010,NA)))))))))))


####################################################################################
#1. Format data for making figures
####################################################################################

names(calcfin)
names(sumtab)

##convert data to long format and merge (NOTE - WATCH COLUMN TITLES)
calcfin_long  <-  gather(calcfin, metric.name, value, count.frost_dor:modSWEdur, factor_key=TRUE)
calcfin_long2 <- separate(data = calcfin_long, col=metric.name, into=c("metric.name", "period"), sep = "_")
sumtab_long <- gather(sumtab, metric.name, value, pval_count.frost_dor:slope_modSWEdur, factor_key=TRUE)
sumtab_long2 <- separate(data = sumtab_long, col = metric.name, into=c("type", "metric.name","period"), sep = "_")

#merge calcfin and sumtab to format data for making plots
plotdata <- merge(calcfin_long2, sumtab_long2, by=c("metric.name", "period", "Site_ID", "LAT", "LON"))

plotdata_split <- spread(plotdata, type, value.y)

#create a column to indicate if there is a positive, significant slope, a negative, significant slope or no significant slope
plotdata_split$bp <- ifelse(plotdata_split$pval<0.05&plotdata_split$slope>0, "pos", 
                          ifelse(plotdata_split$pval<0.05&plotdata_split$slope<0, "neg", 
                                 "ns"))
plotdata_split$slope <- ifelse(plotdata_split$bp=="ns", 0, plotdata_split$slope)

#make the ZONE and bp columns into factors
plotdata_split$ZONE <- as.factor(plotdata_split$ZONE)

plotdata_split <- transform(plotdata_split,
                          bp=factor(bp,levels=c("neg", "ns","pos")))

###read in and format map data
#read in canadian and us data and crop to study extent
can <- readOGR(dsn=".", layer = "gpr_000b11a_e")
can_crop  <-  crop(can, extent(-100, -62, 37.5, 53))

us <- readOGR(dsn=".", layer = "cb_2016_us_state_500k")
us_crop  <-  crop(us, extent(-100, -62, 37.5, 53))

#convert to ggplot object
us_df <- fortify(us_crop)
can_df <- fortify(can_crop)

##################
####FIGURE 1######
##################

#create a dataframe with just start and end dates based on coldness criteria
start.end_all <- plotdata_split[which(plotdata_split$metric.name=="FREEZbeg"|plotdata_split$metric.name=="FROSend"),]
start.end_all$grp  <-  paste(start.end_all$Site_ID,start.end_all$metric.name)

#line plot of trends in start and end of winter based on coldness criteria
start.end_all_plot <- ggplot(data=start.end_all)+aes(WYear, value.x, group=grp, colour=bp)+geom_line(alpha=0.1)+theme(legend.position="none")+
  facet_wrap(~ZONE)+scale_colour_manual(values = c("red", "dark grey", "red"), guide=FALSE)+
  coord_flip()+scale_y_continuous(breaks=c(32,93, 152),
                                  labels=c("Nov", "Jan","Mar"),name = "")+scale_x_reverse()+
  theme_bw()+stat_smooth(geom="line",method="lm", se=F, alpha=0.5, size=1.5)+xlab("Year")+theme(text = element_text(size=20))


#create a dataframe of just the winter duration based on coldness criteria
cold_dur <- plotdata_split[which(plotdata_split$metric.name=="COLDdur"),]

#map of trends in winter duration based on coldness criteria
cold_dur_map <- ggplot()+
  geom_polygon(data=us_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_polygon(data=can_df, alpha=0.5,fill="#d8d0d0", colour="black", size=0.3, aes(x=long, y=lat, group=group))+
  geom_point(data=cold_dur, aes(x = LON, y = LAT, fill=bp, size=abs(slope)*10), color="black", pch=21,alpha=0.5)+
  scale_fill_manual(values = alpha(c( "#d65555", "dark grey","#595bcc"),0.5), guide=F)+
  scale_size(breaks=c(0, 1, 2,3), range = c(2, 8))+
  theme_void()+theme(legend.position ="none")+
  geom_vline(xintercept = -78, color="black", size=1)+geom_vline(xintercept = -87, color="black", size=1)

#combine plots to make figure 1
fig1 <- plot_grid(cold_dur_map, start.end_all_plot, ncol=1, 
                align = 'h', axis = 'b',labels = c("A", "B"), label_size = 30, hjust=0)

jpeg('Figure 1.jpg', width=8, height=9.5, units="in", res=400)
fig1
dev.off()

###################
#####FIGURE 2s#####
###################

#"Flip plots"

#make decade a factor in metfin
metfin$decade <- as.factor(metfin$decade)

#create indices for frostdays and snowcovered days
metfin$frostday <- ifelse(metfin$TMAXfin < 0, 1, 0)
metfin$thawday <- ifelse(metfin$TMAXfin > 0,1,0)
metfin$modSCD <- ifelse(metfin$modSWE > 1, 1, 0)
metfin$modBGD <- ifelse(metfin$modSWE == 0, 1, 0)


#format data for flip plots
flip_plot <- metfin %>% group_by(doy2, decade, ZONE) %>% 
  summarise(thawday=mean(thawday, na.rm=T), modSCD=mean(modSCD, na.rm=T), 
            modBGD=mean(modBGD, na.rm=T), frostday=mean(frostday, na.rm=T),
            TMIN=mean(TMINfin, na.rm=T))

#plot flipplots
frostday_flipplot <- ggplot()+
  annotate("rect",xmin = 30, xmax = 119, ymin = -Inf, ymax = Inf, fill="grey", alpha = 0.5) +
  stat_smooth(aes(x=doy2, y=frostday, group=decade, colour=decade), 
              method = "loess", data=flip_plot, se=F, span = 0.25)+scale_color_brewer(palette = "RdBu", direction = -1)+
  theme_bw()+
  facet_wrap(~ZONE)+scale_x_continuous(breaks=c(32,93,152),
                                       labels=c("Dec", "Feb", "Apr"),name = "")+theme(legend.position = "none")+ylab("Proportion of Frost Days")

modSCD_flipplot <- ggplot()+
  annotate("rect",xmin = 30, xmax = 119, ymin = -Inf, ymax = Inf, fill="grey", alpha = 0.5) +
  stat_smooth(aes(x=doy2, y=modSCD, group=decade, colour=decade), 
              method = "loess", data=flip_plot, se=F, span = 0.25)+scale_color_brewer(palette = "RdBu", direction = -1, name="Decade")+
  theme_bw()+
  facet_wrap(~ZONE)+scale_x_continuous(breaks=c(32,93,152),
                                       labels=c("Dec", "Feb", "Apr"),name = "")+ylab("Proportion of \n Snow Covered Days")

#Kernel density plots

#add columns to identify climate normals
calcfin$Norm = ifelse(calcfin$WYear < 1921 | calcfin$WYear > 2010, NA, ifelse(calcfin$WYear < 1951, "1921 to 1950",
                                                                              ifelse(calcfin$WYear < 1981, "1951 to 1980", "1981 to 2010")))

#subset to only include data from the three climate normals

calcfin_sub = calcfin[calcfin$WYear > 1920 & calcfin$WYear < 2011, ]

#Draw plots for max frost and SCD durations using climate normal x region combinations

fdp  <-  ggplot(calcfin_sub, aes(x = frost.dur.max_dor)) +
  geom_density(aes(fill = Norm), alpha=0.7) +
  labs(y = "Density", fill = "Climate Normal") +
  labs(x = "Maximum Frost Duration (days)") + 
  facet_wrap( ~ ZONE)+
  scale_fill_brewer(palette = "RdBu", direction = -1)+theme_bw()+theme(legend.position = "none")


scdp  <-  ggplot(calcfin_sub, aes(x = modSCD.dur.max_dor)) +
  geom_density(aes(fill = Norm), alpha=0.7) +
  labs(y = "Density", fill = "Climate Normal") +
  labs(x = "Maximum Snow Cover Duration (days)")+
  facet_wrap( ~ ZONE)+
  scale_fill_brewer(palette = "RdBu", direction = -1)+theme_bw()+theme(legend.position = "none")

#align plots with flipplots to make Figure 2
plots_l <- align_plots(frostday_flipplot, fdp, align = 'v', axis = 'l')
plots_r <- align_plots(modSCD_flipplot, scdp, align = 'v', axis = 'l')


jpeg('Figure 2.jpeg', width = 4800, height = 3000, res=360)
plot_grid(plots_l[[2]], plots_r[[2]], plots_l[[1]], plots_r[[1]], rel_widths = c(1, 1.333), labels="AUTO")
dev.off()

#####################
######FIGURE 3#######
#####################

#"Violin plots"

#format data for violin plots
plotdata_violin <- plotdata_split %>% group_by(metric.name, ZONE, Site_ID, period, bp) %>% summarize(slope=mean(slope))

#add column designating different definitions of winter
plotdata_violin$period <- ifelse(plotdata_violin$period=="win", "FGD", 
                             ifelse(plotdata_violin$period=="dor", "DOR",
                                    ifelse(plotdata_violin$period=="cal", "DJF", "AST")))

#divide data into significant and non-significant trends
plotdata_violin_sig <- plotdata_violin[which(plotdata_violin$bp!="ns"),]
plotdata_violin_nonsig <- plotdata_violin[which(plotdata_violin$bp=="ns"),]


####frost days plot

#create new dataframes with just count.frost data
count.frost_sig <- plotdata_violin_sig[which(plotdata_violin_sig$metric.name=="count.frost"),]
count.frost_nonsig <- plotdata_violin_nonsig[which(plotdata_violin_nonsig$metric.name=="count.frost"),]
count.frost_sig$bp <- as.factor(count.frost_sig$bp)

#plot trends in frost days
count.frost.violinplot <- ggplot()+
  geom_violin(data=count.frost_sig, aes(period, slope*10))+theme_bw()+
  geom_jitter(data=count.frost_sig, aes(period, slope*10, colour=bp), width = 0.15)+
  scale_color_manual(values = c("#d65555","#595bcc"))+
  geom_jitter(data=count.frost_nonsig, aes(period, slope*10), colour="grey50", width = 0.15, alpha=0.5)+
  #ggtitle("Frost Days")+
  scale_y_continuous(limits=c(-4.5, 4.5))+
  coord_flip()+theme(legend.position = "none", text=element_text(size=15))+ylab("Rate of Change (# Frost Days days/decade)")+xlab("")

####snowcovered days plot

#create new dataframes with just count.modSCD data
count.modSCD_sig <- plotdata_violin_sig[which(plotdata_violin_sig$metric.name=="count.modSCD"),]
count.modSCD_nonsig <- plotdata_violin_nonsig[which(plotdata_violin_nonsig$metric.name=="count.modSCD"),]
count.modSCD_sig$bp <- as.factor(count.modSCD_sig$bp)

#plot trends in snowcovered days
count.modSCD.violinplot <- ggplot()+
  geom_violin(data=count.modSCD_sig, aes(period, slope*10))+theme_bw()+
  geom_jitter(data=count.modSCD_sig, aes(period, slope*10, colour=bp), width = 0.15)+
  scale_color_manual(values = c("#d65555","#595bcc"))+
  geom_jitter(data=count.modSCD_nonsig, aes(period, slope*10), colour="grey50", width = 0.15, alpha=0.5)+
  #ggtitle("Snow Covered Days")+
  scale_y_continuous(limits=c(-4.5, 4.5))+
  coord_flip()+theme(legend.position = "none", text=element_text(size=15))+ylab("Rate of Change (days/decade)")+xlab("")

####minimum temperature plot

#create new dataframes with just TMIN data
TMIN_sig <- plotdata_violin_sig[which(plotdata_violin_sig$metric.name=="TMIN"),]
TMIN_nonsig <- plotdata_violin_nonsig[which(plotdata_violin_nonsig$metric.name=="TMIN"),]
TMIN_sig$bp <- as.factor(TMIN_sig$bp)

#plot trends in minimum temperatures
TMIN.violinplot <- ggplot()+
  geom_violin(data=TMIN_sig, aes(period, slope*10))+theme_bw()+
  geom_jitter(data=TMIN_sig, aes(period, slope*10, colour=bp), width = 0.15)+
  scale_color_manual(values = c("#595bcc","#d65555"))+
  geom_jitter(data=TMIN_nonsig, aes(period, slope*10), colour="grey50", width = 0.15, alpha=0.5)+
  #ggtitle("Minimum Temperature")+
  scale_y_continuous(limits=c(-0.45, 0.45))+
  coord_flip()+theme(legend.position = "none", text=element_text(size=15))+ylab(expression('Rate of Change ('*~degree*C*' /decade)'))+xlab("")

#join plots together
three_violin_plots <- plot_grid(TMIN.violinplot, count.frost.violinplot, count.modSCD.violinplot, ncol=1, labels = "AUTO")

jpeg('Figure 3.jpeg', width = 3000, height = 3500, res=360)
three_violin_plots
dev.off()

###########################
#####Loess predictions#####
###########################

#These aren't figures but are used in the text to extract values

models  <-  flip_plot %>% 
  group_by(ZONE, decade) %>% 
  do(model = loess(modSCD ~ doy2, data = ., span=0.25)) %>%
  ungroup()
flip_plot_test  <-  left_join(tbl_df(flip_plot ), models, by = c("ZONE", "decade"))

flip_plot_pred <- flip_plot_test %>%
  group_by(ZONE, decade) %>%
  do(modelr::add_predictions(., first(.$model)))




fpw <- flip_plot[which(flip_plot$ZONE=="West"),]
fpc <- flip_plot[which(flip_plot$ZONE=="Central"),]
fpe <- flip_plot[which(flip_plot$ZONE=="East"),]

fpw1910 <- fpw[which(fpw$decade==1910),]
fpw1920 <- fpw[which(fpw$decade==1920),]
fpw1930 <- fpw[which(fpw$decade==1930),]
fpw1940 <- fpw[which(fpw$decade==1940),]
fpw1950 <- fpw[which(fpw$decade==1950),]
fpw1960 <- fpw[which(fpw$decade==1960),]
fpw1970 <- fpw[which(fpw$decade==1970),]
fpw1980 <- fpw[which(fpw$decade==1980),]
fpw1990 <- fpw[which(fpw$decade==1990),]
fpw2000 <- fpw[which(fpw$decade==2000),]
fpw2010 <- fpw[which(fpw$decade==2010),]


fpc1910 <- fpc[which(fpc$decade==1910),]
fpc1920 <- fpc[which(fpc$decade==1920),]
fpc1930 <- fpc[which(fpc$decade==1930),]
fpc1940 <- fpc[which(fpc$decade==1940),]
fpc1950 <- fpc[which(fpc$decade==1950),]
fpc1960 <- fpc[which(fpc$decade==1960),]
fpc1970 <- fpc[which(fpc$decade==1970),]
fpc1980 <- fpc[which(fpc$decade==1980),]
fpc1990 <- fpc[which(fpc$decade==1990),]
fpc2000 <- fpc[which(fpc$decade==2000),]
fpc2010 <- fpc[which(fpc$decade==2010),]

fpe1910 <- fpe[which(fpe$decade==1910),]
fpe1920 <- fpe[which(fpe$decade==1920),]
fpe1930 <- fpe[which(fpe$decade==1930),]
fpe1940 <- fpe[which(fpe$decade==1940),]
fpe1950 <- fpe[which(fpe$decade==1950),]
fpe1960 <- fpe[which(fpe$decade==1960),]
fpe1970 <- fpe[which(fpe$decade==1970),]
fpe1980 <- fpe[which(fpe$decade==1980),]
fpe1990 <- fpe[which(fpe$decade==1990),]
fpe2000 <- fpe[which(fpe$decade==2000),]
fpe2010 <- fpe[which(fpe$decade==2010),]

fpw1910loess_SCD <- predict(loess(fpw1910$modSCD~fpw1910$doy2, span = 0.25), fpw1910$doy2)
fpw1920loess_SCD <- predict(loess(fpw1920$modSCD~fpw1920$doy2, span = 0.25), fpw1920$doy2)
fpw1930loess_SCD <- predict(loess(fpw1930$modSCD~fpw1930$doy2, span = 0.25), fpw1930$doy2)
fpw1940loess_SCD <- predict(loess(fpw1940$modSCD~fpw1940$doy2, span = 0.25), fpw1940$doy2)
fpw1950loess_SCD <- predict(loess(fpw1950$modSCD~fpw1950$doy2, span = 0.25), fpw1950$doy2)
fpw1960loess_SCD <- predict(loess(fpw1960$modSCD~fpw1960$doy2, span = 0.25), fpw1960$doy2)
fpw1970loess_SCD <- predict(loess(fpw1970$modSCD~fpw1970$doy2, span = 0.25), fpw1970$doy2)
fpw1980loess_SCD <- predict(loess(fpw1980$modSCD~fpw1980$doy2, span = 0.25), fpw1980$doy2)
fpw1990loess_SCD <- predict(loess(fpw1990$modSCD~fpw1990$doy2, span = 0.25), fpw1990$doy2)
fpw2000loess_SCD <- predict(loess(fpw2000$modSCD~fpw2000$doy2, span = 0.25), fpw2000$doy2)
fpw2010loess_SCD <- predict(loess(fpw2010$modSCD~fpw2010$doy2, span = 0.25), fpw2010$doy2)

fpc1910loess_SCD <- predict(loess(fpc1910$modSCD~fpc1910$doy2, span = 0.25), fpc1910$doy2)
fpc1920loess_SCD <- predict(loess(fpc1920$modSCD~fpc1920$doy2, span = 0.25), fpc1920$doy2)
fpc1930loess_SCD <- predict(loess(fpc1930$modSCD~fpc1930$doy2, span = 0.25), fpc1930$doy2)
fpc1940loess_SCD <- predict(loess(fpc1940$modSCD~fpc1940$doy2, span = 0.25), fpc1940$doy2)
fpc1950loess_SCD <- predict(loess(fpc1950$modSCD~fpc1950$doy2, span = 0.25), fpc1950$doy2)
fpc1960loess_SCD <- predict(loess(fpc1960$modSCD~fpc1960$doy2, span = 0.25), fpc1960$doy2)
fpc1970loess_SCD <- predict(loess(fpc1970$modSCD~fpc1970$doy2, span = 0.25), fpc1970$doy2)
fpc1980loess_SCD <- predict(loess(fpc1980$modSCD~fpc1980$doy2, span = 0.25), fpc1980$doy2)
fpc1990loess_SCD <- predict(loess(fpc1990$modSCD~fpc1990$doy2, span = 0.25), fpc1990$doy2)
fpc2000loess_SCD <- predict(loess(fpc2000$modSCD~fpc2000$doy2, span = 0.25), fpc2000$doy2)
fpc2010loess_SCD <- predict(loess(fpc2010$modSCD~fpc2010$doy2, span = 0.25), fpc2010$doy2)

fpe1910loess_SCD <- predict(loess(fpe1910$modSCD~fpe1910$doy2, span = 0.25), fpe1910$doy2)
fpe1920loess_SCD <- predict(loess(fpe1920$modSCD~fpe1920$doy2, span = 0.25), fpe1920$doy2)
fpe1930loess_SCD <- predict(loess(fpe1930$modSCD~fpe1930$doy2, span = 0.25), fpe1930$doy2)
fpe1940loess_SCD <- predict(loess(fpe1940$modSCD~fpe1940$doy2, span = 0.25), fpe1940$doy2)
fpe1950loess_SCD <- predict(loess(fpe1950$modSCD~fpe1950$doy2, span = 0.25), fpe1950$doy2)
fpe1960loess_SCD <- predict(loess(fpe1960$modSCD~fpe1960$doy2, span = 0.25), fpe1960$doy2)
fpe1970loess_SCD <- predict(loess(fpe1970$modSCD~fpe1970$doy2, span = 0.25), fpe1970$doy2)
fpe1980loess_SCD <- predict(loess(fpe1980$modSCD~fpe1980$doy2, span = 0.25), fpe1980$doy2)
fpe1990loess_SCD <- predict(loess(fpe1990$modSCD~fpe1990$doy2, span = 0.25), fpe1990$doy2)
fpe2000loess_SCD <- predict(loess(fpe2000$modSCD~fpe2000$doy2, span = 0.25), fpe2000$doy2)
fpe2010loess_SCD <- predict(loess(fpe2010$modSCD~fpe2010$doy2, span = 0.25), fpe2010$doy2)


fpw1910loess_frost <- predict(loess(fpw1910$frostday~fpw1910$doy2, span = 0.25), fpw1910$doy2)
fpw1920loess_frost <- predict(loess(fpw1920$frostday~fpw1920$doy2, span = 0.25), fpw1920$doy2)
fpw1930loess_frost <- predict(loess(fpw1930$frostday~fpw1930$doy2, span = 0.25), fpw1930$doy2)
fpw1940loess_frost <- predict(loess(fpw1940$frostday~fpw1940$doy2, span = 0.25), fpw1940$doy2)
fpw1950loess_frost <- predict(loess(fpw1950$frostday~fpw1950$doy2, span = 0.25), fpw1950$doy2)
fpw1960loess_frost <- predict(loess(fpw1960$frostday~fpw1960$doy2, span = 0.25), fpw1960$doy2)
fpw1970loess_frost <- predict(loess(fpw1970$frostday~fpw1970$doy2, span = 0.25), fpw1970$doy2)
fpw1980loess_frost <- predict(loess(fpw1980$frostday~fpw1980$doy2, span = 0.25), fpw1980$doy2)
fpw1990loess_frost <- predict(loess(fpw1990$frostday~fpw1990$doy2, span = 0.25), fpw1990$doy2)
fpw2000loess_frost <- predict(loess(fpw2000$frostday~fpw2000$doy2, span = 0.25), fpw2000$doy2)
fpw2010loess_frost <- predict(loess(fpw2010$frostday~fpw2010$doy2, span = 0.25), fpw2010$doy2)

fpc1910loess_frost <- predict(loess(fpc1910$frostday~fpc1910$doy2, span = 0.25), fpc1910$doy2)
fpc1920loess_frost <- predict(loess(fpc1920$frostday~fpc1920$doy2, span = 0.25), fpc1920$doy2)
fpc1930loess_frost <- predict(loess(fpc1930$frostday~fpc1930$doy2, span = 0.25), fpc1930$doy2)
fpc1940loess_frost <- predict(loess(fpc1940$frostday~fpc1940$doy2, span = 0.25), fpc1940$doy2)
fpc1950loess_frost <- predict(loess(fpc1950$frostday~fpc1950$doy2, span = 0.25), fpc1950$doy2)
fpc1960loess_frost <- predict(loess(fpc1960$frostday~fpc1960$doy2, span = 0.25), fpc1960$doy2)
fpc1970loess_frost <- predict(loess(fpc1970$frostday~fpc1970$doy2, span = 0.25), fpc1970$doy2)
fpc1980loess_frost <- predict(loess(fpc1980$frostday~fpc1980$doy2, span = 0.25), fpc1980$doy2)
fpc1990loess_frost <- predict(loess(fpc1990$frostday~fpc1990$doy2, span = 0.25), fpc1990$doy2)
fpc2000loess_frost <- predict(loess(fpc2000$frostday~fpc2000$doy2, span = 0.25), fpc2000$doy2)
fpc2010loess_frost <- predict(loess(fpc2010$frostday~fpc2010$doy2, span = 0.25), fpc2010$doy2)

fpe1910loess_frost <- predict(loess(fpe1910$frostday~fpe1910$doy2, span = 0.25), fpe1910$doy2)
fpe1920loess_frost <- predict(loess(fpe1920$frostday~fpe1920$doy2, span = 0.25), fpe1920$doy2)
fpe1930loess_frost <- predict(loess(fpe1930$frostday~fpe1930$doy2, span = 0.25), fpe1930$doy2)
fpe1940loess_frost <- predict(loess(fpe1940$frostday~fpe1940$doy2, span = 0.25), fpe1940$doy2)
fpe1950loess_frost <- predict(loess(fpe1950$frostday~fpe1950$doy2, span = 0.25), fpe1950$doy2)
fpe1960loess_frost <- predict(loess(fpe1960$frostday~fpe1960$doy2, span = 0.25), fpe1960$doy2)
fpe1970loess_frost <- predict(loess(fpe1970$frostday~fpe1970$doy2, span = 0.25), fpe1970$doy2)
fpe1980loess_frost <- predict(loess(fpe1980$frostday~fpe1980$doy2, span = 0.25), fpe1980$doy2)
fpe1990loess_frost <- predict(loess(fpe1990$frostday~fpe1990$doy2, span = 0.25), fpe1990$doy2)
fpe2000loess_frost <- predict(loess(fpe2000$frostday~fpe2000$doy2, span = 0.25), fpe2000$doy2)
fpe2010loess_frost <- predict(loess(fpe2010$frostday~fpe2010$doy2, span = 0.25), fpe2010$doy2)


fpeloess_frost <- cbind.data.frame(fpe1910loess_frost,
                                 fpe1920loess_frost,
                                 fpe1930loess_frost,
                                 fpe1940loess_frost,
                                 fpe1950loess_frost,
                                 fpe1960loess_frost,
                                 fpe1970loess_frost,
                                 fpe1980loess_frost,
                                 fpe1990loess_frost,
                                 fpe2000loess_frost,
                                 fpe2010loess_frost)

fpcloess_frost <- cbind.data.frame(fpc1910loess_frost,
                                 fpc1920loess_frost,
                                 fpc1930loess_frost,
                                 fpc1940loess_frost,
                                 fpc1950loess_frost,
                                 fpc1960loess_frost,
                                 fpc1970loess_frost,
                                 fpc1980loess_frost,
                                 fpc1990loess_frost,
                                 fpc2000loess_frost,
                                 fpc2010loess_frost)


fpwloess_frost <- cbind.data.frame(fpw1910loess_frost,
                                 fpw1920loess_frost,
                                 fpw1930loess_frost,
                                 fpw1940loess_frost,
                                 fpw1950loess_frost,
                                 fpw1960loess_frost,
                                 fpw1970loess_frost,
                                 fpw1980loess_frost,
                                 fpw1990loess_frost,
                                 fpw2000loess_frost,
                                 fpw2010loess_frost)

fpeloess_SCD <- cbind.data.frame(fpe1910loess_SCD,
                               fpe1920loess_SCD,
                               fpe1930loess_SCD,
                               fpe1940loess_SCD,
                               fpe1950loess_SCD,
                               fpe1960loess_SCD,
                               fpe1970loess_SCD,
                               fpe1980loess_SCD,
                               fpe1990loess_SCD,
                               fpe2000loess_SCD,
                               fpe2010loess_SCD)

fpcloess_SCD <- cbind.data.frame(fpc1910loess_SCD,
                               fpc1920loess_SCD,
                               fpc1930loess_SCD,
                               fpc1940loess_SCD,
                               fpc1950loess_SCD,
                               fpc1960loess_SCD,
                               fpc1970loess_SCD,
                               fpc1980loess_SCD,
                               fpc1990loess_SCD,
                               fpc2000loess_SCD,
                               fpc2010loess_SCD)


fpwloess_SCD <- cbind.data.frame(fpw1910loess_SCD,
                               fpw1920loess_SCD,
                               fpw1930loess_SCD,
                               fpw1940loess_SCD,
                               fpw1950loess_SCD,
                               fpw1960loess_SCD,
                               fpw1970loess_SCD,
                               fpw1980loess_SCD,
                               fpw1990loess_SCD,
                               fpw2000loess_SCD,
                               fpw2010loess_SCD)

