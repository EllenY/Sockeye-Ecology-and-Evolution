###################################################
#   Juvenile sockeye salmon project
###################################################

#Author: Ellen Yasumiishi
#Address: 17109 Point Lena Loop Road, Juneau, AK 99801

# Load packages
library(cowplot)
library(googleway)
library(ggrepel)
library(ggspatial)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(splines2)
library(splines)
library(fields)
library(mapdata)
library(maps)
library(akima)
library(stats)
library(car)
library(ggplot2)
library(grid)
library(gridBase)
library(gridExtra)
library(reshape2)
library(Matrix)
library(INLA)
library(effects)

library(VAST)
library(dplyr)
library(magrittr)
library(viridis)
library(plotrix)
library(ggmap)
library(maps)
library(devtools)
library(units)

#Help files
#?FishStatsUtils::make_settings
#?FishStatsUtils::fit_model
#?VAST::make_data
#?FishStatsUtils::plot_results

gc() 
gc(reset=T)
rm(list = ls())

Version = get_latest_version( package="VAST" )
Version
wd<-setwd("C:/Users/ellen.yasumiishi/Work/GitHub/Sockeye-Ecology-and-Evolution")

#Data set to estimate annual indices
#Simply removed missing years
temp<-read.csv("IYSsockeyeCovariates30.csv") 
Year=temp[,1]	      #Survey year
Lat	=temp[,2]       #Latitude
Lon	=temp[,3]       #Longitude
SST	=temp[,4]       #Sea temperature at 20 m depth
A0PollockN =temp[,5]#Catch (numbers) of age-0 pollock	
J_Pink =temp[,6]	  #Catch (kg) of juvenile pink salmon
Calanus	=temp[,7]   #Calanus copepod density
AirT	=temp[,8]     #Air temperature not used
HaulDate	=temp[,9] #Date of surface tow
StationID	=temp[,10]#Station identification
AreaSwept	=temp[,11]#Are swept (km^2) of the surface trawl net tow
kg	=temp[,12]      #Catch (kg) of juvenile sockeye salmon
Sci	=temp[,13]      #Species name
vessel	=temp[,14]  #Vessel ID
Julian	=temp[,15]  #Julian day of the sample
Constant=temp[,16]  #Constant used for Calanus
Lat_rnd=temp[,17]   #Rounded Latitude 
Lon_rnd=temp[,18]   #Rounded Longitude


# Exploratory Plots ========================================================================
do.explore.plot <- TRUE
if(do.explore.plot==TRUE) {
  
  # Catch Distribution by year
  
  g <- ggplot(temp, aes(kg)) +
    theme_linedraw() +
    geom_histogram(fill='blue', color='black') +
    facet_wrap(~factor(Year), scales='fixed') +
    xlab("Weight of Sockeye catch")
  g
  ggsave(file.path(wd, "Wt Distribution.png"), plot=g, height=8, width=9, units='in')

  # Determine proportion of observations with zero catch
  sum.ec <- temp %>% group_by(Year) %>% summarize(zero=sum(SST==0), nonZero=sum(SST>0), 
                                                  n=n(),
                                                  prop.zero=sum(SST==0)/n())
  sum.ec
  write.csv(sum.ec, file=file.path(wd, "Summary of Encounter Prob.csv"))
  
  # Effort Distribution
  g.eff <- ggplot(temp, aes(AreaSwept)) +
    theme_linedraw() +
    geom_histogram(fill='blue', color='black') +
    facet_wrap(~factor(Year), scales='free') +
    xlab("Area Swept")
   g.eff
  ggsave(file.path(wd, "Effort Distribution.png"), plot=g.eff, height=8, width=9, units='in')

    # Plot Map of Catch Rates
  CPU<-log((kg/AreaSwept)+1)
    gmap<-ggplot(temp, aes(x=Lon, y=Lat, color=CPU)) +
    theme_linedraw() +
    scale_color_viridis() +
    geom_point(size=2) +
    facet_wrap(~factor(Year))+
    #labs(title = "Calanus", x = "Longitude", y = "Latitude", color = "log((#/m^2)+1)") +
    labs(title = "Sea temperature (20m depth)", x = "Longitude", y = "Latitude", color = "Celsius") +
    labs(title = "Juvenile sockeye salmon", x = "Longitude", y = "Latitude", color = "log((kg/km^2)+1)") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 20, face = "bold", color = "black"),
          legend.text=element_text(size=12),strip.text.x = element_text(size = 12, face='bold')) 

  ggsave(file.path(wd, "Map of Wt.png"), plot=gmap, height=8, width=9, units='in')
  
  #Mean annual sea temperatures at 20m
  mean.SST <- temp %>% group_by(Year) %>% dplyr::summarize(mean_SST=mean(SST))
  length(SST)
  mean.SST
  
    # Plot Number of Hauls per year
  haul.count <- temp %>% group_by(Year) %>% summarize(n=n())

    g.haul <- ggplot(haul.count, aes(x=Year, y=n, fill=n)) +
    theme_bw() +
    scale_fill_viridis() +
    geom_bar(stat='identity', color='black') +
    ylab('Number of Survey Hauls')
   g.haul
  ggsave(file.path(wd, "Number of Hauls by Year.png"), plot=g.haul, height=4, width=5, units='in')
  
}


########################################
# FIGURE 3   MAP OF STUDY AREA 
########################################

# install.packages("devtools")
devtools::install_github("MattCallahan-NOAA/AKmarineareas", force=TRUE, lib="C:/Program Files/R/R-4.1.3/library")
install.packages("grid", force=TRUE, lib="C:/Program Files/R/R-4.1.3/library")

install.packages("AKmarineareas")
install.packages("remotes")
install.packages("marmap")
install.packages('raster')
remove.packages('marmap')
install.packages('marmap')

library(ggplot2)
library(maps)
library(mapdata)
library(mapproj)
library(grid)
library(AKmarineareas)
library(sf)
library(tidyverse)
library(raster)
library(marmap)
library(dplyr)
library(sp)
library(ggspatial) # Add a scale bar


theme_bw()
#get russia
russia<-rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")%>%
  filter(name=="Russia")
#AK coastlines
AK<-AK_basemap()

#bathymetry data
#resolution ranges from 1 (finest) to 10 (most coarse)
depth <- marmap::as.raster(getNOAA.bathy(lon1=-180,lon2=-129,lat1=50,lat2=72, resolution=10, keep=TRUE, antimeridian=FALSE, path=NULL))
getNOAA.bathy(lon1=-180,lon2=-129,lat1=50,lat2=72, resolution=10, keep=TRUE, antimeridian=FALSE, path=NULL)

#convert to contour
depth_c<-rasterToContour(depth, levels=c(-50,-100,-200))%>%
  st_as_sf()
grid.newpage();grid.draw(roundrectGrob(gp = gpar(lwd = 0)))

pdf(file="FIGURE 3.pdf", height=11, width=8.5) #This generates the figure as a hi-res (pub quality) tiff file. You can change the dimensions of the figure by changing the height and width arguments. You can change the resolution with the 'res' argument.
p<-ggplot()+
  geom_sf(data=russia)+
  geom_sf(data=AK)+
  geom_sf(data=depth_c, aes(color=level), colour=c("grey35","grey50","grey75"))+
  geom_point(data=temp, aes(Lon_rnd, Lat_rnd, group=NULL), shape=4)+
  coord_sf(xlim=c(-178, -150), ylim=c(50,70))+
  theme(panel.background = element_rect(colour = 'black',fill = 'white' ),panel.border = element_rect(colour = "black",fill=NA, size=1))+
  xlab(expression(paste(Longitude)))+
  ylab(expression(paste(Latitude)))+
  geom_text(aes(x=-169.3, y=64, label="St. Lawrence Is.",angle=0), size=4)+
  geom_text(aes(x=-172.9, y=61, label="St. Matthews Is.",angle=0), size=4)+
  geom_text(aes(x=-164.2, y=54.73, label="Unimak Pass.",angle=0), size=4)+
  geom_text(aes(x=-167.1, y=60.6, label="Nunivak Is.",angle=0), size=4)+
  geom_label(aes(x = -169.7, y =60.1, label = "50 m"), fill = "white", label.size = NA, size = 3) +
  geom_label(aes(x = -174.1, y =60.1, label = "100 m"), fill = "white", label.size = NA, size = 3) +
  geom_label(aes(x = -178.1, y =60.1, label = "200 m"), fill = "white", label.size = NA, size = 3) +
  geom_text(aes(x=-168.2, y=62, label="Inner",angle=0), size=5)+
  geom_text(aes(x=-174, y=62, label="Middle",angle=0), size=5)+
  geom_text(aes(x=-177.8, y=62, label="Outer",angle=0), size=5)+
  geom_text(aes(x=-155, y=66, label="Alaska"), size=7)+
  geom_text(aes(x=-176, y=66, label="Russia"), size=7)+
  geom_text(aes(x=-175.1, y=55, label="Bering Sea"), size=7)+
  geom_text(aes(x=-168.6, y=66.25, label="Bering"), size=5)+
  geom_text(aes(x=-169, y=65.9, label="Strait"), size=5)+
  geom_text(aes(x=-168, y=70, label="Chukchi Sea"), size=6)+
  geom_text(aes(x=-155, y=53, label="Gulf of Alaska"), size=6)+
  geom_text(aes(x=-173, y=53, label="Aleutian Islands"), angle=20,size=6)+
  geom_text(aes(x=-165, y=51, label="Pacific Ocean"), size=7)+
  annotation_scale()+
  annotation_north_arrow( height= unit(.75, "cm"), width= unit(.75, "cm"), pad_x = unit(.5, "cm"),
                          pad_y = unit(1, "cm"))

#Alaska Stream
p+theme(legend.position = "none")+
  geom_segment(
    aes(x = -150, y = 56.5, xend = -171, yend = 51.2),
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed", # Describes arrow head (open or closed),
      ends="last"
    ),
    colour = "#999999",
    linewidth = 1,
    angle = 10 # Anything other than 90 or 0 can look unusual
  )+
  geom_text(aes(x=-163.9, y=52.6, label="Alaska Stream"),colour="#999999", angle=25.5,size=5)+
  
  #Alaska Coastal Current
  theme(legend.position = "none")+
  geom_segment(
    aes(x = -156, y = 55.7, xend = -165, yend = 53.4),
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed", # Describes arrow head (open or closed),
      ends="last"
    ),
    colour = "#999999",
    linewidth = 1,
    angle = 10 # Anything other than 90 or 0 can look unusual
  )+
  geom_text(aes(x=-159.8, y=54.4, label="Alaska Coastal Current"),colour="#999999", angle=26,size=5)+
  
  #Aleutian North Slope Current
  theme(legend.position = "none")+
  geom_segment(
    aes(x = -178.5, y = 52.4, xend = -167, yend = 54.75),
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed", # Describes arrow head (open or closed),
      ends="last"
    ),
    colour = "#999999",
    linewidth = 1,
    angle = 10 # Anything other than 90 or 0 can look unusual
  )+
  geom_text(aes(x=-173.9, y=53.8, label="Aleutian North Slope Current"),colour="#999999", angle=22,size=5)+
  
  #Bering Slope Current
  theme(legend.position = "none")+
  geom_segment(
    aes(x = -169.5, y = 55, xend = -178, yend = 59.5),
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed", # Describes arrow head (open or closed),
      ends="last"
    ),
    colour = "#999999",
    linewidth = 1,
    angle = 10 # Anything other than 90 or 0 can look unusual
  )+
  geom_text(aes(x=-174.8, y=57, label="Bering Slope Current"),colour="#999999", angle=315,size=5)+
  
  #Anadyr Current
  theme(legend.position = "none")+
  geom_segment(
    aes(x = -178.6, y = 64.9, xend = -173, yend = 63.5),
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed", # Describes arrow head (open or closed),
      ends="last"
    ),
    colour = "#999999",
    linewidth = 1,
    angle = 10 # Anything other than 90 or 0 can look unusual
  )+
  geom_text(aes(x=-176.5, y=63.8, label="Anadyr Current"),colour="#999999", angle=331,size=5)+
  
  #Bering Strait Current
  theme(legend.position = "none")+
  geom_segment(
    aes(x = -169, y = 64.5, xend = -169, yend = 68),
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed", # Describes arrow head (open or closed),
      ends="last"
    ),
    colour = "#999999",
    linewidth = 1,
    angle = 10 # Anything other than 90 or 0 can look unusual
  )+
  geom_text(aes(x=-168.6, y=66.25, label="Bering"), size=5)+
  geom_text(aes(x=-169, y=65.9, label="Strait"), size=5)+
  
  #Bering  Current
  theme(legend.position = "none")+
  geom_segment(
    aes(x = -163, y = 57.5, xend = -168, yend = 60),
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed", # Describes arrow head (open or closed),
      ends="last"
    ),
    colour = "#999999",
    linewidth = 1,
    angle = 10 # Anything other than 90 or 0 can look unusual
  )+
  
  #Bering Strait Current
  theme(legend.position = "none")+
  geom_segment(
    aes(x = -171.2, y = 51.6, xend = -172.1, yend = 52.8),
    arrow = arrow(
      length = unit(0.03, "npc"), 
      type="closed", # Describes arrow head (open or closed),
      ends="last"
    ),
    colour = "#999999",
    linewidth = 1,
    angle = 10 # Anything other than 90 or 0 can look unusual
  )
dev.off()


######################################################
# FIGURE 4 Base juvenile sockeye Spatiotemporal Model
######################################################
#Data:2002-2018: No sampling in 2013, 2015, 2017, so data are averaged for the VAST model to run.
temp<-read.csv("IYSsockeyeCovariates29.csv") 
Year=temp[,1]	        #Smple year
Lat	=temp[,2]         #Latitude
Lon	=temp[,3]         #Longitude
SST	=temp[,4]         #sea temperature at 20m depth at the station
A0PollockN =temp[,5]	# number of age-0 pollock caught
J_Pink =temp[,6]	    # kg of juvenile pink salmon
Calanus	=temp[,7]     #no./m^2 copepods
AirT	=temp[,8]       #Air temperature
HaulDate	=temp[,9]   #haul ID
StationID	=temp[,10]  #Station ID
AreaSwept	=temp[,11]  #Area swept
kg	=temp[,12]        #juvenile sockeye catch  kg
Sci	=temp[,13]        #species
vessel	=temp[,14]    #vessel
Julian	=temp[,15]    #Julina day of haul
Constant=temp[,16]    # constant 1

settings = make_settings( Version="VAST_v14_0_1",  
                          n_x=500,
                          Region='User', 
                          purpose="index2",
                          bias.correct=TRUE,
                          max_cells = 2000,
                          ObsModel=c(2,1), # catch, Pollock, pink
                          fine_scale=TRUE,
                          treat_nonencounter_as_zero=TRUE,
                          knot_method='grid',
                          use_anisotropy = TRUE)

settings$Options = c( settings$Options, "report_additional_variables"=TRUE )
user_region <- readRDS('user_region.rds')

fit_Spatiotemporal = fit_model( "settings"=settings,
                                "Lat_i"=temp[,'Lat'],
                                "Lon_i"=temp[,'Lon'],
                                "observations_LL" = temp[,c('Lat','Lon')],
                                "t_i"=temp[,'Year'], 
                                "b_i"=(as_units(temp[,'kg'],'kg')), 
                                "a_i"=(as_units(temp[,'AreaSwept'],'km^2')),
                                "getsd"=TRUE,
                                "test_fit"=TRUE,
                                "build_model" = TRUE,
                                input_grid=user_region ,
                                getReportCovariance=TRUE
)

check_fit(parameter_estimates, check_gradients = FALSE, quiet = FALSE)
#Check extrapolation grid total area....
colSums(fit_Spatiotemporal$extrapolation_list$a_el)

#SAVE DATA =================================================================================
x=summary(fit_Spatiotemporal$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"Full_sockeye.csv", row.names = TRUE)
plot_results(fit_Spatiotemporal, n_cells=2000  )


## Below shows to you get the model estimate of density, D_gct,
## for each grid (g), category (c; not used here single
## univariate); and year (t); and link it spatially to a lat/lon
## extrapolation point.  You can do this for any _gct or _gc
## variable in the Report.

years <- unique(temp$Year)
nyrs <- length(years)
#Year=rep(years, each=nrow(ak_map))
## Remake map list locally for recreating plots
mdl <- make_map_info(Region = settings$Region,
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)
## quick dirty AK map
ak_map <- subset(map_data("world"), region=='USA' & subregion=='Alaska')
## Have to duplicate it for each year so can facet below
ak_map <- cbind(ak_map[rep(1:nrow(ak_map), times=nyrs),],
                Year=rep(years, each=nrow(ak_map)))
gmap <- ggplot(ak_map, aes(x = long, y = lat, group = group)) +
  # geom_raster(aes(x=long, y=lat, fill=D))+
  # geom_raster(aes(fill = D), interpolate = TRUE)+
  geom_polygon(fill="black", colour = "white") +
  # scale_colour_gradientn(colors=c('red','orange','yellow','green', 'blue','purple')) +
  #  scale_colour_gradientn(colors=c('darkblue','blue','lightblue','lightgreen','yellow', 'orange','red')) +
  #scale_color_gradientn(colours = rainbow(5))+
  # scale_fill_viridis() + 
  scale_color_viridis_c(option = "D") +  #changed form magma to plasma or virdis
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  #  coord_cartesian(xlim=mdl$Xlim, ylim=mdl$Ylim)
  coord_cartesian(xlim=mdl$Xlim, ylim=mdl$Ylim)+ 
  guides(fill=guide_legend(title='Energy Density'))

#Densities at each extrapolation grid location
fit$Report$D_gct
#Densities at each extrapolation grid location
names(fit$Report)[grepl('_gc|_gct', x=names(fit$Report))]
D_gt <- fit$Report$D_gct[,1,] # drop the category
#dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=c(2003:2019))
dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=years)
## tidy way of doing this, reshape2::melt() does
## it cleanly but is deprecated
#?pivot_longer.sf
#?rownames_to_column
#require(tibble)
library(tidyverse)
#?pivot_longer
D_gt <- D_gt %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  #tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')
D <- merge(D_gt, mdl$PlotDF, by.x='cell', by.y='x2i')
D3<-lapply(((D)), as.numeric)
D4<-data.frame(D3)
head(D4)
summary(exp(D4$D))
summary((D4$D))

#Subet D values greater than .1% of D to plot
D02 <- D4[which((D$Year == "2002")),names(D4) %in% c("Year","D","Lat","Lon")]
D03 <- D4[which(D$Year == "2003"),names(D4) %in% c("Year","D","Lat","Lon")]
D04 <- D4[which(D$Year == "2004"),names(D4) %in% c("Year","D","Lat","Lon")]
D05 <- D4[which(D$Year == "2005"),names(D4) %in% c("Year","D","Lat","Lon")]
D06 <- D4[which(D$Year == "2006"),names(D4) %in% c("Year","D","Lat","Lon")]
D07 <- D4[which(D$Year == "2007"),names(D4) %in% c("Year","D","Lat","Lon")]
D08 <- D4[which(D$Year == "2008"),names(D4) %in% c("Year","D","Lat","Lon")]
D09 <- D4[which(D$Year == "2009"),names(D4) %in% c("Year","D","Lat","Lon")]
D10 <- D4[which(D$Year == "2010"),names(D4) %in% c("Year","D","Lat","Lon")]
D11 <- D4[which(D$Year == "2011"),names(D4) %in% c("Year","D","Lat","Lon")]
D12 <- D4[which(D$Year == "2012"),names(D4) %in% c("Year","D","Lat","Lon")]
D14 <- D4[which(D$Year == "2014"),names(D4) %in% c("Year","D","Lat","Lon")]
D16 <- D4[which(D$Year == "2016"),names(D4) %in% c("Year","D","Lat","Lon")]
D18 <- D4[which(D$Year == "2018"),names(D4) %in% c("Year","D","Lat","Lon")]

D5<-rbind(D02,D03,D04,D05,D06,D07,D08, D09,D10,D11,D12,D14,D16,D18) 
D6<-lapply(D5, as.numeric)
D7<-data.frame(D6)
summary(D7)
write.csv(D7, "SOCKEYE_OUT.csv")
theme_set(theme_bw())
D7$D<=0.001*D  #Plot the top 99.9% of the densities.

g2 <- gmap + 
  theme_bw()+
  #scale_color_gradientn(colours = colorspace::divergingx_hcl(palette="RdGy", n=7, rev=TRUE),oob=scales::squish,limits=c(0,4))+
  #geom_polygon(data = temp, aes(x = Lon, y = Lat, group = NULL), fill=NA,color="lightgray" ,size = 0.005) +
  geom_point(data=D7 , aes(Lon, Lat, color=((D)), group=NULL),  #Use log(D) for GRP to covert back from exp(GRP)
             #geom_point(data=temp2 , aes(Lon, Lat, color=(D), group=NULL),  #Use log(D) for GRP to covert back from exp(GRP)
             na.rm=TRUE, size=2, stroke=0, shape=16) +facet_wrap('Year')+ 
 #labs(color=expression(paste("ln(kg/",km^2,")")))+
  labs(color=paste('log(kg\U00B7km\u00b2+1)'))+
  xlab(expression(paste(Longitude^o,~'W')))+
  ylab(expression(paste(Latitude^o,~'N')))+
  #coord_map(xlim = c(-173, -155),ylim = c(50, 66.5))+
  labs(title="Juvenile sockeye salmon")+
  theme(strip.text.x = element_text(size = 12,face="bold"))  +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major=element_line(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"))
g2
dev.off()

####################################
#  FIGURE 5: DIET PLOT
####################################
library(ggplot2)
library(viridis)
library(hrbrthemes)  
library(tidyverse)

dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dir)
#I converted to annual proportion before reading in
diet <- read.csv("SockDietProp.csv")
#Changed prey item names to appear the way I wanted them to, instead of _ or . separators. 
diet.name <- diet %>% rename("Age-0 Pollock" = Age0pollock, "Arrow Worms" = Arrowworm, "Calanus spp." = Calanus, 
                             "Large Copepods" = Largecopepods, "Other Crustaceans" = Othercrustaceans, "Small Copepods" = Smallcopepods)
#Lon pivot to make ggplot happy. 
diet.Lon <- pivot_longer(diet.name, cols = 2:12, names_to = "PreyItem", values_to = "SCI")
#Quick easy way to remove the gaps on the plot for years with no data.
diet.Lon$Year <- as.factor(diet.Lon$Year)
#Reorder prey items to the order that I want them to appear in. I grouped fish with fish, crustaceans with crustaceans, 
#descending order within each group. 
pdf(file="FIGURE 5.pdf", height=8, width=11) 
p1 <- diet.Lon %>%
  mutate(PreyItem = fct_relevel(PreyItem,
                                "Age-0 Pollock", "Fish", "Euphausiids", 
                                "Amphipods", "Calanus spp.", "Large Copepods", "Small Copepods", 
                                "Other Crustaceans", "Arrow Worms", "Pteropods", "Other")) %>%
  ggplot(aes(Year,PreyItem,fill=SCI))+
  geom_tile(color= "white",size=0.1) + 
  scale_fill_viridis(name="%SCI",option ="C")+
  theme_minimal(base_size = 8)+
  
  theme(text = element_text(size = 12, colour="black"),axis.text=element_text(colour="black"),axis.title = element_text(colour = "black",size = 16))+
  ylab("Prey Item")+
  xlab("Year")
dev.off()


pdf(file="FIGURE 5.pdf", height=8, width=11) 
p1+theme(axis.text.y=element_text(size=16),axis.text.x=element_text(size=16),axis.title.y = element_text(size=24),axis.title.x = element_text(size=24),legend.title = element_text(size = 24),
         legend.text = element_text(size = 16),strip.background = element_rect(colour="white"),axis.ticks=element_blank()) 
dev.off()

####################################
# FIGURE 6 MODEL VALIDATION FIGURE
####################################

# FULL-SAMPLE =======================================================================
#Log transformed Predicted and observed data for Sockeye model ObsModel=c(2,3) full-sample
#log(Predicted)-exp(logSigmaM)^2/2 
#LogSima -0.10784334 for data w/o 2013, 15, 17
#LogSimg w means for 13, 15, 17 -0.2715638 
#Logsigma 27 (@,3)

settings = make_settings( Version="VAST_v14_0_1",  
                          n_x=500,
                          Region='User', 
                          purpose="index2",
                          bias.correct=TRUE,
                          max_cells = 2000,
                          ObsModel=c(2,1), # catch, Pollock, pink
                          fine_scale=TRUE,
                          treat_nonencounter_as_zero=TRUE,
                          knot_method='grid',
                          use_anisotropy = TRUE)

settings$Options = c( settings$Options, "report_additional_variables"=TRUE )
user_region <- readRDS('user_region.rds')

fit = fit_model( "settings"=settings,
                 "Lat_i"=temp[,'Lat'],
                 "Lon_i"=temp[,'Lon'],
                 "observations_LL" = temp[,c('Lat','Lon')],
                 "t_i"=temp[,'Year'], 
                 "b_i"=(as_units(temp[,'kg'],'kg')), 
                 "a_i"=(as_units(temp[,'AreaSwept'],'km^2')),
                 "getsd"=TRUE,
                 "test_fit"=TRUE,
                 "build_model" = TRUE,
                 input_grid=user_region ,
                 getReportCovariance=TRUE
)


Sockeye_Pred_in<-log(fit$Report$R2_i+1)-exp(-0.23875180)^2/2 	#Predicted in-bag-sample log, old -0.06686186
Sockeye_Obs_in<-log(temp[,'kg']+1)	  				#Observed in-bag-sample log
plot(Sockeye_Pred_in,Sockeye_Obs_in)					#Plot
summary(lm(Sockeye_Obs_in~Sockeye_Pred_in))			#Linear regresson

#Not log transformed
Sockeye_Pred_in<-fit$Report$R2_i-exp(-0.23875180)^2/2					#Predicted in-bag-sample 
Sockeye_Obs_in<-(temp[,'kg']) 						#Observed in-bag-sample 
plot(Sockeye_Pred_in,Sockeye_Obs_in)				#Plot
summary(lm(Sockeye_Obs_in~0+Sockeye_Pred_in))		#R2=0.67

#R2=0.8765
#Compare object to fit of out-bag sample
ParHat = fit$ParHat
ParHat
# Generate partitions in data
n_fold = 10
Partition_i = sample( 1:n_fold, size=nrow(temp), replace=TRUE )
prednll_f = rep(NA, n_fold )
# Loop through partitions, refitting each time with a different PredTF_i
for( fI in 1:n_fold ){
  PredTF_i = ifelse( Partition_i==fI, TRUE, FALSE )
  
  settings = make_settings( Version="VAST_v14_0_1",  
                            n_x=500,
                            Region='User', 
                            purpose="index2",
                            bias.correct=TRUE,
                            max_cells = 2000,
                            ObsModel=c(2,1), # catch, Pollock, pink
                            fine_scale=TRUE,
                            treat_nonencounter_as_zero=TRUE,
                            knot_method='grid',
                            use_anisotropy = TRUE)
 # settings$RhoConfig["Beta2"] = 3
  
  # Refit, starting at MLE, without calculating standard errors (to save time)
  fit_new = fit_model("settings"=settings,
                      "Lat_i"=temp[,'Lat'],
                      "Lon_i"=temp[,'Lon'],
                      "observations_LL" = temp[,c('Lat','Lon')],
                      "t_i"=temp[,'Year'], 
                      "b_i"=(as_units(temp[,'kg'],'kg')), 
                      "a_i"=(as_units(temp[,'AreaSwept'],"km^2")),
                      "getsd"=TRUE,
                      "test_fit"=TRUE,
                      "build_model" = TRUE,
                      input_grid=user_region,
                      getReportCovariance=TRUE,"PredTF_i"=PredTF_i, "Parameters"=ParHat)
  
  # Save fit to out-of-bag data
  prednll_f[fI] = fit_new$Report$pred_jnll
}

# Check fit to all out=of-bag data and use as metric of out-of-bag performance
sum( prednll_f )

#1832

#OUT SAMPLE ===============================================================================

#Not log transformed
Sockeye_Pred_out<-fit_new$Report$R2_i
Sockeye_Obs_out<-(temp[,'kg'])
plot(Sockeye_Pred_out, Sockeye_Obs_out)
summary(lm(Sockeye_Obs_out~Sockeye_Pred_out))			#Not logged 0.737 w fake data.

#log transformed data
#log(Predicted)-exp(logSigmaM)^2/2
Sockeye_Pred_out<-log(fit_new$Report$R2_i+1)-exp(-0.2387518 )^2/2
Sockeye_Obs_out<-log(temp[,'kg']+1)
plot(Sockeye_Pred_out, Sockeye_Obs_out)
summary(lm(Sockeye_Obs_out~Sockeye_Pred_out))			#Not logged 0.56 w fake data.

data<-data.frame(Sockeye_Obs_in,Sockeye_Pred_in,Sockeye_Obs_out,Sockeye_Pred_out) 
write.csv(data, file = "data2.csv")

data2<-read.csv("data2.csv") #Not logged in and out samples
Obs=data2[,1]
Sockeye_Obs_in=data2[,2]	
Sockeye_Pred_in=data2[,3]	
Sockeye_Obs_out=data2[,4]	
Sockeye_Pred_out=data2[,5]	

m2<-lm( Sockeye_Pred_in~Sockeye_Obs_in)
summary(m2)

m4<-lm( Sockeye_Pred_out~Sockeye_Obs_out)
summary(m4)

data3<-read.csv("data3.csv") #Not logged in and out samples
Sample=data3[,1]
Sockeye_Obs=data3[,2]	
Sockeye_Pred=data3[,3]	
summary(Sockeye_Pred)
summary(Sockeye_Obs)

pdf(file="FIGURE 6.pdf", height=8.5, width=11) 
pg1 <- ggplot(data3, aes(Sockeye_Pred, Sockeye_Obs))+ylab("")+xlab("")+
  ggtitle(expression(paste("in-sample ", R^2, "=0.87, p<0.0001                                                      out-sample ", R^2, "=0.74, p<0.0001"))) + geom_point(shape = 1) +
  facet_grid(~Sample)+ geom_abline(slope = 1) + facet_wrap(~Sample, nrow = 1) +  theme(aspect.ratio = 1) +
  xlim(-1,200)+
  ylim(-1,200)+
  #  geom_point(          )
  geom_smooth(method = "lm", se = T, fullrange = T, colour = "blue", cex = 1)
grid.arrange(pg1,bottom=textGrob(expression(paste("Predicted")),
                                 gp = gpar(col = "black", fontsize = 20)),left=textGrob(expression(paste("Observed")),rot=90,
                                                                                        gp = gpar(col = "black", fontsize = 20)),nrow = 1, ncol=1)
dev.off()


#################
# FIGURE 7 and 8
#################

###############################################
# Generate annual indices for each covariate
###############################################
settings = make_settings( Version="VAST_v14_0_1",  
                          n_x=500,
                          Region='User', 
                          purpose="index2",
                          bias.correct=TRUE,
                          max_cells = 2000,
                          ObsModel=c(2,1), # catch, Pollock, pink
                          # ObsModel=c(2,4), # calanus
                          fine_scale=TRUE,
                          treat_nonencounter_as_zero=TRUE,
                          knot_method='grid',
                          use_anisotropy = TRUE)

settings$Options = c( settings$Options, "report_additional_variables"=TRUE )
user_region <- readRDS('user_region.rds')
fit_SST_TS = fit_model( "settings"=settings,
                 "Lat_i"=temp[,'Lat'],
                 "Lon_i"=temp[,'Lon'],
                 "observations_LL" = temp[,c('Lat','Lon')],
                 "t_i"=temp[,'Year'], 
                 "b_i"=(as_units(SST,"none")), 
                 "a_i"=temp[,'Constant'],
                 "getsd"=TRUE,
                 "test_fit"=TRUE,
                 "build_model" = TRUE,
                 input_grid=user_region ,
                 getReportCovariance=TRUE
)

x=summary(fit_SST_TS$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"SST_TS.csv", row.names = TRUE)
plot_results(fit_SST_TS, n_cells=2000  )

fit_Cal_TS = fit_model( "settings"=settings,
                      "Lat_i"=temp[,'Lat'],
                      "Lon_i"=temp[,'Lon'],
                      "observations_LL" = temp[,c('Lat','Lon')],
                      "t_i"=temp[,'Year'], 
                      "b_i"=(as_units(log(temp[,'Calanus']+1),"count")), 
                      "a_i"=temp[,'Constant'],
                      "getsd"=TRUE,
                      "test_fit"=TRUE,
                      "build_model" = TRUE,
                      input_grid=user_region ,
                      getReportCovariance=TRUE
)
x=summary(fit_Cal_TS$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"Cal_TS.csv", row.names = TRUE)
plot_results(fit_Cal_TS, n_cells=2000  )

fit_Poll_TS = fit_model( "settings"=settings,
                      "Lat_i"=temp[,'Lat'],
                      "Lon_i"=temp[,'Lon'],
                      "observations_LL" = temp[,c('Lat','Lon')],
                      "t_i"=temp[,'Year'], 
                       "b_i"=(as_units(log(temp[,'A0PollockN']+1),"count")), 
                      "a_i"=(as_units(temp[,'AreaSwept'],'km^2')), 
                      "getsd"=TRUE,
                      "test_fit"=TRUE,
                      "build_model" = TRUE,
                      input_grid=user_region ,
                      getReportCovariance=TRUE
)
x=summary(fit_Poll_TS$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"Poll_TS.csv", row.names = TRUE)
plot_results(fit_Poll_TS, n_cells=2000  )

fit_Pink_TS = fit_model( "settings"=settings,
                      "Lat_i"=temp[,'Lat'],
                      "Lon_i"=temp[,'Lon'],
                      "observations_LL" = temp[,c('Lat','Lon')],
                      "t_i"=temp[,'Year'], 
                      "b_i"=(as_units(log(temp[,'J_Pink']+1),"kg")), 
                      "a_i"=(as_units(temp[,'AreaSwept'],'km^2')), 
                      "getsd"=TRUE,
                      "test_fit"=TRUE,
                      "build_model" = TRUE,
                      input_grid=user_region ,
                      getReportCovariance=TRUE
)

x=summary(fit_Pink_TS$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"Pink_TS.csv", row.names = TRUE)
plot_results(fit_Pink_TS, n_cells=2000  )


########################
#  TIME SERIES PLOTS
########################
temp2<-read.csv("TimeseriesPlots.csv") 

Year=temp2[,1]	
J_Sockeye	=temp2[,2]
J_sockeye_SE=temp2[,3]	
Calanus	=temp2[,4]
Calanus_SE=temp2[,5]	
Pollock	=temp2[,6]	
Pollock_SE=temp2[,7]		
J_Pink	=temp2[,8]	
J_Pink_SE=temp2[,9]		
SST	=temp2[,10]	
SST_SE=temp2[,11]	
Northing=temp2[,12]		
Northing_SE	=temp2[,13]	
EAO	=temp2[,14]	
EAO_SE=temp2[,15]	
lnCalanus=temp2[,16]		
lnCalanus_SE=temp2[,17]		
lnPollock	=temp2[,18]	
lnPollock_SE=temp2[,19]		
lnPink	=temp2[,20]	
lnPink_SE=temp2[,21]
Easting=temp2[,22]
Easting_SE=temp2[,23]

require(cowplot)
require(gridExtra)
require(ggplot2)
require(reprex)
install.packages("reprex")
geom_ribbon(aes(x=Year,
                ymin=Tmin,
                ymax=Tmax),
            fill='blue',alpha=0.2)
theme_bw()

#################################################################################
#  FIGURE 7: Juvenile sockeye salmon annual indices of distribution and abundance
#################################################################################
pdf(file="FIGURE 7.pdf", height=11, width=8.5, paper='special', family ='serif')   

p1<-ggplot(temp2, aes(x=Year, y=J_Sockeye)) +
  theme_bw()+ 
  geom_line()+
  geom_pointrange(aes(ymin=J_Sockeye-J_Sockeye_SE, ymax=J_Sockeye+J_Sockeye_SE))+
  labs(title="", y="Juvenile sockeye salmon (kg)", x="", cex=3)+
  theme(text=element_text(family="sans"),axis.title = element_text(size = 12,colour = "black"),
        axis.text = element_text(size = 12,colour = "black"),
        panel.background = element_rect(colour = "black", linewidth=1.5))+
  geom_point(size=4) +
  geom_segment(aes(x=2001.5,y=0,xend=2005.5,yend=0), color="red",linewidth=2)+
  geom_segment(aes(x=2005.5,y=0,xend=2012.5,yend=0), color="blue",linewidth=2)+
  geom_segment(aes(x=2012.5,y=0,xend=2018.5,yend=0), color="red",linewidth=2)

p2<-ggplot(temp2, aes(x=Year, y=Northing)) +
  theme_bw()+ 
  geom_line()+
  geom_pointrange(aes(ymin=Northing-Northing_SE, ymax=Northing+Northing_SE))+
  labs(title="", y="Northing (km from Equator)", x="", cex=3)+
  theme(text=element_text(family="sans"),axis.title = element_text(size = 12,colour = "black"),
        axis.text = element_text(size = 12,colour = "black"),
        panel.background = element_rect(colour = "black", linewidth=1.5))+
  geom_point(size=4)+
  geom_segment(aes(x=2001.5,y=6300,xend=2005.5,yend=6300), color="red",linewidth=2)+
  geom_segment(aes(x=2005.5,y=6300,xend=2012.5,yend=6300), color="blue",linewidth=2)+
  geom_segment(aes(x=2012.5,y=6300,xend=2018.5,yend=6300), color="red",linewidth=2) 
p3<-ggplot(temp2, aes(x=Year, y=EAO)) +
  theme_bw()+ 
  geom_line()+
  geom_pointrange(aes(ymin=EAO-EAO_SE, ymax=EAO+EAO_SE))+
  labs(title="", y="", x="", cex=3)+
  ylab(bquote('Area occupied '(km^2)))+
  theme(text=element_text(family="sans"),axis.title = element_text(size = 12,colour = "black"),
        axis.text = element_text(size = 12,colour = "black"),
        panel.background = element_rect(colour = "black", linewidth=1.5))+
  geom_point(size=4) +
  geom_segment(aes(x=2001.5,y=0,xend=2005.5,yend=0), color="red",linewidth=2)+
  geom_segment(aes(x=2005.5,y=0,xend=2012.5,yend=0), color="blue",linewidth=2)+
  geom_segment(aes(x=2012.5,y=0,xend=2018.5,yend=0), color="red",linewidth=2)
p4<-ggplot(temp2, aes(x=Year, y=Easting)) +
  theme_bw()+ 
  geom_line()+
  geom_pointrange(aes(ymin=Easting-Easting_SE, ymax=Easting+Easting_SE))+
  labs(title="", y="Easting (km from 180)", x="", cex=3)+
  theme(text=element_text(family="sans"),axis.title = element_text(size = 12,colour = "black"),
        axis.text = element_text(size = 12,colour = "black"),
        panel.background = element_rect(colour = "black", linewidth=1.5))+
  geom_point(size=4)+
  geom_segment(aes(x=2001.5,y=395,xend=2005.5,yend=395), color="red",linewidth=2)+
  geom_segment(aes(x=2005.5,y=395,xend=2012.5,yend=395), color="blue",linewidth=2)+
  geom_segment(aes(x=2012.5,y=395,xend=2018.5,yend=395), color="red",linewidth=2)     

grid.arrange(p1,p2, p3,p4,bottom="Year",top="",left="",nrow = 2, ncol=2)
dev.off() 

##############
#  FIGURE 8
##############


pdf(file="FIGURE 8.pdf", height=11, width=8.5, paper='special')   

p4<-ggplot(temp2, aes(x=Year, y=SST)) +
  theme_bw()+ 
  geom_line()+
  geom_pointrange(aes(ymin=SST-SST_SE, ymax=SST+SST_SE))+
  labs(title="", y="Temperature (Celsius)", x="", cex=3)+
  theme(axis.title = element_text(size = 16,colour = "black"),
        axis.text = element_text(size = 16,colour = "black"),
        panel.background = element_rect(colour = "black", linewidth=1.5))+
  geom_point(size=4)     

p5<-ggplot(temp2, aes(x=Year, y=lnCalanus)) +
  theme_bw()+ 
  geom_line()+
  geom_pointrange(aes(ymin=lnCalanus-lnCalanus_SE, ymax=lnCalanus+lnCalanus_SE))+
  ylab(bquote('Calanus  (ln '(N ~m^-2+1)))+
  # ylab(reprex::bquote('Calanus (ln(#~m^-2+1))'))+
  theme(axis.title = element_text(size = 16,colour = "black"),
        axis.text = element_text(size = 16,colour = "black"),
        panel.background = element_rect(colour = "black", linewidth=1.5))+
  geom_point(size=4)  
p6<-ggplot(temp2, aes(x=Year, y=lnPollock)) +
  theme_bw()+ 
  geom_line()+
  geom_pointrange(aes(ymin=lnPollock-lnPollock_SE, ymax=lnPollock+lnPollock_SE))+
  labs(title="", y="Age-0 Pollock (ln(#+1))", x="", cex=3)+
  theme(axis.title = element_text(size = 16,colour = "black"),
        axis.text = element_text(size = 16,colour = "black"),
        panel.background = element_rect(colour = "black", linewidth=1.5))+
  geom_point(size=4)  
p7<-ggplot(temp2, aes(x=Year, y=lnPink)) +
  theme_bw()+ 
  geom_line()+
  geom_pointrange(aes(ymin=lnPink-lnPink_SE, ymax=lnPink+lnPink_SE))+
  labs(title="", y="Juvenile pink salmon (ln(kg+1))", x="", cex=3)+
  theme(axis.title = element_text(size = 16,colour = "black"),
        axis.text = element_text(size = 16,colour = "black"),
        panel.background = element_rect(colour = "black", linewidth=1.5))+
  geom_point(size=4)  

grid.arrange(p4,p5,p6,p7,bottom="Year",top="",left="",nrow = 2, ncol=2)
dev.off() 


##############################################
# FIGURE 9: Relationship among annual indices
##############################################

pdf(file="FIGURE 9.pdf", height=11, width=8.5, paper='special')   
p1 <- ggplot(temp2, aes(Temp_20m, J_Sockeye)) + 
  geom_point()+
  xlab(NULL)+
  ylab("J. sockeye salmon (kg)")+
  ggtitle(bquote(R^2~'=0.36, p=0.03'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)+
  theme(plot.title = element_text(size = 12, face = "bold"))
p2 <- ggplot(temp5, aes(Calanus/1000000, J_Sockeye)) +
  geom_point()+
  ylab(NULL)+  
  xlab(NULL)+
  ggtitle(bquote('NS'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)
p3 <- ggplot(temp5, aes(lnPollock, J_Sockeye)) +
  geom_point()+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle(bquote(R^2~'=0.45, p=0.01'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)
p4 <- ggplot(temp5, aes(lnPink, J_Sockeye)) +
  geom_point()+
  ylab(NULL)+
  ggtitle(bquote(R^2~'=0.37, p=0.03'))+
  xlab(NULL)+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)
p5 <- ggplot(temp5, aes(Temp_20m, Northing)) + 
  geom_point()+
  xlab(NULL)+ 
  ylab("Northing (km from Equator)")+
  ggtitle(bquote(R^2~'=0.43, p=0.02'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)
p6 <- ggplot(temp5, aes(Calanus/1000000, Northing)) +
  geom_point()+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle(bquote('NS'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)
p7 <- ggplot(temp5, aes(lnPollock, Northing)) +
  geom_point()+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle(bquote(R^2~'=0.41, p=0.02'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)
p8 <- ggplot(temp5, aes(lnPink, Northing)) +
  geom_point()+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle(bquote('NS'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12) 
p9 <- ggplot(temp5, aes(Temp_20m, Easting)) + 
  geom_point()+
  xlab(NULL)+
  ylab("Easting (km from 180)")+
  ggtitle(bquote(R^2~'=0.67, p<0.0001'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)+
  theme(plot.title = element_text(size = 40, face = "bold"))
p10 <- ggplot(temp5, aes(Calanus/1000000, Easting)) +
  geom_point()+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle(bquote('NS'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)
p11 <- ggplot(temp5, aes(lnPollock, Easting)) +
  geom_point()+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle(bquote(R^2~'=0.47, p=0.01'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)
p12 <- ggplot(temp5, aes(lnPink, Easting)) +
  geom_point()+
  xlab(NULL)+
  ylab(NULL)+
  ggtitle(bquote('NS'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)  
p13 <- ggplot(temp5, aes(Temp_20m, EAO/1000)) + 
  geom_point()+
  xlab('Temperature 20m (Celsius)')+
  ylab(expression('Area occupied ' (km^{2}*' x1000'))) +   
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  ggtitle(bquote(R^2~'=0.51, p=0.007'))+
  theme_cowplot(12)
p14 <- ggplot(temp5, aes(Calanus/1000000, EAO/1000)) +
  geom_point()+
  ylab(NULL)+
  xlab("Calanus (# x1000000)")+
  ggtitle(bquote('NS'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)
p15 <- ggplot(temp5, aes(lnPollock, EAO/1000)) +
  geom_point()+
  ylab(NULL)+
  xlab("Age-0 pollock (ln (#+1))")+
  ggtitle(bquote(R^2~'=0.47, p=0.01'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)
p16 <- ggplot(temp5, aes(lnPink, EAO/1000)) +
  geom_point()+
  ylab(NULL)+
  xlab("J. pink salmon (ln(kg+1))")+
  ggtitle(bquote(R^2~'=0.42, p=0.02'))+
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
  theme_cowplot(12)  
plot_grid(p1, p2, p3,p4,p5,p6,p7,p8 ,p9,p10,p11,p12,p13,p14,p15,p16,label_size = 10, cols=4)
dev.off()

#Polynomial models
new<-na.omit(temp5)
summary(lm(new$J_Sockeye~poly(new$Temp_20m, 2))) #R2=0.36, p=0.03
summary(lm(new$J_Sockeye~poly(new$Calanus, 2))) #NS
summary(lm(new$J_Sockeye~poly(new$lnPollock, 2))) #R2=0.45, p=0.01
summary(lm(new$J_Sockeye~poly(new$lnPink, 2))) #R2=0.37, p=0.03
summary(lm(new$Northing~poly(new$Temp_20m, 2))) #R2=0.43, p=0.02
summary(lm(new$Northing~poly(new$Calanus, 2))) #NS
summary(lm(new$Northing~poly(new$lnPollock, 2))) #R20.41, p=0.02
summary(lm(new$Northing~poly(new$lnPink, 2))) #NS
summary(lm(new$Easting~poly(new$Temp_20m, 2))) #R2=0.67, p<0.0001
summary(lm(new$Easting~poly(new$Calanus, 2))) #NS
summary(lm(new$Easting~poly(new$lnPollock, 2))) #R2=0.47, p=0.01
summary(lm(new$Easting~poly(new$lnPink, 2))) #NS
summary(lm(new$EAO~poly(new$Temp_20m, 2))) #R2=0.51, p=0.007
summary(lm(new$EAO~poly(new$Calanus, 2)))
summary(lm(new$EAO~poly(new$lnPollock, 2))) #R2=0.47, p=0.01
summary(lm(new$EAO~poly(new$lnPink, 2))) #R2=0.42, p=0.02

################################################################
#       Covariate effects in spatio-temporal model             #
################################################################
lnSST<-scale(SST)
lnSST<-(SST)

X1<-data.frame(Year, Lat, Lon,lnSST) 
X1_formula = ~ bs( lnSST, degree=2, intercept=FALSE)
X2_formula = ~ bs( lnSST, degree=2, intercept=FALSE)
X1_formula=~lnSST
X2_formula=~lnSST

Calanusnls<-scale(log(Calanus+1))
X1<-data.frame(Lat,Lon, Year, Calanusnls) 
X1_formula=~Calanusnls
X2_formula=~Calanusnls
X1_formula = ~ bs( Calanusnls, degree=2, intercept=FALSE)
X2_formula = ~ bs( Calanusnls, degree=2, intercept=FALSE)

Pinknls<-scale(log(J_Pink+1))
X1<-data.frame(Lat,Lon, Year, Pinknls) 
X1_formula=~Pinknls
X2_formula=~Pinknls
X1_formula = ~ bs( Pinknls, degree=2, intercept=FALSE)
X2_formula = ~ bs( Pinknls, degree=2, intercept=FALSE)

Pollocknls<-scale(log(A0PollockN+1))
X1<-data.frame(Lat,Lon, Year, Pollocknls) 
X1_formula = ~ Pollocknls
X2_formula = ~ Pollocknls
X1_formula = ~ bs( Pollocknls, degree=2, intercept=FALSE)
X2_formula = ~ bs( Pollocknls, degree=2, intercept=FALSE)

X1<-data.frame(Lat,Lon, Year, lnSST,Calanusnls, Pinknls,Pollocknls) 
X1_formula = ~ lnSST+Calanusnls+ Pinknls+Pollocknls
X2_formula = ~ lnSST+Calanusnls+Pinknls+Pollocknls
dat<-data.frame(lnSST,Calanusnls,Pinknls,Pollocknls)


##############################
# SOCKEYE AND COVARIATE MODELS
##############################
settings = make_settings( Version="VAST_v14_0_1",
                           n_x=500,
                           Region='User', 
                           purpose="index2",
                           bias.correct=TRUE,
                           max_cells = 2000,
                           ObsModel=c(2,1), 
                           fine_scale=TRUE,
                           treat_nonencounter_as_zero=TRUE,
                           knot_method='grid',
                           use_anisotropy = TRUE)
user_region <- readRDS('user_region.rds')

fit_SST = fit_model(
  "settings"=settings,
  "Lat_i"=temp[,'Lat'],
  "Lon_i"=temp[,'Lon'],
  "observations_LL" = temp[,c('Lat','Lon')],
  "t_i"=temp[,'Year'],
  "c_iz"=rep(0,nrow(temp)), # for single species 
  "b_i"=(as_units(temp[,'kg'],"kg")),
  "a_i"=(as_units(temp[,'AreaSwept'], 'km^2')),
  "getsd"=TRUE,
  "test_fit"=TRUE,
  "build_model" = TRUE,
  "X1_formula"=X1_formula, 
  "X2_formula"=X2_formula, 
  "vars_to_correct"="Index_cyl",
  "covariate_data"=X1, 
  input_grid=user_region,
  getReportCovariance=TRUE)

x=summary(fit_SST$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"Sock_SST.csv", row.names = TRUE)

# PLOT RESULTs =============================================================================
plot_results( fit_SST, settings=settings, plot_set=c(3,11,12,14,15), n_cells=2000)


###########################################
# Covariate effects on 1st linear predictor
###########################################

#1st Linear predictor for covariates
str(fit$Report) #eta1_gct
fit$Report$eta1_gct

#Covariate effects at each extrapolation grid location
names(fit$Report)[grepl('_gc|_gct', x=names(fit$Report))]
D_gt <- fit$Report$eta1_gct[,1,] # drop the category


####################################
# Density covariate effect eta1_gct
####################################

str(fit$Report) #eta1_gct
fit$Report$eta1_gct
#Densities at each extrapolation grid location
names(fit$Report)[grepl('_gc|_gct', x=names(fit$Report))]
D_gt <- fit$Report$eta1_gct[,1,] # drop the category
D_gt <- fit$Report$eta2_gct[,1,] # drop the category


#dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=c(2003:2019))
dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=years)
## tidy way of doing this, reshape2::melt() does
## it cleanly but is deprecated
#?pivot_longer.sf
#?rownames_to_column
#require(tibble)
library(tidyverse)
#?pivot_longer
D_gt <- D_gt %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  #tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')
D <- merge(D_gt, mdl$PlotDF, by.x='cell', by.y='x2i')
D3<-lapply(D, as.numeric)
D4<-data.frame(D3)
head(D4)
summary(log(D4$D))

#Subet D values greater than .1% of D to plot

D02 <- D4[which(D$Year == "2002"),names(D4) %in% c("Year","D","Lat","Lon")]
D03 <- D4[which(D$Year == "2003"),names(D4) %in% c("Year","D","Lat","Lon")]
D04 <- D4[which(D$Year == "2004"),names(D4) %in% c("Year","D","Lat","Lon")]
D05 <- D4[which(D$Year == "2005"),names(D4) %in% c("Year","D","Lat","Lon")]
D06 <- D4[which(D$Year == "2006"),names(D4) %in% c("Year","D","Lat","Lon")]
D07 <- D4[which(D$Year == "2007"),names(D4) %in% c("Year","D","Lat","Lon")]
D08 <- D4[which(D$Year == "2008"),names(D4) %in% c("Year","D","Lat","Lon")]
D09 <- D4[which(D$Year == "2009"),names(D4) %in% c("Year","D","Lat","Lon")]
D10 <- D4[which(D$Year == "2010"),names(D4) %in% c("Year","D","Lat","Lon")]
D11 <- D4[which(D$Year == "2011"),names(D4) %in% c("Year","D","Lat","Lon")]
D12 <- D4[which(D$Year == "2012"),names(D4) %in% c("Year","D","Lat","Lon")]
D14 <- D4[which(D$Year == "2014"),names(D4) %in% c("Year","D","Lat","Lon")]
D16 <- D4[which(D$Year == "2016"),names(D4) %in% c("Year","D","Lat","Lon")]
D18 <- D4[which(D$Year == "2018"),names(D4) %in% c("Year","D","Lat","Lon")]


D5<-rbind(D02,D03,D04,D05,D06,D07,D08, D09,D10,D11,D12,D14,D16,D18) 
D6<-lapply(D5, as.numeric)
D7<-data.frame(D6)
summary(D7)
write.csv(D7, "Cov_D_OUT.csv")
theme_set(theme_bw())

D9<-D7[D7$D != "NA", ]         # Multiple conditions
#str(D9)
D9<-D7

require("rgdal") 
require("maptools")
require("ggplot2")
require("plyr")

  library(maptools)
  library(rgdal)
  library(sp)
  library(ggplot2)
  library(tidyverse)
  #install.packages("ggtext")
  library(ggtext)

  WGScoor<- readRDS('user_region.rds')
  str(WGScoor)
  summary(WGScoor)
  coordinates(WGScoor)=~Lon+Lat
  proj4string(WGScoor)<- CRS("+proj=longlat +datum=WGS84")
  LLcoor<-spTransform(WGScoor,CRS("+proj=longlat +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
  # raster::shapefile(LLcoor, "MyShapefile.shp")
  #shp <- readOGR(dsn=file.path("MyShapefile.shp"))
  df.WGScoor<-as.data.frame(WGScoor)
  df2.WGScoor<-fortify(df.WGScoor)
  str(df2.WGScoor)

######################################
# APPENDIX COVARIATE EFFECTS FIGURES 1A-7A
######################################

Cal<-log(Calanus+1)
Poll<-log(A0PollockN+1)
Pink<-log(J_Pink+1)
summary(Cal)
summary(Calanus)

pdf(file="FIGURE 7A.pdf", height=8.5, width=11) #This generates the figure as a hi-res (pub quality) tiff file. You can change the dimensions of the figure by changing the height and width arguments. You can change the resolution with the 'res' argument.
#Data are in 
g2 <- gmap + 
  theme_bw()+
  #scale_color_gradientn(colours = colorspace::divergingx_hcl(palette="RdGy", n=7, rev=TRUE),oob=scales::squish,limits=c(0,4))+
  
  #geom_polygon(data = temp, aes(x = Lon, y = Lat, group = NULL), fill=NA,color="lightgray" ,size = 0.005) +
  geom_point(data=D7 , aes(Lon, Lat, color=((D)), group=NULL),  #Use log(D) for GRP to covert back from exp(GRP)
             #geom_point(data=temp2 , aes(Lon, Lat, color=(D), group=NULL),  #Use log(D) for GRP to covert back from exp(GRP)
             
             na.rm=TRUE, size=2, stroke=0, shape=16) +facet_wrap('Year')+ 
  # labs(color = title_exp, shape = title_exp)+
  #labs(color=paste('log((#+1)·km)'))+
  #labs(color=paste('log(kg\U00B7km\u00b2+1)'))+
  # labs(color=paste('log(#+1)\U00B7km\u00b2'))+
  #labs(color=expression(paste("ln(kg/",km^2,")")))+
  #labs(color=paste('ln(kg\U00B7km\u00b2)'))+
  #labs(color=paste('log(#+1)/km\u00b2'))+ #Pink
  # labs(color=paste('log(#+1)/km\u00b2'))+ #Pollock
  # labs(color=paste('ln(#\U00B7m\u00b2+1)'))+
  # labs(color=paste('log(#+1)'))+
  # labs(color=paste('log(kg+1)'))+
labs(color=paste('Covariate\n effect'))+
  ##⋅m-2
  #labs(color=paste('Celsius'))+
  #labs(color=paste('Calanus'))+
  
  #labs(color=paste('log(kg\U00B7km\u00b2+1)'))+
  xlab(expression(paste(Longitude^o,~'W')))+
  ylab(expression(paste(Latitude^o,~'N')))+
  #coord_map(xlim = c(-173, -155),ylim = c(50, 66.5))+
  #labs(y= "Latitude", x = "Longitude", cex=3)+
  #labs(title="Age-0 Pollock")+
  labs(title="Juvenile pink salmon")+
  #labs(title="Calanus")+
  #labs(title="Juvenile sockeye salmon")+
  #labs(title="Temperature 20 m effect")+
  #theme(title = element_text(face = "italic"))+
  theme(strip.text.x = element_text(size = 12,face="bold"))  +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        panel.grid.major=element_line(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = "grey50"))
g2
dev.off()



########################################
# FIGURE 10 LINEAR COVARIATE EFFECT PLOT
########################################
library(mgcv)
library(effects)
library(gridExtra)
#Rerun each respective VAST covariate effect model

pdf(file="FIGURE 10.pdf", height=8.5, width=11)
covariate_data_full = fit_SST$effects$covariate_data_full
catchability_data_full = fit_SST$effects$catchability_data_full
pred3 = Effect.fit_model( fit_SST,
                          focal.predictors = c("SST"),
                          which_formula = "X1", 
                          xlevels = 100)
plot3<-plot(pred3,geom='smooth', xlab="Sea surface temperature (Celsius)", ylab=NULL, main="")

covariate_data_full = fit_Poll$effects$covariate_data_full
catchability_data_full = fit_Poll$effects$catchability_data_full
pred4 = Effect.fit_model( fit_Poll,
                          focal.predictors = c("Pollock"),
                          which_formula = "X1", 
                          xlevels = 100)
plot4<-plot(pred4, xlab="Age-0 Pollock (ln(#+1))", ylab=NULL, main="")
covariate_data_full = fit_Pink$effects$covariate_data_full
catchability_data_full = fit_Pink$effects$catchability_data_full
pred5 = Effect.fit_model( fit_Pink,
                          focal.predictors = c("Pink_salmon"),
                          which_formula = "X1", 
                          xlevels = 100)
plot5<-plot(pred5, xlab="Juvenile pink salmon (ln(kg+1))", ylab=NULL, main="")
grid.arrange(plot3, plot4,plot5,bottom=textGrob(expression(paste("Covariates"))),left=textGrob(expression(paste("1st Linear predictor")),rot=90),nrow = 1, ncol=3)
dev.off()

pdf(file="FIGURE 10.pdf", height=8.5, width=11) 
#plot1<-plot(pred1, ylab="", xlab="Sea temperature (Celsius)", main="Encounter probability")
#plot2<-plot(pred2, ylab="", xlab="Sea temperature", main="Positive catch rate")
plot3<-plot(pred3, ylab="", xlab="Sea temperature",main="")
plot4<-plot(pred4, ylab="", xlab="log(Pink salmon(kg)+1)",main="")
#plot6<-plot(pred6, ylab="", xlab="Pink salmon", main="Positive catch rate")
plot5<-plot(pred5, ylab="", xlab="log(Age-0 pollock(#)+1)",main="")
#plot8<-plot(pred8, ylab="", xlab="log(Age-0 pollock(#)+1)", main="Positive catch rate")
#grid.arrange(plot1, plot5,plot7,plot8,bottom=textGrob(expression(paste("Covariate"), cex=3),gp = gpar(col = "black", fontsize = 16)),left=textGrob(expression(paste("Linear predictor"), cex=3),gp = gpar(col = "black", fontsize = 16),rot=90),nrow = 2, ncol=2)
grid.arrange(plot3,plot4,plot5,bottom=textGrob(expression(paste("Covariate"), cex=3),gp = gpar(col = "black", fontsize = 16)),left=textGrob(expression(paste("1st linear predictor"), cex=3),gp = gpar(col = "black", fontsize = 16),rot=90),nrow = 1, ncol=3)

dev.off()




#QUESTION: How do I only plot certain years/ exclude years modeled with fake data.
# Fix issue in auto-generated plot labels
#https://github.com/James-Thorson-NOAA/VAST/issues/274
fit$years_to_plot = c( 2002:2012,2014,2016,2018)



#########################
# REPLOT covariate effects#
#########################
  library(viridis)
## Below shows to you get the model estimate of density, D_gct,
## for each grid (g), category (c; not used here single
## univariate); and year (t); and link it spatially to a lat/lon
## extrapolation point.  You can do this for any _gct or _gc
## variable in the Report.

years <- unique(temp$Year)
nyrs <- length(years)
#Year=rep(years, each=nrow(ak_map))
## Remake map list locally for recreating plots
mdl <- make_map_info(Region = settings2$Region,
                     spatial_list = fit2$spatial_list,
                     Extrapolation_List = fit2$extrapolation_list)
## quick dirty AK map
ak_map <- subset(map_data("world"), region=='USA' & subregion=='Alaska')
## Have to duplicate it for each year so can facet below
ak_map <- cbind(ak_map[rep(1:nrow(ak_map), times=nyrs),],
                Year=rep(years, each=nrow(ak_map)))
gmap <- ggplot(ak_map, aes(x = long, y = lat, group = group)) +
  # geom_raster(aes(x=long, y=lat, fill=D))+
  # geom_raster(aes(fill = D), interpolate = TRUE)+
  geom_polygon(fill="black", colour = "white") +
  # scale_colour_gradientn(colors=c('red','orange','yellow','green', 'blue','purple')) +
  #  scale_colour_gradientn(colors=c('darkblue','blue','lightblue','lightgreen','yellow', 'orange','red')) +
  #scale_color_gradientn(colours = rainbow(5))+

  # scale_fill_viridis() + 
  scale_color_viridis_c(option = "viridis") +  #changed form magma to plasma or virdis

   theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  #  coord_cartesian(xlim=mdl$Xlim, ylim=mdl$Ylim)
  coord_cartesian(xlim=mdl$Xlim, ylim=mdl$Ylim)+ 
  guides(fill=guide_legend(title='Energy Density'))
#Densities at each extrapolation grid location
names(fit2$Report)[grepl('_gc|_gct', x=names(fit2$Report))]
D_gt <- fit2$Report$eta1_gct[,1,] # Encounter Probability
#D_gt <- fit2$Report$eta2_gct[,1,] # Postivie catch rates

summary(fit2$Report)
print(fit2$Report)
str(fit2$Report)
#dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=c(2003:2019))
dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=years)
## tidy way of doing this, reshape2::melt() does
## it cleanly but is deprecated
library(tidyr)
D_gt <- D_gt %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_longer(-cell, names_to = "Year", values_to='D')
D <- merge(D_gt, mdl$PlotDF, by.x='cell', by.y='x2i')

D3<-lapply(D, as.numeric)
D4<-data.frame(D3)
D02 <- D4[which(D$Year == "2002"),names(D4) %in% c("Year","D","Lat","Lon")]
D03 <- D4[which(D$Year == "2003"),names(D4) %in% c("Year","D","Lat","Lon")]
D04 <- D4[which(D$Year == "2004"),names(D4) %in% c("Year","D","Lat","Lon")]
D05 <- D4[which(D$Year == "2005"),names(D4) %in% c("Year","D","Lat","Lon")]
D06 <- D4[which(D$Year == "2006"),names(D4) %in% c("Year","D","Lat","Lon")]
D07 <- D4[which(D$Year == "2007"),names(D4) %in% c("Year","D","Lat","Lon")]
D08 <- D4[which(D$Year == "2008"),names(D) %in% c("Year","D","Lat","Lon")]
D09 <- D4[which(D$Year == "2009"),names(D4) %in% c("Year","D","Lat","Lon")]
D10 <- D4[which(D$Year == "2010"),names(D4) %in% c("Year","D","Lat","Lon")]
D11 <- D4[which(D$Year == "2011"),names(D4) %in% c("Year","D","Lat","Lon")]
D12 <- D4[which(D$Year == "2012"),names(D4) %in% c("Year","D","Lat","Lon")]
D14 <- D4[which(D$Year == "2014"),names(D4) %in% c("Year","D","Lat","Lon")]
D16 <- D4[which(D$Year == "2016"),names(D4) %in% c("Year","D","Lat","Lon")]
D18 <- D4[which(D$Year == "2018"),names(D4) %in% c("Year","D","Lat","Lon")]

D5<-rbind(D02,D03,D04,D05,D06,D07,D08, D09,D10,D11,D12,D14,D16,D18) 
D6<-lapply(D5, as.numeric)
D7<-data.frame(D6)
summary(D7)
#write.csv(D7, "CATCH_D_OUT.csv")
theme_set(theme_bw())
library(spatstat)
library(maps)
library(maptools)


######################################
# APPRENDIX COVARIATE EFFECTS FIGURES
######################################
g <- gmap +
  geom_point(data=D7, aes(Lon, Lat, color=(as.numeric(exp(D))), group=NULL),
             size=2, stroke=0,shape=16) + facet_wrap('Year')+
             labs(color=paste('Covariate \n effect'))+
  theme(strip.text.x = element_text(size = 6))+
  geom_point(na.rm=TRUE)
  
g
?geom_point




  
  ##################################################### 
  # CREATE USER REGION FOR NEW SURVEYS: DONE for BASIS=================================
  ######################################################
  #DONE
  library(sp) #1.4.5
  library(sf) #0.9.6
  library(rgdal)
  library(here)
  library(splines)
  library(sp)  # vector data
  library(raster)  # raster data
  library(rgdal)  # input/output, projections
  library(rgeos)  # geometry ops
  library(spdep)  # spatial dependence
  ### An example of how to create user-defined extrapolation
  ### regions (extents) for VAST.
  ## Read in the data
  dat <- read.csv("IYSsockeye18.csv")
  
  ### Cecilia O'Leary and Cole Monnahan | December 2020
  
  ### The extrapolation region defines the extent over which the
  ### model predictions are integrated. Any density outside of this
  ### region will not be included in estimates of the index. It is
  ### not used in model fitting. It comes with many built-in
  ### regions but often a user needs to define their own. Here, we
  ### demonstrate two ways of doing this: (1) From a set of points
  ### representing the outer extent of the region; (2) from an
  ### existing shape file.
  
  library(sp) # 1.4.4
  library(sf) # 0.9.6
  
  ### Method 1: use a set of lat/lon coordinates which define the
  ### outer edge of the region. For instance you might want to plot
  ### your data and simply create a region that captures it. The
  ### locator() function can be useful for this as shown
  
  ### Use this to draw points around your data
  plot(dat$Lon, dat$Lat)
  LL <- locator()
  saveRDS(LL, 'extent_LL.rds')
  
  ## Take a data.frame of coordinates in longitude/latitude that
  ## define the outer limits of the region (the extent).
  LL <- readRDS('extent_LL.rds')
  region_extent <- data.frame(long=LL$x, lat=LL$y)
  str(region_extent)
  ## > 'data.frame':	42 obs. of  2 variables:
  ## $ long: num  -166 -166 -165 -165 -164 ...
  ## $ lat : num  53.9 54.1 54.2 54.6 55 ...
  
  #### Turn it into a spatial polygon object
  ## Need to duplicate a point so that it is connected
  region_extent <- rbind(region_extent, region_extent[1,])
  ## https://www.maths.lancs.ac.uk/~rowlings/Teaching/Sheffield2013/cheatsheet.html
  poly <- Polygon(region_extent)
  polys <- Polygons(list(poly), ID='all')
  sps <- SpatialPolygons(list(polys))
  ## I think the F_AREA could be dropped here
  sps <- SpatialPolygonsDataFrame(sps, data.frame(Id=factor('all'), F_AREA=1, row.names='all'))
  proj4string(sps)<- CRS("+proj=longlat +datum=WGS84")
  sps <- spTransform(sps, CRS("+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
  ### Get UTM zone for conversion to UTM projection
  ## retrieves spatial bounding box from spatial data [,1] is
  ## longitude
  lon <- sum(bbox(sps)[1,])/2
  ## convert decimal degrees to utm zone for average longitude, use
  ## for new CRS
  utmzone <- floor((lon + 180)/6)+1
  crs_LL <- CRS('+proj=longlat +ellps=WGS84 +no_defs')
  sps@proj4string <- crs_LL
  
  ### --------------------------------------------------
  ### Create the VAST extroplation grid for method 1 and 2
  ## Convert the final in polygon to UTM
  crs_UTM <- CRS(paste0("+proj=utm +zone=",utmzone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
  region_polygon <- spTransform(sps, crs_UTM)
  
  ### Construct the extroplation grid for VAST using sf package
  ## Size of grid **in meters** (since working in UTM). Controls
  ## the resolution of the grid.
  cell_size <- 2000
  ## This step is slow at high resolutions
  region_grid <- st_make_grid(region_polygon, cellsize = cell_size, what = "centers")
  ## Convert region_grid to Spatial Points to SpatialPointsDataFrame
  region_grid <- as(region_grid, "Spatial")
  region_grid_sp <- as(region_grid, "SpatialPointsDataFrame")
  ## combine shapefile data (region_polygon) with Spatial Points
  ## (region_grid_spatial) & place in SpatialPointsDataFrame data
  ## (this provides you with your strata identifier (here called
  ## Id) in your data frame))
  region_grid_sp@data <- over(region_grid, region_polygon)
  
  ## Convert back to lon/lat coordinates as that is what VAST uses
  region_grid_LL <- as.data.frame(spTransform(region_grid_sp, crs_LL))
  region_df <- with(region_grid_LL,
                    data.frame(Lon=coords.x1,
                               Lat=coords.x2, Id,
                               Area_km2=( (cell_size/1000)^2),
                               row=1:nrow(region_grid_LL)))
  ## Filter out the grid that does not overlap (outside extent)
  region <- subset(region_df, !is.na(Id))
  ## This is the final file needed.
  str(region)
  ## > 'data.frame':	106654 obs. of  5 variables:
  ##  $ Lon     : num  -166 -166 -166 -166 -166 ...
  ##  $ Lat     : num  53.9 53.9 54 53.9 53.9 ...
  ##  $ Id      : Factor w/ 1 level "all": 1 1 1 1 1 1 1 1 1 1 ...
  ##  $ Area_km2: num  4 4 4 4 4 4 4 4 4 4 ...
  ##  $ row     : int  401 402 975 976 977 978 1549 1550 1551 1552 ...
  
  ### Save it to be read in and passed to VAST later.
  saveRDS(region, file = "user_region.rds")
  ### End of creating user extrapolation region object
  ### --------------------------------------------------
  
  ### Quick plots of the process for method 1
  png('user_region.png', width=7, height=7, units='in', res=200)
  par(mfrow=c(2,2))
  with(region_extent, plot(long, lat, main='Extent in points in LL'))
  plot(region_polygon, main='Polygon in UTM', axes=TRUE)
  plot(region_grid, col=ifelse(is.na(region_df$Id), 'red', 'black'),
       axes=TRUE, main='Extrapolation area UTM')
  with(region, plot(Lon, Lat, main='Extrapolation region in LL', pch='.'))
  dev.off()
  
  