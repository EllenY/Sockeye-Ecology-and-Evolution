###################################################
#   Juvenile sockeye salmon project Dec 2020
###################################################
#Ellen updated with fit_model function Dec 30, 2020
#Works DEc 2020
#ISSUE: recreating SST, use an older version with estimates for 2013, 15, 17
.libPaths()
.libPaths()
myPaths <- .libPaths() 
installed.packages(lib.loc = NULL, priority = NULL,
                   noCache = FALSE, fields = NULL,
                   subarch = .Platform$r_arch)

#########################
# Add san serif fonts for PDF files
#########################
install.packages("extrafont")
library(extrafont)
# Find and save information about fonts installed on your system
font_import()
# List the fonts
fonts()
library(extrafont)
# Register the fonts with R
loadfonts()


warnings()
require(devtools)
require(utils)
require(remotes)
library(devtools)
library(utils)
library(remotes)
library("devtools", lib.loc="C:/Program Files/Microsoft/R Open/R-4.0.2/library")
library("utils", lib.loc="C:/Program Files/R/R-4.3.1/library")
library("remotes", lib.loc="C:/Program Files/R/R-4.3.1/library")

find.package("Matrix")
find.package("INLA")
find.package("TMB")
find.package("FishStatsUtils")
find.package("VAST")

#Remove Matrix 1.2.18, install 1.2.17
remove.packages("Matrix",lib="C:/Program Files/R/R-4.3.1/library")
remove.packages("INLA",lib="C:/Program Files/R/R-4.3.1/library")
remove.packages("TMB",lib="C:/Program Files/R/R-4.3.1/library")
remove.packages("FishStatsUtils",lib="C:/Program Files/R/R-4.3.1/library")
remove.packages("VAST",lib="C:/Program Files/R/R-4.3.1/library")
library(remotes)
library(utils)
library(devtools)

#Installing VAST: 
#1. Download TMB, FishStatsUtils-development, and VAST-development.

#Remove Matrix 1.2.18, install 1.2.17
remove.packages("Matrix")
remove.packages("INLA")
remove.packages("TMB")
remove.packages("TMBhelper")

remove.packages("FishStatsUtils")
remove.packages("VAST")
remove.packages("wininet")

#Install Matrix 1.2-17 downloaded from https://cran.r-project.org/src/contrib/Archive/Matrix/
#remotes::install_local("C:\\Users\\ellen.yasumiishi\\Desktop\\Matrix_1.6-1.zip",dep=FALSE, lib="C:/Program Files/R/R-4.3.1/library")

#install TMB Helper & TMB
remotes::install_local( "C:\\Users\\ellen.yasumiishi\\Desktop\\TMB_contrib_R-master\\TMBdebug", lib="C:/Program Files/R/R-4.3.1/library")
remotes::install_local( "C:\\Users\\ellen.yasumiishi\\Desktop\\TMB_contrib_R-master\\TMBhelper", lib="C:/Program Files/R/R-4.3.1/library")
remotes::install_local( "C:\\Users\\ellen.yasumiishi\\Desktop\\TMB_contrib_R-master\\TMBphase", lib="C:/Program Files/R/R-4.3.1/library")
#Failed
# Preinstall INLA before VAST
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("graph")
BiocManager::install("Rgraphviz")
library(graph)
library(Rgraphviz)
tools::write_PACKAGES(type = "win.binary", latestOnly = FALSE)

options(repos=c(
inlabruorg = "https://inlabru-org.r-universe.dev",
INLA = "https://inla.r-inla-download.org/R/testing",
CRAN = "https://cran.rstudio.com"))
install.packages("fmesher")
library(fmesher)

library("devtools")
options(download.file.method = "wininet")

#devtools::install_github(repo = "https://github.com/hrue/r-inla", force=TRUE,ref = "stable", subdir = "rinla", build = FALSE)
remotes::install_local( "C:\\Users\\ellen.yasumiishi\\Desktop\\r-inla-devel\\rinla", dep=FALSE,lib="C:/Program Files/R/R-4.3.1/library")

#install FishStatsUtils: WORKS, download development versions and install https://github.com/James-Thorson-NOAA/FishStatsUtils
remotes::install_local( "C:\\Users\\ellen.yasumiishi\\Desktop\\FishStatsUtils-main.zip", dep=FALSE,lib="C:/Program Files/Microsoft/R Open/R-4.0.2/library")
#install VAST
remotes::install_local("C:\\Users\\ellen.yasumiishi\\Desktop\\VAST-dev.zip",dep=FALSE, lib="C:/Program Files/R/R-4.3.1/library")

# Load package
install.packages('viridis')
install.packages(c("car","gridBase","cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom","akima","splines2","splines","fields", "sf", "rnaturalearth", "rnaturalearthdata"))
library(cowplot)
library(googleway)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(sf)
library(sp)
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

packageVersion("VAST")
?VAST::make_data
?FishStatsUtils::make_settings
?FishStatsUtils::fit_model
?FishStatsUtils::plot_results

example = load_example( data_set="EBS_pollock" )

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 100, 
                          Region = example$Region, 
                          purpose = "index2", 
                          bias.correct = FALSE )

# Run model
fit = fit_model( settings = settings, 
                 Lat_i = example$sampling_data[,'Lat'], 
                 Lon_i = example$sampling_data[,'Lon'], 
                 t_i = example$sampling_data[,'Year'], 
                 b_i = example$sampling_data[,'Catch_KG'], 
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 getReportCovariance=TRUE)

# Plot results
plot( fit )
gc() 
gc(reset=T)
rm(list = ls())

Version = get_latest_version( package="VAST" )
Version
wd<-setwd("C:/Users/ellen.yasumiishi/Work/2023/Sockeye")

#Don't use 2002 kg catch data for pollock. cov29 I use numbers of pollock
#26 file has error with giant outlier in index in 2007: 

#<-read.csv("IYSsockeyeCovariates29SST.csv") #SST estimate and replot 

temp<-read.csv("IYSsockeyeCovariates29.csv") 
Year=temp[,1]	
Lat	=temp[,2]
Lon	=temp[,3]
SST	=temp[,4]
A0PollockN =temp[,5]	
J_Pink =temp[,6]	
Calanus	=temp[,7] # numbers/m2
AirT	=temp[,8]
HaulDate	=temp[,9]
StationID	=temp[,10]
AreaSwept	=temp[,11]
kg	=temp[,12]
Sci	=temp[,13]
vessel	=temp[,14]
Julian	=temp[,15]
Constant=temp[,16]
summary(Year)


#Removed missing years
temp<-read.csv("IYSsockeyeCovariates31.csv") #N=809
#temp<-read.csv("IYSsockeyeCovariates30.csv") 

Year=temp[,1]	
Lat	=temp[,2]
Lon	=temp[,3]
SST	=temp[,4]
A0PollockN =temp[,5]	
J_Pink =temp[,6]	
Calanus	=temp[,7]
AirT	=temp[,8]
HaulDate	=temp[,9]
StationID	=temp[,10]
AreaSwept	=temp[,11]
kg	=temp[,12]
Sci	=temp[,13]
vessel	=temp[,14]
Julian	=temp[,15]
Constant=temp[,16]
Lat_rnd=temp[,17]
Lon_rnd=temp[,18]


#temp<-read.csv("IYSsockeyeCovariates26.csv") 
temp<-read.csv("IYSsockeyeCovariates26base.csv") #no 2013, 2015, 2017 

Year=temp[,1]	
Lat	=temp[,2]
Lon	=temp[,3]
SST	=temp[,4]
Age0_Pollock=temp[,5]	
J_Pink=temp[,6]	
Calanus=temp[,7]	
AirT	=temp[,8]
HaulDate=temp[,9]	
StationID=temp[,10]	
AreaSwept	=temp[,11]
kg	=temp[,12]
Sci	=temp[,13]
vessel=temp[,14]
Date=temp[,15]	
Julian=temp[,16]

# NEW DATASET 2003- don't use 
temp<-read.csv("IYSsockeyeCovariates27.csv") #Has 2013, 15, 17

Year=temp[,1]	
Lat=temp[,2]
Lon	=temp[,3]
SST=temp[,4]
Age0_Pollock=temp[,5]
J_Pink=temp[,6]
Calanus=temp[,7]
AirT	=temp[,8]
HaulDate=temp[,9]
StationID=temp[,10]
AreaSwept=temp[,11]	
kg=temp[,12]	
Sci	=temp[,13]
vessel=temp[,14]
A0PollockN=temp[,15]	
JPinkN	=temp[,16]
Num	=temp[,17]
Leap=temp[,18]
LAT=temp[,19]	
LON=temp[,20]
Date=temp[,21]
Julian=temp[,22]
Constant=temp[,23]

# Exploratory Plots ========================================================================
do.explore.plot <- TRUE
if(do.explore.plot==TRUE) {
  
  # Catch Distribution by year
  
  g <- ggplot(temp, aes(kg)) +
    theme_linedraw() +
    geom_histogram(fill='blue', color='black') +
    # facet_wrap(~factor(StationID))
    facet_wrap(~factor(Year), scales='fixed') +
    xlab("Weight of Sockeye catch")
  g
  ggsave(file.path(wd, "Wt Distribution.png"), plot=g, height=8, width=9, units='in')
  ?facet_wrap
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
  # g.eff
  ggsave(file.path(wd, "Effort Distribution.png"), plot=g.eff, height=8, width=9, units='in')

    # Figure 1APlot Map of Catch Rates
  CPU<-log((kg/AreaSwept)+1)
  
  CPU<-log((kg/AreaSwept)+1)
ggplot(temp, aes(x=Lon, y=Lat, color=CPU)) +
    theme_linedraw() +
    scale_color_viridis() +
    geom_point(size=1) +
    facet_wrap(~factor(Year))+
    #labs(title = "Calanus", x = "Longitude", y = "Latitude", color = "log((#/m^2)+1)") +
    labs(title = "Sea temperature (20m depth)", x = "Longitude", y = "Latitude", color = "Celsius") +
  #labs(title = "Juvenile sockeye salmon", x = "Longitude", y = "Latitude", color = "log((kg/km^2)+1)") +
  theme_bw() +
    theme(axis.text.x = element_text(size = 12), axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 16),
          plot.title = element_text(size = 20, face = "bold", color = "black"),
          legend.text=element_text(size=12),strip.text.x = element_text(size = 12, face='bold')) 

  gsave(file.path(wd, "Map of Wt.png"), plot=gmap, height=8, width=9, units='in')

  
  g <- gmap +
    geom_point(data=D7, aes(Lon, Lat, color=(D), group=NULL),
               size=1.9, stroke=0,shape=16) + facet_wrap('Year')+
    labs(color=paste('Catch rates ()'))+
    theme(strip.text.x = element_text(size = 12))
  
  g
  
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
dev.off()


# SETTINGS ====================================================================================================

#Calculate day
#summary(Julian)
#DOY<-scale(Julian, center=TRUE) #Julian day scaled from 0 to 1
#DOY


######################
# COVARIATE TIME SERIES MODELS 
#####################
settings = make_settings( Version="VAST_v14_0_1",  
                          n_x=500,
                          Region='User', 
                          purpose="index2",
                          bias.correct=TRUE,
                          max_cells = 2000,
                          Options=c('SD_site_logdensity'=TRUE,'Calculate_Range'=TRUE,'Calculate_effective_area'=TRUE),
                          #ObsModel=c(2,1), # catch, Pollock, pink
                          #ObsModel=c(2,4), # Calanus, SST
                          #ObsModel=c(1,4), # SST not matching
                          #ObsModel=c(1,4), # SST not matching n_x 200
                          ObsModel=c(2,4), # SST not matching n_x 200
                          fine_scale=TRUE,
                          treat_nonencounter_as_zero=TRUE,
                          knot_method='grid',
                          use_anisotropy = FALSE)
#settings$FieldConfig["Omega","Component_2"] = 0 #Pollock
#settings$RhoConfig["Beta2"] = 3
#settings$FieldConfig = matrix( c("IID","IID",0,0,"IID","IID"), byrow=TRUE, ncol=2 )
settings$FieldConfig = c(Omega1 = 0, Epsilon1 = 0,Omega2 = "IID",  Epsilon2 = "IID")
#settings$RhoConfig[c("Beta1","Beta2")] = 0
settings$RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0)
?make_settings
?make_data
settings$Options = c( settings$Options, "report_additional_variables"=TRUE )
user_region <- readRDS('user_region.rds')
?as_units
fit_SST = fit_model( "settings"=settings,
                                "Lat_i"=temp[,'Lat'],
                                "Lon_i"=temp[,'Lon'],
                                "observations_LL" = temp[,c('Lat','Lon')],
                                "t_i"=temp[,'Year'], 
                               # "b_i"=temp[,'SST'], 
                               "b_i"=(as_units(temp[,'SST'], NULL)), 
                               #"b_i"=(as_units(log(temp[,'Calanus']+1), 'count')), 
                               #"b_i"=(as_units(log(temp[,'A0PollockN']+1),'count')), 
                               #"b_i"=(as_units(log(temp[,'J_Pink']+1),'kg')), 
                               #"b_i"=(as_units(temp[,'kg'],'kg')), 
                               #"a_i"=rep(1,nrow(temp)),
                               # "a_i"=rep(.000001,nrow(temp)),#Calanus
                               #"a_i"=(as_units(temp[,'AreaSwept'],'km^2')),
                                "a_i"=(as_units(temp[,'Constant'],'km^2' )),
                                "getsd"=TRUE,
                                "test_fit"=TRUE,
                                "build_model" = TRUE,
                                input_grid=user_region ,
                                getReportCovariance=TRUE,
                               newtonsteps = 1
)

#SST 0011 (24)
#Pink 0011
#Pollock 1100
#Calanus 1011 (23)
#Trying cov31 with more 808 data points w areaswept.
#Add areaswept
?valid_udunits()
?units
check_fit(parameter_estimates, check_gradients = FALSE, quiet = FALSE)
#Check extrappolation grid total area....
colSums(fit_SST$extrapolation_list$a_el)
#344,804 km^2
?make_data
#SAVE DATA =================================================================================
x=summary(fit_SST$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"SST240011noaniso.csv", row.names = TRUE)
#v3 is close 0.911 correlation
#try 
plot_results(fit_SST, n_cells=2000  )

?plot_results
######################
# INTERCEPT MODEL 
#####################

# Make settings (turning off bias.correct to save time for example)
settings0 <- settings
settings0$FieldConfig = matrix( c("IID","IID","IID","IID","IID","IID"), byrow=TRUE, ncol=2 )
#settings0$FieldConfig = c(Omega1 = 0, Epsilon1 = 0, Omega2 = 0, Epsilon2 = 0)

settings0$RhoConfig[c("Beta1","Beta2")] = 3

fit_Int = fit_model( settings = settings0,
                      Lat_i = temp[,'Lat'],
                      "Lon_i"=temp[,'Lon'],
                      "observations_LL" = temp[,c('Lat','Lon')],
                      "t_i"=temp[,'Year'], 
                      "b_i"=(as_units(temp[,'kg'],'kg')), 
                      "a_i"=(as_units(temp[,'AreaSwept'],"km^2")),
                      "getsd"=TRUE,
                      "test_fit"=TRUE,
                      "build_model" = TRUE,
                      input_grid=user_region,
                      getReportCovariance=TRUE)

#SAVE DATA =================================================================================
x=summary(fit_Int$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"Sock_Int.csv", row.names = TRUE)

# PLOT RESULTs =============================================================================
plot_results(fit_Int, n_cells=2000)

######################
# INTERCEPT MODEL 
#####################
settings = make_settings( Version="VAST_v14_0_1",  
                          n_x=500,
                          Region='User', 
                          purpose="index2",
                          bias.correct=TRUE,
                          max_cells = 2000,
                          ObsModel=c(2,1), # catch, Pollock, pink
                          # ObsModel=c(2,4), # Calanus
                          fine_scale=TRUE,
                          treat_nonencounter_as_zero=TRUE,
                          knot_method='grid',
                          use_anisotropy = TRUE)
# Make settings (turning off bias.correct to save time for example)
settings0 <- settings
settings0$FieldConfig["Omega","Component_2"] = 0
settings0$FieldConfig["Epsilon","Component_2"] = 0
settings0$FieldConfig["Omega","Component_1"] = 0
settings0$FieldConfig["Epsilon","Component_1"] = 0
settings$RhoConfig["Beta2"] = 3
user_region <- readRDS('user_region.rds')

fit_IntB = fit_model( settings = settings0,
                  Lat_i = temp[,'Lat'],
                  "Lon_i"=temp[,'Lon'],
                  "observations_LL" = temp[,c('Lat','Lon')],
                  "t_i"=temp[,'Year'], 
                  "b_i"=(as_units(temp[,'kg'],'kg')), 
                  "a_i"=(as_units(temp[,'AreaSwept'],"km^2")),
                  "getsd"=TRUE,
                  "test_fit"=TRUE,
                  "build_model" = TRUE,
                  input_grid=user_region,
                  getReportCovariance=TRUE)

#SAVE DATA =================================================================================
x=summary(fit_IntB$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"Sock_IntB.csv", row.names = TRUE)
# PLOT RESULTs =============================================================================
plot_results(fit_IntB, n_cells=2000)


#####################
# Spatial MODEL
#####################
settings1 = make_settings( Version="VAST_v14_0_1",  
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
#settings1 = make_settings(.)
settings1$FieldConfig["Omega","Component_2"] = 0
settings1$FieldConfig["Epsilon","Component_1"] = 0
settings1$FieldConfig["Epsilon","Component_2"] = 0

settings1$RhoConfig["Beta2"] = 3


settings1$Options = c( settings1$Options, "report_additional_variables"=TRUE )
user_region <- readRDS('user_region.rds')

fit_Spatial = fit_model( "settings"=settings1,
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
colSums(fit_Spatial$extrapolation_list$a_el)

#SAVE DATA =================================================================================
x=summary(fit_Spatial$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"1st LP_sockeye_S.csv", row.names = TRUE)

# PLOT RESULTs =============================================================================
plot_results(fit_Spatial, n_cells=2000  )
#####################
# Spatiotemporal MODEL
#####################

settings = make_settings( Version="VAST_v14_0_1",  
                          n_x=500,
                          Region='User', 
                          purpose="index2",
                          bias.correct=FALSE,
                          max_cells = 2000,
                          
                          #ObsModel=c(2,1), # catch, Pollock, pink
                          ObsModel=c(1,4), # SST
                          
                          fine_scale=TRUE,
                          treat_nonencounter_as_zero=TRUE,
                          knot_method='grid',
                          use_anisotropy = TRUE)
settings$FieldConfig["Omega","Component_2"] = 0
settings$FieldConfig["Epsilon","Component_2"] = 0
settings$RhoConfig["Beta2"] = 3

?make_settings
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
                                getReportCovariance=TRUE,
                                newtonsteps = 1
)

check_fit(parameter_estimates, check_gradients = FALSE, quiet = FALSE)
#Check extrapolation grid total area....
colSums(fit_Spatiotemporal$extrapolation_list$a_el)
?make_data
#SAVE DATA =================================================================================
x=summary(fit_Spatiotemporal$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"Full_sockeye.csv", row.names = TRUE)

plot_results(fit_Spatiotemporal, n_cells=2000  )

#####################
# Covariate time series MODEL
#####################
settings = make_settings( Version="VAST_v14_0_1",  
                          n_x=500, #knots
                          Region='User', 
                          purpose="index2",
                          bias.correct=TRUE,
                          max_cells = 2000,
                         # ObsModel=c(2,4), # SST
                          ObsModel=c(2,1), # catch, Pollock, pink
                          #ObsModel=c(2,4), # calanus
                          fine_scale=TRUE,
                          treat_nonencounter_as_zero=TRUE,
                          knot_method='grid',
                          use_anisotropy = TRUE,
                          #use_anisotropy = FALSE, #SST
                          #FieldConfig = c(Omega1 = 0, Epsilon1 = 0, Omega2 = "IID", Epsilon2 = "IID"), #SST
                          FieldConfig = c(Omega1 = 0, Epsilon1 = 0, Omega2 = "IID", Epsilon2 = "IID"),
                          RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0),
)

#SST (24)

settings$Options = c( settings$Options, "report_additional_variables"=TRUE )
user_region <- readRDS('user_region.rds')
?make_data
fit_SST = fit_model( "settings"=settings,
                 "Lat_i"=temp[,'Lat'],
                 "Lon_i"=temp[,'Lon'],
                 "observations_LL" = temp[,c('Lat','Lon')],
                 "t_i"=temp[,'Year'], 
                 "b_i"=(as_units(temp[,'kg'],"kg")), 
                 
                 #"b_i"=(as_units(temp[,'SST'],"kg")), 
                 #"b_i"=(as_units(log(temp[,'A0PollockN']+1),"count")), 
                 #"b_i"=(as_units(log(temp[,'J_Pink']+1),"kg")), 
                 #"b_i"=(as_units(log(temp[,'Calanus']+1),"count")), 
                 #"b_i"=(as_units(temp[,'Calanus'],"count")), 
                 
                 "a_i"=(as_units(temp[,'AreaSwept'],'km^2')), #Calanus
                 
                 #"a_i"=temp[,'Constant'],
                 "getsd"=TRUE,
                 "test_fit"=TRUE,
                 "build_model" = TRUE,
                 input_grid=user_region ,
                 getReportCovariance=TRUE,
                 newtonsteps=1
)
#USE SST O2 E2 only, 24, n_x=500, no aniso.

require(units)
?as_units 
check_fit(parameter_estimates, check_gradients = FALSE, quiet = FALSE)

#Check extrapolation grid total area....
colSums(fit_SST$extrapolation_list$a_el)
344804

#SAVE DATA =================================================================================
x=summary(fit_SST$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"Sock_211111.csv", row.names = TRUE)
#24 not wroking
# PLOT RESULTs =============================================================================
plot_results(fit_SST, n_cells=2000  )

?plot_results
print(plot_results)

################### Covariate effects w S-T model###############
lnSST<-scale(SST)

#X1<-data.frame(Year, Lat, Lon,lnSST) 
#X1_formula=~lnSST
#X2_formula=~lnSST
X1<-data.frame(Year, Lat, Lon,lnSST) 
X1_formula = ~ bs( lnSST, degree=2, intercept=FALSE)
X2_formula = ~ bs( lnSST, degree=2, intercept=FALSE)
#X1_formula=~SST
#X2_formula=~SST
#nonlinear



X1<-data.frame(Year, Lat, Lon,SST) 
X1_formula = ~ bs( SST, degree=2, intercept=FALSE)
X2_formula = ~ bs( SST, degree=2, intercept=FALSE)



Calanusnls<-scale(log(Calanus+1))
X1<-data.frame(Lat,Lon, Year, Calanusnls) 
X1_formula=~Calanusnls
X2_formula=~Calanusnls
X1_formula = ~ bs( Calanusnls, degree=2, intercept=FALSE)
X2_formula = ~ bs( Calanusnls, degree=2, intercept=FALSE)

Pinknls<-scale(J_Pink)

Pinknls<-scale(log(J_Pink+1))
X1<-data.frame(Lat,Lon, Year, Pinknls) 
X1_formula=~Pinknls
X2_formula=~Pinknls
X1_formula = ~ bs( Pinknls, degree=2, intercept=FALSE)
X2_formula = ~ bs( Pinknls, degree=2, intercept=FALSE)
?bs
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
library(corrplot)
M=cor(dat)
M
cor.mtest(dat)
cor.test(lnSST,Calanusnls)
summary(M)
?cor
corrplot(M, method='number')

#################
# COVARIATE MODELS
#################
settings3 = make_settings( Version="VAST_v14_0_1",
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
#settings3$FieldConfig["Omega","Component_2"] = 0
#settings3$FieldConfig["Epsilon","Component_2"] = 0
#settings3$RhoConfig["Beta2"] = 3
#?make_data

fit_SST = fit_model(
  "settings"=settings3,
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
#X1config_cp = array(2,dim=c(1,1)),	#USE FOR annual indices
# X2config_cp = array(2,dim=c(1,1)), 

x=summary(fit_SST$parameter_estimates$SD, select='report')
print(x)
write.csv(x,"Sock_SST.csv", row.names = TRUE)

# PLOT RESULTs =============================================================================
plot_results( fit_SST_Sock, settings=settings3, plot_set=c(3,11,12,14,15,18,19), n_cells=2000)

head(fit$Report)
library(splines)
SST = seq(6,12,length=100)
X1_formula = ~ bs( SST, degree=2, intercept=FALSE)
X = model.matrix( update.formula(X1_formula, ~.+0), data=data.frame(SST=SST) )
gamma_hat = c(7.1, 2.6)
Y = X %*% gamma_hat 
plot( x=SST, y=Y )

1 - (fit_Spatiotemporal$Report$deviance/fit_Spatial$Report$deviance)

#Pink linear
1 - (fit_Pink$Report$deviance/fit_Int$Report$deviance)
1 - (fit_Pink$Report$deviance/fit_IntB$Report$deviance)
1 - (fit_Pink$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_Pink$Report$deviance/fit_Spatiotemporal$Report$deviance)

#Pink nonlinear
1 - (fit_Pink2$Report$deviance/fit_Int$Report$deviance)
1 - (fit_Pink2$Report$deviance/fit_IntB$Report$deviance)
1 - (fit_Pink2$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_Pink2$Report$deviance/fit_Spatiotemporal$Report$deviance)

#SST linear
1 - (fit_SST1$Report$deviance/fit_Int$Report$deviance)
1 - (fit_SST1$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_SST1$Report$deviance/fit_Spatiotemporal$Report$deviance)

#SST nonlinear
1 - (fit_SST$Report$deviance/fit_Int$Report$deviance)
1 - (fit_SST$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_SST$Report$deviance/fit_Spatiotemporal$Report$deviance)

#SST nonlinear
1 - (fit_SST2$Report$deviance/fit_Int$Report$deviance)
1 - (fit_SST2$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_SST2$Report$deviance/fit_Spatiotemporal$Report$deviance)


#Calanus linear
1 - (fit_Calanus1$Report$deviance/fit_Int$Report$deviance)
1 - (fit_Calanus1$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_Calanus1$Report$deviance/fit_Spatiotemporal$Report$deviance)


#Calanus nonlinear
1 - (fit_Calanus2$Report$deviance/fit_Int$Report$deviance)
1 - (fit_Calanus2$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_Calanus2$Report$deviance/fit_Spatiotemporal$Report$deviance)

#Pink linear
1 - (fit_Pink1$Report$deviance/fit_Int$Report$deviance)
1 - (fit_Pink1$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_Pink1$Report$deviance/fit_Spatiotemporal$Report$deviance)

#Pink nonlinear
1 - (fit_Pink2$Report$deviance/fit_Int$Report$deviance)
1 - (fit_Pink2$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_Pink2$Report$deviance/fit_Spatiotemporal$Report$deviance)

#Pollock linear
1 - (fit_Pollock1$Report$deviance/fit_Int$Report$deviance)
1 - (fit_Pollock1$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_Pollock1$Report$deviance/fit_Spatiotemporal$Report$deviance)

#Pollock nonlinear
1 - (fit_Pollock2$Report$deviance/fit_Int$Report$deviance)
1 - (fit_Pollock2$Report$deviance/fit_Spatial$Report$deviance)
1 - (fit_Pollock2$Report$deviance/fit_Spatiotemporal$Report$deviance)




#########################
# REPLOT density by year#
#########################

## Below shows to you get the model estimate of density, D_gct,
## for each grid (g), category (c; not used here single
## univariate); and year (t); and link it spatially to a lat/lon
## extrapolation point.  You can do this for any _gct or _gc
## variable in the Report.

years <- unique(temp$Year)
nyrs <- length(years)
#Year=rep(years, each=nrow(ak_map))
## Remake map list locally for recreating plots
mdl <- make_map_info(Region = settings3$Region,
                     spatial_list = fit_SST_Pink$spatial_list,
                     Extrapolation_List = fit_SST_Pink$extrapolation_list)
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

#############
# APPENDIX Figures Covariate effect PLOTS
#############

#1st Linear predictor for covariates
str(fit_SST_Pink$Report) #eta1_gct
fit_SST_Pink$Report$eta1_gct
#Densities at each extrapolation grid location
names(fit_SST_Sock$Report)[grepl('_gc|_gct', x=names(fit_SST_Pink$Report))]
D_gt <- fit_SST_Sock$Report$eta1_gct[,1,] # drop the category


#Densities
str(fit_SST$Report) #D_gct
fit_SST$Report$D_gct
#Densities at each extrapolation grid location
names(fit_SST$Report)[grepl('_gc|_gct', x=names(fit_SST$Report))]
D_gt <- fit_SST$Report$D_gct[,1,] # drop the category

dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=c(2002:2018))
#dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=years)
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
#D<-exp(D4$D)

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
head(D7)
write.csv(D7, "Pink_Sock_OUT.csv")
theme_set(theme_bw())

####################################
# 1st LP covariate effect eta1_gct
####################################

str(fit_SST_Sock$Report) #eta1_gct
fit_SST_Sock$Report$eta1_gct
#Densities at each extrapolation grid location
names(fit_SST_Sock$Report)[grepl('_gc|_gct', x=names(fit_SST_Sock$Report))]
D_gt <- fit_SST_Sock$Report$eta1_gct[,1,] # drop the category
#D_gt <- fit$Report$eta2_gct[,1,] # drop the category


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

#Subet D values greater than .1% of D to plot
summary(log(D))
summary(D)
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
write.csv(D7, "Cov_PinkLP_OUT.csv")
theme_set(theme_bw())

0.99
summary(D7$D)

0.001*(320.8573)
D7$D[D7$D<=.320857]<-"NA" #Sockeye 99.9%
D7$D[D7$D<=3.20857]<-"NA" #Sockeye 99%
D9<-D7[D7$D != "NA", ]         # Multiple conditions
str(D9)
#D9<-D7



require("rgdal") 
require("maptools")
require("ggplot2")
require("plyr")

  library(maptools)
  library(rgdal)
  library(sp)
  library(ggplot2)


######################################
# APPRENDIX COVARIATE EFFECTS FIGURES
#####################################
#g <- gmap +
#  geom_point(data=D7, aes(Lon, Lat, color=(as.numeric(exp(D))), group=NULL),
#             size=2, stroke=0,shape=16) + facet_wrap('Year')+
#  labs(color=paste('Covariate \n effect'))+
#  theme(strip.text.x = element_text(size = 6))+
#  geom_point(na.rm=TRUE)

#g

######################################
# APPRENDIX COVARIATE VAST FIGURES
#####################################
#Load density data
temp2<-read.csv("lnPollock_D_OUT.csv")


#ggtitle(expression(paste("in-sample ", R^2, "=0.87, p<0.0001                                                      out-sample ", R^2, "=0.74, p<0.0001"))) + geom_point(shape = 1) +
  logD<-log(as.numeric(D9$D))
g <- gmap +
  geom_point(data=temp2, aes(Lon, Lat, color=(as.numeric(((D)))), group=NULL),
             size=2, stroke=0,shape=16) + facet_wrap('Year')+
  labs(color=expression(paste("ln(#/",km^2,")")))+
  #labs(color=expression(paste("Celsius")))+
  theme(strip.text.x = element_text(size = 12))+
  geom_point(na.rm=TRUE)
g

Cal<-log(Calanus+1)
Poll<-log(A0PollockN+1)
Pink<-log(J_Pink+1)
summary(Cal)
summary(Calanus)
#Pollock data are in lnlnPollock_D_OUT.csv
pdf(file="APPENDIX FIGURE 6A.pdf", height=8.5, width=11) #This generates the figure as a hi-res (pub quality) tiff file. You can change the dimensions of the figure by changing the height and width arguments. You can change the resolution with the 'res' argument.
#Data are in 
g2 <- gmap + 
  theme_bw()+
  #scale_color_gradientn(colours = colorspace::divergingx_hcl(palette="RdGy", n=7, rev=TRUE),oob=scales::squish,limits=c(0,4))+
  
  #geom_polygon(data = temp, aes(x = Lon, y = Lat, group = NULL), fill=NA,color="lightgray" ,size = 0.005) +
  geom_point(data=D7 , aes(Lon, Lat, color=(D), group=NULL),  #Use log(D) for GRP to covert back from exp(GRP)
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
  labs(title="Juvenile pink salmon effect")+
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
                     spatial_list = fit_SST_Sock$spatial_list,
                     Extrapolation_List = fit_SST_Sock$extrapolation_list)
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
  coord_cartesian(xlim=temp$Xlim, ylim=temp$Ylim)+ 
  guides(fill=guide_legend(title='Energy Density'))




# MODEL VALIDATION===================================================================

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
settings$FieldConfig["Omega","Component_2"] = 0
settings$FieldConfig["Epsilon","Component_2"] = 0
settings$RhoConfig["Beta2"] = 3

?make_settings
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

#Remove 2013, 2015, 2017
#R2=0.9757


#Not log transformed
Sockeye_Pred_in<-fit$Report$R2_i-exp(-0.23875180)^2/2					#Predicted in-bag-sample 
Sockeye_Obs_in<-(temp[,'kg']) 						#Observed in-bag-sample 
plot(Sockeye_Pred_in,Sockeye_Obs_in)				#Plot
summary(lm(Sockeye_Obs_in~0+Sockeye_Pred_in))		#R2=0.67

R2=0.8765
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
  settings$FieldConfig["Omega","Component_2"] = 0
  settings$FieldConfig["Epsilon","Component_2"] = 0
  settings$RhoConfig["Beta2"] = 3
  
  
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


#ggtitle("Sockeye VAST model: in $R^2$=0.52, p<0.0001: out $R^2$=0.39, p<0.0001")
#= expression(paste("x axis ", ring(A)^2)), y = "y axis")
?pdf
pdf(file="FIGURE 5.pdf", height=8.5, width=11, pointsize=12) 
pg1 <- ggplot(data3, aes(Sockeye_Pred, Sockeye_Obs))+ylab("")+xlab("")+
  ggtitle(expression(paste("in-sample ", R^2, "=0.87, p<0.0001                                               out-sample ", R^2, "=0.74, p<0.0001"))) + 
  geom_point(shape = 1) +
  facet_grid(~Sample)+ geom_abline(slope = 1) +
  facet_wrap(~Sample, nrow = 1) +  
  theme(aspect.ratio = 1,text = element_text(size=12)) +
  xlim(-1,200)+
  ylim(-1,200)+
  #  geom_point(          )
  geom_smooth(method = "lm", se = T, fullrange = T, colour = "blue", cex = 1)
grid.arrange(pg1,bottom=textGrob(expression(paste("Predicted")),
      gp = gpar(col = "black", fontsize = 12)),left=textGrob(expression(paste("Observed")),
      rot=90,  gp = gpar(col = "black", fontsize = 12)),nrow = 1, ncol=1)
dev.off()

?geom_smooth
traceback(max.lines=30)




#---------------------------
# Annual time series data for figures
#---------------------------
Pollock<-log(temp[,'A0PollockN']+1)
hist(temp[,'A0PollockN'])
hist(Pollock)
Pink<-log(temp[,'J_Pink']+1)
hist(temp[,'J_Pink'])
hist(Pink)
Cal<-log(temp[,'Calanus']+1)
summary(Cal)
hist(temp[,'Calanus'])
hist(Cal)
settings = make_settings( Version="VAST_v14_0_1",
                          n_x=500,
                          Region='User', 
                          purpose="index2",
                          bias.correct=TRUE,
                          max_cells = 2000,
                          ObsModel=c(2,4),#Calanus, SST 
                         # ObsModel=c(2,1),#Pollock
                          
                          fine_scale=TRUE,
                          FieldConfig=c(1,1,1,1),
                         #Cal  (    )
                          #(1100) for pollock
                          #(0011) for SST
                          # Pink (1011)
                          treat_nonencounter_as_zero=TRUE,
                          knot_method='grid',
                          use_anisotropy = TRUE)

user_region <- readRDS('user_region.rds')

#?make_settings
?make_data
fit = fit_model(
  "settings"=settings,
  "Lat_i"=temp[,'Lat'],
  "Lon_i"=temp[,'Lon'],
  "observations_LL" = temp[,c('Lat','Lon')],
  "t_i"=temp[,'Year'],
  # "c_i"=as.numeric(temp[,'Sci'])-1, 
  "c_iz"=rep(0,nrow(temp)), # for single species 
 #"b_i"=(as_units(temp[,'SST'])),
  "b_i"=(as_units(Cal, value="count")),
 # "b_i"=(as_units(temp[,'J_Pink'])),
 # "b_i"=(as_units(log(temp[,'A0PollockN']+1))),
 #"b_i"=Pollock, #nl(Pollock+1)
 #"a_i"=(as_units(temp[,'AreaSwept'],"km^2")),
  "a_i"=rep(0.000001,nrow(temp)),
 #"a_i"=temp[,'Constant'],#Error: ‘/ `km2` * `km2`’ is not a unit recognized by udunits or a user-defined unit
 # "a_i"=as_units(rep(1,nrow(temp))),
 #"a_i"=as_units(1, unitless), #For SST
 #"a_i"=(as_units(temp[,'Constant'])),
 # "a_i"=(as_units(temp[,'Constant'], value=unitless)),

  "getsd"=TRUE,
  # "test_fit"=TRUE,
 input_grid=user_region,
  getReportCovariance=TRUE#,
  #"catchability_data" = data.frame(date=DOY), #Day of year for catchability
  #Q2_formula = ~ date
 )
?as_units
?valid_udunits_prefixes
?check_fit
#TMB::sdreport
#SAVE DATA ========================================================================================
x=summary(fit$parameter_estimates$SD, select='report')#Data are now in fit$Report
print(x)
write.csv(x,"SST_24_0011.csv", row.names = TRUE)

#When looking at covariate effects model output_do I fit linear and nonlinear term separately to get effects plots?
#COVARIATE PLOTS ==================================================================================
plot_results( fit,  plot_set=c(3,11,12,14,15,16,17,18,19), n_cells=1000)
plot_results( fit,  n_cells=2000)


###############################################
# PLOT back transformed VAST covariate densities
##############################################

SST_2002<-as.vector(fit$Report$D_gct[,,1])
SST_2003<-as.vector(fit$Report$D_gct[,,2])
SST_2004<-as.vector(fit$Report$D_gct[,,3])
SST_2005<-as.vector(fit$Report$D_gct[,,4])
SST_2006<-as.vector(fit$Report$D_gct[,,5])
SST_2007<-as.vector(fit$Report$D_gct[,,6])
SST_2008<-as.vector(fit$Report$D_gct[,,7])
SST_2009<-as.vector(fit$Report$D_gct[,,8])
SST_2010<-as.vector(fit$Report$D_gct[,,9])
SST_2011<-as.vector(fit$Report$D_gct[,,10])
SST_2012<-as.vector(fit$Report$D_gct[,,11])
SST_2013<-as.vector(fit$Report$D_gct[,,12])
SST_2014<-as.vector(fit$Report$D_gct[,,13])
SST_2015<-as.vector(fit$Report$D_gct[,,14])
SST_2016<-as.vector(fit$Report$D_gct[,,15])
SST_2017<-as.vector(fit$Report$D_gct[,,16])
SST_2018<-as.vector(fit$Report$D_gct[,,17])
library(tidyr)
SST_fitted<-data.frame(SST_2002,SST_2003,SST_2004,SST_2005,SST_2006,SST_2007,SST_2008,SST_2009,SST_2010,SST_2011,SST_2012,SST_2013,SST_2014,SST_2015,SST_2016,SST_2017,SST_2018)
SST_fit2<-data.frame(unlist(SST_fitted))
Year<-rep(c(2002:2018),each=2000)
head(fit$Report)
print(fit$Report)
Lon<-rep(fit$extrapolation_list$Data_Extrap[,1],each=17)
Lat<-rep(fit$extrapolation_list$Data_Extrap[,2],each=17)


lnSST<-scale(SST_fit2)
X1<-data.frame(Year, Lat, Lon,lnSST) 
X1_formula = ~ bs( lnSST, degree=2, intercept=FALSE)
X2_formula = ~ bs( lnSST, degree=2, intercept=FALSE)
X1_formula=~lnSST
X2_formula=~lnSST

##################
# INTERCEPT MODEL
##################

# Make settings (turning off bias.correct to save time for example)
settings0 <- settings2
settings0$FieldConfig = matrix( c("IID","IID",0,0,"IID","IID"), byrow=TRUE, ncol=2 )

settings0$RhoConfig[c("Beta1","Beta2")] = 3 #21
#settings0$RhoConfig[c("Beta1","Beta2")] = 0 #23, 24

# Run model
fit0 = fit_model( settings = settings0,
                  Lat_i = temp[,'Lat'],
                  "Lon_i"=temp[,'Lon'],
                  "observations_LL" = temp[,c('Lat','Lon')],
                  "t_i"=temp[,'Year'], 
                  "b_i"=(as_units(temp[,'kg'],'kg')), 
                  # "b_i"=(as_units(temp[,'Num'],'count')), 
                  "a_i"=(as_units(temp[,'AreaSwept'],"km^2")),
                  "getsd"=TRUE,
                  "test_fit"=TRUE,
                  "build_model" = TRUE,
                  input_grid=user_region,
                  getReportCovariance=TRUE, 
                  "catchability_data" = data.frame(date=DOY), #Day of year for catchability
                  Q2_formula = ~ date)

###### Calculate percent-deviance-explained
1 - (fit2$Report$deviance/fit0$Report$deviance)
#SST, Pink, Pollock w S and SST 0.69
#Base vs int
#VAST lnCalanus & Sockeye vs intercept 77% 0.7656
#lnCalanus & Sockeye vs intercept 
#VAST ln Pink & Sockeye vs intercept
#Poll vs int

#SAVE DATA ========================================================================================
x=summary(fit2$parameter_estimates$SD, select='report')#Data are now in fit$Report
print(x)
write.csv(x,"lnPolll Sock 21.csv", row.names = TRUE)

#When looking at covariate effects model output_do I fit linear and nonlinear term separately to get effects plots?
#COVARIATE PLOTS ==================================================================================
plot_results( fit2,  plot_set=c(3,11,12,14,15,16,17,18,19), n_cells=2000)
#plot_results( fit,  plot_set=c(3,14), n_cells=1000)

#?plot_results
#WORKS!!!!!!!!!!! 

##################################
# FIGURE 10
# COVARIATE EFFECT PLOT
##################################

library(mgcv)
library(effects)
library(gridExtra)
#Rerun each respective model

pdf(file="FIGURE 10.pdf", height=8.5, width=11) 
#Nonlinear SST effect
covariate_data_full = fit_SST$effects$covariate_data_full
catchability_data_full = fit_SST$effects$catchability_data_full
pred3 = Effect.fit_model( fit_SST,
                          focal.predictors = c("SST"),
                          which_formula = "X1", 
                          xlevels = 100)
plot3<-plot(pred3,cex.lab=4, xlab="Sea surface temperature (Celsius)", ylab=NULL, main="")
#Linear pink effect
covariate_data_full = fit_SST_Pink$effects$covariate_data_full
catchability_data_full = fit_SST_Pink$effects$catchability_data_full
pred4 = Effect.fit_model( fit_SST_Pink,
                          focal.predictors = c("Pinknls"),
                          which_formula = "X1", 
                          xlevels = 100)
plot4<-plot(pred4,cex.lab=4, xlab="Juvenile pink salmon (ln(kg+1))", ylab=NULL, main="")
grid.arrange(plot3, plot4,bottom=textGrob(expression(paste("Covariates"))),left=textGrob(expression(paste("Linear predictor")),rot=90),nrow = 1, ncol=2)
dev.off()
?Effect.fit_model
#----------------------------------------------
# Covariate effects 
#----------------------------------------------

pdf(file="FIGURE 9new.pdf", height=8.5, width=11) 
#plot1<-plot(pred1, ylab="", xlab="Sea temperature (Celsius)", main="Encounter probability")
#plot2<-plot(pred2, ylab="", xlab="Sea temperature", main="Positive catch rate")
plot3<-plot(pred3, ylab="", xlab="Sea temperature"^{2}~"", main="Positive catch rate")
plot4<-plot(pred4, ylab="", xlab="log(Pink salmon(kg)+1)", main="Positive catch rate")
#plot6<-plot(pred6, ylab="", xlab="Pink salmon", main="Positive catch rate")
#plot7<-plot(pred7, ylab="", xlab="log(Age-0 pollock(#)+1)", main="Encounter probability")
#plot8<-plot(pred8, ylab="", xlab="log(Age-0 pollock(#)+1)", main="Positive catch rate")
#grid.arrange(plot1, plot5,plot7,plot8,bottom=textGrob(expression(paste("Covariate"), cex=3),gp = gpar(col = "black", fontsize = 16)),left=textGrob(expression(paste("Linear predictor"), cex=3),gp = gpar(col = "black", fontsize = 16),rot=90),nrow = 2, ncol=2)
grid.arrange(plot3,plot4,bottom=textGrob(expression(paste("Covariate"), cex=3),gp = gpar(col = "black", fontsize = 16)),left=textGrob(expression(paste("Linear predictor"), cex=3),gp = gpar(col = "black", fontsize = 16),rot=90),nrow = 2, ncol=2)

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



# EXTRACT DATA AND PLOT EXTERNALLY ===============================================================
#file:///C:/Users/ellen.yasumiishi/Work/R/win-library/3.5/FishStatsUtils/html/00Index.html
# https://github.com/James-Thorson-NOAA/VAST/wiki/Plots-using-ggplot

## A quick demonstration of how to extract map quantities and
### plot them externally. Cole Monnahan | May 2021

library(VAST)                           # 3.8.0
library(ggplot2)                        # 2.10.0
library(dplyr)
library(tidyr)
library(tidyverse)
library(lubridate)
theme_set(theme_bw())


### Run a simple VAST model for a few years
#example <- load_example( data_set="EBS_pollock" )
#dat <- subset(temp$Year, Year>2008)
output <- temp %>%
  filter(Year %in% c(2002,2003,2004,2005, 2006,2007,2008,2009,2010,2011,2012,2014,2016,2018))
dat<-output
list(output$Year)
years <- unique(dat$Year)
nyrs <- length(years)



settings = make_settings( Version="VAST_v13_0_0",  #12
                          n_x=50,
                          Region='User', 
                          purpose="index2",
                          bias.correct=FALSE,
                          max_cells = 2000,
                          ObsModel=c(2,3), 
                          #Options=c('SD_site_logdensity'=1),
                          fine_scale=TRUE,
                          FieldConfig=c(1,1,1,1),	
                          treat_nonencounter_as_zero=TRUE,
                          knot_method='grid')
settings$grid_size_km = 60
user_region <- readRDS('user_region.rds')

fit <- fit_model(settings = settings,
                 Lat_i = dat$Lat,
                 Lon_i = dat$Lon,
                 t_i = dat$Year,
                 getsd = FALSE,
                 b_i = dat$Wt,
                 a_i = dat$AreaSwept,
                 "maximum_distance_from_sample" = 42, #42 in original file, 75 new
                                  input_grid=user_region)

## Remake map list locally for recreating plots
mdl <- make_map_info(Region = 'User',
                     spatial_list = fit$spatial_list,
                     Extrapolation_List = fit$extrapolation_list)
## quick dirty AK map
ak_map <- subset(map_data("world"), region=='USA' & subregion=='Alaska')
## Have to duplicate it for each year so can facet below
ak_map <- cbind(ak_map[rep(1:nrow(ak_map), times=nyrs),],
                Year=rep(years, each=nrow(ak_map)))

gmap <- ggplot(ak_map, aes(x = dat$Lon, y = dat$Lat, group = group)) +
  geom_polygon(fill="black", colour = "white") +
  scale_color_viridis_c(option = "magma") +
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
  coord_cartesian(xlim=mdl$Xlim, ylim=mdl$Ylim)


## Below shows to you get the model estimate of density, D_gct,
## for each grid (g), category (c; not used here single
## univariate); and year (t); and link it spatially to a lat/lon
## extrapolation point.  You can do this for any _gct or _gc

## variable in the Report.
names(fit$Report)[grepl('_gc|_gct', x=names(fit$Report))]
names(fit$Report)

D_gt <- fit$Report$D_gct[,1,] # drop the category, DENSITIES FOR EVERY YEAR

###################
#cAN USE PLOT VARIABLE ON THIS FIT 
#fit$Report$D_gct[,1,]
#200 ROWS, T IS THE RANGE OF YEARS. 
#sUBSET FOR THE 2 COLUMNS i WANT FOR YEARS
########################
1:nrow(D_gt)
dimnames(D_gt) <- list(cell=1:nrow(D_gt), year=years)
## tidy way of doing this, reshape2::melt() does
## it cleanly but is deprecated
D_gt <- D_gt %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>%
  pivot_Loner(-cell, names_to = "Year", values_to='D')
D <- merge(D_gt, mdl$PlotDF, by.x='cell', by.y='x2i')
g <- gmap +
  geom_point(data=D, aes(Lon, Lat, color=log(D), group=NULL),
             ## These settings are necessary to avoid
             ## overlplotting which is a problem here. May need
             ## to be tweaked further.
             size=.3, stroke=0,shape=16) + facet_wrap('Year')
g

g <- gmap +
  geom_point(data=D, aes(Lon, Lat, color=(D), group=NULL),
             ## These settings are necessary to avoid
             ## overlplotting which is a problem here. May need
             ## to be tweaked further.
             size=.3, stroke=0,shape=16) + facet_wrap('Year')
g

## Or you can do it in base R, building your own palette and
## looping through years as needed.


#PLOT....jt

#Plot_map samples from fitted model to calculate SD

plot( fit, plot_set = 3, plot_value = sd )
out = plot( fit, plot_set = 3, plot_value = sd )


# FIGURE 3 USE  MAP OF STUDY AREA ============================================================================
.libPaths()


install.packages("devtools")
require(devtools)
require(marmap)
devtools::install_github("MattCallahan-NOAA/akmarineareas2", force=TRUE, lib="C:/Program Files/R/R-4.2.2/library")

install.packages("grid", force=TRUE, lib="C:/Program Files/R/R-4.2.2/library")
#C:\Users\ellen.yasumiishi\AppData\Local\R\win-library\4.2
#Click on manually in Packages
#remove.packages('marmap')
install.packages("marmap")
install.packages("tidyverse")
install.packages('raster')
install.packages("ggplot2")
install.packages("rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("remotes")
install.packages("urca")
install.packages("rgdal")

library(grid)
library(akmarineareas2)
library(sf)
library(tidyverse)
library(marmap)
library(dplyr)
library(sp)
library(rnaturalearthdata)
library(raster)
library(rgdal)
library(urca)
library(remotes)
library(ggplot2)
library(maps)
library(mapdata)
library(mapproj)

citation("akima")
citation("raster")
citation("grid")
citation("sf")
citation("AKmarineareas")
citation("tidyverse")
citation("marmap")
citation("ggplot2")
citation("rnaturalearth")

#Check package versions
installed.packages()
sessionInfo()

# Set your working directory
wd<-setwd("C:/Users/ellen.yasumiishi/Work/2023/Sockeye")

#Some eastern Bering Sea BASIS data, plot station locations Lon_rnd and Lat_rnd
temp<-read.csv("IYSsockeyeCovariates30.csv") 
Year=temp[,1]	
Lat	=temp[,2]
Lon	=temp[,3]
SST	=temp[,4]
A0PollockN =temp[,5]	
J_Pink =temp[,6]	
Calanus	=temp[,7]
AirT	=temp[,8]
HaulDate	=temp[,9]
StationID	=temp[,10]
AreaSwept	=temp[,11]
kg	=temp[,12]
Sci	=temp[,13]
vessel	=temp[,14]
Julian	=temp[,15]
Constant=temp[,16]
Lat_rnd=temp[,17]
Lon_rnd=temp[,18]
#If you're using decimal degrees rather than an Alaska centered projection try

ak<-ak %>%
st_transform(crs=4326)
use_mit_license()

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
require(raster)
depth_c<-rasterToContour(depth, levels=c(-50,-100,-150))%>%
  st_as_sf()
grid.newpage();grid.draw(roundrectGrob(gp = gpar(lwd = 0)))

#plot
#install.packages("ggspatial")
require(ggspatial) # Add a scale bar

#Add line direction
#https://stackoverflow.com/questions/70249589/how-to-add-multiple-arrows-to-a-path-according-to-line-direction-using-ggplot2
pdf(file="FIGURE 3.pdf", height=11, width=8.5) #This generates the figure as a hi-res (pub quality) tiff file. You can change the dimensions of the figure by changing the height and width arguments. You can change the resolution with the 'res' argument.
#?pdf
#jpeg(file="FIGURE 3.jpg", height=11, width=8.5, res=400, units="in") #This generates the figure as a hi-res (pub quality) tiff file. You can change the dimensions of the figure by changing the height and width arguments. You can change the resolution with the 'res' argument.
p<-ggplot()+
  geom_sf(data=russia)+
  geom_sf(data=ak)+
  geom_sf(data=depth_c, aes(color=level), colour=c("grey35","grey50","grey75"))+
  geom_point(data=temp, aes(Lon_rnd, Lat_rnd, group=NULL), shape=4)+
  coord_sf(xlim=c(-178, -150), ylim=c(50,70))+
  theme(panel.background=element_rect(fill='white', linewidth=1),axis.title = element_text(size = 24),axis.text = element_text(size = 16),
        panel.border = element_rect(colour = "black",fill=NA, linewidth=1))+
   xlab(expression(paste(Longitude)))+
  ylab(expression(paste(Latitude)))+
  geom_text(aes(x=-169.3, y=57.7, label="Pribilof Islands",angle=0), size=4)+
  geom_text(aes(x=-169.3, y=64, label="St. Lawrence Is.",angle=0), size=4)+
  geom_text(aes(x=-172.9, y=61, label="St. Matthews Is.",angle=0), size=4)+
  geom_text(aes(x=-164.2, y=54.73, label="Unimak Pass.",angle=0), size=4)+
  geom_text(aes(x=-167.1, y=60.6, label="Nunivak Is.",angle=0), size=4)+
  #geom_rect(data = temp, aes(xmin = -169 -.4, xmax = -169 + .4, ymin = 60.1 - .4, ymax = 60.1 + .4), fill = "grey80") +
  geom_label(aes(x = -169.7, y =60.1, label = "50 m"), fill = "white", label.size = NA, size = 3) +
  geom_label(aes(x = -174.1, y =60.1, label = "100 m"), fill = "white", label.size = NA, size = 3) +
  geom_label(aes(x = -178.1, y =60.1, label = "200 m"), fill = "white", label.size = NA, size = 3) +
  geom_text(aes(x=-168.2, y=62, label="Inner",angle=0), size=5)+
  geom_text(aes(x=-174, y=62, label="Middle",angle=0), size=5)+
  geom_text(aes(x=-177.8, y=62, label="Outer",angle=0), size=5)+
  geom_text(aes(x=-155, y=66, label="Alaska"), size=7)+
  geom_text(aes(x=-176, y=66, label="Russia"), size=7)+
  geom_text(aes(x=-175.1, y=55, label="Bering Sea"), size=7)+
  geom_text(aes(x=-168, y=70, label="Chukchi Sea"), size=6)+
  geom_text(aes(x=-155, y=53, label="Gulf of Alaska"), size=6)+
  geom_text(aes(x=-173, y=53, label="Aleutian Islands"), angle=20,size=6)+
  geom_text(aes(x=-165, y=51, label="Pacific Ocean"), size=7)+
  annotation_scale()+
  annotation_north_arrow( height= unit(.75, "cm"), width= unit(.75, "cm"), pad_x = unit(.5, "cm"),
                           pad_y = unit(1, "cm"))
#?geom_curve

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

#Help
?arrow
?geom_text
?geom_sf
?annotation_scale
?annotation_north_arrow

#FIGURE 2 AIR TEMPERATURE BAR PLOT ==============================================================
df<-read.csv("AirT.csv") #2002-2018 Made up data for 2013, 2015, 2017, just matches 2012 data
#Includes dummy data for 2013, 2015, 2017 using 2012 data
print(df)
head(df)
tail(df)
yr	=df[,1]	#Survey year
AirTemp=df[,2]#	St Paul airport air temperature
library(ggthemes)
library(tidyverse)

pdf(file="FIGURE 2.pdf", height=8.5, width=11) 
anomaly<-(AirTemp-mean(AirTemp))
p<-ggplot(data=df, aes(x=yr, y=anomaly)) +
geom_bar(stat = 'identity', aes(fill = anomaly>0), position = 'dodge', col = 
'transparent', show.legend=FALSE) + 
theme_bw() + 
  scale_fill_manual(values=c("blue", "red")) +
#scale_fill_discrete(guide = 'none') + 
  xlab("Year") + ylab("Anomaly of summer air temperature (Fahrenheit)")
p+ coord_cartesian(ylim=c(-6,6), xlim=c(2000,2020))
dev.off()


#  FIGURE 4: DIET PLOT =====================================================================
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
pdf(file="FIGURE 4.pdf", height=8, width=11) 
p1 <- diet.Lon %>%
  mutate(PreyItem = fct_relevel(PreyItem,
                                "Age-0 Pollock", "Fish", "Euphausiids", 
                                "Amphipods", "Calanus spp.", "Large Copepods", "Small Copepods", 
                                "Other Crustaceans", "Arrow Worms", "Pteropods", "Other")) %>%
  ggplot(aes(Year,PreyItem,fill=SCI))+
  geom_tile(color= "white",size=0.1) + 
  scale_fill_viridis(name="%SCI",option ="C")+
  theme_minimal(base_size = 8)+

  theme(text = element_text(size = 12))+
  ylab("Prey Item")+
  xlab("Year")
dev.off()

#Clean up plot and export using res that Ellen had previously set, though I changed the aspect ratio to better fit this style of 
#plot. 
pdf(file="FIGURE 6.pdf", height=8, width=11) 
p1+theme(axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),axis.title.y = element_text(size=12),axis.title.x = element_text(size=12),legend.title = element_text(size = 12),
         legend.text = element_text(size = 12),strip.background = element_rect(colour="white"),axis.ticks=element_blank()) 
dev.off()




# APPENDIX FIGURES===============================

#Heatmap code for 12 maps (you can change this, of course!!)

#Step 1.
#Install packages or load libraries
install.packages('fields') #you will need to install each of the libraries below the first time you sue them. After they have been installed, you'll just need to load them using the lines of code below.

library(fields) #for heat map colors
library(maps) #for plotting map region
library(akima) #for interp function
library(mapdata) #for adding landmass overlay
#require(Rcolombos) #for interp function
library(RColorBrewer) #for adding landmass overlay
library(grDevices)
library(ggplot2)
install.packages("interp") #for interp function
require(mapdata)
require(interp)

head(temp)
SalAbund<-log(temp$kg+1)
SalAbund<-(temp$SST)
SalAbund<-log(temp$Calanus+1)
SalAbund<-log(temp$J_Pink+1)
SalAbund<-log(temp$Age0_Pollock+1)
library
#Step 2.
#Create a personalized heatmap color palette. I grabbed 4 colors, but you can use whatever you like. I use page 3 in this pdf: https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf
hmcols<-colorRampPalette(c("darkblue", "blue", "lightblue","lightgreen","yellow","orange","red"))(50)


#Step 3.
#Set up the maximum latitude and Lonitude extent of your data (min and max across all years). The "length=400" argument can be changed to more/less and affects the pixel size of the heatmap interpolation.  Higher length values will make the interpolation look smoother, and vice versa
XO<-seq(min(Lon), max(Lon), length=600)
YO<-seq(min(Lat), max(Lat), length=600)

tiff(file="FIGURE_8_Sockeye.tiff", height=11, width=8.5, units='in', pointsize=16,compression="lzw", res=1200) 

plot.new()
par(mfrow=c(5,3)) #This makes 3 rows and 4 columns....modify as needed
par(mar=c(3,3,1,0.5)) #(5,5,1,0)
par(oma=c(3,3,0,0))

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2003") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2003"], Lat[Year=="2003"], SalAbund[Year=="2003"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2003"],Lat[Year=="2003"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2004") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2004"], Lat[Year=="2004"], SalAbund[Year=="2004"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2004"],Lat[Year=="2004"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2005") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2005"], Lat[Year=="2005"], SalAbund[Year=="2005"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2005"],Lat[Year=="2005"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1, cex.main=1.5,main="2006") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2006"], Lat[Year=="2006"], SalAbund[Year=="2006"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2006"],Lat[Year=="2006"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1, cex.main=1.5,main="2007") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2007"], Lat[Year=="2007"], SalAbund[Year=="2007"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2007"],Lat[Year=="2007"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2008") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2008"], Lat[Year=="2008"], SalAbund[Year=="2008"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2008"],Lat[Year=="2008"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2009") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2009"], Lat[Year=="2009"], SalAbund[Year=="2009"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2009"],Lat[Year=="2009"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2010") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2010"], Lat[Year=="2010"], SalAbund[Year=="2010"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2010"],Lat[Year=="2010"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2011") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2011"], Lat[Year=="2011"], SalAbund[Year=="2011"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2011"],Lat[Year=="2011"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2012") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2012"], Lat[Year=="2012"], SalAbund[Year=="2012"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2012"],Lat[Year=="2012"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2014") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2014"], Lat[Year=="2014"], SalAbund[Year=="2014"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2014"],Lat[Year=="2014"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2016") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2016"], Lat[Year=="2016"], SalAbund[Year=="2016"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2016"],Lat[Year=="2016"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')

plot(c(-173.1,-159),c(54.9,60), type='n',axes=F, xlab='', ylab='', cex.lab=0.001,cex.axis=1,las=1,cex.main=1.5, main="2018") #you can change the lat/Lon limits, as needed
axis(2, tick = TRUE, labels = TRUE, cex.axis=1, las=1) #This is to get latitude labeled on all of the lefthand-most maps
axis(1, tick = TRUE, labels = TRUE, cex.axis=1) #This is to get latitude labeled on all of the lefthand-most maps
image(interp(Lon[Year=="2018"], Lat[Year=="2018"], SalAbund[Year=="2018"], duplicate='mean', xo=XO, yo=YO),extrap=FALSE, col=hmcols, zlim=c(min(SalAbund), max(SalAbund)), add=T)  
maps::map('worldHires',fill=T,xlim=c(-173.1,-159), ylim=c(54.9,60),add=T, col="grey67")
points(Lon[Year=="2018"],Lat[Year=="2018"],pch=16, cex=0.5) #If you want to plot sampling stations
box(lty = 'solid', col = 'black')
plot.new()

image.plot(zlim=c(min(SalAbund), max(SalAbund)),legend.only=TRUE, legend.cex=0.1,legend.shrink=0.9, smallplot=c(0.8,0.99,0.1,0.9), axis.args=c(cex.axis=1), col=hmcols)

mtext(expression(paste("Latitude")), outer=T, side=2, line=1, cex=1.5)
mtext(expression(paste("Longitude")), outer=T, side=1, line=1, cex=1.5)
dev.off()

#Step 6.
#Make a scalebar

tiff(file="Sockeye_scalebar.tiff", height=30, width=15, units='cm', compression="lzw", res=200)
plot.new()
par(mar=c(0,0,0,0))
image.plot(zlim=c(min(SalAbund), max(SalAbund)),legend.only=TRUE, legend.lab="", smallplot=c(0.5,0.6,0.1,0.9), axis.args=c(cex.axis=3), col=hmcols)
dev.off()


#mtext(expression(paste("Latitude")), outer=T, side=2, line=1, cex=1.5)
library(maps)
library(geosphere)
library(dplyr)
maps::map('world',
    col="#b3b3b3", fill=TRUE, bg="white", lwd=0.05,
    mar=rep(0,4),border=0, ylim=c(-80,80) 
)
points(x=cldrd$longitude, y=cldrd$latitude, col=c("magenta"), cex=(.7), pch=16)    
points(x=outlrd$longitude, y=outlrd$latitude, col=c("black"), cex=(.2), pch=16)


################
# Matt Callahan code
#################
require(devtools)
devtools::install_github("MattCallahan-NOAA/AKmarineareas")

###################
# FIGURE 1 MAP OF STUDY AREA
###################

library(ggplot2)
library(maps)
library(mapdata)
install.packages("mapproj")
library(mapproj)
library(sf)
library(marmap)
library(tidyverse)
library(rnaturalearth)
tiff(file="FIGURE 2.tiff", height=28, width=21, units='cm', compression="lzw", res=1200) #This generates the figure as a hi-res (pub quality) tiff file. You can change the dimensions of the figure by changing the height and width arguments. You can change the resolution with the 'res' argument.

ak<-map_data('worldHires','USA:Alaska')
akmap<-ggplot()+geom_polygon(data=ak,aes(long,lat,group=group),fill="gray",color="black")+
  #theme(panel.background=element_rect(fill='white'),
  #      panel.border = element_rect(colour = "black", fill=NA, size=2))+
  xlab(expression(paste(Longitude^o,~'W')))+
  ylab(expression(paste(Latitude^o,~'N')))+
  coord_map(xlim = c(-173, -155),ylim = c(50, 66.5))+
  geom_text(aes(x=-160, y=62, label="Alaska"), size=6)+
  geom_text(aes(x=-168, y=58.5, label="Bering Sea"), size=6)+
  geom_text(aes(x=-165, y=55.5, label="Aleutian Chain"), size=6, angle=35)+
  geom_text(aes(x=-164, y=51, label="Pacific Ocean"), size=7)+
  geom_text(aes(x=-170, y=66, label="Bering"), size=6)+
  geom_text(aes(x=-170, y=65.5, label="Strait"), size=6)
akmap


dev.off()
#warnings
#use linewidth and not element_rect()
lifecycle::last_lifecycle_warnings()

###############
# FIGURE 1 MAP 2
################

#Get data
remove.packages("marmap")
update.packages("marmap")
install.packages("remotes")
library(remotes)
.libPaths()

remotes::install_local( "C:\\Users\\ellen.yasumiishi\\Desktop\\marmap-master.zip")
#Say no to updating packages.
install.packages("rgdal")

  devtools::install_github("ericpante/marmap", lib="C:\\Users\\ellen.yasumiishi\\AppData\\Local\\R\\win-library\\4.2")
install.packages("raster")
C:/Program Files/R/R-4.2.2/library
devtools::install_github("ericpante/marmap", lib="C:\\Program Files\\R\\R-4.2.2\\library")

library(rgdal)
library(raster)
library(ggplot2)
library(mapdata)
library(marmap)
library(sp)
library(maps)
library(sf)
# get bathymetry data
install.packages("githubinstall", lib="C:\\Users\\ellen.yasumiishi\\AppData\\Local\\R\\win-library\\4.2")
library(githubinstall)

#Check library paths
.libPat
sessionInfo()
hs()
#Check packages and R versions.

b = marmap::getNOAA.bathy(lon1=-175,  lon2=-155, lat1=50, lat2=66, resolution = 2, keep=TRUE)

#Install older version of raster try 3.4-13. https://github.com/ericpante/marmap/issues/25#issuecomment-965834490
install.packages("versions")
library(versions)
remove.packages('raster')
versions::install.versions('raster', '3.4-13')
library(raster)


marmap::getNOAA.bathy(lon1, lon2, lat1, lat2, resolution = 4,
              keep = FALSE, antimeridian = FALSE, path = NULL)


LON1<-round(ak$long,2)
LAT1<-round(ak$lat,2)

b <- marmap::getNOAA.bathy(floor(min(round(ak$long,2))
),ceiling(max(round(ak$long,2))),
floor(min(round(ak$lat,2)) ),ceiling(max(round(ak$lat,2)))  ,
resolution = 1, keep = TRUE)

b <- marmap::getNOAA.bathy(floor(min(LON1)
),ceiling(max(LON1)),
floor(min(LAT1) ),ceiling(max(LAT1))  ,
resolution = 1, keep = TRUE)

remove.packages('marmap')
library(devtools)
install_github("ericpante/marmap")

packageVersion("marmap")
remove.packages("marmap")

library(marmap) ; library(mapdata)

# Get bathymetric data
dat <- getNOAA.bathy(-90,-70,-20,-2,res=4, keep=TRUE)

# Create nice color palettes
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

## First option for plotting
plot(dat, land=TRUE, n=100, lwd=0.03)
map("worldHires", res=0, add=TRUE)

# Second option
plot(dat, im=TRUE, land=TRUE, bpal=list(c(min(dat),0,blues),c(0,max(dat),greys)), lwd=.05, las=1 )
map("worldHires", res=0, lwd=0.7, add=TRUE)

# Add -200m and -1000m isobath
plot(dat, deep=-200, shallow=-200, step=0, lwd=0.5, drawlabel=TRUE, add=TRUE)
plot(dat, deep=-1000, shallow=-1000, step=0, lwd=0.3, drawlabel=TRUE, add=TRUE)




## Querying NOAA database ...
## This may take seconds to minutes, depending on grid size
## Building bathy matrix ...
library(marmap)
a <- marmap::getNOAA.bathy(-10,10,-20,0)

# convert bathymetry to data frame
bf = fortify.bathy(b)

# get regional polygons
reg = map_data("world2Hires")
reg = subset(reg, region %in% c('Canada', 'USA'))

# convert lat longs
reg$long = (360 - reg$long)*-1

# set map limits
lons = c(-67.5, -63.5)
lats = c(42, 45)

# make plot
ggplot()+
  
  # add 100m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-100),
               size=c(0.3),
               colour="grey")+
  
  # add 250m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-250),
               size=c(0.6),
               colour="grey")+
  
  # add coastline
  geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
               fill = "darkgrey", color = NA) + 
  
  # add polygon
  geom_polygon(data = ply, aes(x = lon, y = lat),
               color = "black", alpha = 0.3) +
  
  # add line
  geom_path(data = lin, aes(x = lon, y = lat),
            colour = "black", alpha = 1, size=0.3)+
  
  # add points
  geom_point(data = pts, aes(x = lon, y = lat),
             colour = "black", fill = "grey", 
             stroke = .5, size = 2, 
             alpha = 1, shape = 21)+
  
  # configure projection and plot domain
  coord_map(xlim = lons, ylim = lats)+
  
  # formatting
  ylab("")+xlab("")+
  theme_bw()

##########################
#
##########################

library(ggplot2)
library(mapdata)

install.packages("oce")
install.packages("ocedata")
library(oce)
library(ocedata)
install.packages(lib=C:/Program Files/R/R-4.2.2/library)
.libPaths()
library(marmap)

# get bathymetry data

b = marmap::getNOAA.bathy(floor(lon1 = -175), ceiling(lon2 = -155), floor(lat1 = 50), ceiling(lat2 = 70), 
                  resolution = 1)
## Querying NOAA database ...
## This may take seconds to minutes, depending on grid size
## Building bathy matrix ...
#ERROR Error in if (ncol(x) == 3 & !exists("bathy", inherits = FALSE)) { : argument is of length zero
#SOLUTION need to reinstall rgdal and raster

# make a simple track line
lin = data.frame(
  lon = c(-65.17536, -65.37423, -65.64541, -66.06122, -66.15161),  
  lat = c(43.30837, 42.94679, 42.87448, 42.92871, 42.72985)
)

# make a few points
pts = data.frame(
  lon = c(-65.3, -65.7, -64.1),
  lat = c(43.4, 43, 42.9)
)

# build a polygon (in this case the 'Roseway Basin Area To Be Avoided')
ply = data.frame(
  lon = c(-64.916667,-64.983333,-65.516667, -66.083333),
  lat = c(43.266667,  42.783333, 42.65, 42.866667)
)
#-------------------------------
data("coastlineWorldFine")

# convert bathymetry
bathyLon = as.numeric(rownames(b))
bathyLat = as.numeric(colnames(b))
bathyZ = as.numeric(b)
dim(bathyZ) = dim(b)

# define plotting region
mlon = mean(pts$lon)
mlat = mean(pts$lat)
span = 300

# plot coastline (no projection)
plot(coastlineWorldFine, clon = mlon, clat = mlat, span = span)

# plot bathymetry
contour(bathyLon,bathyLat,bathyZ,
        levels = c(-50, -100, -150, -200, -250),
        lwd = c(1, 1, 2, 2, 3),
        lty = c(3, 1, 3, 1, 3),
        drawlabels = F, add = TRUE, col = 'darkgray')

# add depth legend
legend("bottomright", seg.len = 3, cex = 0.8,
       lwd = c(1, 1, 2, 2, 3),
       lty = c(3, 1, 3, 1, 3),
       legend = c("50", "100", "150", "200", "250"),
       col = 'darkgray', title = "Depth [m]", bg= "white")

# add map data
points(pts, pch = 16, col = 'red')
lines(lin, col = 'blue')
polygon(ply, lty = 2)


plot.new()
# convert bathymetry to data frame
  # get bathymetry data
  b = getNOAA.bathy(lon1 = -175, lon2 = -155, lat1 = 50, lat2 = 70, 
                    resolution = 1)
bf = fortify.bathy(b)

# get regional polygons
reg = map_data("world2Hires")
reg = subset(reg, region %in% c('Canada', 'USA'))


# convert lat longs
reg$long = (360 - reg$long)*-1

# set map limits
lons = c(-175, -155)
lats = c(50, 70)

# make plot
ggplot()+
  # add 50m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-200),
               size=c(0.3),
               colour="black")+

  # add coastline
  geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
               fill = "darkgrey", color = "black") + 
  
  # configure projection and plot domain
  coord_map(xlim = lons, ylim = lats)+
  
  # formatting
  ylab("")+xlab("")+
  theme_bw()+

  geom_text(aes(x=-162, y=62, label="Alaska"), size=6)+
  geom_text(aes(x=-168, y=58.5, label="Bering Sea"), size=6)+
  geom_text(aes(x=-165, y=55.5, label="Aleutian Chain"), size=6, angle=35)+
  geom_text(aes(x=-164, y=51, label="Northeastern Pacific Ocean"), size=7)+
  geom_text(aes(x=-170, y=66, label="Bering"), size=6)+
  geom_text(aes(x=-170, y=65.5, label="Strait"), size=6)+
  geom_text(aes(x=-164, y=70, label="Chukchi Sea"), size=6)

legend("bottomright", seg.len = 3, cex = 0.7,
       lwd = c(2),
       lty = c(1),
       legend = c( "200"),
       col = 'darkgray', title = "Depth [m]", bg = "white")
  dev.off()
  
 
  # CREATE USER REGION FOR NEW SURVEYS: DONE for BASIS=======================================================================================
  
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
  saveRDS(region, file = "user_region2.rds")
  ### End of creating user extrapolation region object
  ### --------------------------------------------------
  
  ### Quick plots of the process for method 1
  png('user_region2.png', width=7, height=7, units='in', res=200)
  par(mfrow=c(2,2))
  with(region_extent, plot(long, lat, main='Extent in points in LL'))
  plot(region_polygon, main='Polygon in UTM', axes=TRUE)
  plot(region_grid, col=ifelse(is.na(region_df$Id), 'red', 'black'),
       axes=TRUE, main='Extrapolation area UTM')
  with(region, plot(Lon, Lat, main='Extrapolation region in LL', pch='.'))
  dev.off()
  
  
  
  temp<-read.csv("IYSsockeyeCovariates22.csv") #No data added for missing years.
  ### Making data object
  #Error in FishStatsUtils::make_covariates(formula = X1_formula, covariate_data = covariate_data,  : 
  #Year 2013 not found in `covariate_data` please specify covariate values for all years
  
  #temp<-read.csv("IYSsockeyeCovariates23.csv") 
  #2013, 2015, 2017 are estimated as catch plus 1 from the previous year 2012, 2014, 2016. Zero catch are set as 0.
  
  temp<-read.csv("IYSsockeyeCovariates26.csv") 
  #2013, 2015, 2017 are estimated as the average catch at each station for all years, again plus .01, again minus .01. Zero catch are set as 0.
  
  
  #temp<-read.csv("IYSsockeyeCovariates21.csv") 
  #temp<-read.csv("IYSsockeyeCovariates18.csv") #original data 
  
  #temp<-read.csv("IYSsockeyeCovariates24.csv") 
  #2013, 2015, 2017 are estimated as catch plus 0.1 from the previous year 2012, 2014, 2016. Zero catch are set as 0.
  
  head(temp)
  Year=temp[,1]
  Lat=temp[,2]	
  Lon=temp[,3]	
  SST=temp[,4]	
  Age0_Pollock=temp[,5]	
  J_Pink=temp[,6]	
  Calanus=temp[,7]	
  AirT	=temp[,8]
  HaulDate=temp[,9]	
  StationID=temp[,10]	
  AreaSwept=temp[,11]	
  kg=temp[,12]	
  Sci=temp[,13]
  vessel=temp[,14]
  
  
  #Run 1: IYSsockeyeCovariates23.csv
  #Run 1: IYSsockeyeCovariates18.csv logkappa 1 and 2 failed
  
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
  Temp_20m	=temp2[,10]	
  Temp_20m_SE=temp2[,11]	
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

  pdf(file="FIGURE 7.pdf", height=11, width=8.5, paper='special', family ='serif')   

 p1<-ggplot(temp2, aes(x=Year, y=J_Sockeye)) +
   theme_bw()+ 
   geom_line()+
   geom_pointrange(aes(ymin=J_Sockeye-J_Sockeye_SE, ymax=J_Sockeye+J_Sockeye_SE))+
   labs(title="", y="J. sockeye salmon (kg)", x="", cex=3)+
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
  ?grid.arrange

   pdf(file="FIGURE 8.pdf", height=11, width=8.5, paper='special', family ='serif')   
   
  p4<-ggplot(temp2, aes(x=Year, y=Temp_20m)) +
    theme_bw()+ 
    geom_line()+
    geom_pointrange(aes(ymin=Temp_20m-Temp_20m_SE, ymax=Temp_20m+Temp_20m_SE))+
    labs(title="", y="Temperature (Celsius)", x="", cex=3)+
    theme(text=element_text(family="sans"),axis.title = element_text(size = 12,colour = "black"),
          axis.text = element_text(size = 12,colour = "black"),
          panel.background = element_rect(colour = "black", linewidth=1.5))+
    geom_point(size=4)+
    geom_segment(aes(x=2001.5,y=7.75,xend=2005.5,yend=7.75), color="red",linewidth=2)+
    geom_segment(aes(x=2005.5,y=7.75,xend=2012.5,yend=7.75), color="blue",linewidth=2)+
    geom_segment(aes(x=2012.5,y=7.75,xend=2018.5,yend=7.75), color="red",linewidth=2)       

    p5<-ggplot(temp2, aes(x=Year, y=lnCalanus)) +
    theme_bw()+ 
    geom_line()+
    geom_pointrange(aes(ymin=lnCalanus-lnCalanus_SE, ymax=lnCalanus+lnCalanus_SE))+
    ylab(bquote('Calanus  (ln '(N ~m^-2+1)))+
    #ylab(bquote('Calanus (ln(#~m^-2+1))'))+
    theme(text=element_text(family="sans"),axis.title = element_text(size = 12,colour = "black"),
          axis.text = element_text(size = 12,colour = "black"),
          panel.background = element_rect(colour = "black", linewidth=1.5))+
    geom_point(size=4) +
    geom_segment(aes(x=2001.5,y=400000, xend=2005.5,yend=400000), color="red",linewidth=2)+
    geom_segment(aes(x=2005.5,y=400000,xend=2012.5,yend=400000), color="blue",linewidth=2)+
    geom_segment(aes(x=2012.5,y=400000,xend=2018.5,yend=400000), color="red",linewidth=2)   
  p6<-ggplot(temp2, aes(x=Year, y=lnPollock)) +
    theme_bw()+ 
    geom_line()+
    geom_pointrange(aes(ymin=lnPollock-lnPollock_SE, ymax=lnPollock+lnPollock_SE))+
    labs(title="", y="Age-0 Pollock (ln(#+1))", x="", cex=3)+
    theme(text=element_text(family="sans"),axis.title = element_text(size = 12,colour = "black"),
          axis.text = element_text(size = 12,colour = "black"),
          panel.background = element_rect(colour = "black", linewidth=1.5))+
    geom_point(size=4) +
    geom_segment(aes(x=2001.5,y=700000,xend=2005.5,yend=700000), color="red",linewidth=2)+
    geom_segment(aes(x=2005.5,y=700000,xend=2012.5,yend=700000), color="blue",linewidth=2)+
    geom_segment(aes(x=2012.5,y=700000,xend=2018.5,yend=700000), color="red",linewidth=2)   

    p7<-ggplot(temp2, aes(x=Year, y=lnPink)) +
    theme_bw()+ 
    geom_line()+
    geom_pointrange(aes(ymin=lnPink-lnPink_SE, ymax=lnPink+lnPink_SE))+
    labs(title="", y="Juvenile pink salmon (ln(kg+1))", x="", cex=3)+
    theme(text=element_text(family="sans"),axis.title = element_text(size = 12,colour = "black"),
          axis.text = element_text(size = 12,colour = "black"),
          panel.background = element_rect(colour = "black", linewidth=1.5))+
    geom_point(size=4)+
    geom_segment(aes(x=2001.5,y=-100000,xend=2005.5,yend=-100000), color="red",linewidth=2)+
    geom_segment(aes(x=2005.5,y=-100000,xend=2012.5,yend=-100000), color="blue",linewidth=2)+
    geom_segment(aes(x=2012.5,y=-100000,xend=2018.5,yend=-100000), color="red",linewidth=2)    
  
  grid.arrange(p4,p5,p6,p7,bottom="Year",top="",left="",nrow = 2, ncol=2)
  dev.off() 
  
  #Style 2
  
  p1 <- ggplot(SST) +
    geom_ribbon(aes(x=Year,
                    ymin=SST-SST_SE,
                    ymax=SST+SST_SE),
                fill='blue',alpha=0.2) +
    geom_line(aes(x=Year,
                  y=Mean),
              color='blue') +
    labs(x="",y="",title = "Ensemble mean with min/max ribbon") 
  print(p1)
 
  
  plot_grid(p1, p2, labels = c('A', 'B'))
  ?grid.arrange
  
  
  plot_variable(
    fit2$Report$eta1_gct[,1,],
    map_list = NULL,
    panel_labels = NULL,
    projargs = "+proj=longlat",
    map_resolution = "medium",
    file_name = "density",
    working_dir = paste0(getwd(), "/"),
    Format = "png",
    Res = 200,
    add = FALSE,
    outermargintext = c("Eastings", "Northings"),
    zlim = NULL,
    col = viridisLite::viridis,
    mar = c(0, 0, 2, 0),
    oma = c(4, 4, 0, 0),
    legend_x = c(0, 0.05),
    legend_y = c(0.05, 0.45),
    cex.legend = 1,
    mfrow = NULL,
    land_color = "grey",
    n_cells = NULL,
    xlim = NULL,
    ylim = NULL,
    country = NULL,
    contour_nlevels = 0,
    fun = mean,
    format = "sf",
    cex.points = 1,
    legend_digits = 1
  )
  
  
  
  dev.off()
  
  ####################
  # Timeseries Plots 2
  ####################
  temp3<-read.csv("TimeseriesPlotsLong.csv") 
  Year=temp3[,1]		
  Mean=temp3[,2]		
  SE=temp3[,3]
  Min=temp3[,4]
  Max=temp3[,5]
  library(ggplot2)
  library(cowplot)
  library(raster)
  library(ggplot2)
  library(xts)
  
  p1 <- ggplot(temp3, group=Index) +
    geom_ribbon(aes(x=Year,
                    ymin=Min,
                    ymax=Max),
                fill='blue',alpha=0.2) +
    geom_line(aes(x=Year,
                  y=Mean),
              color='blue') +
    labs(x="",y="",title = "Ensemble mean with min/max ribbon") 
  print(p1)
  p1 <- ggplot(mtcars, aes(disp, mpg)) + 
    geom_point()
  p2 <- ggplot(mtcars, aes(qsec, mpg)) +
    geom_point()
  
  plot_grid(p1, p2, labels = c('A', 'B'))
  
  
  #================== CORRELATION MATRIX ====================
  #  FIGURE 


  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(lattice)
  library(Hmisc)
  library("PerformanceAnalytics")
  library(corrplot)  #Package corrplot version 0.84

  library(terra)
  library(ggplot2)
  library(cowplot)
  
  
  #if(!require(devtools)) install.packages("devtools")
  #devtools::install_github("kassambara/ggpubr")
  #install.packages("ggpubr")
  
  # Correlation matrix
  temp1<-data.frame(J_Sockeye, Northing, Easting, EAO, Temp_20m, lnCalanus,	lnPollock,	lnPink)

  
  #Correlation panel
  panel.cor<-function(x,y){
    usr<-par("usr"); on.exit(par(usr))
    par(usr=c(0,1,0,1))
    r<-round(cor(x,y), digits=2)
    txt<-paste0("R=", r)
    cex.cor<-0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex=cex.cor*r)
  }
  pairs(temp1, lower.panel=panel.cor)
  
  library(psych)
    library(ggpubr)
pairs.panels(temp1, method="pearson", hist.col="#00AFBB", density=TRUE)
  
  res <- cor(temp1)
  round(res, 2)
   #temp2<-data.frame(Year,		SST,		CPUE_RAW, CPUE_VAST,	WT_RAW,		WT_VAST,		ED_RAW,	ED_VAST,		GRP_RAW,GRP_VAST)
  temp2<-data.frame(Year, J_Sockeye, Northing, Easting, EAO, Temp_20m_SE, lnCalanus,	lnPollock,	lnPink)
  
  p1 <- ggplot(temp2, aes(SST, CPUE_VAST))+ 
    geom_text(label = (Year)) +
    geom_smooth(method = "lm", se = FALSE)  
  
  p2 <- ggplot(temp2, aes(SST, Pollock))+
    geom_smooth(method = "lm", se = FALSE)+ 
    geom_text(label = (Year)) 
  
  p3 <- ggplot(temp2, aes(SST, ED_VAST)) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(label = (Year))
  
  p4 <- ggplot(temp2, aes(SST, ED_RAW)) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text(label = (Year))
  
  ggarrange(p1, p2,p3,p4 + rremove("x.text"), 
            ncol = 2, nrow = 2)
  
  # ++++++++++++++++++++++++++++
  # flattenCorrMatrix
  # ++++++++++++++++++++++++++++
  # cormat : matrix of the correlation coefficients
  # pmat : matrix of the correlation p-values
 temp1<-data.frame(J_Sockeye, Northing, Easting, EAO, Temp_20m, lnCalanus,	lnPollock,	lnPink)
  
  flattenCorrMatrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  =(cormat)[ut],
      p = pmat[ut]
    )
  }
  
  res1 <- rcorr(as.matrix(temp1))
  res1
  flattenCorrMatrix(res1$r, res1$P)
  
  my_data <- temp1
  chart.Correlation(my_data, histogram=TRUE, pch=19)
  res <- cor(my_data)
  round(res, 2)
  
  # Insignificant correlation are crossed
  col<-scale_fill_viridis(option = "D")
  ?scale_fill_viridis
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  col <- colorRampPalette(c("#24868EFF", "#20A486FF", "#47C16EFF", "#8FD744FF", "#E3E418FF","#FDE725FF"))
  #col <- colorRampPalette(c("#22A884FF", "#7AD151FF","#FDE725FF"))
  
  corrplot(res1$r, type="upper", diag=FALSE, order='original',
           p.mat = res1$P, addCoef.col = "black", method = "circle", col = COL2('RdYlBu', 10),sig.level = 0.05,
           tl.col = "black", tl.srt = 45) #Text label color and rotation
 ?corrplot
  plot(J_Sockeye,SST)
  .05/7
  ?corrplot
  
  #########
  #  FIGURE MATRIX
  #########
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot(res1$r, method = "color", col = col(200),  
           type = "upper", order = "hclust", 
           addCoef.col = "black", # Add coefficient of correlation
           tl.col = "darkblue", tl.srt = 45, #Text label color and rotation
           # Combine with significance level
           p.mat = p_mat, sig.level = 0.01,  
           # hide correlation coefficient on the principal diagonal
           diag = FALSE 
  )
  
  #NONlinear
  X3_formula=bs(SST, degree=2, intercept=FALSE)
  
  ?bs
  bs(SST,degree=2)
  summary(X3_formula<-lm(J_Sockeye~splines::bs(Temp_20m, df=3), data=temp1))
  plot(temp1)
  ht<-seq(57,73,length.out=17)
  lines(ht,predict(X3_formula, data.frame(J_Sockeye=Temp_20m)))
  
  ?panel.smooth
  
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use = "complete.obs"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex =  cex.cor * (1 + r) / 2)
  }
  
  panel.hist <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
  }
  
  panel.lm <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                        cex = 1, col.smooth = "black", ...) {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    #abline(stats::lm(y ~ x),  col = col.smooth, ...)
    abline(stats::loess(y ~ x),  col = col.smooth, ...)
    
  } 
  panel.smooth(x, y, col = par("col"), bg = NA, pch = par("pch"),
               cex = 1, col.smooth = 2, span = .8, iter = 3
               )
  
  ?loess
  splines::bs(J_Sockeye~SST, degree=3, intercept=FALSE)
  lo <- loess(y~x)
  plot(x,y)
  lines(predict(lo), col='red', lwd=2)
  
  ?pairs
  
  temp4<-data.frame(J_Sockeye, Temp_20m, Calanus, lnPollock, lnPink)
  temp5<-data.frame(J_Sockeye, Temp_20m, Calanus, lnPollock, lnPink, Northing, Easting, EAO)
 
   pdf(file="FIGURE matrix.pdf", height=8.5, width=11, paper='special')   
   
  pairs( 
    temp5,
    lower.panel = panel.cor,
    diag.panel  = panel.hist,
    upper.panel = panel.smooth,
    pch = 20,
    size=1,
     cex.labels=2,font.labels=1.5,  cex = 2, cex.axis = 1.5
  )

  dev.off()
   ?pairs
 ggplot(temp5, aes(Temp_20m,J_Sockeye)) + geom_point() +
    geom_smooth(method = "gam", formula = y ~ poly(x, 2))
  
  geom_smooth(method = "gam", formula = y ~ poly(x, 2))

  
  summary(temp5,formula = J_Sockeye ~ poly(Temp_20m, 2))
  summary(lm(J_Sockeye~Temp_20m) )
  ?lm
  data(kyphosis) gam(Kyphosis~s(Age,4)+Number,family=binomial,data=kyphosis, trace=TRUE) data(airquality) gam(Ozone^(1/3)~lo(Solar.R)+lo(Wind,Temp),data=airquality,na=na.gam.replace) gam(Kyphosis~poly(Age,2)+s(Start),data=kyphosis,family=binomial,subset=Number>2) data(gam.data) 
  Gam.object<-gam(J_Sockeye~s(Temp_20m,2),data=temp5) 
  summary(Gam.object) 
  plot(Gam.object,se=TRUE)
  library(plotly)
  library(GGally)
  install.packages("gam")
  library(gam)
  data(flea)
  
  p <- ggpairs(temp5, columns = 1:5,
               lower = list(continuous = "cor", combo = "box_no_facet", discrete = "count", na = "na"),
               upper = list(continuous = "points", combo = "facethist", discrete = "facetbar", na =
                              "na"),
               diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag"),) 
  
  ggplotly(p)
  ?ggpairs
  
  # load libraries ggplot2 and ggally 
  library("ggplot2", lib.loc="C:/Program Files/R/R-4.2.2/library")
  library("GGally", lib.loc="C:/Program Files/R/R-4.2.2/library")
  install.packages("GGally")
  
  install.packages("ggplot2")
  
  install.packages(ggplot2 )
  ?install.packages
  library(GGally) 
  library(ggplot2)

  # create pairs plot 
  ggpairs( temp4 )
  
  
?pairs
  
  plot(Temp_20m,J_Sockeye)
  
  ######################
  # COW PLOT
  #######################
  ggtitle(bquote('Number VS'~Number^2))
  
  
  pdf(file="FIGURE matrix2.pdf", height=11, width=8.5, paper='special')   
 p1 <- ggplot(temp5, aes(Temp_20m, J_Sockeye)) + 
    geom_point()+
   xlab(NULL)+
   ylab("J. sockeye salmon (kg)")+
   ggtitle(bquote(R^2~'=0.36, p=0.03'))+
    geom_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE)+
   theme_cowplot(12)
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
    theme_cowplot(12)
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
 
 
  
  #############
  # Fields krig
  #############
  
  # Make a raster layer
  library(raster)
  
  # projection
  utm.prj = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"   
  
  # create a SpatialPointsDataFrame
  coordinates(SST) = ~Lon+Lat   
  proj4string(SST) <-CRS(utm.prj)

  # create an empty raster object to the extent of the points
  rast <- raster(ext=extent(divetemps),crs = CRS(utm.prj), resolution = 500) # 500 m x 500 m
  rast
  
  # rasterize your irregular points 
  rasOut<-raster::rasterize(divetemps, rast, divetemps$depthbin1, fun = mean) # we use a mean function here to regularly grid the irregular input points
  plot(rasOut)
  
  library(fields)
  # Function to Krig 
  krigR <- function(rast){ 
    xy <- data.frame(raster::xyFromCell(rast, 1:ncell(rast)))
    v <- getValues(rast)
    krg <- fields::Krig(xy, v) 
    ras.int <- raster::interpolate(rast, krg)
    proj4string(ras.int) <- proj4string(rast)
    return(ras.int)
  }
  
  surface = krigR(rasOut) 
  plot(surface)
  
  
  ############
  #Kriging####
  ############
  #https://stackoverflow.com/questions/50594542/how-to-extract-specific-values-with-point-coordinates-from-kriging-interpolation
  temp<-read.csv("IYSsockeyeCovariates29.csv") 
  
  PG<-read.csv("IYSsockeyeCovariates29.csv", header=T, stringsAsFactors=FALSE)
  install.packages("geoR")
  library("geoR")
  x<-(PG$Lat) #Latitude
  y<-(PG$Lon) #Longitude
  
  ?expand.grid
  ?seq
  #Grid
  summary(temp$Lat)
  summary(temp$Lon)
  loci<-expand.grid(x=seq(54.9, 60, length=10), y=seq(-173, -159, length=10))
  names(loci)<-c("x", "y")
  mix<-cbind(rep(100), loci$x, loci$y, loci$x*loci$y)
  temp$SST
  #Model
  pH1.mod<-lm(SST~y*x, data=PG, x=T)
  pH1.kg<-cbind(pH1.mod$x[,3], pH1.mod$x[,2], pH1.mod$residuals)
  #Transform to geographic data
  pH1.geo<-as.geodata(pH1.kg,)
  #?as.geodata
  #Variogram
  pH1.vario<-variog(pH1.geo, max.dist=35)
  pH1.vario.mod<-eyefit(pH1.vario)
  #Cross validation
  pH1.valcruz<-xvalid(pH1.geo, model=pH1.vario.mod)
  #Kriging
  pH1.krig<-krige.conv(pH1.geo, loc=loci, krige=krige.control(obj.model=pH1.vario.mod[[1]]))
  #Predictive model
  pH1a.yhat<-mix %*% pH1.mod$coefficients + pH1.krig$predict
  #Exchange Kriging prediction values
  pH1.krig$predict<-pH1.yhat
  #Image
  image(pH1.krig2)
  contour(pH1.krig2, add=TRUE)
  
  #Tree matrix####:
  
  CoA<-read.csv("CoAr.csv", header=T)
  #Data
  X<-temp$Lon
  Y<-temp$Lat
  xa<-(CoA$X)
  ya<-(CoA$Y)
  points(xa,ya, col=4)
  
  TreeDF<-(cbind.data.frame(xa, ya, CoA$Species, CoA$DBH, CoA$Height, stringsAsFactors = TRUE))
  m<-(cbind(xa, ya, 1:305)) 
  as.matrix(m)
  
  
  
  #EXtract data point
  library(raster)
  
  r <- SpatialPointsDataFrame(loci, data.frame(predict = pH1.krig$predict))
  gridded(r) <- T
  r <- as(r,'RasterLayer')
  
  COA2<-cbind(CoA, pH1val=pH1.arb)
  
  pts <- SpatialPointsDataFrame(CoA[,c('X','Y')],CoA)
  
  pH1.arb <-extract(r, pts)
  
  ################
  #gstat
  ################
  library(sp); library(raster); library(gstat)
  v.utm # SpatialPointsDataFrame with >30,000 points
  
  # Remove points with identical positons
  zd = zerodist(v.utm)
  nzd = v.utm[-zd[,1],] # Layer with no identical positions
  
  # Make a raster layer covering point layer
  resolution=1e4
  e = extent(as.matrix(v.utm@coords))+resolution
  r = raster(e,resolution=resolution) 
  proj4string(r) = proj4string(v.utm)
  
  # r is a 181x157 raster
  
  # Fit variogram
  fv = fit.variogram(variogram(AVGDEPTH~1, nzd),model=vgm(6000,"Exp",1,5e5,1))
  
  # Krige on random sample of 500 points - works fine
  size=500
  ss=nzd[sample.int(nrow(nzd),size),]
  depth.krig = krige(AVGDEPTH~1,ss,as(r,"SpatialPixelsDataFrame"),
                     model=depth.fit)
  
  # Krige on random sample of 5000 points - never seems to end
  size=5000
  ss=nzd[sample.int(nrow(nzd),size),]
  depth.krig = krige(AVGDEPTH~1,ss,as(r,"SpatialPixelsDataFrame"),
                     model=depth.fit)
  
  
  ################
  # gstat chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/gstat/vignettes/gstat.pdf
  ###################
  library(sp) 
  temp<-read.csv("IYSsockeyeCovariates29.csv") 
  
  
  data(temp) 
    class(temp)
    coordinates(temp) = ~Lat+Lon
    coordinates(temp)=~Lat_rnd+Lon_rnd
    class(temp)
    summary(temp)
    coordinates(temp)[1:5,]
    bubble(temp, "SST", col=c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC"), main = "SST (C)")
    
    ?col
    gridded(temp) = TRUE
    class(temp)
    image(meuse.grid["dist"]) 
     title("SST") 
     library(gstat) 
     SST.idw = idw(SST~1, temp, temp) 
     #[inverse distance weighted interpolation] 
     class(SST.idw) 
     "SpatialPixelsDataFrame" 
      spplot(SST.idw["var1.pred"], main = "SST weighted interpolations")
      plot(log(SST)~sqrt(dist), meuse) 
      abline(lm(log(zinc)~sqrt(dist), meuse))
    
      lzn.kriged = krige(log(SST)~1, temp, temp, model = lzn.fit)
      spplot(lzn.kriged["var1.pred"])
      
      #Data on a regular grid
      lzn.condsim = krige(log(SST)~1, temp, meuse, model = lzn.fit,  nmax=30,nsim=4) 
                            spplot(lzn.condsim, main = "four conditional simulations") 10
      
      
    data(meuse.grid)
    summary(meuse.grid)
    

    library(sp)
    data(temp)
    coordinates(temp) = ~Lat_rnd+Lon_rnd
    data(temp.grid)
    gridded(meuse.grid) = ~x+y
    m <- vgm(.59, "Sph", 874, .04)
    # ordinary kriging:
    x <- krige(log(SST)~1, temp, temp, model = m)
    spplot(x["var1.pred"], main = "ordinary kriging predictions")
    spplot(x["var1.var"], main = "ordinary kriging variance")
    
    temp<-read.csv("IYSsockeyeCovariates29.csv") 
    library(gstat)
    library(sp)
    data(temp)
    coordinates(temp) = ~Lat+Lon
    data(meuse.grid)
    gridded(temp) = ~Lat_rnd+Lon_rnd
    m <- vgm(.59, "Sph", 874, .04)
    # ordinary kriging:
    x <- krige(log(SST)~1, temp, temp, model = m, grouping=Year)
    spplot(x["var1.pred"], main = "ordinary kriging predictions")
    spplot(x["var1.var"], main = "ordinary kriging variance")
    # simple kriging:
    x <- krige(log(zinc)~1, meuse, meuse.grid, model = m, beta = 5.9)
    # residual variogram:
    m <- vgm(.4, "Sph", 954, .06)
    # universal block kriging:
    x <- krige(log(zinc)~x+y, meuse, meuse.grid, model = m, block = c(40,40))
    spplot(x["var1.pred"], main = "universal kriging predictions")
    # krige0, using user-defined covariance function and multiple responses in y:
    # exponential variogram with range 500, defined as covariance function:
    v = function(x, y = x) { exp(-spDists(coordinates(x),coordinates(y))/500) }
    # krige two variables in a single pass (using 1 covariance model):
    y = cbind(meuse$zinc,meuse$copper,meuse$lead,meuse$cadmium)
    x <- krige0(zinc~1, meuse, meuse.grid, v, y = y)
    meuse.grid$zinc = x[,1]
    spplot(meuse.grid["zinc"], main = "zinc")
    meuse.grid$copper = x[,2]
    spplot(meuse.grid["copper"], main = "copper")
    # the following has NOTHING to do with kriging, but --
    # return the median of the nearest 11 observations:
    x = krige(zinc~1, meuse, meuse.grid, set = list(method = "med"), nmax = 11)
    # get 25%- and 75%-percentiles of nearest 11 obs, as prediction and variance:
    x = krige(zinc~1, meuse, meuse.grid, nmax = 11,
              set = list(method = "med", quantile = 0.25))
    # get diversity (# of different values) and mode from 11 nearest observations:
    x = krige(zinc~1, meuse, meuse.grid, nmax = 11, set = list(method = "div"))
    
    
    ###############
    https://rdrr.io/cran/marmap/man/griddify.html
    
    ###########
    # Sockeye model with VAST covariate densities outputs 
    ###########
    
    BioCov_2002<-as.vector(fit$Report$D_gct[,,1])
    BioCov_2003<-as.vector(fit$Report$D_gct[,,2])
    BioCov_2004<-as.vector(fit$Report$D_gct[,,3])
    BioCov_2005<-as.vector(fit$Report$D_gct[,,4])
    BioCov_2006<-as.vector(fit$Report$D_gct[,,5])
    BioCov_2007<-as.vector(fit$Report$D_gct[,,6])
    BioCov_2008<-as.vector(fit$Report$D_gct[,,7])
    BioCov_2009<-as.vector(fit$Report$D_gct[,,8])
    BioCov_2010<-as.vector(fit$Report$D_gct[,,9])
    BioCov_2011<-as.vector(fit$Report$D_gct[,,10])
    BioCov_2012<-as.vector(fit$Report$D_gct[,,11])
    BioCov_2013<-as.vector(fit$Report$D_gct[,,12])
    BioCov_2014<-as.vector(fit$Report$D_gct[,,13])
    BioCov_2015<-as.vector(fit$Report$D_gct[,,14])
    BioCov_2016<-as.vector(fit$Report$D_gct[,,15])
    BioCov_2017<-as.vector(fit$Report$D_gct[,,16])
    BioCov_2018<-as.vector(fit$Report$D_gct[,,17])
    
    library(tidyr)
    BioCov_fitted<-data.frame(BioCov_2002,BioCov_2003,BioCov_2004,BioCov_2005,BioCov_2006,BioCov_2007,BioCov_2008,BioCov_2009,BioCov_2010,BioCov_2011,BioCov_2012,BioCov_2013,BioCov_2014,BioCov_2015,BioCov_2016,BioCov_2017,BioCov_2018)
    BioCov_fit2<-data.frame(unlist(BioCov_fitted))
    Year<-rep(c(2002:2018),each=2000)
    head(fit$Report)
    print(fit$Report)
    Lon<-rep(fit$extrapolation_list$Data[,"Lon"],each=17)
    Lat<-rep(fit$extrapolation_list$Data[,"Lat"],each=17)
    
    #1. run VAST on pollock data
    #2.run density estimates as covariates in the juvenile sockeye VAST model
    
    #lnBioCov<-scale(log(BioCov_fit2+1))
    lnBioCov<-scale(BioCov_fit2)
    head(BioCov_fit2)
    summary(BioCov_fit2)
    plot(lnBioCov)
    dev.off()
    X1<-data.frame(Year, Lat, Lon,lnBioCov) 
    X1_formula = ~ bs( lnBioCov, degree=2, intercept=FALSE)
    X2_formula = ~ bs( lnBioCov, degree=2, intercept=FALSE)
    X1_formula=~lnBioCov
    X2_formula=~lnBioCov
    
    settings2 = make_settings( Version="VAST_v14_0_1",
                               n_x=500,
                               Region='User', 
                               purpose="index2",
                               bias.correct=TRUE,
                               max_cells = 2000,
                               ObsModel=c(2,1), 
                               fine_scale=TRUE,
                               FieldConfig = matrix( c("IID","IID","IID","IID","IID","IID"), byrow=TRUE, ncol=2 ),
                               treat_nonencounter_as_zero=TRUE,
                               knot_method='grid',
                               use_anisotropy = TRUE)
    user_region <- readRDS('user_region.rds')
    
    fit2 = fit_model(
      "settings"=settings2,
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
      getReportCovariance=TRUE,
      "catchability_data" = data.frame(date=DOY), #Day of year
      Q2_formula = ~ date)
    
    
    check_fit(fit2$parameter_estimates, check_gradients = FALSE, quiet = FALSE)
    #Check extrapolation grid total area....
    colSums(fit2$extrapolation_list$a_el)
    
    #SAVE DATA =================================================================================
    x=summary(fit2$parameter_estimates$SD, select='report')
    print(x)
    write.csv(x,"Sockeye21SST_nl.csv", row.names = TRUE)
    
    # PLOT RESULTs =============================================================================
    plot_results(fit2, n_cells=2000  )
    plot_results( fit2, settings=settings, plot_set=c(3,11,12,14,15), n_cells=2000)
    
    
    