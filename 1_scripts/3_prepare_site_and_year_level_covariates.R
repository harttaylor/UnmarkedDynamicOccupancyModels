# Prepare all your covariates to run in unmarked 
# need to structure data as site-level, yearly-level, and observation-level covariates (detection covariates)
# data scaling, standardization, or any other necessary manipulation
# Basically this is how each type of covariate should be prepared: 
# 1. Site-level covariate: constant for each site and do not change over years or visits (e.g., habitat type, elevation, lat). A dataframe where each row represents a site (219 rows in my case)
# Ensure the order of sites in this data frame matches the order in the occupancy and observation-level data.
# 2. yearly-level covariates: change from year-to-year but are constant across all visits within a year for a site (i.e. a dataframe with 219 sites x 25 years = 5475 rows).
# Each row corresponds to a specific site-year combination. Columns are the yearly covariates (e.g., annual rainfall, temperature).
# 3. Observation-level covariates (detection covariates): These are covariates that can change with each visit (e.g., time of day, weather conditions).
# need a dataframe with 219 sites × 25 years × 4 visits = 21900 rows. Each row corresponds to a specific site-year-visit combination.


# site-level covariates preparation 
# Load Data
load("0_data/raw/Calling Lake/CallingLakeVegData.RData")  # stand variables. These are data from the Beaudoin satellite layers
names(CL)
nlevels(as.factor(CL$SS))#292 sites, not all of which will be used in models

points.sitedata <- read.csv("0_data/raw/Calling Lake/CLpoints_covariates.csv", header = TRUE)

#Formatting/cleaning data because it is messy 
points.sitedata.u<-unique(points.sitedata[,c("SS","PKEY","YYYY",
                                             "latitude",
                                             "longitude",
                                             "EASTING",
                                             "NORTHING",
                                             "Patch.Size..ha.",
                                             "Treatment")])
points.sitedata.u[points.sitedata.u$Treatment=="Riparian Strip",]
#some of the "extra" points within buffer sites lack latitude and longitude data

points.sitedata.u<-points.sitedata.u[!is.na(points.sitedata.u$latitude),]
#'data.frame':	219 sites,  9 variables; stations should be in same order as other files

# Merge file with stand variables (age, % conifer) with file with correct years 
StandVariables <- merge(CL, points.sitedata.u, by = c("PKEY","SS"))
str(StandVariables)
nlevels(as.factor(StandVariables$SS))#219

#get unique values in the columns "SS", "YYYY", and veg variables
standvariables.u<-StandVariables
standvariables.u$PKEY<-NULL#drop PKEY
standvariables.u <- unique(standvariables.u)
nrow(standvariables.u)#2814

#filter for 1993 only-affecting initial occupancy 
#there is only CL data for 2001 and more recent, so use year 2001
#for the vegetation in 1993
standvariables.2001 <- filter(standvariables.u, YYYY == "2001")
write.csv(standvariables.2001, "0_data/raw/Calling Lake/StandVariables2001.csv")

#read back in. This file gives us vegetation data that could influence
#probability of initial occupancy
standvariables.init <- read.csv("0_data/raw/Calling Lake/StandVariables2001.csv")
standvariables.init$WhiteSpruce <- standvariables.init$Species_Pice_Gla_v1
standvariables.init$whitespruce.s<-standvariables.init$WhiteSpruce/max(standvariables.init$WhiteSpruce)
standvariables.init$whitespruce.s.2<-(standvariables.init$whitespruce.s)^2
range(standvariables.init$whitespruce.s.2)#0 1

standvariables.init$Age <- standvariables.init$Structure_Stand_Age_v1
#these ages are from 2001, subtract 8 to get ages in 1993
standvariables.init$Age <- standvariables.init$Structure_Stand_Age_v1 -8
range(standvariables.init$Age)#20.25667 94.50667
standvariables.init$age.s<-standvariables.init$Age/max(standvariables.init$Age)
standvariables.init$age.s.2<- (standvariables.init$age.s-0.5)^2


#new site-year data for Calling Lake
CL.93.14<-read.csv("0_data/raw/Calling Lake/new footprint/all-footprints-areadist.age-allsurveys-justCallingLake-1993-2014.csv", header=TRUE)
names(CL.93.14)
CL.15.18<-read.csv("0_data/raw/Calling Lake/new footprint/all-footprints-areadist.age-CL.2015-2018.csv", header=TRUE)
names(CL.15.18)

CL.93.18<-bind_rows(CL.93.14,
                    CL.15.18)

CL.93.18.SS<-CL.93.18#33801
CL.93.18.SS$PKEY<-NULL#33801
CL.93.18.SS.u<-unique(CL.93.18.SS)#7965
#We'll use this below

#Footprint data from 1993. We'll use this with vegetation data from 1993/2002 to
#model initial occupancy
CL.93only.SS.u<-CL.93.18.SS.u[CL.93.18.SS.u$YEAR==1993,]#328
names(CL.93only.SS.u)
names(which(colSums(is.na(CL.93only.SS.u))>0))

#add new footprint data for 1993 to vegetation data for 1993
veghf.1993<-merge(standvariables.init, CL.93only.SS.u, by=c("SS"))



# Standardize environmental covariates and scale from 0-1
sitecovs<-points.sitedata.u[!is.na(points.sitedata.u$latitude),]
#'data.frame':	219 obs. of  11 variables; stations should be in same order as other files
#TREATMENT variable
sitecovs$Treatment
#scale/standardize your variables before you add site cover data to 
#unmarked-data-frame
#harvest variables

sitecovs$HARVEST_NEAR_DIST<-ifelse(sitecovs$HARVEST_NEAR_DIST==0,1,sitecovs$HARVEST_NEAR_DIST)
sitecovs$harvest.dist.s<-sitecovs$HARVEST_NEAR_DIST/max(sitecovs$HARVEST_NEAR_DIST)
range(sitecovs$harvest.dist.s)
sitecovs$harvest.dist.cen<-(sitecovs$harvest.dist.s-0.5)
sitecovs$harvest.dist2<-(sitecovs$harvest.dist.s-0.5)^2
sitecovs$harvest.dist05<-sqrt(sitecovs$harvest.dist.s)
sitecovs$harvest.dist.i<-1/(sitecovs$harvest.dist.s)
#road variables 
sitecovs$ROADS_NEAR_DIST<-ifelse(sitecovs$ROADS_NEAR_DIST==0,1,sitecovs$ROADS_NEAR_DIST)
sitecovs$road.dist.s<-sitecovs$ROADS_NEAR_DIST/max(sitecovs$ROADS_NEAR_DIST)
range(sitecovs$road.dist.s)
sitecovs$road.dist.cen<-(sitecovs$road.dist.s-0.5)
sitecovs$road.dist2<-(sitecovs$road.dist.s-0.5)^2
sitecovs$road.dist05<-sqrt(sitecovs$road.dist.s)
sitecovs$road.dist.i<-1/(sitecovs$road.dist.s)
# pipeline variables 
sitecovs$PIPELINE_NEAR_DIST<-ifelse(sitecovs$PIPELINE_NEAR_DIST==0,1,sitecovs$PIPELINE_NEAR_DIST)
sitecovs$pipeline.dist.s<-sitecovs$PIPELINE_NEAR_DIST/max(sitecovs$PIPELINE_NEAR_DIST)
range(sitecovs$pipeline.dist.s)
sitecovs$pipeline.dist.cen<-(sitecovs$pipeline.dist.s-0.5)
sitecovs$pipeline.dist2<-(sitecovs$pipeline.dist.s-0.5)^2
sitecovs$pipeline.dist05<-sqrt(sitecovs$pipeline.dist.s)
sitecovs$pipeline.dist.i<-1/(sitecovs$pipeline.dist.s)
#seismicline variables 
sitecovs$SEISMIC_NEAR_DIST <-ifelse(sitecovs$SEISMIC_NEAR_DIST==0,1,sitecovs$SEISMIC_NEAR_DIST)
sitecovs$seismic.dist.s<-sitecovs$SEISMIC_NEAR_DIST/max(sitecovs$SEISMIC_NEAR_DIST)
range(sitecovs$seismic.dist.s)
sitecovs$seismic.dist.cen<-(sitecovs$seismic.dist.s-0.5)
sitecovs$seismic.dist2<-(sitecovs$seismic.dist.s-0.5)^2
sitecovs$seismic.dist05<-sqrt(sitecovs$seismic.dist.s)
sitecovs$seismic.dist.i<-1/(sitecovs$seismic.dist.s)
site.cov<-read.csv(paste0("0_data/processed/species weights/",i,"/",i,"site.cov.check.csv"), header=TRUE)



# This is species weight matrices to account for spatial autocorrelation in my model because my sites are so tightly clustered together 
names = c("BTNW")

for (i in names){
wt4600<-read.csv(paste0("0_data/processed/species weights/",i,"/",i,"_weights_4600.csv"),header=TRUE)
weights4600<-as.data.frame(wt4600[,2])
#these files contain "weights" for the observed presences at each station.
#the weights account for spatial autocorrelation in detections from stations at
#the same sites.

#weights are calculated as the proportion of stations within 4600 m
#of each station in either year 1 each year (not including the central station) where a given
#species has been detected at least once in that year.
#We expect that the weights should be (+) correlated with probability of initial
#occupancy (year 1) at each station, (+) correlated with probability of colonization
#of stations unoccupied in previous year, and (-) correlated with probability of extinction
#(failure to return to station occupied in previous year). These predictions
#were supported by results from occupancy analysis.
site.cov<-c(sitecovs,weights4600)
site.cov$weights4600<-site.cov$`wt4600[, 2]`
site.cov$`wt4600[, 2]`<-NULL
write.csv(site.cov, file=paste0("0_data/processed/species weights/",i,"/",i,"site.cov.check.csv"))
#data frame check. Right now, site.cov is a list
#these files contain "weights" for the observed presences at each station.
#the weights account for spatial autocorrelation in detections from stations at
#the same sites.

#weights are calculated as the proportion of stations within 4600 m
#of each station in either year 1 each year (not including the central station) where a given
#species has been detected at least once in that year.
#We expect that the weights should be (+) correlated with probability of initial
#occupancy (year 1) at each station, (+) correlated with probability of colonization
#of stations unoccupied in previous year, and (-) correlated with probability of extinction
#(failure to return to station occupied in previous year). These predictions
#were supported by results from occupancy analysis.

site.cov<-c(veghf.1993,weights4600)
site.cov$weights4600<-site.cov$`wt4600[, 2]`
site.cov$`wt4600[, 2]`<-NULL
write.csv(site.cov, file="0_data/processed/Calling Lake/species weights/BTNW/BTNWsite.cov.check.csv")
#data frame check. Right now, site.cov is a list




site.cov<-read.csv("0_data/processed/Calling Lake/species weights/BTNW/BTNWsite.cov.check.csv", header=TRUE)

wts4600<-read.csv("0_data/processed/Calling Lake/species weights/BTNW/BTNW_weights_4600.stacked.csv", header=TRUE) 
autolog.4600<-wts4600$values

#Get yearly covariates 
#old yearly covariates: gives us the list of sites we want for filtering the new data
yearlycov<-read.csv("0_data/raw/Calling Lake/yearlycovariatesfromanotherproject.csv", header=TRUE)
}

for (i in names){
  ifelse(!dir.exists(file.path(paste0("2_outputs/",i,"/"))), dir.create(file.path(paste0("2_outputs/",i,"/"))), FALSE)
  site.cov<-read.csv(paste0("0_data/processed/species weights/",i,"/",i,"site.cov.check.csv"), header=TRUE)
}#read site.cov back in to get it as a data frame.


#merge stand covariates with site covariates 
site.cov <- merge(standvariables, site.cov, by = "SS")
site.cov <- select(site.cov,- X)
site.cov <- select(site.cov, - X.x)
write.csv(site.cov, "0_data/processed/site.cov.csv") # here is the final product 


wts4600<-read.csv(paste0("0_data/processed/species weights/",i,"/",i,"_weights_4600.stacked.csv"), header=TRUE)
autolog.4600<-wts4600$values





# Get yearly covariates 
yearlycov<-read.csv("0_data/raw/yearlycovariatesfromanotherproject.csv", header=TRUE)
str(yearlycov)
yearlycov.u<-unique(yearlycov[,c("SS","YYYY",
                                 "Logged100","Logged200","Logged400",
                                 "Logged600","CutblockEdgeDistY",
                                 "PIPELINE_DIST","ROAD_DIST",
                                 "CUTLINE_DIST","TRAIL_DIST","WATER_DIST",
                                 #"L100","L400","RoadEdge","SeismicEdge","TrailEdge",
                                 #"PipelineEdge","WaterEdge",
                                 "YearOfNearestCutblock","AgeOfNearestCutblock")])
nrow(yearlycov.u)#5256, still too many observations
nlevels(as.factor(yearlycov.u$SS))#219, right number of sites
range(yearlycov.u$YYYY)#1993 2018
#now get rid of 2004 (year 11) because that year was
#been excluded from the visits used in these models

yearlycov.u<-yearlycov.u[!yearlycov.u$YYYY==2004,]
nrow(yearlycov.u)#there should be 5475, ordered by station then by year

logged100.s<-yearlycov.u$Logged100/max(yearlycov.u$Logged100)
logged200.s<-yearlycov.u$Logged100/max(yearlycov.u$Logged200)
logged400.s<-yearlycov.u$Logged100/max(yearlycov.u$Logged400)
logged600.s<-yearlycov.u$Logged100/max(yearlycov.u$Logged600)
harvest.dist.y.s<-yearlycov.u$CutblockEdgeDistY/max(yearlycov.u$CutblockEdgeDistY)
harvest.age.y.s<-yearlycov.u$AgeOfNearestCutblock/max(yearlycov.u$AgeOfNearestCutblock)
y.seismic.dist<-yearlycov.u$CUTLINE_DIST/max(yearlycov.u$CUTLINE_DIST)
y.pipeline.dist<-yearlycov.u$PIPELINE_DIST/max(yearlycov.u$PIPELINE_DIST)
y.road.dist<-yearlycov.u$ROAD_DIST/max(yearlycov.u$ROAD_DIST)


yearly.cov<-data.frame(logged100.s,logged200.s,
                       logged400.s,logged600.s,harvest.dist.y.s,
                       harvest.age.y.s, y.seismic.dist, y.road.dist, y.pipeline.dist)
write.csv(yearly.cov, "0_data/processed/Calling Lake/yearlycovariates.csv")
