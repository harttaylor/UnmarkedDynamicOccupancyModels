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

library(dplyr)
# site-level covariates preparation 
# Load Data
load("0_data/raw/CallingLakeVegData.RData")  # stand variables. These are data from the Beaudoin satellite layers
names(CL)
nlevels(as.factor(CL$SS))#292 sites, not all of which will be used in models

points.sitedata <- read.csv("0_data/raw/CLpoints_covariates.csv", header = TRUE)

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
write.csv(standvariables.2001, "0_data/raw/StandVariables2001.csv")

#read back in. This file gives us vegetation data that could influence
#probability of initial occupancy
standvariables.init <- read.csv("0_data/raw/StandVariables2001.csv")
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
CL.93.14<-read.csv("0_data/raw/all-footprints-areadist.age-allsurveys-justCallingLake-1993-2014.csv", header=TRUE)
names(CL.93.14)
CL.15.18<-read.csv("0_data/raw/all-footprints-areadist.age-CL.2015-2018.csv", header=TRUE)
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



#add new footprint data for 1993 to vegetation data for 1993
veghf.1993<-merge(standvariables.init, CL.93only.SS.u, by=c("SS"))

#standardize environmental covariates and scale from 0-1
#'data.frame':	219 obs. of  402 variables; stations should be in same order as other files

#scale/standardize your variables before you add site cover data to 
#unmarked-data-frame
#harvest variables
veghf.1993$NEAR.DIST.harvest<-ifelse(veghf.1993$NEAR.DIST.harvest==0,1,veghf.1993$NEAR.DIST.harvest)
veghf.1993$harvest.dist.s<-veghf.1993$NEAR.DIST.harvest/max(veghf.1993$NEAR.DIST.harvest)
range(veghf.1993$harvest.dist.s)
veghf.1993$harvest.dist.cen<-(veghf.1993$harvest.dist.s-0.5)
veghf.1993$harvest.dist2<-(veghf.1993$harvest.dist.s-0.5)^2
veghf.1993$harvest.dist05<-sqrt(veghf.1993$harvest.dist.s)
veghf.1993$harvest.dist.i<-1/(veghf.1993$harvest.dist.s)

#road variables: changed "paved" to "unimproved" 
veghf.1993$NEAR.DIST.unimproved.road<-ifelse(veghf.1993$NEAR.DIST.unimproved.road==0,1,veghf.1993$NEAR.DIST.unimproved.road)
veghf.1993$road.dist.s<-veghf.1993$NEAR.DIST.unimproved.road/max(veghf.1993$NEAR.DIST.unimproved.road)
range(veghf.1993$road.dist.s)
veghf.1993$road.dist.cen<-(veghf.1993$road.dist.s-0.5)
veghf.1993$road.dist2<-(veghf.1993$road.dist.s-0.5)^2
veghf.1993$road.dist05<-sqrt(veghf.1993$road.dist.s)
veghf.1993$road.dist.i<-1/(veghf.1993$road.dist.s)

#pipeline variables
veghf.1993$NEAR.DIST.pipeline<-ifelse(veghf.1993$NEAR.DIST.pipeline==0,1,veghf.1993$NEAR.DIST.pipeline)
veghf.1993$pipeline.dist.s<-veghf.1993$NEAR.DIST.pipeline/max(veghf.1993$NEAR.DIST.pipeline)
range(veghf.1993$pipeline.dist.s)
veghf.1993$pipeline.dist.cen<-(veghf.1993$pipeline.dist.s-0.5)
veghf.1993$pipeline.dist2<-(veghf.1993$pipeline.dist.s-0.5)^2
veghf.1993$pipeline.dist05<-sqrt(veghf.1993$pipeline.dist.s)
veghf.1993$pipeline.dist.i<-1/(veghf.1993$pipeline.dist.s)

#seismicline variables 
veghf.1993$NEAR.DIST.conventional.seismic <-ifelse(veghf.1993$NEAR.DIST.conventional.seismic==0,1,veghf.1993$NEAR.DIST.conventional.seismic)
veghf.1993$seismic.dist.s<-veghf.1993$NEAR.DIST.conventional.seismic/max(veghf.1993$NEAR.DIST.conventional.seismic)
range(veghf.1993$seismic.dist.s)
veghf.1993$seismic.dist.cen<-(veghf.1993$seismic.dist.s-0.5)
veghf.1993$seismic.dist2<-(veghf.1993$seismic.dist.s-0.5)^2
veghf.1993$seismic.dist05<-sqrt(veghf.1993$seismic.dist.s)
veghf.1993$seismic.dist.i<-1/(veghf.1993$seismic.dist.s)

write.csv(veghf.1993, "0_data/processed/site.cov.csv")


#Spatial weights for Black-throated Green Warbler to account for 
#spatial correlation in occupancy models
wt250<-read.csv("0_data/processed/species weights/BTNW/BTNW_weights_250.csv",header=TRUE)
weights250<-as.data.frame(wt250[,2])
#these files contain "weights" for the observed presences at each station.
#the weights account for spatial autocorrelation in detections from stations at
#the same sites.

#weights are calculated as the proportion of stations within 250 m
#of each station in either year 1 each year (not including the central station) where a given
#species has been detected at least once in that year.
#We expect that the weights should be (+) correlated with probability of initial
#occupancy (year 1) at each station, (+) correlated with probability of colonization
#of stations unoccupied in previous year, and (-) correlated with probability of extinction
#(failure to return to station occupied in previous year). These predictions
#were supported by results from occupancy analysis.

site.cov<-c(veghf.1993,weights250)
site.cov$weights250<-site.cov$`wt250[, 2]`
site.cov$`wt250[, 2]`<-NULL
write.csv(site.cov, file="0_data/processed/species weights/BTNW/BTNWsite.cov.check.csv")
#data frame check. Right now, site.cov is a list


site.cov<-read.csv("0_data/processed/species weights/BTNW/BTNWsite.cov.check.csv", header=TRUE)

wts250<-read.csv("0_data/processed/species weights/BTNW/BTNW_weights_250.stacked.csv", header=TRUE) 
autolog.250<-wts250$values



#Get yearly covariates 
#old yearly covariates: gives us the list of sites we want for filtering the new data
yearlycov<-read.csv("0_data/raw/yearlycovariatesfromanotherproject.csv", header=TRUE)


yearlycov.u<-CL.93.18.SS.u[CL.93.18.SS.u$SS %in% levels(as.factor(yearlycov$SS)),]
#now get rid of 2004 (year 11) because that year was
#been excluded from the visits used in these models

yearlycov.u<-yearlycov.u[!yearlycov.u$YEAR==2004,]
nrow(yearlycov.u)#there should be 5475, ordered by station then by year

logged150.s<-yearlycov.u$PROP150.harvest/max(yearlycov.u$PROP150.harvest)
logged_1sqk.s<-yearlycov.u$PROP565.harvest/max(yearlycov.u$PROP565.harvest)
y.prop565.harvest <- yearlycov.u$PROP565.harvest
harvest.dist.y.s<-yearlycov.u$NEAR.DIST.harvest/max(yearlycov.u$NEAR.DIST.harvest)
#How to deal with harvest age
yearlycov.u$MEANAGE.150.harvest<-ifelse(is.na(yearlycov.u$MEANAGE.150.harvest),
                                        max(site.cov$Structure_Stand_Age_v1),
                                        #102.5067, age of oldest forests according to Beaudoin stand age layers
                                        yearlycov.u$MEANAGE.150.harvest)

harvestage150.s<-yearlycov.u$MEANAGE.150.harvest/max(yearlycov.u$MEANAGE.150.harvest)
yearlycov.u$MEANAGE.565.harvest<-ifelse(is.na(yearlycov.u$MEANAGE.565.harvest),
                                        max(site.cov$Structure_Stand_Age_v1),
                                        #102.5067, age of oldest forests according to Beaudoin stand age layers
                                        yearlycov.u$MEANAGE.565.harvest)

harvestage_1sqk.s<-yearlycov.u$MEANAGE.565.harvest/max(yearlycov.u$MEANAGE.565.harvest)

seismic150.s<-yearlycov.u$PROP150.conventional.seismic/max(yearlycov.u$PROP150.conventional.seismic)
seismic_1sqk.s<-yearlycov.u$PROP565.conventional.seismic/max(yearlycov.u$PROP565.conventional.seismic)
y.seismic.dist<-yearlycov.u$NEAR.DIST.conventional.seismic/max(yearlycov.u$NEAR.DIST.conventional.seismic)
#How to deal with conventional.seismic age
yearlycov.u$MEANAGE.150.conventional.seismic<-ifelse(is.na(yearlycov.u$MEANAGE.150.conventional.seismic),
                                                     max(site.cov$Structure_Stand_Age_v1),
                                                     #102.5067, age of oldest forests according to Beaudoin stand age layers
                                                     yearlycov.u$MEANAGE.150.conventional.seismic)

seismicage150.s<-yearlycov.u$MEANAGE.150.conventional.seismic/max(yearlycov.u$MEANAGE.150.conventional.seismic)
yearlycov.u$MEANAGE.565.conventional.seismic<-ifelse(is.na(yearlycov.u$MEANAGE.565.conventional.seismic),
                                                     max(site.cov$Structure_Stand_Age_v1),
                                                     #102.5067, age of oldest forests according to Beaudoin stand age layers
                                                     yearlycov.u$MEANAGE.565.conventional.seismic)

seismicage_1sqk.s<-yearlycov.u$MEANAGE.565.conventional.seismic/max(yearlycov.u$MEANAGE.565.conventional.seismic)


y.pipeline.dist<-yearlycov.u$NEAR.DIST.pipeline/max(yearlycov.u$NEAR.DIST.pipeline)
y.road.dist<-yearlycov.u$NEAR.DIST.unimproved.road/max(yearlycov.u$NEAR.DIST.unimproved.road)

# proportion of footprints within 150m and 565m of point count station (yearly covariate)

#conventional seismic
y.seismic150.s<-yearlycov.u$PROP150.conventional.seismic/max(yearlycov.u$PROP150.conventional.seismic)

#unimproved road
y.prop150.unimprovedroad <- yearlycov.u$PROP150.unimproved.road/max(yearlycov.u$PROP150.unimproved.road)
y.prop565.unimprovedroad <- yearlycov.u$PROP565.unimproved.road/max(yearlycov.u$PROP565.unimproved.road)

#pipeline
y.prop150.pipeline <- yearlycov.u$PROP150.pipeline/max(yearlycov.u$PROP150.pipeline)
y.prop565.pipeline <- yearlycov.u$PROP565.pipeline/max(yearlycov.u$PROP565.pipeline)

#low impact seismic- not doing anymore because there isnt any 
#prop150.lowimpactseis <- yearlycov.u$PROP150.low.impact.seismic/max(yearlycov.u$PROP150.low.impact.seismic)
#prop565.lowimpactseis <- yearlycov.u$PROP565.low.impact.seismic/max(yearlycov.u$PROP565.low.impact.seismic)


yearly.cov<-data.frame(y.prop565.harvest, logged150.s,logged_1sqk.s,seismic150.s,seismic_1sqk.s,
                       harvest.dist.y.s, harvestage150.s, harvestage_1sqk.s, 
                       y.seismic.dist, seismicage150.s, seismicage_1sqk.s, 
                       y.road.dist, y.prop150.unimprovedroad, y.pipeline.dist,
                       y.prop565.unimprovedroad, y.prop150.pipeline, y.prop565.pipeline)

write.csv(yearly.cov, "0_data/processed/yearlycovariates.csv")
