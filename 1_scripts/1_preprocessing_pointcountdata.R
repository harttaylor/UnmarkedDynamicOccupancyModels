#library(XLConnect)
library(dplyr)
library(tidyr)
setwd("C:/Users/hartt/Documents/Chapter 1")

#1. Get point counts from uncut forests

spp<-read.csv(file="0_data/raw/RQST_LIONEL_93_18_CL_PTCOUNT50m1ha.csv", header=TRUE)   
#pulls Calling Lake point counts into R  
#212676 obs. of  9 variables

#For this analysis, we want to use all the years of data, but not all 
#of the sites. Also, we are using birds within 100 m of point counts, but at the 1-ha
#sites, we are only using birds occurring within the original boundaries of those
#sites (50 m)

spp_CLCON<-spp[which(substring(spp$PKEY,1,5)=="CL:C-"),]
#filters out control point counts from Calling Lake point counts 
spp_CLISO<-spp[which(substring(spp$PKEY,1,5)=="CL:F-"),]
#filters out forest fragment point counts from Calling Lake point counts 
spp_CLRIP<-spp[which(substring(spp$PKEY,1,5)=="CL:R-"),]
#filters out riparian site point counts from Calling Lake point counts 

spp_cc<-rbind(spp_CLCON,spp_CLISO,spp_CLRIP)
#163917 obs. of  9 variables

spp_cc<-spp_cc[which(spp_cc$DURATION==1),]
#filters out species detected during 5-minute point counts
#163899 obs. of  9 variables



#Each observation in this data frame is the number of a species 
#within a given distance interval (0-50 m, 50-100 m, 100-Inf) and
#a given duration interval within a point count visit to each station.

#A station visit can have multiple observations.

#The column WTABUND (weighted abundance) upweights the number counted 
#according to behavioral observations suggestive of local breeding. The
#highest weight would be an abundance of 2 for an observation (1 male+
#1 female). Observations can also be downweighted if they do not represent
#significant evidence of use.
#So birds that were observed before or after point counts or as flyovers
#or flying through without "using" the point count have already been 
#treated as having "zero" abundance. Any waterfowl or waders within the processed
#data files for this study should probably then be treated as "present"
#since flyovers will have been eliminated.

#2. remove birds > 100 m distant
spp_cc<-spp_cc[!spp_cc$DISTANCE>2,]
#removes observations >100 m away because it's likely such individuals have
#been counted at an adjacent station
#DISTANCE==1 = 0-50 m
#DISTANCE==2 = 50-100 m
#DISTANCE==3 = 100 m - unlimited distance
nrow(spp_cc)
#151346 obs.


#3. merge point count data with variables affecting detectability.
pk<-read.csv(file="0_data/raw/RQST_LIONEL_93_18_CL_PKEY.csv", header=TRUE)
#pulls Calling Lake point count data (date, time of morning, observer) into R  

merge1<-merge(pk, spp_cc, by = c("PCODE","PKEY"))
#a left-handed join: by default, the only observations that will be merged
#are those observations with PKEY values in both data frames 
# Filter the merge1 dataframe for only the sites in sites_to_keep
merge1_filtered <- merge1[merge1$SS %in% sites_to_keep, ]

# Check the structure of the filtered dataframe
str(merge1_filtered)
#variables from "pk" in left hand columns followed by variables from "spp_cc"
#the same variables from both data frames, e.g. PCODE, are renamed, e.g. PCODE.x and PCODE.y
#unless they are the variables used in merging, e.g. PKEY

merge1_filtered$VISIT<-paste0(merge1_filtered$SS,"_",merge1_filtered$YYYY,"_",merge1_filtered$ROUND)
#creates a variable for station-visit
write.csv(merge1_filtered, file = "0_data/processed/1_merge1.csv")

#4. calculate total abundance per species per site visit in each cutblock 
merge1<-read.csv("0_data/processed/1_merge1.csv", header=TRUE)
merge1$ABUND<-is.atomic(merge1$WT_ABUND)
tapply.spp<-tapply(merge1$ABUND, list(merge1$VISIT, merge1$SPECIES), sum, na.rm=TRUE)
tapply.spp<-data.frame(tapply.spp)
tapply.spp$VISIT<-row.names(tapply.spp)
write.csv(tapply.spp, file = "0_data/processed/2_tapply.spp.csv")  
#in the past, I manually created "3_spp_count" from "2_tapply.spp". Automate this
#with R.
spp_count<-tapply.spp%>%
  mutate_at(vars(ALFL:YWAR), ~replace(., is.na(.), 0))
write.csv(spp_count, file = "0_data/processed/3_spp_count.csv")  

#5. merge summarized bird count data per visit and visit-specific environmental 
#covariates in "merge1" by the variable VISIT
merge1<-read.csv("0_data/processed/1_merge1.csv", header=TRUE)
birdcount<-read.csv("0_data/processed/3_spp_count.csv", header=TRUE)
VISIT<-merge1$VISIT
SITE<-merge1$SITE
SS<-merge1$SS
YYYY<-merge1$YYYY
YEAR<-merge1$YYYY-1993
ROUND<-merge1$ROUND
JULIAN<-merge1$JULIAN
HR<-merge1$HR
MIN<-merge1$MIN
TIME<-merge1$HR+(merge1$MIN/60)

newdf<-data.frame(VISIT,SITE,SS,YYYY,YEAR,ROUND,JULIAN,HR,MIN,TIME)
#create "newdf" to have a smaller number of covariates to add
#to bird counts. Right now, "merge1" and "newdf" have multiple rows
#per site visit because the bird data have not been summarized
#in these files.
nrow(newdf)#151295 observations
#you need to whittle the data frame down to 1 row per visit,
#otherwise you will have many false duplicate observations when 
#you add the summarized bird data.
newdf<-unique(newdf)
nrow(newdf)#25498 observations


merge2<-merge(newdf, birdcount, by = c("VISIT"), all.x=TRUE)  
#a left outer join. Visits in "newdf" missing data from "birdcount"
#will still have observations in the merged file, but missing bird data
#(if there are any) will be "NA"
str(merge2)
nrow(merge2)#25498 observations
write.csv(merge2, "0_data/processed/4_merge2.csv") #output check

