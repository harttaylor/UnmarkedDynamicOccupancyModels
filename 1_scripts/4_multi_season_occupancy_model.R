# Load required data, create an unmarkedMultiFrame, and fit a multi-season occupancy model using the colext function from the unmarked package

setwd("C:/Users/hartt/Documents/Chapter 1/UnmarkedDynamicOccupancyModels")

# Load required libraries
library(unmarked)
library(scales)
library(MuMIn)
library(jpeg)
library(ggplot2)
library(dplyr)

# Load Data
site.cov<-read.csv("0_data/processed/species weights/BTNW/BTNWsite.cov.check.csv", header=TRUE)
yearlycov<-read.csv("0_data/raw/Calling Lake/yearlycovariatesfromanotherproject.csv", header=TRUE)


"'logLik.try-error'"<-function(object,...) structure(-Inf, nobs = 1L, df = 1L, class = "logLik")
                      
#purpose of this line is to assign a ridiculously small log-likelihood to a model if the model fails to run, so that
#the model can still be compared to other models within the model.sel function as part of automated model selection.

#Single Species MSOM

data<-read.csv(paste0("0_data/processed/6_rearranged.BTNW.csv"),header=TRUE)
y<-data[,c(3:ncol(data))]#colums 2 through 101
S <- nrow(data) # number of sites  (219)
J <- 4 # number of secondary sampling occasions    (4 per year)
T <- (ncol(data)-2+1)/J # number of primary periods =(ncol(data)-2)/J   (25 years of observations over 26 years)


# Create an unmarkedMultiFrame, which is the data format required for unmarked models
bird<-unmarkedMultFrame(y=y,
                        siteCovs=site.cov,
                        numPrimary=25,
                        yearlySiteCovs=yearly.cov,
                        obsCovs=list(Julian=recast.julian[,c(3:ncol(recast.julian))], 
                                     Time=recast.time[,c(3:ncol(recast.time))],
                                     jcen=recast.jcen[,c(3:ncol(recast.jcen))], 
                                     jcen2=recast.jcen2[,c(3:ncol(recast.jcen2))], 
                                     tcen=recast.tcen[,c(3:ncol(recast.tcen))], 
                                     tcen2=recast.tcen2[,c(3:ncol(recast.tcen2))]))  # I have my obscovs like this, but yours could be in a single dataframe like the site and yearly covs


# I always get this warning message, its ok:
#Warning messages:
#1: siteCovs contains characters. Converting them to factors. 
#2: yearlySiteCovs contains characters. Converting them to factors. 
#cannot use dredge function because there are missing values
#psiformula=probability of initial occupancy in Year 1
#gammaformula=probability of colonizing previously unoccupied location
#epsilonformula=probability of going extinct at previously occupied location
#pformula=probability of detection given occupancy of location



# Define the MSOM model using colext function from the unmarked package
# psiformula: Probability of initial occupancy
# gammaformula: Probability of colonization
# epsilonformula: Probability of local extinction
# pformula: Probability of detection given presence
# I tested covarites in different model stages, but depending on how many covariates you are using you could do it in one model stage 

# Fit the models 
#Stage 1
#test for evidence of spatial autocorrelation and effect  
null<-try(colext(psiformula = ~1, gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))    #null model
autolog250<-try(colext(psiformula = ~weights250, gammaformula = ~autolog.250, epsilonformula = ~autolog.250, pformula = ~1, data=bird))    #autocorrelation model
#positive nonzero effect of autocorrelation variable on psi and gamma but not epsilon
#in other words, stations with more points within 250 m that contain BTNW are more likely
#to be occupied by BTNW in year 1, and unoccupied stations with more BTNW points within a given year
#are more likely to have BTNW in the next year.
#test for effect of age/ percent conifer on initial occupancy and the interaction and quadratics of age 
age<-try(colext(psiformula = ~age.s +weights250
                , gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))    #stand age model 
age.2 <-try(colext(psiformula = ~ age.s+age.s.2 +weights250
                   , gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))
#poly(age.s, 2) 
#Age n.s. Maybe makes sense since in 1993 all of the forest around points was pretty old
whitespruce<-try(colext(psiformula = ~whitespruce.s +weights250
                        , gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))
whitespruce.2 <-try(colext(psiformula = ~whitespruce.s+whitespruce.s.2 +weights250
                           , gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))
#poly(whtespruce.s, 2) 
#Check interaction and additive age/conifer terms 
age.conifer<-try(colext(psiformula = ~whitespruce.s + age.s +weights250
                        , gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))
ageXconifer<-try(colext(psiformula = ~whitespruce.s * age.s +weights250
                        , gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))
#aggregation hypothesis (Mattson et al. 2013)
#model1 <- colext(psiformula = ~weights250, gammaformula = ~autolog.250, 
# epsilonformula = ~autolog.250, pformula = ~1, data = bird)
results.stage1 <- model.sel('Additiveage.conifer'= age.conifer, 
                            'Interactionage.conifer'= ageXconifer, 
                            'Age'=age, 
                            'Age2'=age.2, 
                            'WSpruce2'=whitespruce.2, 
                            'WSpruce' = whitespruce, 
                            'Null'=null, 
                           'Autologistic'=autolog250, 
                            rank=AIC)           
bestmodel.stage1<-get.models(results.stage1, subset = 1)[[1]]


fms <- fitList('Additiveage.conifer'= age.conifer, 
               'Interactionage.conifer'= ageXconifer, 
               'Age'=age, 
               'Age2'=age.2, 
               'WSpruce2'=whitespruce.2, 
               'WSpruce' = whitespruce, 
               'Null'=null, 
               'Autologistic'=autolog250)          
ms<-modSel(fms)
ms.coef<-coef(ms)
ms.SE<-SE(ms)
output<-as(ms,"data.frame")
write.csv(output, file=paste0("2_outputs/BTNW.stage1modeltable_Mar2023.csv"))

#Stage 2- Look at distance to different linear features 
harvest.psi<-try(colext(psiformula =update(bestmodel.stage1@psiformula, ~.+harvest.dist.s), 
                              gammaformula= bestmodel.stage1@gamformula, 
                              epsilonformula= bestmodel.stage1@epsformula, 
                              pformula= ~1, data=bird))

distharvest565.gam<-try(colext(psiformula =bestmodel.stage1@psiformula, 
                               gammaformula= update(bestmodel.stage1@gamformula, ~.+harvest.dist.y.s), 
                               epsilonformula= bestmodel.stage1@epsformula, 
                               pformula= ~1, data=bird))
distharvest565.eps<-try(colext(psiformula =bestmodel.stage1@psiformula,
                               gammaformula= bestmodel.stage1@gamformula, 
                               epsilonformula= update(bestmodel.stage1@epsformula, ~.+harvest.dist.y.s), 
                               pformula= ~1, data=bird))
#Try seeing how distance to roads affects colonization and extinction
gam.road<-try(colext(psiformula =update(bestmodel.stage1@psiformula, ~.+road.dist.s), 
gammaformula= update(bestmodel.stage1@gamformula, ~.+y.road.dist), 
epsilonformula= bestmodel.stage1@epsformula, 
pformula= ~1, data=bird))

eps.road<-try(colext(psiformula =update(bestmodel.stage1@psiformula, ~.+road.dist.s), 
gammaformula= bestmodel.stage1@gamformula, 
epsilonformula= update(bestmodel.stage1@epsformula, ~.+y.road.dist), 
pformula= ~1, data=bird))

distpipeline565.psi<-try(colext(psiformula =update(bestmodel.stage1@psiformula,~.+pipeline.dist.s), 
                                gammaformula= bestmodel.stage1@gamformula, 
                                epsilonformula= bestmodel.stage1@epsformula, 
                                pformula= ~1, data=bird))
distpipeline565.eps<-try(colext(psiformula =bestmodel.stage1@psiformula,
                                gammaformula= bestmodel.stage1@gamformula, 
                                epsilonformula= update(bestmodel.stage1@epsformula, ~.+y.pipeline.dist), 
                                pformula= ~1, data=bird))
distpipeline565.gam <- try(colext(psiformula =bestmodel.stage1@psiformula,
                                  gammaformula= update(bestmodel.stage1@gamformula,~.+y.pipeline.dist), 
                                  epsilonformula= bestmodel.stage1@epsformula, 
                                  pformula= ~1, data=bird))

results.stage2 <- model.sel('PsiHarvestdist'=harvest.psi, 'GamDistHarv'=distharvest565.gam, 
                            'EpsDistHarv'=distharvest565.eps,'GamDistRoad'=gam.road, 'EpsDistRoad'=eps.road, 
                            'EpsDistPipe'=distpipeline565.eps, 'PsiDistPipe'=distpipeline565.psi,
                            'GamDistPipe'=distharvest565.gam,'Additiveage.conifer'= age.conifer, 
                            'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 
                            'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, 
                            'Autologistic'=autolog250, rank=AIC)

bestmodel.stage2<-get.models(results.stage2, subset = 1)[[1]]

fms <- fitList('PsiHarvestdist'=harvest.psi, 'GamDistHarv'=distharvest565.gam, 
               'EpsDistHarv'=distharvest565.eps,'GamDistRoad'=gam.road, 'EpsDistRoad'=eps.road, 
               'EpsDistPipe'=distpipeline565.eps, 'PsiDistPipe'=distpipeline565.psi,
               'GamDistPipe'=distharvest565.gam,'Additiveage.conifer'= age.conifer, 
               'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 
               'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, 
               'Autologistic'=autolog250)    

ms<-modSel(fms)
ms.coef<-coef(ms)
ms.SE<-SE(ms)
output<-as(ms,"data.frame")
write.csv(output, file=paste0("2_outputs/BTNW.stage2modeltable_Mar2023.csv"))


# Stage 3 test if model is improved by using distances/proportions to different footprints as covariates
#HARVEST
# Try using 565m, especially for harvest because there is higher variability when icnlduing 565m, intead of just 150m 
propharvest565.gam<-try(colext(psiformula =bestmodel.stage2@psiformula, 
                             gammaformula= update(bestmodel.stage2@gamformula, ~.+logged_1sqk.s), 
                            epsilonformula= bestmodel.stage2@epsformula, 
                             pformula= ~1, data=bird))

propharvest565.eps<-try(colext(psiformula =bestmodel.stage2@psiformula,
                           gammaformula= bestmodel.stage2@gamformula, 
                            epsilonformula= update(bestmodel.stage2@epsformula, ~.+logged_1sqk.s), 
                            pformula= ~1, data=bird))

#CONVENTIONAL SEISMIC
#convseis.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+ seismic.dist.s), 
                      #  gammaformula= bestmodel.stage2@gamformula, 
                      #  epsilonformula= bestmodel.stage2@epsformula, 
                      #  pformula= ~1, data=bird))#AIC: 16992.7 
#Try seeing how distance to seismic lines affects colonization and extinction
#gam.seismic<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+seismic.dist.s), 
                       # gammaformula= update(bestmodel.stage2@gamformula, ~.+y.seismic.dist), 
                       # epsilonformula= bestmodel.stage2@epsformula, 
                       # pformula= ~1, data=bird))
#eps.seismic<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+seismic.dist.s), 
                       # gammaformula= bestmodel.stage2@gamformula, 
                       # epsilonformula= update(bestmodel.stage2@epsformula, ~.+y.seismic.dist), 
                       # pformula= ~1, data=bird))
#Try seeing how proportion seismic lines affects colonization and extinction
#Try looking at just 565m
#propconvseis150.gam<-try(colext(psiformula =bestmodel.stage2@psiformula, 
                             #  gammaformula= update(bestmodel.stage2@gamformula, ~.+seismic150.s), 
                              # epsilonformula= bestmodel.stage2@epsformula, 
                              # pformula= ~1, data=bird))
#propconvseis565.gam<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+ seismic.dist.s), 
                      #    gammaformula= update(bestmodel.stage2@gamformula, ~.+seismic_1sqk.s), 
                      #    epsilonformula= bestmodel.stage2@epsformula, 
                       #   pformula= ~1, data=bird))
#propconvseis150.eps<-try(colext(psiformula =bestmodel.stage2@psiformula, 
                              # gammaformula= bestmodel.stage2@gamformula, 
                              # epsilonformula= update(bestmodel.stage2@epsformula, ~.+seismic150.s), 
                              # pformula= ~1, data=bird))
#propconvseis565.eps<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+ seismic.dist.s), 
                         #   gammaformula= bestmodel.stage2@gamformula, 
                         #   epsilonformula= update(bestmodel.stage2@epsformula, ~.+seismic_1sqk.s), 
                           # pformula= ~1, data=bird))

#UNIMRPOVED ROADS
#road.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+road.dist.s), 
                   #  gammaformula= bestmodel.stage2@gamformula, 
                    # epsilonformula= bestmodel.stage2@epsformula, 
                     #pformula= ~1, data=bird)) 


#Try seeing how proportion of roads affects colonization and extinction
#look just at 565m for now 
#proproad150.gam<-try(colext(psiformula =bestmodel.stage2@psiformula, 
                    # gammaformula= update(bestmodel.stage2@gamformula, ~.+prop150.unimprovedroad), 
                    # epsilonformula= bestmodel.stage2@epsformula, 
                    # pformula= ~1, data=bird))
proproad565.gam<-try(colext(psiformula =bestmodel.stage2@psiformula, 
                        gammaformula= update(bestmodel.stage2@gamformula, ~.+y.prop565.unimprovedroad), 
                        epsilonformula= bestmodel.stage2@epsformula, 
                         pformula= ~1, data=bird))
#distroad565.gam<-try(colext(psiformula =bestmodel.stage2@psiformula, 
                          #  gammaformula= update(bestmodel.stage2@gamformula, ~.+y.road.dist), 
                          #  epsilonformula= bestmodel.stage2@epsformula, 
                           # pformula= ~1, data=bird))


#proproad150.eps<-try(colext(psiformula =bestmodel.stage2@psiformula, 
                    # gammaformula= bestmodel.stage2@gamformula, 
                    # epsilonformula= update(bestmodel.stage2@epsformula, ~.+prop150.road), 
                     #pformula= ~1, data=bird))
proproad565.eps<-try(colext(psiformula =bestmodel.stage2@psiformula,
                         gammaformula= bestmodel.stage2@gamformula, 
                         epsilonformula= update(bestmodel.stage2@epsformula, ~.+y.prop565.unimprovedroad), 
                         pformula= ~1, data=bird))
#distroad565.eps<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+road.dist.s),
                        #    gammaformula= bestmodel.stage2@gamformula, 
                        #    epsilonformula= update(bestmodel.stage2@epsformula, ~.+y.road.dist), 
                         #   pformula= ~1, data=bird))


#PIPELINE
proppipeline565.gam<-try(colext(psiformula =update(bestmodel.stage2@psiformula,~.+pipeline.dist.s), 
                                gammaformula= bestmodel.stage2@gamformula, 
                                epsilonformula= bestmodel.stage2@epsformula, 
                                pformula= ~1, data=bird))


#proppipeline565.eps<-try(colext(psiformula =bestmodel.stage2@psiformula,
                            #    gammaformula= bestmodel.stage2@gamformula, 
                            #    epsilonformula= update(bestmodel.stage2@epsformula, ~.+y.prop565.pipeline), 
                            #    pformula= ~1, data=bird))

#'GamDistSeismic'=gam.seismic, 'EpsDistSeis'=eps.seismic, 'PsiDistRoad'= road.psi, 'GamDistRoad'=gam.road, 'EpsRoad'=eps.road,

results.stage3 <- model.sel('PsiHarvestdist'=harvest.psi, 'GamDistHarv'=distharvest565.gam, 
                            'EpsDistHarv'=distharvest565.eps,'GamDistRoad'=gam.road, 'EpsDistRoad'=eps.road, 
                           'EpsDistPipe'=distpipeline565.eps, 
                            'PsiDistPipe'=distpipeline565.psi,
                            'GamDistPipe'=distharvest565.gam,'Additiveage.conifer'= age.conifer, 
                            'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 
                            'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, 
                            'Autologistic'=autolog250, 'GamProp565Harv'=propharvest565.gam,
                            'EpsProp565Harv'=propharvest565.eps, 'GamProp565Road'=proproad565.gam,
                            'EpsPropRoad565'=proproad565.eps,
                            'GamProp565Pipe'=proppipeline565.gam, rank=AIC)
bestmodel.stage3<-get.models(results.stage3, subset = 1)[[1]]

#plogis(2.3)
fms <- fitList('PsiHarvestdist'=harvest.psi, 'GamDistHarv'=distharvest565.gam, 
               'EpsDistHarv'=distharvest565.eps,'GamDistRoad'=gam.road, 'EpsDistRoad'=eps.road, 
               'EpsDistPipe'=distpipeline565.eps, 'PsiDistPipe'=distpipeline565.psi,
               'GamDistPipe'=distharvest565.gam,'Additiveage.conifer'= age.conifer, 
               'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 
               'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, 
               'Autologistic'=autolog250, 'GamProp565Harv'=propharvest565.gam,
               'EpsProp565Harv'=propharvest565.eps, 'GamProp565Road'=proproad565.gam,
               'EpsProp565Pipe' = proppipeline565.eps, 'EpsPropRoad565'=proproad565.eps,
               'GamProp565Pipe'=proppipeline565.gam)    

ms<-modSel(fms)#changed modsel to modSel
ms.coef<-coef(ms)
ms.SE<-SE(ms)
output<-as(ms,"data.frame")
write.csv(output, file=paste0("2_outputs/BTNW.stage3modeltable_Mar2023.csv"))


#Stage 4 - detection covariates
det.julian<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~Julian, data=bird))     
det.julian.sq<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~jcen+jcen2, data=bird))     
det.time<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~Time, data=bird))     
det.time.sq<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~tcen+tcen2, data=bird))     
det.juliantime<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~jcen+jcen2+Time, data=bird))     
det.juliantime.sq<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~Julian+tcen+tcen2, data=bird))     


results.stage4 <- model.sel('PsiHarvestdist'=harvest.psi, 'GamDistHarv'=distharvest565.gam, 
                            'EpsDistHarv'=distharvest565.eps,'GamDistRoad'=gam.road, 'EpsDistRoad'=eps.road, 
                            'EpsDistPipe'=distpipeline565.eps, 'PsiDistPipe'=distpipeline565.psi,
                            'GamDistPipe'=distharvest565.gam,'Additiveage.conifer'= age.conifer, 
                            'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 
                            'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, 
                            'Autologistic'=autolog250, 'GamProp565Harv'=propharvest565.gam,
                            'EpsProp565Harv'=propharvest565.eps, 'GamProp565Road'=proproad565.gam,
                            'EpsPropRoad565'=proproad565.eps,
                            'GamProp565Pipe'=proppipeline565.gam,
                           'Julian' = det.julian, 
                            'Time' = det.time, 'Julian.sq' = det.julian.sq, 'Time.sq' = det.time.sq, 
                            'Julian+Time' = det.juliantime, rank=AIC)                 

bestmodel.stage4<-get.models(results.stage4, subset = 1)[[1]]
plogis(-1.07)
fms <- fitList('PsiHarvestdist'=harvest.psi, 'GamDistHarv'=distharvest565.gam, 
               'EpsDistHarv'=distharvest565.eps,'GamDistRoad'=gam.road, 'EpsDistRoad'=eps.road, 
               'EpsDistPipe'=distpipeline565.eps, 'PsiDistPipe'=distpipeline565.psi,
               'GamDistPipe'=distharvest565.gam,'Additiveage.conifer'= age.conifer, 
               'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 
               'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, 
               'Autologistic'=autolog250, 'GamProp565Harv'=propharvest565.gam,
               'EpsProp565Harv'=propharvest565.eps, 'GamProp565Road'=proproad565.gam,
               'EpsPropRoad565'=proproad565.eps,
               'GamProp565Pipe'=proppipeline565.gam,
               'Julian' = det.julian, 
               'Time' = det.time, 'Julian.sq' = det.julian.sq, 'Time.sq' = det.time.sq, 
               'Julian+Time' = det.juliantime)           
ms<-modSel(fms)
ms.coef<-coef(ms)
ms.SE<-SE(ms)
output<-as(ms,"data.frame")
write.csv(output, file=paste0("2_outputs/BTNW.stage4modeltable_Mar2023.csv"))








########################
#use bootstrap to get estimates of psi (occupancy) from derived parameter estimates
bestmodel.stage4.boot<-nonparboot(bestmodel.stage4, B=100)
str(bestmodel.stage4.boot@projected)
#a multi-dimensional array: num [1:2, 1:25, 1:219]

smooth.occ <- smoothed(bestmodel.stage4.boot)[2,] 
Occupancy <- data.frame(smooth.occ)
#this is the average probability of occupancy
#(proportion of sites occupied)
#in each year across all stations
Occupancy$year <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
Occupancy$smooth.occSE<-bestmodel.stage4.boot@projected.mean.bsse[1,]
Occupancy$smooth.occ95UCL<-Occupancy$smooth.occ+1.96*Occupancy$smooth.occSE
Occupancy$smooth.occ95LCL<-Occupancy$smooth.occ-1.96*Occupancy$smooth.occSE
Occupancy$smooth.unocc<-1-Occupancy$smooth.occ

ggplot(Occupancy, aes(x = year, y = smooth.occ)) +
  geom_errorbar(aes(ymin=smooth.occ95LCL, ymax=smooth.occ95UCL), colour="black", width=.1) +
  geom_line() +
  geom_point() + 
  #ggtitle("Average Site Occupancy vs. Time") +
  ylab("Average Occupancy Probability") +
  xlab("Year")+my.theme

#Average Extinction Rate
(0.6121221-0.2857127)/25
#0.0086002  *100 = just under 1% of sites per year

Occupancy$linearpredUnocc<-Occupancy[1,]$smooth.unocc+(Occupancy$year-1)*((0.8196931-0.6046926)/25)
ggplot(Occupancy, aes(x = year, y = smooth.unocc)) +
  geom_line() +
  geom_point() + 
  ggtitle("Average Site Extinction vs. Time") +
  ylab("Extinction Probability") + ylim(0,1)+
  xlab("Year")+my.theme

ggplot(Occupancy, aes(x = year, y = linearpredUnocc)) +
  geom_line() +
  geom_point() + 
 # ggtitle("Average Site Extinction vs. Time") +
  ylab("Extinction Probability") + ylim(0.2,0.6)+
  xlab("Year")+my.theme

#use bootstrap to get estimates of psi (occupancy) from derived parameter estimates
bestmodel.stage4.boot<-nonparboot(bestmodel.stage4, B=100)
str(bestmodel.stage4.boot@projected)
#a multi-dimensional array: num [1:2, 1:25, 1:219]
#3d array whose number is equal to 219(#sites)*25(#years)*2(1st vector for prob. absence,
#2nd vector for prob. occupancy)
#http://www.statmethods.net/advstats/bootstrapping.html for a description of the bootstrap object structure
#boot( ) calls the statistic function R times. Each time, it generates a set of 
#random indices, with replacement, from the integers 1:nrow(data). These indices 
#are used within the statistic function to select a sample. The statistics are 
#calculated on the sample and the results are accumulated in the bootobject. The 
#bootobject structure includes

#basic bootstrap object element description 
#t0           The observed values of k statistics applied to the orginal data.  
#t            An R x k matrix where each row is a bootstrap replicate of the k statistics. 

#You can access these as bootobject$t0 and bootobject$t.

#Once you generate the bootstrap samples, print(bootobject) and plot(bootobject) can be used to examine the 
#results. If the results look reasonable, you can use boot.ci( ) function to obtain confidence intervals for 
#the statistic(s). 
bestmodel.stage4.boot@projected[2,25,1:219]
#this is how you call up predicted probability of 
#occupancy based on the bootstrap samples for site 1
#(which would be CL:C-001-1-1)
#in the colext unmarked Frame (one of the unlogged 
#"control" sites) in each of the 25 years

bestmodel.stage4.boot@projected[2,,9]
#this is how you call up predicted probability of 
#occupancy based on the bootstrap samples for site 9
#(CL:C-010-3:2)
bestmodel.stage4.boot@projected[2,,]
#this is how you call up predicted probability of occupancy 
#based on the bootstrap samples 
#in the colext unmarked Frame, for all sites in each of 
#the 25 years. An 25x219 matrix

bestmodel.stage4.boot@projected[1,,]
#this is how you call up predicted probability that site 
#is unoccupied based on the original data for site 1
#in the colext unmarked Frame
smooth.occ <- smoothed(bestmodel.stage4.boot)[2,] 
Occupancy <- data.frame(smooth.occ)
#I THINK this is the average probability of occupancy
#(proportion of sites occupied)
#in each year across all stations
smooth.occ$year <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)
bestmodel.stage4.boot@projected.mean.bsse
#the standard errors for projected occupancy. 
#Note that SE is same for probability of
#presence and probability of absence
bestmodel.stage4.boot@projected.mean.bsse[2,]
#Also note that there is 1 SE calculated for each year 
#but not for each site
#unlike probability of occupancy

#Making some plots 
library(ggplot2)
model4 <- colext(psiformula = ~weights250 + pipeline.dist.s, 
                 gammaformula = ~autolog.250,
                 epsilonformula = ~autolog.250 + harvest.dist.y.s, 
                 pformula = ~jcen + jcen2, data = bird)
#make graph showing harvest distance and extinction probability 
max(site.cov$NEAR.DIST.harvest)#1000

#Look at exticntion probability and amount logged within 1 sq km 
nd <- data.frame(weights250=c(rep(0.5,11)),
                 autolog.250=c(rep(0.5,11)),
                 pipeline.dist.s=c(rep(0.5,11)),
                 harvest.dist.y.s=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
      
P.Ext.HarvDist <- predict(bestmodel.stage4, type="ext", newdata=nd)

Actual.HarvestDistance<-nd$harvest.dist.y.s*1000
nd.preds<-cbind(nd, P.Ext.HarvDist, Actual.HarvestDistance)

ggplot(nd.preds, aes(x = Actual.HarvestDistance, y = P.Ext.HarvDist$Predicted)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=0) +
  geom_line() +
  geom_point() + 
  ylab("Extinction Probability") +
  xlab("Distance to harvest (m)") +theme_classic()
max(P.Ext.Logged1sqkm$Predicted)
min(P.Ext.Logged1sqkm$Predicted)


#Look at initial occupancy probability and distance to pipeline 
nd <- data.frame(weights250=c(rep(0.5,11)),
                 autolog.250=c(rep(0.5,11)),
                 harvest.dist.y.s=c(rep(0.5,11)),
                 pipeline.dist.s=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

P.Occ.PipeDist <- predict(bestmodel.stage4, type="ext", newdata=nd)
max(sitecovs$NEAR.DIST.pipeline)#1000
Actual.PipeDistance<-nd$pipeline.dist.s*1000
nd.preds<-cbind(nd, P.Occ.PipeDist, Actual.PipeDistance)

ggplot(nd.preds, aes(x = Actual.PipeDistance, y = P.Occ.PipeDist$Predicted)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=0) +
  geom_line() +
  geom_point() + 
  ylab("Initial Occupancy Probability") +
  xlab("Distance to pipeline (m)") +theme_classic()
max(P.Ext.Logged1sqkm$Predicted)
min(P.Ext.Logged1sqkm$Predicted)
#Make graph showing colnization probability with and without spatial effect 
#Colonization probability and proportion logged within 1km2 of point count station 
model1 <- colext(psiformula = ~weights250 + PROP565.harvest, 
                 gammaformula = ~autolog.250 + logged_1sqk.s,
                 epsilonformula = ~autolog.250 + harvest.dist.y.s, 
                 pformula = ~jcen + jcen2, data = bird)

nd <- data.frame(weights250=c(rep(0.5,11)),
                autolog.250=c(rep(0.5,11)),
                PROP565.harvest=c(rep(0.5,11)),
                Treatment=c(rep("Control",11)),
                harvest.dist.y.s=c(rep(0.5,11)),
                logged_1sqk.s=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) #data.frame(harvestage150.s)
P.Spat.Harvest <- predict(model1, type="col", newdata=nd)

nd.preds<-cbind(nd, P.Col.Logged1km,Actual.HarvestAmount1km)

ggplot(nd.preds, aes(x = Actual.HarvestAmount1km, y = P.Col.Logged1km$Predicted)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.05) +
  geom_line() +
  geom_point() + 
  ylab("Colonization Probability") +
  xlab("Proportion Harvested (within 565m of point count station)")+theme_classic()

#Colonization probability and proportion logged within 1km2 of point count station without spatial autocovariates
model2 <- colext(psiformula = ~weights250 + PROP565.harvest, 
                 gammaformula = ~ logged_1sqk.s,
                 epsilonformula = ~autolog.250 + harvest.dist.y.s, 
                 pformula = ~jcen + jcen2, data = bird)

nd <- data.frame(weights250=c(rep(0.5,11)),
                 autolog.250=c(rep(0.5,11)),
                 PROP565.harvest=c(rep(0.5,11)),
                 Treatment=c(rep("Control",11)),
                 harvest.dist.y.s=c(rep(0.5,11)),
                 logged_1sqk.s=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) #data.frame(harvestage150.s)
P.Col.Logged1km <- predict(model2, type="col", newdata=nd)

nd.preds<-cbind(nd, P.Col.Logged1km,Actual.HarvestAmount1km)

ggplot(nd.preds, aes(x = Actual.HarvestAmount1km, y = P.Col.Logged1km$Predicted)) +
  #geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.05) +
  geom_line() +
  geom_point() + 
  ylab("Colonization Probability") +
  xlab("Proportion Harvested (within 565m of point count station)")+theme_classic()

#Show spatial vs non-spatial model on same graph 
nd$P.Spatial.Weights <- predict(model1, type="col", newdata=nd)
nd$P.NoSpatial.Weights <- predict(model2, type="col", newdata=nd)
#P.Prop.Harvest <- predict(model4, type="psi", newdata=nd)

nd$Actual.Prop.Harvest <- nd$logged_1sqk.s*0.6506915

write.csv(nd, "0_data/processed/nd.csv")
nddd <- read.csv("0_data/processed/nd.csv")
dataTH <- read.csv("0_data/processed/ndTH.csv", header = T)
dataTAE <- read.csv("0_data/processed/ndTH.csv", header = T)

ggplot(nddd, aes(x = Actual.Prop.Harvest, y = MeanPred, group=Model, color=Model)) +
  geom_errorbar(aes(ymin=LCLPred, ymax=UCLPred), colour="gray", width=.05) +
  geom_line() +
  scale_fill_manual(values=c("#999999", "#E69F00"), 
                    name="Model",
                    breaks=c("Spatial", "NonSpatial"),
                    labels=c("Spatial", "Non-Spatial")) +
  geom_point() + 
  ylab("Predicted Colonization Probability") +
  xlab("Proportion Harvested (within 565m of point count site)")+my.theme
#Initial occupancy probability and proportion harvested within 565m 
model4 <- colext(psiformula = ~weights250, gammaformula = ~autolog.250 + logged_1sqk.s, 
                 epsilonformula = ~autolog.250 + harvestage150.s, pformula = ~jcen + jcen2, data = bird)

nd <- data.frame(weights250=c(rep(0.5,11)),
                 autolog.250=c(rep(0.5,11)),
                 PROP565.harvest=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                 logged_1sqk.s=c(rep(0.5,11)),
                 harvest.dist.y.s=c(rep(0.5,11)),
                 Treatment=c(rep("Control",11)))

bestmodel.stage4
P.Prop.Harvest <- predict(model4, type="psi", newdata=nd)
max(prop565.harvest)

Actual.Prop.Harvest <- nd$PROP565.harvest*0.6506915

nd.preds<-cbind(nd, P.Prop.Harvest,Actual.Prop.Harvest)

ggplot(nd.preds, aes(x = Actual.Prop.Harvest, y = P.Prop.Harvest$Predicted)) +
  #geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.05) +
  geom_line() +
  geom_point() + 
  ylab("Initial Occupancy Probability") +
  xlab("Proportion harvested (within 565m of point count station)")+my.theme

#Initial occupancy probability and autologistic variable 



ndLL <- data.frame(weights250=c(rep(0.5,11)),
                 autolog.250=c(rep(0.5,11)),
                 PROP565.harvest=c(rep(0.5,11)),
                 logged_1sqk.s=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),#ALIGN WITH PROP565
                 harvest.dist.y.s=c(rep(0.5,11)),#ALIGN WITH PROP565
                 Treatment=c(rep("Control",11)))

bestmodel.stage4

model4 <- colext(psiformula = ~weights250+ Treatment, gammaformula = ~autolog.250 +Treatment +logged_1sqk.s, 
                 epsilonformula = ~autolog.250 +Treatment + logged_1sqk.s, pformula = ~jcen + jcen2, data = bird)

model1 <- colext(psiformula = ~Treatment, gammaformula = ~ Treatment +logged_1sqk.s, 
                 epsilonformula = ~Treatment + logged_1sqk.s, pformula = ~jcen + jcen2, data = bird)

nd$P.Spatial.Weights <- predict(model1, type="col", newdata=nd)
nd$P.NoSpatial.Weights <- predict(model2, type="col", newdata=nd)
#P.Prop.Harvest <- predict(model4, type="psi", newdata=nd)

ndLL$Actual.Prop.Harvest <- ndLL$logged_1sqk.s*0.6506915

write.csv(ndLL, "0_data/processed/ndTH.csv")
dataTH <- read.csv("0_data/processed/ndTH.csv", header = T)
dataTAE <- read.csv("0_data/processed/ndTH.csv", header = T)

ggplot(dataTAE, aes(x = Actual.Prop.Harvest, y = MeanPred, group=Model, color=Model)) +
 geom_errorbar(aes(ymin=LCLPred, ymax=UCLPred), colour="gray", width=.05) +
  geom_line() +
  scale_fill_manual(values=c("#999999", "#E69F00"), 
                    name="Model",
                    breaks=c("Spatial", "NonSpatial"),
                    labels=c("Spatial", "Non-Spatial")) +
  geom_point() + 
  ylab("Colonization Probability") +
  xlab("Proportion Harvested (within 565m of point count site)")+theme_classic()

# Expected DETECTION over range Julian day
# julian
summary(recast.jcen)
summary(recast.jcen2)
max(recast.jcen)
jcen <- recast.jcen


E.det <- predict(model4, type="det", newdata=newDataDetvobs, appendData=TRUE)
head(E.det)
Det.jcen <- ggplot(E.det, aes(x = jcen, y = E.det$Predicted)) +
  geom_errorbar(aes(ymin = E.det$Predicted-SE,
                    ymax = E.det$Predicted+SE),
                width = 0) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Julian Date", y = "Predicted Detection") +
  my.theme
Det.jcen

plot(Predicted ~ VOb, E.det, type="l", ylim=c(0,1),
     xlab="visual obstruction (standardized)",
     ylab="Expected detection probability")

lines(lower ~ VObs, E.det, type="l", col=gray(0.5))
lines(upper ~ VObs, E.det, type="l", col=gray(0.5))

#test for a confounding effect between harvest amount and treatment 
install.packages("asbio")
boxplot(site.cov$PROP565.harvest ~ site.cov$Treatment)
library(asbio)
??asbio
TukeyHSD(aov(site.cov$PROP565.harvest ~ site.cov$Treatment))
m<-lm(site.cov$PROP565.harvest ~ site.cov$Treatment)
summary(m)
anova(m)
pairw.anova(y=prop565.harvest, x=site.cov$Treatment, method="tukey")
z <- pairw.anova(y=site.cov$PIPELINE_NEAR_DIST, x=site.cov$Treatment, method="tukey")
plot.pairw(z, type = 1, las=1)









Ext.pipeliendist <- ggplot(P.Ext.pipeline, 
              aes(x = P.Ext.pipeline$Pipeline_dist, y = P.Ext.pipeline$Predicted)) +
  geom_errorbar(aes(ymin = P.Ext.pipeline$Predicted-SE,
                    ymax = P.Ext.pipeline$Predicted+SE),
                width = 0) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Distance to Pipeline (m)", y = "Predicted Extinction") +
  my.theme
Ext.pipeliendist
ggsave(filename=paste0("2_outputs/",i,"/",i,"P.ExtPipeline.png"), plot= Ext.pipeliendist)
model4
#Gamma and pipeline dist
fm <- colext(~weights250, ~PIPELINE_NEAR_DIST + Treatment, ~PIPELINE_NEAR_DIST + Treatment, ~PIPELINE_NEAR_DIST, data = bird)
dist_pipeline <- read.csv("0_data/processed/site.covDistpipeline.csv")
#y.pipeline.dist <- site.cov$PIPELINE_NEAR_DIST
#treatment <- site.cov$Treatment
nd <- data.frame(dist_pipeline)
P.Psi.Pipeline <- predict(fm, type="psi", newdata=nd)
P.Psi.Pipeline$Pipeline_dist <- dist_pipeline$PIPELINE_NEAR_DIST


Psi.pipeliendist <- ggplot(P.Psi.Pipeline, 
                           aes(x = Pipeline_dist, y = P.Psi.Pipeline$Predicted)) +
  geom_errorbar(aes(ymin = P.Psi.Pipeline$Predicted-SE,
                    ymax = P.Psi.Pipeline$Predicted+SE),
                width = 0) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Distance to Pipeline (m)", y = "Predicted Occupancy") +
  my.theme
Psi.pipeliendist
ggsave(filename=paste0("2_outputs/",i,"/",i,"P.Psi.Pipeline.png"), plot= Psi.pipeliendist)

#Now do Distance to harvest
fm <- colext(~PROP565.harvest, ~PROP565.harvest, ~PROP565.harvest, ~PROP565.harvest, data = bird)
prop_harvest <- read.csv("0_data/processed/site.cov.harvest.csv")
prop_harvest$PROP565.harvest
nd <- data.frame(prop_harvest)
P.Psi.Harvest <- predict(fm, type="psi", newdata=nd)
#predict(fm, type="col", newdata=nd)
P.Ext.Harvest <- predict(fm, type="ext", newdata=nd)
#predict(fm, type="det", newdata=nd)
#Extinction and pipeline dist
plot(prop_harvest$PROP565.harvest, P.Psi.Harvest$Predicted)
P.Psi.Harvest$Harvest_dist <- dist_harvest$HARVEST_NEAR_DIST
#Psi and pipeline dist

Psi.harvestdist <- ggplot(P.Psi.Harvest, 
                           aes(x = Harvest_dist, y = P.Psi.Harvest$Predicted)) +
  geom_errorbar(aes(ymin = P.Psi.Harvest$Predicted-SE,
                    ymax = P.Psi.Harvest$Predicted+SE),
                width = 0) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Distance to Harvest (m)", y = "Predicted Occupancy") +
  my.theme
Psi.harvestdist
ggsave(filename=paste0("2_outputs/",i,"/",i,"P.Psi.Harvest.png"), plot= Psi.harvestdist)

# Now do distance to roads 
fm <- colext(~ROADS_NEAR_DIST, ~ROADS_NEAR_DIST, ~ROADS_NEAR_DIST, ~ROADS_NEAR_DIST, data = bird)
dist_road <- read.csv("0_data/processed/site.covDistroad.csv")
#y.pipeline.dist <- site.cov$PIPELINE_NEAR_DIST
#treatment <- site.cov$Treatment
nd <- data.frame(dist_road)
P.Psi.Road <- predict(fm, type="psi", newdata=nd)
#predict(fm, type="col", newdata=nd)
P.Ext.Road <- predict(fm, type="ext", newdata=nd)
#predict(fm, type="det", newdata=nd)

P.Psi.Road$Road_dist <- dist_road$ROADS_NEAR_DIST
#Psi and pipeline dist

Psi.roaddist <- ggplot(P.Psi.Road, 
                          aes(x = Road_dist, y = P.Psi.Road$Predicted)) +
  geom_errorbar(aes(ymin = P.Psi.Road$Predicted-SE,
                    ymax = P.Psi.Road$Predicted+SE),
                width = 0) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Distance to Road (m)", y = "Predicted Occupancy") +
  my.theme
Psi.roaddist
ggsave(filename=paste0("2_outputs/",i,"/",i,"P.Psi.Road.png"), plot= Psi.roaddist)

# Now do distance to seismic lines 
fm <- colext(~SEISMIC_NEAR_DIST, ~SEISMIC_NEAR_DIST, ~SEISMIC_NEAR_DIST, ~SEISMIC_NEAR_DIST, data = bird)
dist_seismic <- read.csv("0_data/processed/site.covDistseismic.csv")
#y.pipeline.dist <- site.cov$PIPELINE_NEAR_DIST
#treatment <- site.cov$Treatment
nd <- data.frame(dist_seismic)
P.Psi.Seismic <- predict(fm, type="psi", newdata=nd)
#predict(fm, type="col", newdata=nd)
P.Ext.Seismic <- predict(fm, type="ext", newdata=nd)
#predict(fm, type="det", newdata=nd)

P.Psi.Seismic$seismic_dist <- dist_seismic$SEISMIC_NEAR_DIST
#Psi and pipeline dist

Psi.seismicdist <- ggplot(P.Psi.Seismic, 
                       aes(x = seismic_dist, y = P.Psi.Seismic$Predicted)) +
  geom_errorbar(aes(ymin = P.Psi.Seismic$Predicted-SE,
                    ymax = P.Psi.Seismic$Predicted+SE),
                width = 0) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Distance to Seismic Lines (m)", y = "Predicted Occupancy") +
  my.theme
Psi.seismicdist
ggsave(filename=paste0("2_outputs/",i,"/",i,"P.Psi.Seismic.png"), plot= Psi.seismicdist)
#Estimates of occupancy probability in years T > 1 must be derived from the estimates
#of first-year occupancy and the two parameters governing the dynamics, extinction/survival
#and colonization. unmarked does this automatically in two ways. First, the population-level 
#estimates of occupancy probability  t =  t????1t????1 + (1 ???? t????1)
#are calculated and stored in the slot named projected. Slots can be accessed using the @
#operator, e.g. fm@projected. In some cases, interest may lie in making inference about the
#proportion of the sampled sites that are occupied, rather than the entire population of sites.
#These estimates are contained in the smoothed slot of the fitted model. Thus, the projected
#values are estimates of population parameters, and the smoothed estimates are of the
#finite-sample quantities. Discussions of the differences can be found in Weir et al. (2009).

#e.g.smoothed=smoothed(bestmodel.stage2.boot)[2,], SE=bestmodel.stage2.boot@smoothed.mean.bsse[2,]
#Note that the smoothed estimates of occupancy probability are overall calculations for all sites, not differing 
#among individual sites

#below is how you get occupancies for specific sites
psi.unlogged.c_1_1_1<-bestmodel.stage4.boot@projected[2,,1]


#Get average psi values across control sites 
psi.Control<-bestmodel.stage4.boot@projected[2,,1:93]
Psi.Control <- data.frame(psi.Control)
write.csv(Psi.Control, "0_data/processed/Predicted.Psi.Occupancy.csv")
Psi.Treatments.avgs <- read.csv("0_data/processed/Predicted.Psi.Occupancy.csv")
Psi.all.trt <- ggplot(Psi.Treatments.avgs, aes(x = Year, y = Predicted.Occupancy, color = Treatment)) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = "Year", y = "Smoothed occupancy (derived)") +
  my.theme
count(year, genus) %>%
  ggplot(mapping = aes(x = year, y = n, color = genus)) +
  geom_line()
ControlOccu
ggsave(filename=paste0("2_outputs/",i,"/",i,"_psicontrolstations.png"), plot= ControlOccu)

#Get average psi values across fragmented sites 
psi.Fragment<-bestmodel.stage4.boot@projected[2,,94:186]
psi.Fragment <- data.frame(psi.Fragment)
write.csv(psi.Fragment, "0_data/processed/Predicted.Psi.OccupancyFragments.csv")
Psi.fragment.avgs <- read.csv("0_data/processed/Predicted.Psi.OccupancyFragments.csv")
#Plot fragment sites occupancy
FragmentOccu <- ggplot(Psi.fragment.avgs, aes(x = Year, y = Average.Occupancy)) +
  #geom_errorbar(aes(ymin = smoothed_occ-SE,
  #   ymax = smoothed_occ+SE),
  # width = 0) +
  geom_point(size = 3, col = "darkgreen") +
  geom_line(col= "darkgreen") +
  labs(x = "Year", y = "Smoothed occupancy (derived)") +
  my.theme
FragmentOccu
ggsave(filename=paste0("2_outputs/",i,"/",i,"_psifragmentstations.png"), plot= FragmentOccu)

#Get average psi values across riparian sites 
psi.riparian<-bestmodel.stage4.boot@projected[2,,187:219]
psi.riparian <- data.frame(psi.riparian)
write.csv(psi.riparian, "0_data/processed/Predicted.Psi.OccupancyRiparian.csv")
Psi.riparian.avgs <- read.csv("0_data/processed/Predicted.Psi.OccupancyRiparian.csv")
RiparianOccu <- ggplot(Psi.riparian.avgs, aes(x = Year, y = Predicted.Occupancy)) +
  #geom_errorbar(aes(ymin = smoothed_occ-SE,
  #   ymax = smoothed_occ+SE),
  # width = 0) +
  geom_point(size = 3, col = "black") +
  geom_line(col= "black") +
  labs(x = "Year", y = "Smoothed occupancy (derived)") +
  my.theme
RiparianOccu
ggsave(filename=paste0("2_outputs/",i,"/",i,"_psiriparianstations.png"), plot= RiparianOccu)
psi.unlogged.c_1_3_1<-bestmodel.stage4.boot@projected[2,,3]
psi.unlogged.c_10_1_1<-bestmodel.stage4.boot@projected[2,,4]
psi.unlogged.c_10_1_2<-bestmodel.stage4.boot@projected[2,,5]
psi.unlogged.c_10_2_1<-bestmodel.stage4.boot@projected[2,,6]
psi.unlogged.c_10_2_2<-bestmodel.stage4.boot@projected[2,,7]
psi.unlogged.c_10_3_1<-bestmodel.stage4.boot@projected[2,,8]
psi.unlogged.c_10_3_2<-bestmodel.stage4.boot@projected[2,,9]

#get occupancy probabilities, SE, LCL, UCL for all sites
bestmodel.stage4.boot@projected[2,,]#change 25x219 to 5475x1

library(reshape2)
occ.df<-as.data.frame(bestmodel.stage4.boot@projected[2,,])
occ.stacked<-stack(occ.df)
P.occ<-occ.stacked$values
P.occ.SE<-rep((bestmodel.stage4.boot@projected.mean.bsse[2,]),219)
P.occ.LCL<-P.occ-1.96*P.occ.SE
P.occ.UCL<-P.occ+1.96*P.occ.SE
preds<-cbind(yearlycov.u, P.occ, P.occ.SE, P.occ.LCL, P.occ.UCL)

#At this point, you should now be able to use plot, ggplot2 or 
#matplot to visualize how predicted occupancy varies with  
#covariates

plot(preds$YYYY, preds$P.occ)
plot(preds$CutblockEdgeDistY, preds$P.occ)
plot(preds$PIPELINE_DIST, preds$P.occ) 
library(mgcv)
x <- gam(preds$P.occ~s(preds$PIPELINE_DIST), family = "betar")
plot(x)
summary(x)
install.packages("gratia")
library(ggeffects)
install.packages("effects")
library(effects)
library(gratia)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
plot_model(x, type = "eff")
draw(x)
occu_pipeline_pred_plot <- ggplot(preds, aes(PIPELINE_DIST , P.occ)) +
  geom_smooth(aes(ymin = P.occ.LCL, ymax = P.occ.UCL), alpha = 0.5, linetype = "dashed") +
  geom_path(size = 1) 
occu_pipeline_pred_plot

ggplot(data = preds) +
  geom_line(preds$CutblockEdgeDistY, aes(preds$YYYY, preds$P.occ), colour="steelblue") +
  geom_line(preds$PIPELINE_DIST, aes(preds$YYYY, preds$P.occ), colour="red") +
  my.theme
 
my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

library(ggplot2)#to be done
# Plot of mean predicted occupancy across sites for each year (primary period)
ggplot(data=Occupancy) +
  geom_line(aes(x= year, y= smooth.occ), colour="steelblue") +
  
  my.theme
 
Occu

ggsave(filename=paste0("2_outputs/",i,"/",i,"_psi.allstations.png"), plot=Occu)

tiff(paste0("2_outputs/",i,"/",i,"_psi.allstations.tiff"), units="in", width=12, height=8, res=300)


matplot#to be done

save(bestmodel.stage4.boot, 
     bestmodel.stage4,file=paste0("2_outputs/",i,"/",i,"_bestmodelandbootstraps.RData"))

########################


