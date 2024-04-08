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
site.cov<-read.csv("0_data/processed/site.cov.csv")
yearlycov<-read.csv("0_data/processed/Calling Lake/yearlycovariates.csv")

#Detection covariates - these were prepared in script 2
recast.julian<-read.csv("0_data/processed/8_rearranged.julian.csv",header=TRUE)
recast.jcen<-read.csv("0_data/processed/10_rearranged.jcen.csv",header=TRUE)
recast.jcen2<-read.csv("0_data/processed/12_rearranged.jcen2.csv",header=TRUE)
recast.time<-read.csv("0_data/processed/14_rearranged.time.csv",header=TRUE)
recast.tcen<-read.csv("0_data/processed/16_rearranged.tcen.csv",header=TRUE)
recast.tcen2<-read.csv("0_data/processed/18_rearranged.tcen2.csv",header=TRUE)


"'logLik.try-error'"<-function(object,...) structure(-Inf, nobs = 1L, df = 1L, class = "logLik")
                      
#purpose of this line is to assign a ridiculously small log-likelihood to a model if the model fails to run, so that
#the model can still be compared to other models within the model.sel function as part of automated model selection.

#Example: Single Species

data<-read.csv(paste0("0_data/processed/6_rearranged.",i,".csv"),header=TRUE) # this is bird data that was prepared in scripts 1 & 2
y<-data[,c(3:ncol(data))]#colums 2 through 101
S <- nrow(data) # number of sites  (219)
J <- 4 # number of secondary sampling occasions    (4 per year)
T <- (ncol(data)-2)/J # number of primary periods =(ncol(data)-2)/J   (25 years of observations over 26 years)

bird<-unmarkedMultFrame(y=y,
                        siteCovs=site.cov,
                        numPrimary=25,
                        yearlySiteCovs=yearly.cov,
                        obsCovs=list(Julian=recast.julian[,c(3:ncol(recast.julian))], 
                                     Time=recast.time[,c(3:ncol(recast.time))],
                                     jcen=recast.jcen[,c(3:ncol(recast.jcen))], 
                                     jcen2=recast.jcen2[,c(3:ncol(recast.jcen2))], 
                                     tcen=recast.tcen[,c(3:ncol(recast.tcen))], 
                                     tcen2=recast.tcen2[,c(3:ncol(recast.tcen2))]))  

#Warning messages:
#1: siteCovs contains characters. Converting them to factors. 
#2: yearlySiteCovs contains characters. Converting them to factors. 
#cannot use dredge function because there are missing values

#psiformula=probability of initial occupancy in Year 1
#gammaformula=probability of colonizing previously unoccupied location
#epsilonformula=probability of going extinct at previously occupied location
#pformula=probability of detection given occupancy of location

#Stage 1
#test for evidence of spatial autocorrelation and effect of stand age and %conifer 
null<-try(colext(psiformula = ~1, gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))    #null model
#autolog.4600<-try(colext(psiformula = ~weights4600, gammaformula = ~autolog.4600, epsilonformula = ~autolog.4600, pformula = ~1, data=bird))    #autocorrelation model
#test for effect of age/ percent conifer on initial occupancy and the interaction and quadratics of age 
age<-try(colext(psiformula = ~age.s, gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))    #stand age model 
age.2 <-try(colext(psiformula = ~ poly(age.s, 2), gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))
whitespruce<-try(colext(psiformula = ~whtespruce.s, gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))
whitespruce.2 <-try(colext(psiformula = ~ poly(whtespruce.s, 2), gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))
#Check interaction and additive age/conifer terms 
age.conifer<-try(colext(psiformula = ~whtespruce.s + age.s, gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))
ageXconifer<-try(colext(psiformula = ~whtespruce.s * age.s, gammaformula = ~1, epsilonformula = ~1, pformula = ~1, data=bird))
#aggregation hypothesis (Mattson et al. 2013)
#model1 <- colext(psiformula = ~weights4600, gammaformula = ~autolog.4600, 
                   # epsilonformula = ~autolog.4600, pformula = ~1, data = bird)
results.stage1 <- model.sel('Additiveage.conifer'= age.conifer, 'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, rank=AIC)           
bestmodel.stage1<-get.models(results.stage1, subset = 1)[[1]]

fms <- fitList('Additiveage.conifer'= age.conifer, 'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null)          
ms<-modSel(fms)
ms.coef<-coef(ms)
ms.SE<-SE(ms)
output<-as(ms,"data.frame")
write.csv(output, file=paste0("2_outputs/",i,".stage1modeltable.csv"))

#Call:
#  colext(psiformula = ~weights4600, gammaformula = ~autolog.4600, 
#         epsilonformula = ~autolog.4600, pformula = ~1, data = bird)

#Initial:
#  Estimate   SE     z P(>|z|)
#(Intercept)    -10.3 3.75 -2.75 0.00599
#weights4600     18.1 6.25  2.89 0.00385

#As the proportion of points within 4600 m with BTNW detections in year 1 increases,
#probability of initial occupancy by BTNW increases.

#Colonization:
#  Estimate    SE      z  P(>|z|)
#(Intercept)     -2.74 0.231 -11.84 2.36e-32
#autolog.4600     2.45 0.538   4.56 5.21e-06

#As the proportion of points within 4600 m with BTNW detections in year N increases,
#probability of occupancy by BTNW in year N at previously unoccupied stations increases.

#Extinction:
#  Estimate    SE     z  P(>|z|)
#(Intercept)    -0.876 0.265 -3.30 0.000964
#autolog.4600   -1.726 0.561 -3.08 0.002079

#As the proportion of points within 4600 m with BTNW detections in year N increases,
#probability of occupancy by BTNW being absent in year N at previously occupied stations decreases.

#Detection:
#  Estimate     SE     z  P(>|z|)
#-0.462 0.0261 -17.7 2.52e-70

#AIC: 16994.38 
#Stage 2- Colonization, exticntion and initial occupancy vary with treatments 
gameps.trt<-try(colext(psiformula =bestmodel.stage1@psiformula, 
                        gammaformula= update(bestmodel.stage1@gamformula, ~.+ Treatment),
                        epsilonformula= update(bestmodel.stage1@epsformula, ~.+ Treatment),
                        pformula= ~1, data=bird))
#colonization and extinction vary with treatment
eps.trt<-try(colext(psiformula= bestmodel.stage1@psiformula, 
                     gammaformula= bestmodel.stage1@gamformula,
                     epsilonformula= update(bestmodel.stage1@epsformula, ~.+ Treatment),
                    pformula= ~1, data=bird))
#initial site occupancy, colonization and extinction vary with treatment 
psigameps.trt<-try(colext(psiformula= update(bestmodel.stage1@psiformula, ~.+ Treatment), 
                     gammaformula= update(bestmodel.stage1@gamformula, ~.+ Treatment), 
                     epsilonformula= update(bestmodel.stage1@epsformula, ~.+ Treatment),
                     pformula= ~1, data=bird))

results.stage2 <- model.sel('GamEpsTrt'= gameps.trt, 'EpsTrt' = eps.trt, 'PsiGamEpsTrt' = psigameps.trt, 'Additiveage.conifer'= age.conifer, 'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, rank=AIC)           
bestmodel.stage2<-get.models(results.stage2, subset = 1)[[1]]
#plogis(2.3)
fms <- fitList('GamEpsTrt'= gameps.trt, 'EpsTrt' = eps.trt, 'PsiGamEpsTrt' = psigameps.trt, 'Additiveage.conifer'= age.conifer, 'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null)          
ms<-modSel(fms)
ms.coef<-coef(ms)
ms.SE<-SE(ms)
output<-as(ms,"data.frame")
write.csv(output, file=paste0("2_outputs/",i,".stage2modeltable.csv"))


#STAGE 3 - test if model is improved by using distances to different footprints as covariates
#HARVEST
harvest.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+harvest.dist.s), 
                        gammaformula= bestmodel.stage2@gamformula, 
                        epsilonformula= bestmodel.stage2@epsformula, 
                        pformula= ~1, data=bird))#AIC: 16992.7 
#harvest05.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+harvest.dist05), 
                         # gammaformula= bestmodel.stage2@gamformula, 
                         # epsilonformula= bestmodel.stage2@epsformula, 
                         # pformula= ~1, data=bird))#AIC: 16993.23
#harvest.i.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+harvest.dist.i), 
                         # gammaformula= bestmodel.stage2@gamformula, 
                         # epsilonformula= bestmodel.stage2@epsformula, 
                         # pformula= ~1, data=bird))#AIC: 16995.46
#Black-throated Green Warbler initial occupancy increases slightly
#with distance from harvest.
#A linear function of distance to harvest has the lowest AIC, lower
#than the best stage 1 model but not by much. You COULD use either
#the best stage 1 model or harvest.psi model.

#Try seeing how distance to harvest affects colonization and extinction
gam.harvest<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+harvest.dist.s), 
                        gammaformula= update(bestmodel.stage2@gamformula, ~.+harvest.dist.y.s), 
                        epsilonformula= bestmodel.stage2@epsformula, 
                        pformula= ~1, data=bird))#AIC: 16992.91
#probability of colonizing previously unoccupied sites increases 
#with distance from nearest harvest, but not significantly so

eps.harvest<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+harvest.dist.s), 
                        gammaformula= bestmodel.stage2@gamformula, 
                        epsilonformula= update(bestmodel.stage2@epsformula, ~.+harvest.dist.y.s), 
                        pformula= ~1, data=bird))

#ROAD AND PIPELINE
road.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+road.dist.s), 
                        gammaformula= bestmodel.stage2@gamformula, 
                        epsilonformula= bestmodel.stage2@epsformula, 
                        pformula= ~1, data=bird)) 

#Try seeing how distance to roads affects colonization and extinction
gam.road<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+road.dist.s), 
                        gammaformula= update(bestmodel.stage2@gamformula, ~.+y.road.dist+y.pipeline.dist), 
                        epsilonformula= bestmodel.stage2@epsformula, 
                        pformula= ~1, data=bird))

eps.road<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+road.dist.s), 
                        gammaformula= bestmodel.stage2@gamformula, 
                        epsilonformula= update(bestmodel.stage2@epsformula, ~.+y.road.dist+y.pipeline.dist), 
                        pformula= ~1, data=bird))

#PIPELINES
pipeline.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+pipeline.dist.s), 
                    gammaformula= bestmodel.stage2@gamformula, 
                    epsilonformula= bestmodel.stage2@epsformula, 
                     pformula= ~1, data=bird))#AIC: 16992.7 
#pipeline05.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+pipeline.dist05), 
                     #  gammaformula= bestmodel.stage2@gamformula, 
                      # epsilonformula= bestmodel.stage2@epsformula, 
                       #pformula= ~1, data=bird))#AIC: 16993.23
#pipeline.i.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+pipeline.dist.i), 
                    #   gammaformula= bestmodel.stage2@gamformula, 
                     #  epsilonformula= bestmodel.stage2@epsformula, 
                     #   pformula= ~1, data=bird))#AIC: 16995.46

#Try seeing how distance to pipelines affects colonization and extinction
gam.pipeline<-try(colext(psiformula =bestmodel.stage2@psiformula, 
                   gammaformula= update(bestmodel.stage2@gamformula, ~.+y.pipeline.dist), 
                    epsilonformula= bestmodel.stage2@epsformula, 
                    pformula= ~1, data=bird))
eps.pipeline<-try(colext(psiformula =bestmodel.stage2@psiformula, 
                    gammaformula= bestmodel.stage2@gamformula, 
                    epsilonformula= update(bestmodel.stage2@epsformula, ~.+y.pipeline.dist), 
                    pformula= ~1, data=bird))
#SEISMIC LINES
seismic.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+seismic.dist.s), 
                     gammaformula= bestmodel.stage2@gamformula, 
                     epsilonformula= bestmodel.stage2@epsformula, 
                     pformula= ~1, data=bird))#AIC: 16992.7 
#seismic05.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+seismic.dist05), 
                      # gammaformula= bestmodel.stage2@gamformula, 
                      # epsilonformula= bestmodel.stage2@epsformula, 
                      # pformula= ~1, data=bird))#AIC: 16993.23
#seismic.i.psi<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+seismic.dist.i), 
                      # gammaformula= bestmodel.stage2@gamformula, 
                      # epsilonformula= bestmodel.stage2@epsformula, 
                      # pformula= ~1, data=bird))#AIC: 16995.46

#Try seeing how distance to seismic lines affects colonization and extinction
gam.seismic<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+seismic.dist.s), 
                     gammaformula= update(bestmodel.stage2@gamformula, ~.+y.seismic.dist), 
                     epsilonformula= bestmodel.stage2@epsformula, 
                     pformula= ~1, data=bird))
eps.seismic<-try(colext(psiformula =update(bestmodel.stage2@psiformula, ~.+seismic.dist.s), 
                     gammaformula= bestmodel.stage2@gamformula, 
                     epsilonformula= update(bestmodel.stage2@epsformula, ~.+y.seismic.dist), 
                     pformula= ~1, data=bird))


results.stage3 <- model.sel('OccuHarvestdistS'=harvest.psi, 'OccuRoaddistS'=road.psi, 'OccuseismicdistS'=seismic.psi, 'OccuPipelinedist' = pipeline.psi,
                            'ColHarvest'=gam.harvest, 'ColRoad'=gam.road, 'ColPipeline' = gam.pipeline, 'ColSeismic'=gam.seismic, 'ExtHarvest' = eps.harvest, 'ExtPipeline' = eps.pipeline, 'ExtRoad'=eps.road, 'ExtSeismic'=eps.seismic, 'GamEpsTrt'= gameps.trt, 'Additiveage.conifer'= age.conifer, 'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, rank=AIC)       
bestmodel.stage3<-get.models(results.stage3, subset = 1)[[1]]
#Initial:
#Estimate    SE     z P(>|z|)
#(Intercept)       -10.09 3.938 -2.56 0.01037
#weights4600        18.24 6.614  2.76 0.00582
#pipeline.dist.s    -1.15 0.866 -1.33 0.18474

#Colonization:
 # Estimate    SE      z  P(>|z|)
#(Intercept)     -2.71 0.229 -11.86 1.83e-32
#autolog.4600     2.47 0.531   4.65 3.32e-06

#Extinction:
 # Estimate    SE     z  P(>|z|)
#(Intercept)        -1.15 0.274 -4.19 2.76e-05
#autolog.4600       -1.93 0.574 -3.36 7.68e-04
#pipeline.dist.s     1.49 0.316  4.72 2.34e-06

#Detection:
 # Estimate     SE     z  P(>|z|)
#-0.454 0.0261 -17.4 5.37e-68

#AIC: 16974.86 
#There is strong spatial autocorrelation, so the proportion ofneighbouring sites within 
#the buffer that had btnw detections is a strong predictor of btnw occupancy at focal site
# AND extinction risk is positively related to pipeline distance- occupied points are less
#likely to lose warblers if they are close to a pipeline 

fms <- fitList('OccuHarvestdistS'=harvest.psi, 'OccuRoaddistS'=road.psi, 'OccuseismicdistS'=seismic.psi,
               'ColHarvest'=gam.harvest, 'ColRoad'=gam.road, 'ColSeismic'=gam.seismic, 'ExtHarvest' = eps.harvest, 'ExtRoad'=eps.road, 'ExtSeismic'=eps.seismic, 'GamEpsTrt'= gameps.trt, 'Additiveage.conifer'= age.conifer, 'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, 'Autologistic'=autolog.4600)
ms<-modSel(fms)
ms.coef<-coef(ms)
ms.SE<-SE(ms)
output<-as(ms,"data.frame")
write.csv(output, file=paste0("2_outputs/",i,"/",i,".stage3modeltable.csv"))
#probability of extinction at previously occupied sites decreases 
#with distance from nearest harvest

#Stage 4 - detection covariates
det.julian<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~Julian, data=bird))     
det.julian.sq<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~jcen+jcen2, data=bird))     
det.time<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~Time, data=bird))     
det.time.sq<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~tcen+tcen2, data=bird))     
det.juliantime<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~jcen+jcen2+Time, data=bird))     
det.juliantime.sq<-try(colext(psiformula = bestmodel.stage3@psiformula, gammaformula = bestmodel.stage3@gamformula, epsilonformula = bestmodel.stage3@epsformula, pformula = ~Julian+tcen+tcen2, data=bird))     

results.stage4 <- model.sel('OccuHarvestdistS'=harvest.psi, 'OccuRoaddistS'=road.psi, 'OccuseismicdistS'=seismic.psi,
                            'ColHarvest'=gam.harvest, 'ColRoad'=gam.road, 'ColSeismic'=gam.seismic, 'ExtHarvest' = eps.harvest, 'ExtRoad'=eps.road, 'ExtSeismic'=eps.seismic, 'GamEpsTrt'= gameps.trt, 'Additiveage.conifer'= age.conifer, 'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, 'Julian' = det.julian, 'Time' = det.time, 'Julian.sq' = det.julian.sq, 'Time.sq' = det.time.sq, 'Julian+Time' = det.juliantime, rank=AIC)           
bestmodel.stage4<-get.models(results.stage4, subset = 1)[[1]] #'Julian.sq+Time.sq' = det.juliantime.sq, 

fms <- fitList('OccuHarvestdistS'=harvest.psi, 'OccuRoaddistS'=road.psi, 'OccuseismicdistS'=seismic.psi,
               'ColHarvest'=gam.harvest, 'ColRoad'=gam.road, 'ColSeismic'=gam.seismic, 'ExtHarvest' = eps.harvest, 'ExtRoad'=eps.road, 'ExtSeismic'=eps.seismic, 'GamEpsTrt'= gameps.trt, 'Additiveage.conifer'= age.conifer, 'Interactionage.conifer'= ageXconifer, 'Age'=age, 'Age2'=age.2, 'WSpruce2'=whitespruce.2, 'WSpruce' = whitespruce, 'Null'=null, 'Julian' = det.julian, 'Time' = det.time, 'Julian.sq' = det.julian.sq, 'Time.sq' = det.time.sq, 'Julian+Time' = det.juliantime)           
ms<-modSel(fms)#, 'Julian.sq+Time.sq' = det.juliantime.sq 
ms.coef<-coef(ms)
ms.SE<-SE(ms)
output<-as(ms,"data.frame")
write.csv(output, file=paste0("2_outputs/",i,"/",i,".stage4modeltable.csv"))

model4 <- colext(psiformula = ~weights4600, gammaformula = ~autolog.4600 + Treatment, 
                 epsilonformula = ~autolog.4600 + y.pipeline.dist + Treatment, pformula = ~jcen + jcen2, data = bird)


#test for a confounding effect between harvest amount and treatment 
install.packages("asbio")
boxplot(site.cov$PIPELINE_NEAR_DIST ~ site.cov$Treatment)
library(asbio)
??asbio
TukeyHSD(aov(site.cov$PIPELINE_NEAR_DIST ~ site.cov$Treatment))
m<-lm(site.cov$PIPELINE_NEAR_DIST ~ site.cov$Treatment)
summary(m)
anova(m)
pairw.anova(y=PIPELINE_NEAR_DIST, x=site.cov$Treatment, method="tukey")
z <- pairw.anova(y=site.cov$PIPELINE_NEAR_DIST, x=site.cov$Treatment, method="tukey")
plot.pairw(z, type = 1, las=1)
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

bestmodel.stage4.boot@projected[2,,1]
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
#plot julian day 
model4 <- colext(psiformula = ~weights4600 + road.dist.s, gammaformula = ~autolog.4600 + Treatment, 
                 epsilonformula = ~autolog.4600 + y.pipeline.dist + Treatment + y.road.dist, pformula = ~I(jcen)+I(jcen^2), data = bird)

#Initial:
#Estimate   SE     z P(>|z|)
#(Intercept)    -9.04 4.35 -2.08  0.0374
#weights4600    16.07 7.22  2.23  0.0260

#Colonization:
 # Estimate    SE     z  P(>|z|)
#(Intercept)               -2.193 0.254 -8.62 6.42e-18
#autolog.4600               1.892 0.546  3.47 5.29e-04
#TreatmentFragment         -0.368 0.154 -2.38 1.71e-02
#TreatmentRiparian Strip   -0.721 0.225 -3.21 1.33e-03

#Extinction:
#  Estimate    SE      z  P(>|z|)
#(Intercept)               -2.114 0.327 -6.454 1.09e-10
#autolog.4600              -0.144 0.613 -0.235 8.14e-01
#y.pipeline.dist            5.334 0.527 10.126 4.25e-24
#TreatmentFragment         -1.477 0.214 -6.897 5.31e-12
#TreatmentRiparian Strip   -1.605 0.316 -5.073 3.92e-07

#Detection:
 # Estimate       SE     z  P(>|z|)
#(Intercept) -0.31161 0.032002 -9.74 2.10e-22
#jcen        -0.00820 0.002186 -3.75 1.76e-04
#jcen2       -0.00123 0.000152 -8.09 6.14e-16

#AIC: 16801.62 
#nd <- data.frame(dist_pipeline_cov)
#P.Psi.Julian <- predict(fm, type="det", newdata=nd)

#Makin gplots 
library(ggplot2)
# Expected DETECTION over range Julian day
# julian
summary(recast.jcen)
summary(recast.jcen2)
max(recast.jcen)
jcen <- recast.jcen
jcen<-select(recast.jcen, -SS)
jcen<-as.data.frame(recast.jcen)
jcen.stacked<- stack(jcen)
jcen.stacked<-na.omit(jcen.stacked)
range(jcen.stacked$values)
jcen2<-select(recast.jcen2, -SS)
jcen2.stacked<-stack(jcen2)
jcen2.stacked<-na.omit(jcen2.stacked)
range(jcen2.stacked$values)


nd <- data.frame(weights4600=c(rep(0.5,12)),
                 autolog.4600=c(rep(0.5,12)),
                 Treatment=c(rep("Control",12)),
                 y.pipeline.dist=c(rep(0.5,12)),
                 jcen=seq(-24.644, 26.355, by = 4.6))
Actual.Juliandate<-nd.preds$jcen+160.644                 
P.Det.Julian <- predict(model4, type="det", newdata=nd)
#Actual.Juliandate<-nd$jcen *26.35574
nd.preds<-cbind(nd, P.Det.Julian,Actual.Juliandate)
#Actual.Juliandate2<-nd$jcen2*694.6247
ggplot(nd.preds, aes(x = Actual.Juliandate, y = Predicted)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1) +
  geom_line() +
  geom_point() + 
  # ggtitle("Extinction Probability vs. Pipeline Distance") +
  ylab("Detection Probability") +
  xlab("Julian Date")+my.theme

#Making awesome ggplot of extinction probability and distance to pipeline 
model4 <- colext(psiformula = ~weights4600, gammaformula = ~autolog.4600 + Treatment, 
                 epsilonformula = ~autolog.4600 + y.pipeline.dist + Treatment, pformula = ~jcen + jcen2, data = bird)


nd <- data.frame(weights4600=c(rep(0.5,11)),
                 autolog.4600=c(rep(0.5,11)),
                 Treatment=c(rep("Control",11)),
                 y.pipeline.dist=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))#data.frame(dist_pipeline_cov)
P.Ext.PipelineDist <- predict(model4, type="ext", newdata=nd)

Actual.PipelineDist<-nd$y.pipeline.dist*2853.62749
nd.preds<-cbind(nd, P.Ext.PipelineDist,Actual.PipelineDist)

ggplot(nd.preds, aes(x = Actual.PipelineDist, y = Predicted)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1) +
  geom_line() +
  geom_point() + 
 # ggtitle("Extinction Probability vs. Pipeline Distance") +
  ylab("Extinction Probability") +
  xlab("Distance to Pipeline (m)")+my.theme

#Distance to road and exticntion probability plot 
nd <- data.frame(weights4600=c(rep(0.5,11)),
                 autolog.4600=c(rep(0.5,11)),
                 Treatment=c(rep("Control",11)),
                 road.dist.s=c(rep(0.5,11)),
                 y.road.dist=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                 y.pipeline.dist=c(rep(0.5,11))) #data.frame(dist_pipeline_cov)
P.Ext.RoadDist <- predict(model4, type="ext", newdata=nd)
range(site.cov$ROADS_NEAR_DIST)
Actual.RoadDist<-nd$y.road.dist*3542.76851
nd.preds<-cbind(nd, P.Ext.RoadDist,Actual.RoadDist)

ggplot(nd.preds, aes(x = Actual.RoadDist, y = Predicted)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1) +
  geom_line() +
  geom_point() + 
  # ggtitle("Extinction Probability vs. Road Distance") +
  ylab("Extinction Probability") +
  xlab("Distance to Road (m)")+my.theme

# Distance to road Psi plot
nd <- data.frame(weights4600=c(rep(0.5,11)),
                 autolog.4600=c(rep(0.5,11)),
                 Treatment=c(rep("Control",11)),
                 road.dist.s=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                 y.road.dist=c(rep(0.5,11)),
                 y.pipeline.dist=c(rep(0.5,11))) #data.frame(dist_pipeline_cov)
P.Psi.RoadDist <- predict(model4, type="psi", newdata=nd)
range(site.cov$ROADS_NEAR_DIST)
Actual.RoadDist<-nd$road.dist.s*3542.76851
nd.preds<-cbind(nd, P.Psi.RoadDist,Actual.RoadDist)

ggplot(nd.preds, aes(x = Actual.RoadDist, y = Predicted)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=.1) +
  geom_line() +
  geom_point() + 
  # ggtitle("Extinction Probability vs. Road Distance") +
  ylab("Initial Occupancy Probability") +
  xlab("Distance to Road (m)")+my.theme

#predict(fm, type="col", newdata=nd)
P.Ext.pipeline <- predict(fm, type="ext", newdata=nd)
#Plot distance to pipeline and exticntion probability
write.csv(site.cov, "0_data/processed/site.cov.csv")
dist_pipeline_cov <- read.csv("0_data/processed/site.covDisttopipe.csv")

fm <- colext(~PIPELINE_NEAR_DIST, ~PIPELINE_NEAR_DIST, ~PIPELINE_NEAR_DIST, ~PIPELINE_NEAR_DIST, data = bird)

nd <- data.frame(dist_pipeline_cov)
P.Psi.Pipeline <- predict(fm, type="psi", newdata=nd)
#predict(fm, type="col", newdata=nd)
P.Ext.pipeline <- predict(fm, type="ext", newdata=nd)
#predict(fm, type="det", newdata=nd)
#Extinction and pipeline dist
plot(dist_pipeline_cov$PIPELINE_NEAR_DIST, P.Ext.pipeline$Predicted)
P.Ext.pipeline$Pipeline_dist <- dist_pipeline_cov$PIPELINE_NEAR_DIST
#Psi and pipeline dist
plogis(-0.31161)
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
fm <- colext(~weights4600, ~PIPELINE_NEAR_DIST + Treatment, ~PIPELINE_NEAR_DIST + Treatment, ~PIPELINE_NEAR_DIST, data = bird)
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
fm <- colext(~HARVEST_NEAR_DIST, ~HARVEST_NEAR_DIST, ~HARVEST_NEAR_DIST, ~HARVEST_NEAR_DIST, data = bird)
dist_harvest <- read.csv("0_data/processed/site.covDistharvest.csv")
#y.pipeline.dist <- site.cov$PIPELINE_NEAR_DIST
#treatment <- site.cov$Treatment
nd <- data.frame(dist_harvest)
P.Psi.Harvest <- predict(fm, type="psi", newdata=nd)
#predict(fm, type="col", newdata=nd)
P.Ext.Harvest <- predict(fm, type="ext", newdata=nd)
#predict(fm, type="det", newdata=nd)
#Extinction and pipeline dist
plot(dist_harvest$HARVEST_NEAR_DIST, P.Psi.Harvest$Predicted)
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

#smoothed probability of occupancy over time
smooth.df<-cbind(Year=seq(1:25), MeanOccProb=smoothed(bestmodel.stage4.boot)[2,], MeanOccSE=bestmodel.stage4.boot@projected.mean.bsse[2,])
smooth.df<-data.frame(smooth.df)

ggplot(smooth.df, aes(x = Year, y = MeanOccProb)) +
  geom_errorbar(aes(ymin=MeanOccProb-MeanOccSE, ymax=MeanOccProb+MeanOccSE), colour="black", width=.1) +
  geom_line() +
  geom_point() + 
  ggtitle("Occupancy Probability Over Time") +
  ylab("Average Occupancy Probability") +
  xlab("Year")+my.theme

#Occupancy probability over time differences between treatments with SE bars 
Psi.Treatments.avgs <- read.csv("0_data/processed/Predicted.Psi.OccupancyTreatments.csv")
Psi.all.trt <- ggplot(Psi.Treatments.avgs, aes(x = Year, y = Occ.Fr, color = Treatment)) +
  geom_errorbar(aes(ymin=Occ.Fr-SE, ymax=Occ.Fr+SE)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", y = "Average Occupancy Probability") +
  my.theme
count(year, genus) %>%
  ggplot(mapping = aes(x = year, y = n, color = genus)) +
  geom_line()
ggplot.df.long<-matplot.df%>%
  pivot_longer(cols=OccProb.V1:OccProb.V219, 
               names_to="SiteID",
               values_to="OccProb")
str(matplot.df.long)
matplot.df.long$Treatment<-ifelse(matplot.df.long$SiteID %in% colnames(matplot.df[2:94]), "Control",
                                  ifelse(matplot.df.long$SiteID %in% colnames(matplot.df[95:187]),"Fragment","Riparian"))

ggplot(matplot.df.long, aes(x = Year, y = OccProb, group = SiteID, color = Treatment)) +
  geom_line(lwd = 0.5, show.legend = TRUE, alpha=0.25) + 
  #scale_color_manual(values = country_colors) +
  #theme_bw() + theme(strip.text = element_text(size = rel(1.1))) +
  ggtitle("Occupancy Probability Over Time") +
  ylab("Occupancy Probability") +
  xlab("Year")+my.theme

#matrix plot
matplot.df<-cbind(Year=c(1993,1994,1995,1996,1997,1998,
                         1999,2000,2001,2002,2003,2005,
                         2006,2007,2008,2009,2010,2011,
                         2012,2013,2014,2015,2016,2017,
                         2018), OccProb=as.data.frame(bestmodel.stage4.boot@projected[2,,]))
matplot.df<-data.frame(matplot.df)
library(graphics)
matplot(matplot.df$Year, 
        y=as.matrix(matplot.df[2:220]), 
        type = "l", 
        xlab="Year",
        ylab="Predicted Occupancy")

#colnames(matplot.df[2:94])<-paste0("Control",colnames(matplot.df[2:94]))
library(tidyr)
library(dplyr)
matplot.df.long<-matplot.df%>%
  pivot_longer(cols=OccProb.V1:OccProb.V219, 
               names_to="SiteID",
               values_to="OccProb")
str(matplot.df.long)
matplot.df.long$Treatment<-ifelse(matplot.df.long$SiteID %in% colnames(matplot.df[2:94]), "Control",
                                  ifelse(matplot.df.long$SiteID %in% colnames(matplot.df[95:187]),"Fragment","Riparian"))

ggplot(matplot.df.long, aes(x = Year, y = OccProb, group = SiteID, color = Treatment)) +
  geom_line(lwd = 0.5, show.legend = TRUE, alpha=0.25) + 
  #scale_color_manual(values = country_colors) +
  #theme_bw() + theme(strip.text = element_text(size = rel(1.1))) +
  ggtitle("Occupancy Probability Over Time") +
  ylab("Occupancy Probability") +
  xlab("Year")+my.theme


 
my.theme <- theme_classic() +
  theme(text=element_text(size=16, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

save(bestmodel.stage4.boot, 
     bestmodel.stage4,file=paste0("2_outputs/",i,"/",i,"_bestmodelandbootstraps.RData"))

########################
#Test for correlation between pipeline and other linear features 
cor(site.cov$PIPELINE_NEAR_DIST, site.cov$ROADS_NEAR_DIST)
# 0.6484541
cor(site.cov$PIPELINE_NEAR_DIST, site.cov$HARVEST_NEAR_DIST)
# -0.2361143
cor(site.cov$PIPELINE_NEAR_DIST, site.cov$SEISMIC_NEAR_DIST)
# 0.07366551
#Test for correlation between harvest and other linear features 
cor(site.cov$HARVEST_NEAR_DIST, site.cov$ROADS_NEAR_DIST)
# -0.2858552
cor(site.cov$HARVEST_NEAR_DIST, site.cov$SEISMIC_NEAR_DIST)
# -0.03081123
cor(site.cov$SEISMIC_NEAR_DIST, site.cov$ROADS_NEAR_DIST)
# 0.06715956





