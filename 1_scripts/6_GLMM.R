
# set session working directory

# reading in the data
dat <- read.csv("0_data/raw/CLpoints_covariates.csv", header = T)

# getting new footrpint data and formatting it for analyses
#1993-2014
footprint <- read.csv("0_data/processed/all-footprints-areadist.age-allsurveys-justCallingLake-1993-2014.csv")
footprint_u <- footprint %>% distinct(SS, YEAR, .keep_all = TRUE)

#2015-2018
footprint2015 <- read.csv("0_data/processed/all-footprints-areadist.age-CL.2015-2018.csv")
footprint2015_u <- footprint2015 %>% distinct(SS, YEAR, .keep_all = TRUE)

# merge them together 
footprint <- rbind(footprint2015, footprint)
# change column name YEAR to YYYY so it matches with dat 
footprint <- rename(footprint, YYYY = YEAR)
#merge dat with footprint to get BTNW presence, treatment, pacth size and new footprint covariates together 
data <- merge(footprint, dat, by = c("SS", "YYYY"))
#5529, still too many rows, delete year 2011
#footprint <- footprint %>% filter(!(YEAR=="2011"))

# factor variables - random effects
year <- as.factor(data$YYYY)
site <- as.factor(data$SITE)
cluster <- as.factor(data$SS)

# explanatory variables - categorical
# abundance/occupancy 
trt <- as.factor(data$Treatment)
area <- as.factor(data$Patch.Size..ha.)

# explanatory variables - covariates
# abundance/occupancy
harv <- data$HARVEST_NEAR_DIST
pipe <- data$PIPELINE_NEAR_DIST
road <- data$ROADS_NEAR_DIST
seis <- data$SEISMIC_NEAR_DIST

prop_harv <- data$PROP150.harvest
prop_pipe <- data$PROP150.pipeline
prop_seismic <- data$PROP150.conventional.seismic
prop_road <- data$PROP565.gravel.road
prop_seismic565 <- data$PROP565.conventional.seismic

# detection
jul <- data$JULIAN
sun <- data$MIN_SUN

# spatial
lat <- data$latitude
long <- data$longitude

# response variables
abu <- data$BTNW
occ <- ifelse(abu > 0, 1, 0)

# modelling
library(glmmTMB)
install.packages("glmmTMB")
install.packages('TMB', type = 'source')
library(TMB)
library(mgcv)
library(gratia)
library(DHARMa)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(gstat)
library(raster)
library(sp)
library(dplyr)

# glmmTMB with distance to pipeline as predictor 
mod_glmm_pipe <- glmmTMB(occ ~ trt * pipe +  
                      + I(jul) + I(jul^2) +
                      (1 | year) + 
                      (1 | site / cluster),
                    family = binomial, REML = T)
summary(mod_glmm_pipe)
res <- residuals(mod_glmm_pipe)

plot(gam(res~s(lat), family = "gaussian"), rug = T, shade=TRUE, seWithMean=TRUE)
plot(gam(res~s(long), family = "gaussian"), rug = T, shade=TRUE, seWithMean=TRUE)

plot_model(mod_glmm_pipe, type = "eff")
plot_model(mod_glmm_pipe, type = "int")
plot_model(mod_glmm_pipe, type = "re")
boxplot(res~year)
lines(lowess(res~year), col="red")

# looking at spatial autocorrelation
res <-  residuals(mod_glmm_pipe, type="pearson")
?coordinates
# with package gstat
E <- res
mydata <- data.frame(E, long, lat)
mydata <- na.omit(mydata)
coordinates(mydata) <- c("long", "lat")
bubble(mydata, "E", main="Model residuals", xlab="Long.", ylab="Lat.", col=c("1","8"))
# to reestablsh color values to default, col=c("1","8"))

# semivariogram
vario <- variogram(E~1,mydata)
plot(vario, main="Semivariogram", xlab="Distance", ylab="Residuals semivariance")
# HIGH SPATIAL AUTOCORRELATIO

# para esta funcion lat y long van separadas, y K lo mimso, objeto con los residuos
install.packages("ncf")
library(ncf)
corrspline <- spline.correlog(long, lat, E, resamp = 99, type = "boot")
plot(corrspline)

# glmmTMB with proportion pipeline as predictor 
mod_glmm_pipep <- glmmTMB(occ ~ trt * prop_pipe +  
                           + I(jul) + I(jul^2) +
                           (1 | year) + 
                           (1 | site / cluster), data = data,
                         family = binomial, REML = T)
summary(mod_glmm_pipep)
res <- residuals(mod_glmm_pipep)
plot(gam(res~s(lat), family = "gaussian"), rug = T, shade=TRUE, seWithMean=TRUE)
plot(gam(res~s(long), family = "gaussian"), rug = T, shade=TRUE, seWithMean=TRUE)

plot_model(mod_glmm_pipep, type = "eff")
plot_model(mod_glmm_pipep, type = "int")
plot_model(mod_glmm_pipep, type = "re")
boxplot(res~year)
lines(lowess(res~year), col="red")

# glmmTMB with proportion seismic line in 565m buffer as predictor 
mod_glmm_seismicp565 <- glmmTMB(occ ~ trt * prop_seismic565 +  
                            + I(jul) + I(jul^2) +
                            (1 | year) + 
                            (1 | site / cluster), data = data,
                          family = binomial, REML = T)
summary(mod_glmm_seismicp565)

plot_model(mod_glmm_seismicp565, type = "eff")
plot_model(mod_glmm_seismicp565, type = "int")
plot_model(mod_glmm_pipep565, type = "re")
boxplot(res~year)
lines(lowess(res~year), col="red")

# glmmTMB with seismic line as predictor 
mod_glmm_seis <- glmmTMB(occ ~ trt * seis +  
                           + I(jul) + I(jul^2) +
                           (1 | year) + 
                           (1 | site / cluster),
                         family = binomial, REML = T)
summary(mod_glmm_seis)
res <- residuals(mod_glmm_seis)
plot(gam(res~s(lat), family = "gaussian"), rug = T, shade=TRUE, seWithMean=TRUE)
plot(gam(res~s(long), family = "gaussian"), rug = T, shade=TRUE, seWithMean=TRUE)

plot_model(mod_glmm_seis, type = "eff")
plot_model(mod_glmm_seis, type = "int")
plot_model(mod_glmm_seis, type = "re")
boxplot(res~year)
lines(lowess(res~year), col="red")

# glmmTMB with distance to harvest as predictor 
mod_glmm_harv <- glmmTMB(occ ~ trt * harv +  
                           + I(jul) + I(jul^2) +
                           (1 | year) + 
                           (1 | site / cluster),
                         family = binomial, REML = T)
summary(mod_glmm_harv)
res <- residuals(mod_glmm_harv)
plot(gam(res~s(lat), family = "gaussian"), rug = T, shade=TRUE, seWithMean=TRUE)
plot(gam(res~s(long), family = "gaussian"), rug = T, shade=TRUE, seWithMean=TRUE)

plot_model(mod_glmm_harv, type = "eff")
plot_model(mod_glmm_harv, type = "int")
plot_model(mod_glmm_harv, type = "re")

# glmmTMB with proportion harvest in 150m as predictor 
mod_glmm_pharv <- glmmTMB(occ ~ trt * prop_harv +  
                           + I(jul) + I(jul^2) +
                           (1 | year) + 
                           (1 | site / cluster),
                         family = binomial, REML = T)
summary(mod_glmm_pharv)
res <- residuals(mod_glmm_pharv)
plot(gam(res~s(lat), family = "gaussian"), rug = T, shade=TRUE, seWithMean=TRUE)
plot(gam(res~s(long), family = "gaussian"), rug = T, shade=TRUE, seWithMean=TRUE)

plot_model(mod_glmm_pharv, type = "eff")
plot_model(mod_glmm_pharv, type = "int")
plot_model(mod_glmm_pharv, type = "re")

# glmmTMB with road as predictor 
mod_glmm_road <- glmmTMB(occ ~ trt * road +  
                           + I(jul) + I(jul^2) +
                           (1 | year) + 
                           (1 | site / cluster),
                         family = binomial, REML = T)
summary(mod_glmm_road)

plot_model(mod_glmm_road, type = "eff")
plot_model(mod_glmm_road, type = "int")
plot_model(mod_glmm_road, type = "re")

# glmmTMB with proportion road in 150m buffer as predictor 
mod_glmm_proad <- glmmTMB(occ ~ trt * prop_road +  
                           I(jul) + I(jul^2) +
                           (1 | year) + 
                           (1 | site / cluster),
                         family = binomial, REML = T)
summary(mod_glmm_proad)

plot_model(mod_glmm_proad, type = "eff")
plot_model(mod_glmm_proad, type = "int")

# interaction model trt * seis


plot(gam(res~s(long), family = "gaussian"))
plot(mod)

mod <- gam(occ ~ trt * seis + 
             s(jul, k=3) + 
             s(year,bs="re") +
             s(lat),
           family = "binomial")
plot_model(mod, type = "eff")




res <- residuals.gam(mod)
plot(gam(res ~ s(lat), family = "gaussian"))



#draw(mod)
summary(mod)
# draw(mod)

plot_model(mod, transform = NULL, type = "re")

sim <- simulateResiduals(mod)

plot(simulateResiduals(mod))
summary(mod)

mod

rec <- recalculateResiduals(mod)
library(mgcViz)
testSpatialAutocorrelation(rec, x = long, y = lat,
                           distMat = NULL, 
                           plot = T)


## First specify the packages of interest
packages = c("rgdal","spatstat","pgirmess","ncf","spdep", "geoR","gstat","raster", "waveslim","fields","vegan","reshape2","adespatial","ggplot2")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
