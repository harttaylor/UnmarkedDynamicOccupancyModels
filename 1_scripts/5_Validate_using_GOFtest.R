library(AICcmodavg) # For GOF tests

# Mackenzie-Bailey GOF test
# Simulate capture history data (if model correct). Compare obs X2 to sim X2
# Must simulate 1000-5000 times for good distribution estimates
# Likely to take a Very Long Time on large data sets
# with nsim = 10, takes my Mackbook ~ 60 seconds
mb.boot <- AICcmodavg::mb.gof.test(model4, nsim = 50) # Must be much higher than five to be useful
print(mb.boot, digit.vals = 4, digits.chisq = 4)

#get stack of predicted occuancies and actual occupancies stack of same lenggth
#function AUC
##model summary
summary(bestmodel.stage4.boot)

## ----sitepredictions-------------------
# Get the smoothed occupancy probabilities for just site 25 (at year 1:5) 
# [2, , 25] = 2 for occupied, " " for all 5 surveys, 25 for site 25; remove [] to print whole array
bestmodel.stage4.boot@smoothed[2, , 25]

## ----smoothedmeanpredictions---------------------------------------------
# Mean smoothed occupancy probabilities can be accessed using:
# dynamic_occ_m1@smoothed.mean, or
# smoothed(dynamic_occ_m1)

# Calculate SE for derived occupancy predictions using bootstrap, larger B requires longer time
m1 <- nonparboot(bestmodel.stage4.boot, 
                 B = 10)

# Predicted occupancy in each year (with SE)
# the "[2,]" calls the occupied estimates,
# "[1,]" for unoccupied estimates
predicted_occupancy <- data.frame(year = c(1:25),
                                  smoothed_occ = smoothed(bestmodel.stage4.boot)[2,],
                                  SE = m1@smoothed.mean.bsse[2,])

## ----plotderivedoccupancy----
occu <- ggplot(predicted_occupancy, 
       aes(x = year, y = smoothed_occ)) +
  geom_errorbar(aes(ymin = smoothed_occ-SE,
                    ymax = smoothed_occ+SE),
                width = 0) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = "Year", y = "Smoothed occupancy (derived)") +
  my.theme
occu
ggsave(filename=paste0("2_outputs/",i,"/",i,"_psiallstations.png"), plot= occu)


