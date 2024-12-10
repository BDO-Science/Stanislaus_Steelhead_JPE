#Script for using SSN2 package to fit SSN models to trout density dataset & obtain block kriging populatation estimates as described in Isaak, D.J., Ver Hoef, J.M., Peterson, E.E., Horan, D.L. and Nagel, D.E., 2017. Canadian Journal of Fisheries and Aquatic Sciences 74: 147-156.

#Script runs using helper files in TroutDensity_BlockKrige.ssn object which is available at the SSN/STARS website as an example dataset

## Load SSN2 package into R
library("SSN2")

#check package versions
sessionInfo() 

## Import the observed and prediction point data from the .ssn directory. Prediction points represent midpoints of all reaches throughout the Salt River network plus prediction points spaced at 100m intervals for the stream blocks where population estimates will be made 
SaltFish <- ssn_import("C:\\Users\\disaak\\Desktop\\SaltBlockKrigPopEstimator\\TroutDensity_BlockKrige.ssn", predpts = c("preds", "cottonwood", "crow", "dry", "jackknife", "salt", "spring", "strawberry", "stump", "swift", "tincup", "willow"), overwrite = TRUE)

## Describe the dataset and variable names in the .ssn object, see Isaak et al. 2017 for definitions
summary(SaltFish)

library(ggplot2)

#view the stream network with overlay of observed sites & trout density values
ggplot() +
geom_sf(data = SaltFish$edges) +
geom_sf(data = SaltFish$obs, aes(color = trout_100m), size = 4) +
scale_color_viridis_c(limits = c(0, 132), option = "H") +
theme_bw()


## Create distance matrices for observed sites
ssn_create_distmat(SaltFish, overwrite = TRUE)
## Create distance matrices for prediction sites
ssn_create_distmat(SaltFish, predpts = "preds", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "cottonwood", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "crow", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "dry", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "jackknife", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "salt", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "spring", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "strawberry", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "stump", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "swift", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "tincup", only_predpts = TRUE, overwrite = TRUE)
ssn_create_distmat(SaltFish, predpts = "willow", only_predpts = TRUE, overwrite = TRUE)


#empirical Torgegrams to describe spatial structure among observed values
tg_int <- Torgegram(trout_100m ~ 1, SaltFish, cutoff = 25000)
plot(tg_int)

##Covariates: CANOPY SLOPE S1_93_11
##Response variables: trout_100m

#Use ml estimation for intial fits and model comparisons based on AIC
#After top model selected, use REML for final parameter estimates
##Fit dataset as nonspatial GLM multiple linear regression (MLR) described in Table 2 of Isaak et al. 2017
ssn_GLMmod <- ssn_glm(formula = trout_100m ~ CANOPY + SLOPE + S1_93_11, family = "gaussian", ssn.object = SaltFish, estmethod = "ml")
summary(ssn_GLMmod)
varcomp(ssn_GLMmod)

#Fit dataset as SSN model 1 in Table 2.
ssn_SSNmod1 <- ssn_glm(formula = trout_100m ~ CANOPY + SLOPE + S1_93_11, family = "gaussian", ssn.object = SaltFish, taildown_type = "exponential", tailup_type = "exponential", additive = "afvArea", estmethod = "ml")
summary(ssn_SSNmod1)
varcomp(ssn_SSNmod1)

#Fit dataset as SSN model 2 in Table 2.
ssn_SSNmod2 <- ssn_glm(formula = trout_100m ~ CANOPY + SLOPE + S1_93_11, family = "gaussian", ssn.object = SaltFish, taildown_type = "exponential", tailup_type = "exponential", euclid_type = "exponential", additive = "afvArea", estmethod = "ml")
summary(ssn_SSNmod2)
varcomp(ssn_SSNmod2)

#Fit dataset as SSN model 3 in Table 2.
ssn_SSNmod3 <- ssn_glm(formula = trout_100m ~ S1_93_11, family = "gaussian", ssn.object = SaltFish, taildown_type = "exponential", tailup_type = "exponential", additive = "afvArea", estmethod = "ml")
summary(ssn_SSNmod3)
varcomp(ssn_SSNmod3)

#Fit dataset as SSN model 4 in Table 2.
ssn_SSNmod4 <- ssn_glm(formula = trout_100m ~ 1, family = "gaussian", ssn.object = SaltFish, taildown_type = "exponential", tailup_type = "exponential", additive = "afvArea", estmethod = "ml")
summary(ssn_SSNmod4)
varcomp(ssn_SSNmod4)

#Fit dataset as SSN model 5 in Table 2.
ssn_SSNmod5 <- ssn_glm(formula = trout_100m ~ 1, family = "gaussian", ssn.object = SaltFish, taildown_type = "exponential", tailup_type = "exponential", euclid_type = "exponential", additive = "afvArea", estmethod = "ml")
summary(ssn_SSNmod5)
varcomp(ssn_SSNmod5)


#Provides AIC, deviance, etc. metrics for model comparisons
glances(ssn_GLMmod, ssn_SSNmod1, ssn_SSNmod2, ssn_SSNmod3, ssn_SSNmod4, ssn_SSNmod5)

#Provides loocv prediction averages
loocv_GLMmod <- loocv(ssn_GLMmod)
loocv_GLMmod$MSPE

loocv_SSNmod1 <- loocv(ssn_SSNmod1)
loocv_SSNmod1$MSPE

loocv_SSNmod2 <- loocv(ssn_SSNmod2)
loocv_SSNmod2$MSPE

loocv_SSNmod3 <- loocv(ssn_SSNmod3)
loocv_SSNmod3$MSPE

loocv_SSNmod4 <- loocv(ssn_SSNmod4)
loocv_SSNmod4$MSPE
 
loocv_SSNmod5 <- loocv(ssn_SSNmod5)
loocv_SSNmod5$MSPE


#augment observed data for GLM & SSN with loocv predictions & graph scatterplots
aug_GLMmod <- augment(ssn_GLMmod)
loocv_GLMmod <- loocv(ssn_GLMmod, cv_predict = TRUE)
aug_GLMmod$loocv <- loocv_GLMmod$cv_predict
print (aug_GLMmod)
ggplot(aug_GLMmod, aes(x=trout_100m, y=loocv)) + geom_point()

aug_SSNmod3 <- augment(ssn_SSNmod3)
loocv_SSNmod3 <- loocv(ssn_SSNmod3, cv_predict = TRUE)
aug_SSNmod3$loocv <- loocv_SSNmod3$cv_predict
print (aug_SSNmod3)
ggplot(aug_SSNmod3, aes(x=trout_100m, y=loocv)) + geom_point()


#augment network prediction points with kriged prediction values from SSN model
#then map predictions from SSN model as small circles 
aug_SSNpreds <- augment(ssn_SSNmod3, newdata = "preds", type = "response", se_fit = TRUE)
ggplot() +
geom_sf(data = SaltFish$edges) +
geom_sf(data = aug_SSNpreds, aes(color = .fitted), size = 2) +
scale_color_viridis_c(limits = c(0, 132), option = "H") +
theme_bw()

#then map predictions from SSN model as small circles overlaid on large circles representing the observed trout density values
aug_SSNpreds <- augment(ssn_SSNmod3, newdata = "preds", type = "response", se_fit = TRUE)
ggplot() +
geom_sf(data = SaltFish$edges) +
geom_sf(data = SaltFish$obs, aes(color = trout_100m), size = 4) + 
geom_sf(data = aug_SSNpreds, aes(color = .fitted), size = 2) +
scale_color_viridis_c(limits = c(0, 132), option = "H") +
theme_bw()

#map prediction SEs from SSN model as small circles overlaid on large circles representing observed trout density sites to create a spatial uncertainty map
ggplot() +
geom_sf(data = SaltFish$edges) +
geom_sf(data = aug_SSNpreds, aes(color = .se.fit), size = 2) +
scale_color_viridis_c(limits = c(10, 34), option = "H") +
theme_bw()


#Apply SSN to make predictions at midpoints of dense 100m reaches within stream blocks
dry <- augment(ssn_SSNmod3, newdata = "dry", type = "response", se_fit = TRUE)
jackknife <- augment(ssn_SSNmod3, newdata = "jackknife", type = "response", se_fit = TRUE)
salt <- augment(ssn_SSNmod3, newdata = "salt", type = "response", se_fit = TRUE)
spring <- augment(ssn_SSNmod3, newdata = "spring", type = "response", se_fit = TRUE)
stump <- augment(ssn_SSNmod3, newdata = "stump", type = "response", se_fit = TRUE)
willow <- augment(ssn_SSNmod3, newdata = "willow", type = "response", se_fit = TRUE)
#cottonwood <- augment(ssn_SSNmod3, newdata = "cottonwood", type = "response", se_fit = TRUE)
#crow <- augment(ssn_SSNmod3, newdata = "crow", type = "response", se_fit = TRUE)
#strawberry <- augment(ssn_SSNmod3, newdata = "strawberry", type = "response", se_fit = TRUE)
#swift <- augment(ssn_SSNmod3, newdata = "swift", type = "response", se_fit = TRUE)
#tincup <- augment(ssn_SSNmod3, newdata = "tincup", type = "response", se_fit = TRUE)

#Plot trout density predictions for 100m reaches within stream blocks
ggplot() +
geom_sf(data = SaltFish$edges) +
geom_sf(data = dry, aes(color = .fitted), size = 2) + geom_sf(data = jackknife, aes(color = .fitted), size = 2) + geom_sf(data = salt, aes(color = .fitted), size = 2) + geom_sf(data = spring, aes(color = .fitted), size = 2) + geom_sf(data = stump, aes(color = .fitted), size = 2) + geom_sf(data = willow, aes(color = .fitted), size = 2) +
scale_color_viridis_c(limits = c(0, 132), option = "H") +
theme_bw()

#Obtain estimates of mean trout density & SEs for each stream block 
dryblock <- predict(ssn_SSNmod3, "dry", block = TRUE, se.fit = TRUE)
print(dryblock)
jackknifeblock <- predict(ssn_SSNmod3, "jackknife", block = TRUE, se.fit = TRUE)
print(jackknifeblock)
saltblock <- predict(ssn_SSNmod3, "salt", block = TRUE, se.fit = TRUE)
print(saltblock)
springblock <- predict(ssn_SSNmod3, "spring", block = TRUE, se.fit = TRUE)
print(springblock)
stumpblock <- predict(ssn_SSNmod3, "stump", block = TRUE, se.fit = TRUE)
print(stumpblock)
willowblock <- predict(ssn_SSNmod3, "willow", block = TRUE, se.fit = TRUE)
print(willowblock)

#Final step, expand the mean density estimates to population estimates by multiplying with the number of 100m reaches in each block. Calculate 95% CIs around population estimates using the SE estimates.
