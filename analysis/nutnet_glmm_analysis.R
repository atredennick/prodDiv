##  Script to reproduce findings from Adler et al. 2011, Science
##
##  Global-extent Generalized Linear Mixed Effects Model version
##
##
##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 10-03-2015


# Clear the workspace
rm(list=ls())

####
####  Load necessary libraries --------------------------------------------------
####
library(lme4)
source("rsquaredglmm.R") #fxn for Rsqaures (https://github.com/jslefche/rsquared.glmm)



####
####  Source the data read-in/clean script -------------------------------------
####
good.data <- read.csv("../nutnet_data/nutnet_proddiv_data.csv")
good.data <- good.data[complete.cases(good.data$live.biomass),]
good.data$logbiomass <- log(good.data$live.biomass)


####
####  Fit global relationship with GLMM ----------------------------------------
####

### Run GLMM with block nested in site random effects
### Using Poisson 

### Note that you MAY get convergence warnings (depending on LME4 version).
### I ran this with all difference combos of optimizers like
### Ben Bolker suggests here: http://stackoverflow.com/a/21370041
### I got the same results each time (estimates and significance tests)
### and these results correspond with that of Fraser et al., so I
### think it is fine and has to do with sensitive convergence tests in
### LME4 (see here: https://github.com/lme4/lme4/issues/120)



### Make explicit nesting covariate for random effects efficiency
good.data$siteblock <-paste(good.data$site,good.data$block,sep="")

##  Global quadratic fit
# Write model formula
mod.formula <- as.formula(richness~logbiomass+I(logbiomass^2)+(1|site/siteblock))
# mod.formula <- as.formula(richness~logbiomass+I(logbiomass^2)+(1|site) +(1|siteblock)) # same result
# Fit the model
global.quadratic.pois <- glmer(mod.formula, data=good.data, family="poisson",
                      control=glmerControl(optCtrl=list(maxfun=2e4)))
relgrad <- with(global.quadratic.pois@optinfo$derivs, solve(Hessian, gradient))
if(max(abs(relgrad)) > 0.002) { stop("relative gradient too small") }

# Get summary
summary(global.quadratic.pois)
# Make predictions w/o random effects for RMSE
pred.quad <- predict(global.quadratic.pois, type = "response",re.form=NA)
rmse.quad <- sqrt(mean((good.data$richness - pred.quad)^2))


##  Global linear fit
# Write model formula
mod.formula <- as.formula(richness~logbiomass+(1|block/site))
# Fit the model
global.linear.pois <- glmer(mod.formula, data=good.data, family="poisson",
                               control=glmerControl(optCtrl=list(maxfun=2e4)))
relgrad <- with(global.linear.pois@optinfo$derivs, solve(Hessian, gradient))
if(max(abs(relgrad)) > 0.002) { stop("relative gradient too small") }

# Make predictions for RMSE
pred.linear <- predict(global.linear.pois, type = "response",re.form=NA)
rmse.linear <- sqrt(mean((good.data$richness - pred.linear)^2))

# Get Rsquares
r2.lin <- rsquared.glmm(global.linear.pois)
r2.quad <- rsquared.glmm(global.quadratic.pois)

write.csv(data.frame(model=c("quadratic", "linear"),
                     rmse=c(rmse.quad, rmse.linear),
                     r2_marg=c(r2.quad$Marginal, r2.lin$Marginal),
                     r2_cond=c(r2.quad$Conditional, r2.lin$Conditional)),
          "../results/nutnet_global_rmses_varexpl.csv")

### Compare models
global_model_anova <- anova(global.quadratic.pois, global.linear.pois, test = "chi")
capture.output(global_model_anova,file="../results/nutnet_global_model_comparison_ANOVA.txt")

