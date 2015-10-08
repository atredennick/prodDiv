##  Script to reproduce findings from Fraser et al. 2015, Science
##
##  Mainly reproduces Fig. 2A and relevant information for Table 1
##
##
##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 9-27-2015


# Clear the workspace
rm(list=ls())

####
####  Load necessary libraries --------------------------------------------------
####
library(lme4)
library(rdryad) #for reading in Fraser data from Dryad
source("rsquaredglmm.R") #fxn for R2 estimates (https://github.com/jslefche/rsquared.glmm)



####
####  Source the data read-in/clean script -------------------------------------
####
source("read_clean_Fraser_data.R")



####
####  Fit global relationship with GLMM ----------------------------------------
####

### Run GLMM with grid nested in site random effects
### Using Poisson 

### Note that you MAY get convergence warnings (depending on LME4 version).
### I ran this with all different combos of optimizers like
### Ben Bolker suggests here: http://stackoverflow.com/a/21370041
### I got the same results each time (estimates and significance tests)
### and these results correspond with that of Fraser et al., so I
### think it is fine and has to do with sensitive convergence tests in
### LME4 (see here: https://github.com/lme4/lme4/issues/120). I also
### ran the site-level hierarchical Bayesian models using rstan to make
### sure we get consistent results (see /prodDiv/results/*bayes and
### /prodDiv/analysis/bayes/). We do get consistent results. Thus, given
### the information at hand, feel safe to ignore the warning messages.



### Make explicit nesting covariate for random effects efficiency
good.data <- within(good.data, site.grid <- (pi:grid)[drop = TRUE])



####
####  Fit global-extent GLMMs --------------------------------------------------
####
##  Global quadratic fit
# Write model formula
mod.formula <- as.formula(sr~log10.tot.bio+I(log10.tot.bio^2)+(1|pi)+(1|site.grid))
# Fit the model
global.quadratic.pois <- glmer(mod.formula, data=good.data, family=quasipoisson(link=log),
                               control=glmerControl(optCtrl=list(maxfun=2e4), 
                                                    optimizer = c("Nelder_Mead", "bobyqa")))
relgrad <- with(global.quadratic.pois@optinfo$derivs, solve(Hessian, gradient))
if(max(abs(relgrad)) > 0.002) { stop("relative gradient too large") }


good.data$obs <- seq(1:nrow(good.data))
mod.formula <- as.formula(sr~log10.tot.bio+I(log10.tot.bio^2)+(1|pi)+(1|site.grid)+(1|obs))
global.quadratic.quasipois <- glmer(mod.formula, data=good.data, family="poisson",
                                    control=glmerControl(optCtrl=list(maxfun=2e4), 
                                                         optimizer = c("Nelder_Mead", "bobyqa")))
summary(global.quadratic.quasipois)