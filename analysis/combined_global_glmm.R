##  Script to run across study global-extent analysis
##  for the bivariate relationship between productivity and diversity.
##  Uses data from:
##      1. Adler et al. (2011) Science
##      2. Fraser et al. (2015) Science

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Date created: 10-09-2015



####
####  Load libraries -----------------------------------------------------------
####
library(lme4)
library(rdryad)
source("rsquaredglmm.R") #fxn for R2 estimates (https://github.com/jslefche/rsquared.glmm)



####
####  Read in and combine data sets --------------------------------------------
####
### Fraser data
source("read_clean_Fraser_data.R")
fraser.data <- good.data
fraser.data <- fraser.data[complete.cases(fraser.data$biomass),]
fraser.data <- fraser.data[which(fraser.data$biomass>0),] # exclude 0s
fraser.data$loglivebiomass <- log(fraser.data$biomass)

### Adler data
adler.data <- read.csv("../nutnet_data/nutnet_proddiv_data.csv")
adler.data <- adler.data[complete.cases(adler.data$live.biomass),]
adler.data$loglivebiomass <- log(adler.data$live.biomass)

### Combine data sets
all.data <- fraser.data[,c("pi", "grid", "sr", "loglivebiomass")]
all.data$study <- "fraser"
adler.tmp <- adler.data[,c("site", "block", "richness", "loglivebiomass")]
adler.tmp$study <- "adler"
colnames(adler.tmp) <- colnames(all.data)
all.data <- rbind(all.data, adler.tmp)
colnames(all.data) <- c("site", "block", "richness", "loglivebiomass", "study")



####
####  Fit one big GLMM ---------------------------------------------------------
####
### Make explicit nesting for plots in sites
all.data$siteblock <- paste0(all.data$site,all.data$block)
# all.data$loglivebiomass <- scale(all.data$loglivebiomass) #helps glmer converge, but gets same answer as non-scaled

### Write the quadratic version GLMM
### Includes study and block nested within site random effects
quad.mod.form <- as.formula(richness ~ loglivebiomass + I(loglivebiomass^2) +
                              (1|study) + (1|site) + (1|siteblock))

quad.mod <- glmer(quad.mod.form, data=all.data, family="poisson",
                  control=glmerControl(optCtrl=list(maxfun=2e4), 
                                       optimizer = c("bobyqa","Nelder_Mead")))
### Check model diagnostics since warnings thrown
relgrad <- with(quad.mod@optinfo$derivs, solve(Hessian, gradient))
if(max(abs(relgrad)) > 0.001) { stop("relative gradient too large") }



lin.mod.form <- as.formula(richness ~ loglivebiomass +
                              (1|study) + (1|site) + (1|siteblock))
lin.mod <- glmer(lin.mod.form, data=all.data, family="poisson",
                  control=glmerControl(optCtrl=list(maxfun=2e4), 
                                       optimizer = c("bobyqa","Nelder_Mead")))
### Check model diagnostics since warnings thrown
relgrad <- with(lin.mod@optinfo$derivs, solve(Hessian, gradient))
if(max(abs(relgrad)) > 0.001) { stop("relative gradient too large") }


### Look at summary and R2
summary(quad.mod)
summary(lin.mod)
r2.quad <- rsquared.glmm(quad.mod)
r2.lin <- rsquared.glmm(lin.mod)

### Make predicitons; get RMSE
pred.quad <- predict(quad.mod, type = "response",re.form=NA)
rmse.quad <- sqrt(mean((all.data$richness - pred.quad)^2))
pred.lin <- predict(lin.mod, type = "response",re.form=NA)
rmse.lin <- sqrt(mean((all.data$richness - pred.lin)^2))


### Write out the results
write.csv(data.frame(model=c("quadratic", "linear"),
                     rmse=c(rmse.quad, rmse.lin),
                     r2_marg=c(r2.quad$Marginal, r2.lin$Marginal),
                     r2_cond=c(r2.quad$Conditional, r2.lin$Conditional)),
          "../results/combined_global_rmses_varexpl.csv")

