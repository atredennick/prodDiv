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
library(ggplot2)
source("../rsquaredglmm.R") #fxn for R2 estimates (https://github.com/jslefche/rsquared.glmm)



####
####  Source the data read-in/clean script -------------------------------------
####
source("../read_clean_Fraser_data.R")



####
####  Plot all sites -----------------------------------------------------------
####
ggplot(good.data, aes(x=log10.tot.bio, y=sr))+
  geom_point()+
  facet_wrap("pi")



####
####  Pull out "Bartha" data -----------------------------------------------------
####
bartha.data <- subset(good.data, pi=="Bartha")
reg.mod <- glm(sr~log10.tot.bio+I(log10.tot.bio^2), family="poisson",
               data=bartha.data)
summary(reg.mod)

mixed.mod <- glmer(sr~log10.tot.bio+I(log10.tot.bio^2) + (1|grid), 
                   family="poisson", data=bartha.data,
                   control=glmerControl(optCtrl=list(maxfun=2e4), 
                                        optimizer = c("Nelder_Mead", "bobyqa")))
summary(mixed.mod)

tempdata <- bartha.data
tempdata <- tempdata[order(tempdata$log10.tot.bio),]
pred.reg <- predict(reg.mod, newdata = tempdata, type = "response", re.form=NA)
pred.mixed <- predict(mixed.mod, newdata = tempdata, type = "response", re.form=NA)



####
####  Plot Bartha example ------------------------------------------------------
####
tempdata$predreg <- pred.reg
tempdata$predmix <- pred.mixed

ggplot(tempdata)+
  geom_point(aes(x=log10.tot.bio, y=sr, color=grid))+
  geom_line(aes(x=log10.tot.bio, y=predreg), linetype=1, lwd=2)+
  geom_line(aes(x=log10.tot.bio, y=predmix), linetype=2, lwd=2)+
  xlab("log(Live+litter biomass)")+
  ylab("Species richness")+
  guides(color=FALSE)+
  theme_bw()


