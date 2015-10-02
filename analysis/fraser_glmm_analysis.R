f##  Script to reproduce findings from Fraser et al. 2015, Science
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
### I ran this with all difference combos of optimizers like
### Ben Bolker suggests here: http://stackoverflow.com/a/21370041
### I got the same results each time (estimates and significance tests)
### and these results correspond with that of Fraser et al., so I
### think it is fine and has to do with sensitive convergence tests in
### LME4 (see here: https://github.com/lme4/lme4/issues/120)



### Make explicit nesting covariate for random effects efficiency
good.data <- within(good.data, site.grid <- (pi:grid)[drop = TRUE])

##  Global quadratic fit
# Write model formula
mod.formula <- as.formula(sr~log10.tot.bio+I(log10.tot.bio^2)+(1|pi)+(1|site.grid))
# Fit the model
global.quadratic.pois <- glmer(mod.formula, data=good.data, family="poisson",
                      control=glmerControl(optCtrl=list(maxfun=2e4)))
relgrad <- with(global.quadratic.pois@optinfo$derivs, solve(Hessian, gradient))
if(max(abs(relgrad)) > 0.002) { stop("relative gradient too small") }

# Get summary
summary(global.quadratic.pois)
# Make predictions w/o random effects for RMSE
pred.quad <- predict(global.quadratic.pois, type = "response",re.form=NA)
rmse.quad <- sqrt(mean((good.data$sr - pred.quad)^2))


##  Global linear fit
# Write model formula
mod.formula <- as.formula(sr~log10.tot.bio+(1|pi)+(1|site.grid))
# Fit the model
global.linear.pois <- glmer(mod.formula, data=good.data, family="poisson",
                               control=glmerControl(optCtrl=list(maxfun=2e4)))
relgrad <- with(global.linear.pois@optinfo$derivs, solve(Hessian, gradient))
if(max(abs(relgrad)) > 0.002) { stop("relative gradient too small") }

# Make predictions for RMSE
pred.linear <- predict(global.linear.pois, type = "response",re.form=NA)
rmse.linear <- sqrt(mean((good.data$sr - pred.linear)^2))

write.csv(data.frame(model=c("quadratic", "linear"),
                     rmse=c(rmse.quad, rmse.linear)),
          "../results/global_rmses.csv")

### Compare models
global_model_anova <- anova(global.quadratic.pois, global.linear.pois, test = "chi")
capture.output(global_model_anova,file="../results/global_model_comparison_ANOVA.txt")


####
####  Fit site-level (pi) GLMMS ------------------------------------------------
####
p.value  =  0.1   ## this is the alpha level used by Adler et al. 2011

# Get vector of investigator names (= sites)
pi.names <- as.character(sort(unique(good.data$pi)))

# set up output dataframe for poisson
site.poisson.log10.glmm.results <- data.frame(
  investigator  =  pi.names,
  sampsize  =  NA,
  term1.coef = NA,
  term1.se = NA,
  term1.pval = NA,
  term2.coef = NA,
  term2.se = NA,
  term2.pval = NA,
  type = NA,
  min.biomass = NA,
  max.biomass = NA,
  intercept = NA)

site.final.log10.glmm.results <- data.frame(
  investigator  =  pi.names,
  sampsize  =  NA,
  term1.coef = NA,
  term1.se = NA,
  term1.pval = NA,
  term2.coef = NA,
  term2.se = NA,
  term2.pval = NA,
  type = NA,
  min.biomass = NA,
  max.biomass = NA,
  intercept = NA)

####
#### Begin loop through 28 investigators (sites) -------------------------------
####
for (i in 1:length(pi.names)){
  # get tmp data
  tmp.data <- good.data[good.data$pi == pi.names[i], ]
  # get sample size
  tmpsampsize <- nrow(na.omit(tmp.data[c("sr","log10.tot.bio")]))
  tmprow <- which(site.poisson.log10.glmm.results$investigator==pi.names[i])
  site.poisson.log10.glmm.results[tmprow,"sampsize"] <- tmpsampsize
  
  ###
  ### Run poisson generalized linear mixed effects model (GLMM)
  ###
  if(pi.names[i] != "Wilson"){
  tmp.glmm <- glmer(sr ~ log10.tot.bio + I(log10.tot.bio^2) + (1|grid),
                    family = "poisson", data = tmp.data)
  ## Test to see if warnings matter (https://github.com/lme4/lme4/issues/120)
  relgrad <- with(tmp.glmm@optinfo$derivs, solve(Hessian, gradient))
  if(max(abs(relgrad)) > 0.002) { stop("relative gradient too small") }
  }
  
  if(pi.names[i] == "Wilson"){
    tmp.glmm <- glm(sr ~ log10.tot.bio + I(log10.tot.bio^2),
                      family = "poisson", data = tmp.data)
  }
  
  glmm.summ <- c(t(summary(tmp.glmm)$coef)[-3,2:3])
  glmm.int <- c(t(summary(tmp.glmm)$coef)[1])
  
  # fill poisson output dataframe with quadratic model output
  mod.columns <- c("term1.coef","term1.se","term1.pval",
                   "term2.coef","term2.se","term2.pval")
  site.poisson.log10.glmm.results[tmprow,mod.columns] <- glmm.summ
  
  # get intercept estimate
  site.poisson.log10.glmm.results[tmprow,"intercept"] <-  glmm.int
  
  # note type of error distribution used for current model (poisson)
  site.poisson.log10.glmm.results[tmprow,"type"] <- "poisson"
  
  # get min and max biomass values
  tmp.minmax.biomass <- range(tmp.data$log10.tot.bio,na.rm = T)
  biocols <- c("min.biomass","max.biomass")
  site.poisson.log10.glmm.results[tmprow,biocols] <- tmp.minmax.biomass
  
  
  ###
  ### Check significance of quadratic term; fit linear if not
  ###
  tmp.lin.glmm <- glmer(sr~log10.tot.bio + (1|grid),
                           family = "poisson", data = tmp.data)
  tmplin.summ <- t(summary(tmp.lin.glmm)$coef)[,2][-3]
  tmplin.int <- t(summary(tmp.lin.glmm)$coef)[1]
  mod.lin.columns <- c("term1.coef","term1.se","term1.pval")
  mod.quad.columns <- c("term2.coef","term2.se","term2.pval")
  
  # Reset summaries IF linear model is better fit
  if(anova(tmp.glmm,tmp.lin.glmm,test = "Chisq")$Pr[2] > p.value) {
    # re-assign linear terms
    site.poisson.log10.glmm.results[tmprow,mod.lin.columns] <- tmplin.summ
    site.poisson.log10.glmm.results[tmprow,"intercept"] <- tmplin.int
    
    # set quadratic terms to NA
    site.poisson.log10.glmm.results[tmprow,mod.quad.columns] <- NA
  }
  
  
  # now put the appropriate results in the "final" output table; 
  site.final.log10.glmm.results[tmprow,] <- site.poisson.log10.glmm.results[tmprow,]
  print(paste("Done with", pi.names[i], "site-level regression."))
  
} # end site-level regression loop

# remove the intermediate output table
rm(site.poisson.log10.glmm.results)



####
####  Categorize relationships by significant form -----------------------------
####
site.final.log10.glmm.results$form <- ifelse(
  site.final.log10.glmm.results$term2.coef<0,"CD",
  ifelse(site.final.log10.glmm.results$term2.coef>0,"CU",NA));

site.final.log10.glmm.results[is.na(site.final.log10.glmm.results$term2.coef) & !is.na(site.final.log10.glmm.results$term1.coef),"form"] <- 
  ifelse(
    site.final.log10.glmm.results[is.na(site.final.log10.glmm.results$term2.coef) & !is.na(site.final.log10.glmm.results$term1.coef),"term1.pval"]>p.value,"NS",
    ifelse(
      site.final.log10.glmm.results[is.na(site.final.log10.glmm.results$term2.coef)&!is.na(site.final.log10.glmm.results$term1.coef),"term1.coef"]>0,"POS","NEG"))

site.final.log10.glmm.results[is.na(site.final.log10.glmm.results$form),"form"] <- "NS";




# if GLM is not significant, place NAs
site.final.log10.glmm.results[site.final.log10.glmm.results$form == "NS",c("term1.coef","term1.se","term1.pval","term2.coef","term2.se","term2.pval")] <- NA

# order the form variable:
site.final.log10.glmm.results$form <- ordered(site.final.log10.glmm.results$form,levels = c("NS","POS","NEG","CU","CD"));



####
####  Write results to CSV -----------------------------------------------------
####
write.csv(site.final.log10.glmm.results,"../results/site_scale_results.csv");
# table(site.final.log10.glmm.results)



####
####  Make plot a la Fraser et al. code ----------------------------------------
####

# get the results from the site extent analyses, identify  significant ones:
sig.results <- site.final.log10.glmm.results[site.final.log10.glmm.results$form!= "NS",]

# set default par
par.default <- par(no.readonly = TRUE)


# --------------
# Begin plotting function (this is not generalized; boo.) 
# --------------
scatterhist.lines.log10  =  function(x, y, xlab = "", ylab = ""){
  zones = matrix(c(2,0,1,3), ncol = 2, byrow = TRUE)
  layout(zones, widths = c(4/5,1/5), heights = c(1/5,4/5))
  xhist  =  hist(log10(x), plot = FALSE,nclass = 25)
  yhist  =  hist(y, plot = FALSE,nclass = 25)
  top  =  max(c(xhist$counts, yhist$counts))
  par(mar = c(4.2,4.8,0.8,0.8))
  plot(x,y,log = "x",las = 1,col = "grey",cex = 0.75,cex.lab = 1.5,cex.axis = 1.3,
       ylab  =  expression("Number of species " ~ (m^{-2})),
       xlab  =  expression("Total biomass " ~ (g~m^{-2})))
  
  
  # add GLM regression lines
  for (i in 1:nrow(sig.results)){
    tempdata <- good.data[good.data$pi == sig.results$investigator[i],]
    tempdata <- tempdata[order(tempdata$log10.tot.bio),]
    
    # get regression terms:
    if(!is.na(sig.results[i,"term2.coef"])){
      lines(tempdata$tot.bio+1,
            predict(glmer(sr~log10.tot.bio+I(log10.tot.bio^2)+(1|grid),tempdata,family = "poisson"),newdata = tempdata,type = "response", re.form=NA),
            col = ifelse(sig.results[i,"form"] == "CD","red","purple"),lwd = 2.3)
    }
    else {
      lines(tempdata$tot.bio+1,
            predict(glmer(sr~log10.tot.bio+(1|grid),tempdata,family = "poisson"),newdata = tempdata,type = "response",re.form=NA),
            col = ifelse(sig.results[i,"form"] == "NEG","green","darkgrey"),lwd = 2.3)
    }
  }
  
  # add poisson global prediction
  tempdata <- good.data
  tempdata <- tempdata[order(tempdata$log10.tot.bio),]
  global.pred <- predict(global.quadratic.pois, tempdata, type = "response",re.form=NA)
  lines(tempdata$tot.bio, global.pred, col="black", lwd=3)
  
  par(mar = c(0,4.8,1,1))
  barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0)
  par(mar = c(4.2,0,1,1))
  barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0, horiz = TRUE);
  
  par(new = TRUE) # overlay existing plot
  par(mar = c(0,1,0,0),mgp = c(3,0.4,0),tcl = -0.3) # strip out the margins for the inset plot
  par(fig = c(0.13,0.36,0.6,0.77)) # fig shrinks and places relative to figure region
  barplot(table(site.final.log10.glmm.results$form),
          names = c("NS","Linear+","Linear-","Concave+","Concave-"),
          las = 2,col = c("white","darkgrey","green","purple","red"),
          cex.axis = 1,cex.lab = 1,cex.names = 1.1,ylim = c(0,20),space = 0.1);
  mtext("Frequency",side = 2,line = 1.5,cex = 1)
  
}



# plot the figure
png("../results/fig2A_ATT.png", width = 9, height=8, units = "in", res = 100)
par(par.default)
scatterhist.lines.log10(good.data$tot.bio+1,good.data$sr);
par(par.default)
dev.off()

