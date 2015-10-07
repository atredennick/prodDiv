##  Script to reproduce:
# Fraser, L.H., Pither, J., Jentsch, A., Sternberg, M., Zobel, M., and the 
# HerbDivNet global network.  2015. Worldwide evidence of a unimodal 
# relationship between productivity and plant species richness.  Science.

##  Here I model the coefficients hierarchically, rather than just error terms.
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date:   7-20-2015

####
####  Load libraries
####
library(rstan)
library(ggmcmc)
library(reshape2)
library(plyr)
library(parallel)
library(rdryad)


####
####  Set some global options for iters and chains
####
num_iters <- 200
num_chains <- 3


####
####  Read in data and format as Fraser et al.
####
setwd("/Users/atredenn/Desktop/")
alldata <- dryad_getfile("http://datadryad.org/bitstream/handle/10255/dryad.90398/fraser_plotdata.csv?sequence=1")

# first 10 columns are summary data
alldata.summary <- alldata[,c("grid","original.name","community.type","country","pi","plot","biomass","litter","tot.bio","sr")]  

# column headings are self-explanatory; but note that "pi" is what we equate with "site" throughout the paper

# exclude grids that had not litter measurements conducted in all 64 quadrats
complete.grid.numbers <- 
  setdiff(unique(alldata.summary$grid),names(table(alldata.summary[is.na(alldata.summary$litter),"grid"])[table(alldata.summary[is.na(alldata.summary$litter),"grid"]) == 64]))

num.complete.grids <- length(complete.grid.numbers)

# now we have 151 grids with complete or mostly complete measurements
# (10 have between 1 and 4 quadrats missing litter)

# now subset data to only those grids
alldata.litter <- alldata.summary[!is.na(match(alldata.summary$grid,as.numeric(complete.grid.numbers))),]

# there are some quadrats with zero species, and two (7023, 7030)
# in grid 110 that have valid litter biomass (zero) but no live species;
# these remain as VALID data points

# put "NA" in all other quadrats
na.rownames <- row.names(alldata.litter[alldata.litter$sr  ==  0,][setdiff(row.names(alldata.litter[alldata.litter$sr  ==  0,]),c("7023","7030")),])

# remaining quadrats with neither biomass nor species can be "NA"
alldata.litter[na.rownames,"biomass"] <- NA
alldata.litter[na.rownames,"litter"] <- NA
alldata.litter[na.rownames,"tot.bio"] <- NA
alldata.litter[na.rownames,"sr"] <- NA

# create a log10 total biomass (+1) column
alldata.litter$log10.tot.bio <- log10(alldata.litter$tot.bio+1) 

# ---------------------------------------------------------
# DONE data input / cleaning
# ---------------------------------------------------------

y = alldata.litter$sr
x = alldata.litter$log10.tot.bio
groups = as.numeric(alldata.litter$pi)
N.groups = max(groups)
grid = alldata.litter$grid
N.grids = max(grid)
plot(x,y,col=groups)

#8 missing values in x, 2 in y; omit
dat = na.exclude(data.frame(x,y,groups,grid))

#######################################################
#######################################################
##############  QUADRATIC MODEL FIT ###################
#######################################################
#######################################################
####
####  Write out the STAN model
####
model_string <- "
data{
  int<lower=0> N; // observations
  int<lower=0> grids[N]; // grid id
  int<lower=0> G; // number of grids per study
  int<lower=0> y[N]; // observation vector
  vector[N] x; // biomass vector
}
parameters{
  real b0;
  real b1;
  real b2;
  vector[G] bgrids;
  real<lower=0.00001> sig_grid;
}
transformed parameters{
  real mu[N];
  for(n in 1:N)
    mu[n] <- b0 + b1*x[n] + b2*(x[n]^2) + bgrids[grids[n]];
}
model{
// Priors
  b0 ~ normal(0,100);
  b1 ~ normal(0,100);
  b2 ~ normal(0,100);
  bgrids ~ normal(0, sig_grid); 
  sig_grid ~ cauchy(0,5);

// Likelihood
  y ~ poisson_log(mu);
}
"


####
####  Format data for STAN model
####
dat$study <- as.numeric(as.factor(dat$groups)) # makes sure the ids are continuous
dat$newgrids <- as.numeric(as.factor(dat$grid))
nstudy <- length(unique(dat$study))
ngrids <- length(unique(dat$newgrids))


####
####  Fit STAN models
####
dat.now <- dat[dat$study==1, ]
datalist <- list(N=nrow(dat.now), grids=dat.now$newgrids, 
                 G=length(unique(dat.now$newgrids)),
                 y=dat.now$y, x=dat.now$x)
pars=c("b0", "b1", "b2")
model <- stan(model_code=model_string, data=datalist, chains=0, pars=pars)

out_fits <- data.frame(site_id=NA, b0=NA, b1=NA, b2=NA)
for(i in 1:length(unique(dat$study))){
  dat.now <- dat[dat$study==i, ]
  dat.now$newgrids <- dat.now$newgrids - min(dat.now$newgrids) + 1
  if(i==11) dat.now[which(dat.now$newgrids==151),"newgrids"] <- 14
  datalist <- list(N=nrow(dat.now), grids=dat.now$newgrids, 
                   G=length(unique(dat.now$newgrids)),
                   y=dat.now$y, x=dat.now$x)
  pars=c("b0", "b1", "b2")
  fit <- stan(fit=model, data=datalist, chains=2, 
              iter = 2000, warmup = 1000, pars=pars)
  outfit <- summary(fit)$summary[,"mean"]
  out_fits[i,pars] <- outfit[1:3]
  out_fits[i,"site_id"] <- i
}
out_fits$pi <- sort(unique(alldata.litter$pi))
write.csv(out_fits, "../../results/stan_bayes_quadratic_estimates.csv")

### Compare to LME4 results
library(lme4)
out_fits_lme <- data.frame(site_id=NA, b0=NA, b1=NA, b2=NA)
for(i in 1:length(unique(dat$study))){
  dat.now <- dat[dat$study==i, ]
  dat.now$newgrids <- dat.now$newgrids - min(dat.now$newgrids) + 1
  if(i==11) dat.now[which(dat.now$newgrids==151),"newgrids"] <- max(dat.now[which(dat.now$newgrids!=151),"newgrids"])+1
  num.grids <- length(unique(dat.now$newgrids))
  if(num.grids > 2){
    mod <- as.formula(y~x+I(x^2)+(1|newgrids))
    fit <- glmer(mod, family = "poisson", data=dat.now,
                 control=glmerControl(optCtrl=list(maxfun=2e4)))
    relgrad <- with(fit@optinfo$derivs, solve(Hessian, gradient))
    if(max(abs(relgrad)) > 0.002) { stop("relative gradient too small") }
  }
  if(num.grids < 3){
    mod <- as.formula(y~x+I(x^2))
    fit <- glm(mod, family = "poisson", data=dat.now)
  }
  
  out_fits_lme[i,pars] <- summary(fit)$coef[,1]
  out_fits_lme[i,"site_id"] <- i
  print(i)
}
out_fits_lme$pi <- sort(unique(alldata.litter$pi))
write.csv(out_fits_lme, "../../results/lme4_quadratic_estimates.csv")

# ggplot(alldata.litter, aes(x=log10.tot.bio, y=sr))+
#   geom_point(alpha=0.5, size=0.95)+
#   stat_smooth(formula = y~x, se=FALSE, method = "lm", col="steelblue")+
#   stat_smooth(formula = y~x+I(x^2), se=FALSE, method = "lm", col="darkred")+
#   facet_wrap("pi")

