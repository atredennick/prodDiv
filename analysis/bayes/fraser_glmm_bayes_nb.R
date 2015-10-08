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
plot(x,y,col=pi)

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
  int<lower=0> study[N]; // study id
  int<lower=0> grids[N]; // grid id
  int<lower=0> S; // number of studies
  int<lower=0> G; // number of grids per study
  int<lower=0> y[N]; // observation vector
  vector[N] x; // biomass vector
}
parameters{
  real b0_mu;
  real b1_mu;
  real b2_mu;
  vector[S] b0;
  vector[S] b1;
  vector[S] b2;
  vector[G] bgrids;
  real<lower=0.00001> sig_b0_study;
  real<lower=0.00001> sig_b1_study;
  real<lower=0.00001> sig_b2_study;
  real<lower=0.00001> sig_grid;
}
transformed parameters{
  real mu[N];
  for(n in 1:N)
    mu[n] <- b0[study[n]] + b1[study[n]]*x[n] + b2[study[n]]*(x[n]^2) + bgrids[grids[n]];
}
model{
// Priors
  b0_mu ~ normal(0,100);
  b1_mu ~ normal(0,100);
  b2_mu ~ normal(0,100);
  b0 ~ normal(b0_mu, sig_b0_study);
  b1 ~ normal(b1_mu, sig_b1_study);
  b2 ~ normal(b2_mu, sig_b2_study);
  bgrids ~ normal(0, sig_grid); 
  sig_b0_study ~ cauchy(0,5);
  sig_b1_study ~ cauchy(0,5);
  sig_b2_study ~ cauchy(0,5);
  sig_grid ~ cauchy(0,5);
  phi ~ cauchy(0,3);
// Likelihood
  y ~ neg_binomial_2_log(mu, phi);
}
"


####
####  Format data for STAN model
####
dat$study <- as.numeric(as.factor(dat$groups)) # makes sure the ids are continuous
dat$newgrids <- as.numeric(as.factor(dat$grid))
nstudy <- length(unique(dat$study))
ngrids <- length(unique(dat$newgrids))
inits <- list()
inits[[1]] <- list(b0_mu=0, b1_mu=0, b2_mu=0,
                   b0=rep(0,nstudy), b1=rep(0,nstudy), b2=rep(0,nstudy),
                   bgrids=rep(0.1,ngrids), sig_b0_study=0.1,
                   sig_b1_study=0.1, sig_b2_study=0.1,
                   sig_grid=0.1, phi=1)
inits[[2]] <- list(b0_mu=0.5, b1_mu=0.5, b2_mu=0.5,
                   b0=rep(0.5,nstudy), b1=rep(0.5,nstudy), b2=rep(0.5,nstudy),
                   bgrids=rep(0.2,ngrids), sig_b0_study=0.2,
                   sig_b1_study=0.2, sig_b2_study=0.2,
                   sig_grid=0.2, phi=0.5)
inits[[3]] <- list(b0_mu=-1, b1_mu=-1, b2_mu=-1,
                   b0=rep(-1,nstudy), b1=rep(-1,nstudy), b2=rep(-1,nstudy),
                   bgrids=rep(0.15,ngrids), sig_b0_study=0.15,
                   sig_b1_study=0.15, sig_b2_study=0.15,
                   sig_grid=0.15, phi=0.15)


####
####  Fit STAN model
####
datalist <- list(N=nrow(dat), study=dat$study, grids=dat$newgrids, 
                 S=length(unique(dat$study)), G=length(unique(dat$newgrids)),
                 y=dat$y, x=dat$x)
pars=c("b0", "b1", "b2", "b0_mu", "b1_mu", "b2_mu", "bgrids", "phi")
mcmc_samples <- stan(model_code=model_string, data=datalist, pars=pars, chains=0)
test <- stan(model_code=model_string, data=datalist, pars=pars, chains=1,iter = 500, warmup = 100)



rng_seed <- 123
sflist <-
  mclapply(1:num_chains, mc.cores=num_chains,
           function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                            seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                            iter=num_iters, warmup=num_iters/2, init=list(inits[[i]])))
fit <- sflist2stanfit(sflist)
long <- ggs(fit)
# ggs_caterpillar(long, family = "b2")


####
####  Plot mean response on data
####
##  Global response
newx <- seq(min(dat$x), max(dat$x), 0.01)
mean_b0 <- mean(unlist(long[long$Parameter=="b0_mu","value"]))
mean_b1 <- mean(unlist(long[long$Parameter=="b1_mu","value"]))
mean_b2 <- mean(unlist(long[long$Parameter=="b2_mu","value"]))
newy <- exp(mean_b0 + mean_b1*newx + mean_b2*newx^2)
par(mar = c(4.2,4.8,0.8,0.8))
plot(alldata.litter[!is.na(alldata.litter$log10.tot.bio),"tot.bio"],alldata.litter[!is.na(alldata.litter$log10.tot.bio),"sr"], ylab = expression("Number of species " ~ (m^{-2})),
     col="lightgrey",log="x",las=1,
     xlab=expression("Total biomass" ~ (g ~ m^{-2})),
     ylim=c(-2,54),cex.lab=1.5,cex.axis=1.3, cex=0.8)

##  Study responses
long <- as.data.frame(long)
long$study <- substr(long$Parameter, 4, length(long$Parameter))
long$study <- as.numeric(unlist(strsplit(long$study, split=']'))) 
long$param <- substr(long$Parameter, 1, 2)
for(do_study in 1:nstudy){
  tmp_long <- subset(long, study==do_study)
  b0 <- mean(tmp_long[tmp_long$param=="b0", "value"])
  b1 <- mean(tmp_long[tmp_long$param=="b1", "value"])
  b2 <- mean(tmp_long[tmp_long$param=="b2", "value"])
  
  tmp_dat <- subset(dat, study==do_study)
  tmp_newx <- seq(min(tmp_dat$x), max(tmp_dat$x), 0.01)
  tmp_newy <- exp(b0 + b1*tmp_newx + b2*tmp_newx^2)
  lines(10^tmp_newx, tmp_newy, lwd=1.3, col="red")
}
lines(10^newx, newy, lwd=3) # global prediction line

####
####  Make predictions on data for RMSE
####
bgrids <- long[grep("bgrids", long$Parameter), ]
bgrids$grid <- substr(bgrids$Parameter, 8, length(bgrids$Parameter))
bgrids$grid <- as.numeric(unlist(strsplit(bgrids$grid, split=']'))) 
b0 <- subset(long, param=="b0")
b0 <- subset(b0, Parameter!="b0_mu")
b1 <- subset(long, param=="b1")
b1 <- subset(b1, Parameter!="b1_mu")
b2 <- subset(long, param=="b2")
b2 <- subset(b2, Parameter!="b2_mu")
sim_model <- function(x, params){
  out <- exp(params[1] + params[2]*x + params[3]*x^2)
  return(out)
}

all_preds <- data.frame(newgrids=NA, study=NA, x=NA, y=NA, groups=NA, grid=NA,
                        b0=NA, b1=NA, b2=NA, bgrids=NA, pred=NA, iter=NA, chain=NA)
for(c in 1:num_chains){
  for(i in 1:num_iters){
    tmp_b0 <- subset(b0, Iteration==i & Chain==c)
    tmp_b0 <- tmp_b0[,c("value","study")]
    colnames(tmp_b0) <- c("b0", "study")
    tmp_b1 <- subset(b1, Iteration==i & Chain==c)
    tmp_b1 <- tmp_b1[,c("value","study")]
    colnames(tmp_b1) <- c("b1", "study")
    tmp_b2 <- subset(b2, Iteration==i & Chain==c)
    tmp_b2 <- tmp_b2[,c("value","study")]
    colnames(tmp_b2) <- c("b2", "study")
    tmp_bgrids <- subset(bgrids, Iteration==i & Chain==c)
    tmp_bgrids <- tmp_bgrids[,c("value","grid")]
    colnames(tmp_bgrids) <- c("bgrids", "grid")
    tmp_dat <- merge(dat, tmp_b0, by="study")
    tmp_dat <- merge(tmp_dat, tmp_b1, by="study")
    tmp_dat <- merge(tmp_dat, tmp_b2, by="study")
    tmp_dat <- merge(tmp_dat, tmp_bgrids, by.x="newgrids", by.y="grid")
    tmp_dat$pred <- with(tmp_dat, b0+b1*x+b2*x^2+bgrids)
    tmp_dat$iter <- i
    tmp_dat$chain <- c
    all_preds <- rbind(all_preds, tmp_dat)
  }
}
all_preds_quad <- all_preds[2:nrow(all_preds),]

####
####  Aggregate and calculate RMSEs
####
all_preds_quad$sqerror <- with(all_preds_quad, (exp(pred)-y)^2)
quad_rmse <- ddply(all_preds_quad, .(iter, chain), summarise,
                   rmse = sqrt(mean(sqerror)))
plot(density(quad_rmse$rmse), main="")




#######################################################
#######################################################
################  LINEAR MODEL FIT ####################
#######################################################
#######################################################
####
####  Write out the STAN model
####
# model_string <- "
# data{
#   int<lower=0> N; // observations
#   int<lower=0> study[N]; // study id
#   int<lower=0> grids[N]; // grid id
#   int<lower=0> S; // number of studies
#   int<lower=0> G; // number of grids per study
#   int<lower=0> y[N]; // observation vector
#   vector[N] x; // biomass vector
# }
# parameters{
#   real b0_mu;
#   real b1_mu;
#   vector[S] b0;
#   vector[S] b1;
#   vector[G] bgrids;
#   real<lower=0.00001> sig_b0_study;
#   real<lower=0.00001> sig_b1_study;
#   real<lower=0.00001> sig_grid;
# }
# transformed parameters{
#   real mu[N];
#   for(n in 1:N)
#     mu[n] <- b0[study[n]] + b1[study[n]]*x[n] + bgrids[grids[n]];
# }
# model{
#   // Priors
#   b0_mu ~ uniform(-300,300);
#   b1_mu ~ uniform(-300,300);
#   b0 ~ normal(b0_mu, sig_b0_study);
#   b1 ~ normal(b1_mu, sig_b1_study);
#   bgrids ~ normal(0, sig_grid); 
#   sig_b0_study ~ cauchy(0,5);
#   sig_b1_study ~ cauchy(0,5);
#   sig_grid ~ cauchy(0,5);
#   
#   // Likelihood
#   y ~ poisson_log(mu);
# }
# "
# 
# 
# ####
# ####  Format data for STAN model
# ####
# dat$study <- as.numeric(as.factor(dat$groups)) # makes sure the ids are continuous
# dat$newgrids <- as.numeric(as.factor(dat$grid))
# nstudy <- length(unique(dat$study))
# ngrids <- length(unique(dat$newgrids))
# inits <- list()
# inits[[1]] <- list(b0_mu=0, b1_mu=0,
#                    b0=rep(0,nstudy), b1=rep(0,nstudy),
#                    bgrids=rep(0.1,ngrids), sig_b0_study=0.1,
#                    sig_b1_study=0.1, sig_grid=0.1)
# inits[[2]] <- list(b0_mu=0.5, b1_mu=0.5,
#                    b0=rep(0.5,nstudy), b1=rep(0.5,nstudy),
#                    bgrids=rep(0.5,ngrids), sig_b0_study=0.5,
#                    sig_b1_study=0.5, sig_grid=0.5)
# inits[[3]] <- list(b0_mu=-1, b1_mu=-1,
#                    b0=rep(-1,nstudy), b1=rep(-1,nstudy),
#                    bgrids=rep(0.15,ngrids), sig_b0_study=0.15,
#                    sig_b1_study=0.15, sig_grid=0.15)
# 
# 
# ####
# ####  Fit STAN model
# ####
# datalist <- list(N=nrow(dat), study=dat$study, grids=dat$newgrids, 
#                  S=length(unique(dat$study)), G=length(unique(dat$newgrids)),
#                  y=dat$y, x=dat$x)
# pars=c("b0", "b1", "b0_mu", "b1_mu", "bgrids")
# mcmc_samples <- stan(model_code=model_string, data=datalist, pars=pars, 
#                      chains=0)
# 
# rng_seed <- 1234
# sflist <-
#   mclapply(1:num_chains, mc.cores=num_chains,
#            function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
#                             seed=rng_seed, chains=1, chain_id=i, refresh=-1,
#                             iter=num_iters, warmup=num_iters/2, init=list(inits[[i]])))
# fit <- sflist2stanfit(sflist)
# long <- ggs(fit)
# ggs_caterpillar(long, family="b2")
# 
# 
# ####
# ####  Plot mean response on data
# ####
# ##  Global response
# newx <- seq(min(dat$x), max(dat$x), 0.01)
# mean_b0 <- mean(unlist(long[long$Parameter=="b0_mu","value"]))
# mean_b1 <- mean(unlist(long[long$Parameter=="b1_mu","value"]))
# newy <- exp(mean_b0 + mean_b1*newx)
# par(mar = c(4.2,4.8,0.8,0.8))
# plot(alldata.litter[!is.na(alldata.litter$log10.tot.bio),"tot.bio"],alldata.litter[!is.na(alldata.litter$log10.tot.bio),"sr"], ylab = expression("Number of species " ~ (m^{-2})),
#      col="lightgrey",log="x",las=1,
#      xlab=expression("Total biomass" ~ (g ~ m^{-2})),
#      ylim=c(-2,54),cex.lab=1.5,cex.axis=1.3, cex=0.8)
# 
# ##  Study responses
# long <- as.data.frame(long)
# long$study <- substr(long$Parameter, 4, length(long$Parameter))
# long$study <- as.numeric(unlist(strsplit(long$study, split=']'))) 
# long$param <- substr(long$Parameter, 1, 2)
# for(do_study in 1:nstudy){
#   tmp_long <- subset(long, study==do_study)
#   b0 <- mean(tmp_long[tmp_long$param=="b0", "value"])
#   b1 <- mean(tmp_long[tmp_long$param=="b1", "value"])
#   
#   tmp_dat <- subset(dat, study==do_study)
#   tmp_newx <- seq(min(tmp_dat$x), max(tmp_dat$x), 0.01)
#   tmp_newy <- exp(b0 + b1*tmp_newx)
#   lines(10^tmp_newx, tmp_newy, lwd=1.3, col="red")
# }
# lines(10^newx, newy, lwd=3) # global prediction line
# 
# 
# ####
# ####  Make predictions on data for RMSE
# ####
# bgrids <- long[grep("bgrids", long$Parameter), ]
# bgrids$grid <- substr(bgrids$Parameter, 8, length(bgrids$Parameter))
# bgrids$grid <- as.numeric(unlist(strsplit(bgrids$grid, split=']'))) 
# b0 <- subset(long, param=="b0")
# b0 <- subset(b0, Parameter!="b0_mu")
# b1 <- subset(long, param=="b1")
# b1 <- subset(b1, Parameter!="b1_mu")
# 
# all_preds <- data.frame(newgrids=NA, study=NA, x=NA, y=NA, groups=NA, grid=NA,
#                         b0=NA, b1=NA, bgrids=NA, pred=NA, iter=NA, chain=NA)
# for(c in 1:num_chains){
#   for(i in 1:num_iters){
#     tmp_b0 <- subset(b0, Iteration==i & Chain==c)
#     tmp_b0 <- tmp_b0[,c("value","study")]
#     colnames(tmp_b0) <- c("b0", "study")
#     tmp_b1 <- subset(b1, Iteration==i & Chain==c)
#     tmp_b1 <- tmp_b1[,c("value","study")]
#     colnames(tmp_b1) <- c("b1", "study")
#     tmp_bgrids <- subset(bgrids, Iteration==i & Chain==c)
#     tmp_bgrids <- tmp_bgrids[,c("value","grid")]
#     colnames(tmp_bgrids) <- c("bgrids", "grid")
#     tmp_dat <- merge(dat, tmp_b0, by="study")
#     tmp_dat <- merge(tmp_dat, tmp_b1, by="study")
#     tmp_dat <- merge(tmp_dat, tmp_bgrids, by.x="newgrids", by.y="grid")
#     tmp_dat$pred <- with(tmp_dat, b0+b1*x+bgrids)
#     tmp_dat$iter <- i
#     tmp_dat$chain <- c
#     all_preds <- rbind(all_preds, tmp_dat)
#   }
# }
# all_preds_lin <- all_preds[2:nrow(all_preds),]
# 
# ####
# ####  Aggregate and calculate RMSEs
# ####
# all_preds_lin$sqerror <- with(all_preds_lin, (exp(pred)-y)^2)
# lin_rmse <- ddply(all_preds_lin, .(iter, chain), summarise,
#                    rmse = sqrt(mean(sqerror)))
# 
# 
# ####
# ####  Plot quadratic and linear RMSE distributions
# ####
# commbw <- density(lin_rmse$rmse, adjust=2)[["bw"]]
# plot(density(lin_rmse$rmse, bw=commbw), main="", lwd=3, xlim=c(2.5,2.6), las=1,
#      xlab("Root mean square error (RMSE)"), cex.lab=1.5,cex.axis=1.3, cex=0.8)
# lines(density(quad_rmse$rmse, bw=commbw), col="steelblue", lwd=3)
# 
