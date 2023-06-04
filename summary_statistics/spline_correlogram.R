## Exploratory script to investigate spline correlogram and
## its use in summary statistics
## cleaned up version of spline_corr.R

## load packages, function, data
###########################################################
load("./data/pine2.RData")
source("./general/functions.R")
library(ncf)
library(ggplot2)
library(dplyr)
library(Matrix)
library(directlabels)

## Data prep
############################################################
## Seeds ##
# split seed data by plot
seeds_by_plot <- split(seeds, f=seeds$plot)

# get distances within plots 
seeds_u <- lapply(seeds_by_plot, function(z)
{(z
  %>% dplyr::select(x,y)
  %>% dist()
)})
## maximum distance within plots
max_dist_seeds <- max(unlist(seeds_u))

## Seedlings ##
# split seedling data by plot
seedlings_by_plot <- split(seedlings, f=seedlings$plot)

# get distances within plots 
seedlings_u <- lapply(seedlings_by_plot, function(z)
{(z
  %>% dplyr::select(x,y)
  %>% dist()
)})
## maximum distance within plots
max_dist_seedlings <- max(unlist(seedlings_u))



## validating degrees of freedom
############################################################

## default degrees of freedom (sqrt(number of data points))
############################################################
## observed data
sp_seed <- spline.correlog(seeds$x,
                           seeds$y,
                           seeds$count,
                           xmax=60, #approx max_dist_seeds
                           resamp=0,
                           save=FALSE
)
sp_seedling <- spline.correlog(seedlings$x,
                               seedlings$y,
                               seedlings$count,
                               xmax=70, #approx max_dist_seedlings
                               resamp=0,
                               save=FALSE
)
obs_seed <- data.frame(x=c(sp_seed$real$predicted$x),
                       y=c(sp_seed$real$predicted$y))

obs_seedling <- data.frame(x=c(sp_seedling$real$predicted$x),
                           y=c(sp_seedling$real$predicted$y))

# perhaps df too high, df \approx 10
ggplot(obs_seed, aes(x,y))+
  geom_line()
# definitely too large, df \approx 27
ggplot(obs_seedling, aes(x,y))+
  geom_line()


# Testing other df
############################################################
deg_free = 8
# observed data
sp_seed <- spline.correlog(seeds$x,
                      seeds$y,
                      seeds$count,
                      xmax=60,
                      resamp=0,
                      save=FALSE,
                      df=deg_free
)
sp_seedling <- spline.correlog(seedlings$x,
                               seedlings$y,
                               seedlings$count,
                               xmax=70,
                               resamp=0,
                               save=FALSE,
                               df=deg_free
)

obs_seed <- data.frame(x=c(sp_seed$real$predicted$x),
                        y=c(sp_seed$real$predicted$y))

obs_seedling <- data.frame(x=c(sp_seedling$real$predicted$x),
                                  y=c(sp_seedling$real$predicted$y))


## generate simulated data
seed_base <- 111
sims_seed_traps <- list()
sims_seed_quadrats <- list()
sims_seedling <- list()
params <- list()

for (j in 1:10){
  
  # set seed
  set.seed(seed_base + j)
  
  # sample from priors
  sp <- seed_priors()
  ep <- establishment_priors()
  
  # simulate seed and seedling data
  sims_seed_traps[[j]] <- simulate_seedfall(sp$sc,sp$sh,sp$nugget_var,sp$signal_var,sp$mu,seeds_u)
  sims_seed_quadrats[[j]] <- simulate_seedfall(sp$sc,sp$sh,sp$nugget_var,sp$signal_var,sp$mu,seedlings_u)
  sims_seedling[[j]] <- simulate_seedlings(ep$sc, ep$sh, ep$mu, ep$signal_var, ep$nugget_var, seedlings_u, sims_seed_quadrats[[j]])
  params[[j]] <- cbind(sp$sc,sp$sh,sp$total_cv, sp$prop_nug, sp$nugget_var,sp$signal_var,sp$mu)
  names(params[[j]]) <- paste0("sim_",j)
  
}


# compute spline correlograms for simulated data
spcorr_sim_seeds <- lapply(sims_seed_traps, function(z){
  sp<-spline.correlog(x=seeds$x,
                      y=seeds$y,
                      z=z,
                      xmax=60,
                      resamp=0,
                      save=FALSE
                      ,df=deg_free
  )
  data.frame(x=c(sp$real$predicted$x),
             y=c(sp$real$predicted$y))}
)  %>% bind_rows(.id="sim")

spcorr_sim_seedlings <- lapply(sims_seedling, function(z){
  sp<-spline.correlog(x=seedlings$x,
                      y=seedlings$y,
                      z=z,
                      xmax=70,
                      resamp=0,
                      save=FALSE
                      ,df=deg_free
  )
  data.frame(x=c(sp$real$predicted$x),
             y=c(sp$real$predicted$y))}
)  %>% bind_rows(.id="sim")


# seeds - simulated and observed(black)
gg_seeds <- ggplot(spcorr_sim_seeds, aes(x,y,colour=sim))+
  geom_line()+ 
  geom_line(data=obs_seed, aes(x,y),colour="black")
direct.label(gg_seeds,'first.points')

# seedlings - simulated and observed(black)
gg_seedlings <- ggplot(spcorr_sim_seedlings, aes(x,y,colour=sim))+
  geom_line()+ 
  geom_line(data=obs_seedling, aes(x,y),colour="black")
direct.label(gg_seedlings,'first.points')


## considering log/quantile breaks on spline correlogram for summary statistics
######################################################################
hist(unlist(seeds_u),breaks=50)
hist(log(unlist(seeds_u)),breaks=50)
quantile(unlist(seeds_u),probs=seq(0.2,0.9,0.2))
x <-exp(seq(log(5),log(50),length.out=6))
x
y <- 1:6
plot(x,y)

hist(unlist(seedlings_u),breaks=50)
x <-exp(seq(log(5),log(60),length.out=6))
x

seeds_u_dist <- seeds_u %>% unlist() %>% data.frame() %>% setNames("distance")
seedlings_u_dist <- seedlings_u %>% unlist() %>% data.frame() %>% setNames("distance")
## distances between seed traps and seedling quadrats
(ggplot(seeds_u_dist, aes(x=distance,fill="seed", after_stat(density)))
  + geom_histogram(breaks=seq(0,60,by=1),alpha=0.5)
  + geom_histogram(data=seedlings_u_dist, aes(x=distance,fill="seedling"),breaks=seq(0,70,by=1),alpha=0.5)
  + theme_bw()
  + ggtitle("Within Plot Distances"))

# testing how to get the correlation value that correspond to approximate log
# distance breaks
ind <-floor(exp(seq(log(5),log(300),length.out=6)))
test <- spcorr_sim_seeds %>% filter(sim==1)
test$x[ind]
sapply(test$x, function(a, x) {x[which.min(abs(a-x))]}, x)
x <-exp(seq(log(5),log(50),length.out=6))
#sapply(x, function(a=test$x, x) {a[which.min(abs(a-x))]})
# indices of test$x that are closest to x
sapply(x, function(x){which.min(abs(x-test$x))})
which.min(abs(x[1]-test$x))


## testing spline correlograms for seeds and correlation length
##############################################################
sims <- list()
params <- list()

for (j in 1:30){
  
  # set seed
  set.seed(seed_base + j)
  
  # sample from seed priors
  sp <- seed_priors()
  
  # simulate seeds
  sims[[j]] <- simulate_seedfall(sp$sc,sp$sh,sp$nugget_var,sp$signal_var,sp$mu,seeds_u)
  params[[j]] <- cbind(sp$sc,sp$sh,sp$total_cv, sp$prop_nug, sp$nugget_var,sp$signal_var,sp$mu)
  names(params[[j]]) <- paste0("sim_",j)
  
}

params <- params %>% bind_rows(.id="sim")
params <- params %>% setNames(c("sim","scale","shape","total_cv","prop_nug","nug_var","sig_var","mean"))

# simulated seed spline correlograms
spcorr_sims <- lapply(sims, function(z){
  sp<-spline.correlog(x=seeds$x,
                      y=seeds$y,
                      z=z,
                      xmax=60,
                      save=FALSE,
                      resamp=0)
  data.frame(x=c(sp$real$predicted$x),
             y=c(sp$real$predicted$y))}
)  %>% bind_rows(.id="sim")

# simulated seed correlation length
spcorr_xint <- lapply(sims, function(z){
  sp<-spline.correlog(x=seeds$x,
                      y=seeds$y,
                      z=z,
                      xmax=60,
                      save=FALSE,
                      resamp=0)
  data.frame(xint=c(sp$real$x.intercept))}
)  %>% bind_rows(.id="sim")

# observed seed spline correlogram
sp <- spline.correlog(seeds$x,
                      seeds$y,
                      seeds$count,
                      xmax=60,
                      save=FALSE,
                      resamp=0
)
# observed data
orig_data_xint <- data.frame(xint=c(sp$real$x.intercept))
orig_data <- data.frame(x=c(sp$real$predicted$x),
                        y=c(sp$real$predicted$y))

# simulated data sets vs. observed(black)
gg <- ggplot(spcorr_sims, aes(x,y,colour=sim))+
  geom_line()+ 
  geom_line(data=orig_data, aes(x,y),colour="black")
gg
direct.label(gg,'first.points')

# correlation length of simulated data vs. observed(black)
# a lot are NA
length(which(is.na(spcorr_xint$xint)))
gg <- ggplot(spcorr_xint, aes(x=xint,y=0,colour=sim))+
  geom_point(alpha=0.5)+ 
  geom_point(data=orig_data_xint, aes(x=xint,y=0),colour="black")
gg
direct.label(gg,'first.points')
