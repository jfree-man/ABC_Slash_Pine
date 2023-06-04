## Frequentist parameter estimation
######################################################################

## Load packages
library(dplyr) 
library(tidyr)
library(purrr)
library(abc)
library(ggplot2)
library(HDInterval)
library(cowplot)
library(stringr)
library(Hmisc)
library(directlabels)
library(ncf)
library(RColorBrewer)
library(gdata)
library(readr)
library(Matrix)
library(forcats)
library(stringr) #str_replace
library(fdrtool) #rhalfnorm
library(directlabels)
library(boot)
library(glmmTMB)
library(TMB)
library(bbmle)
#library(gtools) #inv.logit
load("./data/pine2.RData")
source("./general/functions.R")
source("./general/abc_local.R")

## Data prep
######################################################################
# make coordinates factors
seed_data <- seeds
seed_data$pos <- numFactor(seed_data$x, seed_data$y)
# plot grouping variable (to consider covariance only within plots)
seed_data$plot_group <- factor(seed_data$plot)
# observation level grouping variable (for nugget effect)
seed_data$obs <- factor(1:nrow(seed_data))

## MLE
######################################################################

## fit Poisson log-link model
fit.mat <- glmmTMB(count ~ mat(pos + 0 | plot_group),family=poisson(link="log"),data=seed_data)
## add  nugget effect
fit.mat.nugget <- update(fit.mat,.~.+(1|obs)) 

# not sure if AIC comparison is appropriate, first model is better fit?
AICtab(fit.mat,fit.mat.nugget)

# in order, mean, sd, range and smoothness
exp(fit.mat$sdr$par.fixed)

# in order, mean, sd, range, smoothness, and nugget
exp(fit.mat.nugget$sdr$par.fixed)

# exponentiating CIs
# mean
exp(2.05309714 + c(-1,1)*1.96*0.2895127)
# sd
exp(-0.01288669 + c(-1,1)*1.96*0.1963667)
# range
exp(5.87611832 + c(-1,1)*1.96*1.8595592)
# smoothness
exp(-1.53081845 + c(-1,1)*1.96*0.4780609)
# nugget
exp(-7.70129616 + c(-1,1)*1.96*2246.5250095)#huge standard error


# is this the distance matrix used in the model?
dist_M <- fit.mat$modelInfo$reStruc$condReStruc$`pos + 0 | plot_group`$dist

# all distances
dist_obs <- as.matrix(dist(seeds[c("x","y")]))


# compare full distance matrix to that from glmmTMB fit
fit_dist <- sort(c(dist_M))
obs_dist <- sort(c(dist_obs))
all.equal(fit_dist, obs_dist)
plot(fit_dist, obs_dist) 

min(fit_dist-obs_dist)#minor differences in rounding only

## Profile Likelihood
######################################################################

# slicing = fix all other parameters and vary range, see negative log likelihood
# profiling = for each value of range parameter, optimize the other parameters

tt <- tmbprofile(fit.mat.nugget$obj, 3)
plot(value~theta, tt) #- plots negative log-likelihood (there is a global min but perhaps still relatively flat)


tt <- tmbprofile(fit.mat.nugget$obj, 5)
plot(value~theta, tt, log="y") #- weird
ggplot(tt, aes(theta, value))+ geom_point()+ylim(c(348,349))#scale_y_log10()



# what is the log-likelihood for these values of the four parameters (on the log scale)
# fit.mat$obj$fn(c(2.1, -0.01, 1, -1.53))
# confidence interval for the range on the log scale (very wide)
# 5.9 + c(-1,1)*2*sqrt(3.45)
# correlation between 2 parameters
# range and smoothness are somewhat negatively correlated but model is sophisticated enough to handle this
cov2cor(vcov(fit.mat, full = TRUE))



## Example from glmmTMB
###########################################################
d <- data.frame(z = as.vector(volcano),
                x = as.vector(row(volcano)),
                y = as.vector(col(volcano)))
set.seed(1)
d$z <- d$z + rnorm(length(volcano), sd=15) #(plus nugget ?)
d <- d[sample(nrow(d), 100), ] # spatial data set, want to be able to construc
#t original image
volcano.data <- array(NA, dim(volcano))
volcano.data[cbind(d$x, d$y)] <- d$z
image(volcano.data, main="Spatial data", useRaster=TRUE)

d$pos <- numFactor(d$x, d$y) #concatenates coordinates and makes them factors?
d$group <- factor(rep(1, nrow(d))) #creates factor variable with one level, why?

# mean of 1, not sure why?
# thing the grouping variable doesn't matter in this case because all the same group
# why are we adding 0 to pos...
f <- glmmTMB(z ~ 1 + exp(pos + 0 | group), data=d)
confint(f, "sigma")
###########################################################
