## Script to generate summary statistics for observed data
################################################################

## Load data, packages, functions
load("./data/pine2.RData")
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(gdata)
library(geoR)
library(readr)
library(Matrix) #bdiag
library(forcats)
library(stringr) #str_replace
library(fdrtool) #rhalfnorm
library(directlabels)
library(boot)
library(abc)
library(ncf)
source("./general/functions.R")

## Seeds ##
## split seed data by plot
seeds_by_plot <- split(seeds, f=seeds$plot)

## Compute within-plot distances
seeds_u <- lapply(seeds_by_plot, function(z)
{(z
  %>% dplyr::select(x,y)
  %>% dist()
)})
## seed locations
seed_locations <- seeds %>% dplyr::select(x,y)
## maximum distance within plots
max_dist_seeds <- max(unlist(seeds_u))
## create log distance breaks for seed spline correlograms
dist_vec_seeds <- exp(seq(log(5),log(50),length.out=6))


## Seedlings ##
## split seedling data by plot
seedlings_by_plot <- split(seedlings, f=seedlings$plot)

## Compute within-plot distances
seedlings_u <- lapply(seedlings_by_plot, function(z)
{(z
  %>% dplyr::select(x,y)
  %>% dist()
)})
## seedling locations
seedling_locations <- seedlings %>% dplyr::select(x,y)
## maximum distance within plots
max_dist_seedlings <- max(unlist(seedlings_u))
## create log distance breaks for seedling spline correlograms
dist_vec_seedlings <- exp(seq(log(5),log(60),length.out=6))


## number of summary statistics (8 = mean and sd of seeds and seedlings, e.intercept, x.intercept for seeds and seedlings)
num_stats <- 8 + length(dist_vec_seeds)+length(dist_vec_seeds)


## seed spline correlogram
obs_seeds_splinecorr <- spline.correlog(x=seed_locations$x,
                                        y=seed_locations$y,
                                        z=seeds$count,
                                        xmax=60,
                                        resamp=0,
                                        save=FALSE,
                                        df=8)
## indices of seed spline correlogram at dist_vec_seeds
seed_ind <- sapply(dist_vec_seeds, function(x){which.min(abs(x-obs_seeds_splinecorr$real$predicted$x))})

## observed spline correlograms for seedlings
obs_seedlings_splinecorr <- spline.correlog(x=seedling_locations$x,
                                            y=seedling_locations$y,
                                            z=seedlings$count,
                                            xmax=70,
                                            resamp=0,
                                            save=FALSE,
                                            df=8)
## indices of seedling spline correlogram at dist_vec_seeds
seedling_ind <- sapply(dist_vec_seedlings, function(x){which.min(abs(x-obs_seedlings_splinecorr$real$predicted$x))})


## save observed statistics
obs_stats <- (matrix(nrow=1,ncol=num_stats)
              %>% data.frame()
              %>% setNames(c("seeds_mean_pts","seeds_sd_pts",
                             sapply(seq_along(dist_vec_seeds),
                                    function(i){paste0("seeds_pt_",i)}),
                             "seeds_xint","seeds_eint",
                             "seedlings_mean_pts","seedlings_sd_pts",
                             sapply(seq_along(dist_vec_seedlings),
                                    function(i){paste0("seedlings_pt_",i)}),
                             "seedlings_xint", "seedlings_eint"))
)
obs_stats[1,]<- c(mean(seeds$count),
                  sd(seeds$count),
                  obs_seeds_splinecorr$real$predicted$y[seed_ind],
                  obs_seeds_splinecorr$real$x.intercept,
                  obs_seeds_splinecorr$real$e.intercept,
                  mean(seedlings$count),
                  sd(seedlings$count),
                  obs_seedlings_splinecorr$real$predicted$y[seedling_ind],
                  obs_seedlings_splinecorr$real$x.intercept,
                  obs_seedlings_splinecorr$real$e.intercept)


saveRDS(obs_stats, paste0("./summary_statistics/output/observed_stats.RDS"))
