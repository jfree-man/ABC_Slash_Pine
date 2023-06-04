## Exploratory script for testing simulations using
## R parallel packages: foreach, doParallel
##
## This uses an older ABC version in which variogram values
## were used as summary statistics instead of spline correlogram 
## values
######################################################################

## Load packages
load("./data/pine2.RData")
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(gdata)
library(geoR)
library(readr)
library(Matrix)
library(forcats)
library(stringr) #str_replace
library(fdrtool) #rhalfnorm
library(directlabels)
library(boot)
library(abc)
library(foreach)
library(doParallel)
source("./general/functions.R")

## Data prep
######################################################################
## Seeds ##

# split seed data by plot
seeds_by_plot <- split(seeds, f=seeds$plot)
# get distances within plots 
seeds_u <- lapply(seeds_by_plot, function(z)
{(z
  %>% dplyr::select(x,y)
  %>% dist()
)})
# seed locations
seed_locations <- seeds %>% dplyr::select(x,y)
# maximum distance within plots
max_dist_seeds <- max(unlist(seeds_u))

# create 8 bins, equally spaced
# should it be more equal number of points in each bin?
bin_lims_seeds <- seq(0,max_dist_seeds, l=9)

# observed variogram bins
obs_seed_bins <- variog(coords = seed_locations,
                        data = seeds$count,
                        option = "bin",
                        breaks=bin_lims_seeds)

## Seedlings ##
# split seedling data by plot
seedlings_by_plot <- split(seedlings, f=seedlings$plot)
# get distances within plots 
seedlings_u <- lapply(seedlings_by_plot, function(z)
{(z
  %>% dplyr::select(x,y)
  %>% dist()
)})
# seedling locations
seedling_locations <- seedlings %>% dplyr::select(x,y)
# maximum distance within plots
max_dist_seedlings <- max(unlist(seedlings_u))
# create 8 bins, equally spaced
bin_lims_seedlings <- seq(0,max_dist_seedlings, l=9)
# observed variogram bins
obs_seedling_bins <- variog(coords = seedling_locations,
                  data = seedlings$count,
                  option = "bin",
                  breaks=bin_lims_seedlings)

# Initialize data frames
# observed summary statistics (both seeds and seedlings)
obs_stats <- (matrix(c(mean(seeds$count),
                       sd(seeds$count),
                       obs_seed_bins$v,
                       mean(seedlings$count),
                       sd(seedlings$count),
                       obs_seedling_bins$v),nrow=1)
            %>% data.frame()
            %>% setNames(c("seeds_mean","seeds_sd",
                           sapply(seq_along(obs_seed_bins$u),
                                  function(i){paste0("seeds_bin_",i)}),
                           "seedlings_mean","seedlings_sd", 
                           sapply(seq_along(obs_seedling_bins$u),
                                  function(i){paste0("seedlings_bin_",i)})))
            )
# number of simulations to compute
num_sims <- 200

# data frame to store all summary statistics and parameters from simulations
all_sims <- (matrix(nrow=num_sims,ncol=ncol(obs_stats)+10) 
             %>% data.frame()
             %>% setNames(c("seeds_mean","seeds_sd",
                            sapply(seq_along(obs_seed_bins$u),
                                   function(i){paste0("seeds_bin_",i)}),
                            "seedlings_mean","seedlings_sd", 
                            sapply(seq_along(obs_seedling_bins$u),
                                   function(i){paste0("seedlings_bin_",i)}),
                            "seed_Matern Scale","seed_Matern Shape","seed_Total CV",
                            "seed_Nugget Proportion","seed_Mean","est_Matern Scale",
                            "est_Matern Shape","est_Total CV","est_Nugget Proportion",
                            "est_Mean"))
)

## doParallel and foreach test
######################################################################

# seed base
seed_base <- 111


cl <- makePSOCKcluster(2) #cluster of 2 cores
registerDoParallel(cl)
tic_p <- Sys.time()

x <- foreach(i=1:num_sims, 
             .export=c("seed_priors",
                       "simulate_seedfall",
                       "establishment_priors",
                       "simulate_seedlings",
                       "corr_matern",
                       "get_gnorm"),
             .packages=c("dplyr",
                         "Matrix",
                         "RandomFields",
                         "geoR"
                         )) %dopar% {
  abc_simulations_vario(seed_base = seed_base,
                        i = i,
                        num_params = 10,
                        num_stats_vario = ncol(obs_stats),
                        bin_lims_seeds = bin_lims_seeds,
                        bin_lims_seedlings = bin_lims_seedlings, 
                        seeds_u = seeds_u, 
                        seed_locations = seed_locations, 
                        seedlings_u = seedlings_u, 
                        seedling_locations = seedling_locations,
                        seed_prior_function = seed_priors,
                        est_prior_function = establishment_priors)
                         }
stopCluster(cl)
toc_p <- Sys.time()
parallel_time <- toc_p-tic_p


## sequential test
######################################################################
tic_s <- Sys.time()
for (i in 1:num_sims){
  results_i <- abc_simulations_vario(seed_base = seed_base,
                        i = i,
                        num_params = 10,
                        num_stats_vario = ncol(obs_stats),
                        bin_lims_seeds = bin_lims_seeds,
                        bin_lims_seedlings = bin_lims_seedlings, 
                        seeds_u = seeds_u, 
                        seed_locations = seed_locations, 
                        seedlings_u = seedlings_u, 
                        seedling_locations = seedling_locations,
                        seed_prior_function = seed_priors,
                        est_prior_function = establishment_priors)
}
toc_s <- Sys.time()
sequential_time <- toc_s-tic_s

sequential_time
parallel_time
