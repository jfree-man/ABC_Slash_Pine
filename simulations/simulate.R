## Script to generate simulated data using the META-Farm package
## Ran on Digital Research Alliance of Canada Graham HPC cluster
################################################################

## Load data, packages, functions
load("~/ABC_Slash_Pine/data/pine2.RData")
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
source("~/ABC_Slash_Pine/general/functions.R")

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

## number of parameters 
num_params <- 10

## get arguments from table.dat
args <- commandArgs(trailingOnly=TRUE)
## number of simulations (loop size)
chunk_size <- as.numeric(args[2])
## unique file number for output
file_number <- as.numeric(args[1])
## ensure random seed is different for each script run
seed_base <- as.numeric(args[1])* chunk_size


## data frame to store all summary statistics and parameters from simulations
all_sims <- (matrix(nrow=chunk_size,ncol=num_stats+num_params) 
             %>% data.frame()
             %>% setNames(c("seeds_mean_pts","seeds_sd_pts",
                            sapply(seq_along(dist_vec_seeds),
                                   function(i){paste0("seeds_pt_",i)}),
                            "seeds_xint","seeds_eint",
                            "seedlings_mean_pts","seedlings_sd_pts", 
                            sapply(seq_along(dist_vec_seedlings),
                                   function(i){paste0("seedlings_pt_",i)}),
                            "seedlings_xint", "seedlings_eint",
                            "seed_Matern Scale","seed_Matern Shape","seed_Total CV",
                            "seed_Nugget Proportion","seed_Mean","est_Matern Scale",
                            "est_Matern Shape","est_Total CV","est_Nugget Proportion",
                            "est_Mean"))
)

## run simulations
sims <- list()
for (i in 1:chunk_size){
  sims[[i]] <- abc_simulations_splinecorr(seed_base = seed_base,
                                         i = i,
                                         num_params = num_params,
                                         num_stats_sp = num_stats,
                                         dist_vec_seeds = dist_vec_seeds,
                                         dist_vec_seedlings = dist_vec_seedlings,
                                         seeds_u = seeds_u, 
                                         seed_locations = seed_locations, 
                                         seedlings_u = seedlings_u, 
                                         seedling_locations = seedling_locations,
                                         seed_prior_function = seed_priors,
                                         est_prior_function = establishment_priors
  )
  
}

## combine simulations
all_sims <- bind_rows(sims)
params <- all_sims[,21:30]
sum_stats <- all_sims[,1:20]

# save parameters and summary stats to file
saveRDS(params, paste0("~/ABC_Slash_Pine/simulations/output/params_",chunk_size,"_",file_number,".RDS"))
saveRDS(sum_stats, paste0("~/ABC_Slash_Pine/simulations/output/summary_statistics_",chunk_size,"_",file_number,".RDS"))
