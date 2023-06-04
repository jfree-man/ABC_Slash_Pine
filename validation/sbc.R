## Script run Simulation Based Calibration using the META-Farm package
## Ran on Digital Research Alliance of Canada Graham HPC cluster
######################################################################

## SBC Steps
## 1. Draw one parameter set from the joint prior.
## 2. Generate one simulated data set using the parameter set from 1.
# 1&2 taken care of by abc_simulations_splinecorr functions in functions.R

## 3. Draw M parameters from posterior after ABC
## 4. Compute rank, count number of posterior draws that are less than prior draw.

## ------------------------------------------------------
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
library(ncf)
load("~/ABC_Slash_Pine/data/pine2.RData")
source("~/ABC_Slash_Pine/general/functions.R")
source("~/ABC_Slash_Pine/general/abc_local.R")


## Get inputs from table.dat
args <- commandArgs(trailingOnly=TRUE)
## number of rank statistics to compute for each script run
num_rank_stat <- as.numeric(args[1])
## number of posterior sample draws (L)
L <- as.numeric(args[2])
## file number to save files uniquely for each run
file_number <- as.numeric(args[3])
## random seed (ensure it unique for each script run)
seed_base <- num_rank_stat * file_number


## Read in simulations
sim_params <- readRDS("~/ABC_Slash_Pine/simulations/output/params_meta_splinecorr.RDS")
sim_stats <- readRDS("~/ABC_Slash_Pine/simulations/output/summary_statistics_meta_splinecorr.RDS")

## remove summary statistics that we don't use
sim_stats <- sim_stats %>% dplyr::select(-seeds_xint,-seeds_eint, -seedlings_xint, -seedlings_eint)


## Optionally try substracting mean and dividing by sd
# sim_params <- (sim_params
#                %>% mutate(across(.cols=everything(), ~ scale.default(.x)[,]))
#                  )

# sim_stats <- (sim_stats
#                %>% mutate(across(.cols=everything(), ~ scale.default(.x)[,]))
#               )


## scale only but don't center (centering gets cancelled out in Euclidean distance computation in ABC)
sim_stats <- (sim_stats
              %>% mutate(across(.cols=everything(), ~drop(scale(.x,center=FALSE,scale=TRUE))))
)

## data frame to store rank statistics
rank_stat<- (matrix(nrow=num_rank_stat,ncol=10)
             %>% data.frame()
             %>% setNames(c("seed_Matern Scale","seed_Matern Shape","seed_Total CV",
                            "seed_Nugget Proportion","seed_Mean","est_Matern Scale",
                            "est_Matern Shape","est_Total CV","est_Nugget Proportion",
                            "est_Mean"))
)

for (j in 1:num_rank_stat){
  set.seed(seed_base+j)
  
  ## 1. Draw one parameter set from the joint prior.
  ## 2. Generate one simulated data set using the parameter set from 1.
  obs_stats <- (sim_stats
                %>% mutate(row_num=row_number())
                %>% drop_na()
                %>% sample_n(1)
  )
  obs_params <- (sim_params
                 %>% filter(row_number()==obs_stats$row_num)
  )
  sim_params_j <- sim_params %>% anti_join(obs_params, by=c("seed_Matern Scale", "seed_Matern Shape", "seed_Total CV",
                                           "seed_Nugget Proportion", "seed_Mean", "est_Matern Scale", "est_Matern Shape",
                                           "est_Total CV", "est_Nugget Proportion", "est_Mean"))
  sim_stats_j <- sim_stats %>% anti_join(obs_stats, by=c("seeds_mean_pts", "seeds_sd_pts", "seeds_pt_1", "seeds_pt_2",
                                                         "seeds_pt_3", "seeds_pt_4", "seeds_pt_5", "seeds_pt_6", "seedlings_mean_pts",
                                                         "seedlings_sd_pts", "seedlings_pt_1", "seedlings_pt_2", "seedlings_pt_3",
                                                         "seedlings_pt_4", "seedlings_pt_5", "seedlings_pt_6"))
  obs_stats <- obs_stats %>% dplyr::select(-row_num)
  
  # Run ABC
  # abc_results <- abc(target=obs_stats,
  #                    param=sim_params_j,
  #                    sumstat=sim_stats_j,
  #                    tol=0.005,
  #                    method = "rejection"
  # )
  
  # Run ABC local function (no scaling)
  abc_results <- abc_local(target=obs_stats,
                     param=sim_params_j,
                     sumstat=sim_stats_j,
                     tol=0.01,
                     method = "rejection"
  )
  
  ## 3. Draw M parameters from posterior after ABC
  post_samp <- (abc_results$unadj.values
                %>% data.frame()
                %>% sample_n(L)
  )
  
  ## 4. Compute rank, count number of posterior draws that are less than prior draw.
  for (i in 1:ncol(post_samp)){
    rank_stat[j,i] <- length(which(post_samp[,i] < obs_params[,i]))
  }
  
  if (j %% 100 ==0){
    cat(j,"\n")
  }
}

saveRDS(rank_stat, paste0("~/ABC_Slash_Pine/validation/output/sbc_scaleonly/rankstat_",file_number,".RDS"))
