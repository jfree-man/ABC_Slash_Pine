## Script to aggregate simulated data files using the META-Farm package
## Ran on Digital Research Alliance of Canada Graham HPC cluster
#######################################################################

## load packages
library(dplyr)
library(tidyverse)

## combine parameter files
param_data <- (list.files(path = "~/ABC_Slash_Pine/simulations/output/",
                         pattern="^params_[0-9]+_[0-9]+\\.RDS$",
                         full.names = TRUE)
             %>% map_dfr(readRDS)
)

## combine summary statistic files
sum_data <- (list.files(path = "~/ABC_Slash_Pine/simulations/output/",
                        pattern="^summary_statistics_[0-9]+_[0-9]+\\.RDS$",
                        full.names = TRUE)
             %>% map_dfr(readRDS)
)

## save aggregated data (all 1 million simulations)
saveRDS(param_data,file="~/ABC_Slash_Pine/simulations/output/params_meta_splinecorr.RDS")
saveRDS(sum_data, file="~/ABC_Slash_Pine/simulations/output/summary_statistics_meta_splinecorr.RDS")
