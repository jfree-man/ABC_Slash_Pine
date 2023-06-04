## Script to test ABC using a known parameter set
#######################################################################

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
library(boot)
library(ggpubr) #ggarrange
load("./data/pine2.RData")
source("./general/functions.R")

## select a "good" parameter set from within priors
#######################################################################
seed_shape <- 1.5
seed_alpha <- 30 
seed_scale <- seed_alpha/(2*sqrt(seed_shape))
seed_nugget_prop <- 0.24 
seed_total_cv <- 0.5 
seed_signal_var <- (seed_total_cv)^2*(1-seed_nugget_prop)
seed_nugget_var <- (seed_total_cv)^2*seed_nugget_prop
seed_mu = log(30) 

est_shape <- 3
est_alpha <- 5
est_scale <- est_alpha/(2*sqrt(est_shape))
est_nugget_prop <- 0.5
est_total_cv <- 1.5
est_mu = -1
est_signal_var <- (est_total_cv)^2*(1-est_nugget_prop)
est_nugget_var <- (est_total_cv)^2*est_nugget_prop

known_params <- matrix(c(seed_scale,
                         seed_shape,
                         seed_nugget_prop,
                         seed_total_cv,
                         seed_mu,
                         est_scale,
                         est_shape,
                         est_nugget_prop,
                         est_total_cv,
                         est_mu),nrow=1) %>% data.frame()
names(known_params) <- c("seed_Matern.Scale",
                         "seed_Matern.Shape",
                         "seed_Nugget.Proportion",
                         "seed_Total.CV",
                         "seed_Mean",
                         "est_Matern.Scale",
                         "est_Matern.Shape",
                         "est_Nugget.Proportion",
                         "est_Total.CV",
                         "est_Mean") 

## Data prep
#######################################################################
# split seed data by plot
seeds_by_plot <- split(seeds, f=seeds$plot)

# get distances within plots 
# for now, ignoring distances between plots
seeds_u <- lapply(seeds_by_plot, function(z)
{(z
  %>% dplyr::select(x,y)
  %>% dist()
)})
# seed locations
seed_locations <- seeds %>% dplyr::select(x,y)
# maximum distance within plots
max_dist_seeds <- max(unlist(seeds_u))
# create log breaks along spline correlogram
dist_vec_seeds <- exp(seq(log(5),log(50),length.out=6))


## Seedlings ##
# split seedling data by plot
seedlings_by_plot <- split(seedlings, f=seedlings$plot)

# get distances within plots 
# for now, ignoring distances between plots
seedlings_u <- lapply(seedlings_by_plot, function(z)
{(z
  %>% dplyr::select(x,y)
  %>% dist()
)})
# seedling locations
seedling_locations <- seedlings %>% dplyr::select(x,y)
# maximum distance within plots
max_dist_seedlings <- max(unlist(seedlings_u))
# create log breaks along spline correlogram
dist_vec_seedlings <- exp(seq(log(5),log(60),length.out=6))

## "Observed" data set using known parameters
#######################################################################
set.seed(111)
# generate one "observed" data set

sim_seeds_traps <- simulate_seedfall(seed_scale,
                                     seed_shape,
                                     seed_nugget_var,
                                     seed_signal_var,
                                     seed_mu,
                                    seeds_u)
# observed spline correlograms for seeds
sim_seeds_pts <- spline.correlog(x=seed_locations$x,
                                 y=seed_locations$y,
                                 z=sim_seeds_traps,
                                 xmax=60,
                                 resamp=0,
                                 save=FALSE,
                                 df=8)
# get indices closest to dist_vec_seeds
seed_ind <- sapply(dist_vec_seeds, function(x){which.min(abs(x-sim_seeds_pts$real$predicted$x))})

# simulate seedfall (at seedling locations)
sim_seeds_land <- simulate_seedfall(seed_scale,
                                    seed_shape,
                                    seed_nugget_var,
                                    seed_signal_var,
                                    seed_mu,
                                    seedlings_u)

# simulate seedlings
sim_seedlings <- simulate_seedlings(est_scale,
                                    est_shape,
                                    est_mu,
                                    est_signal_var,
                                    est_nugget_var,
                                    seedlings_u,
                                    sim_seeds_land)

# compute spline correlograms
sim_seedling_pts <- spline.correlog(x=seedling_locations$x,
                                    y=seedling_locations$y,
                                    z=sim_seedlings,
                                    xmax=70,
                                    resamp=0,
                                    save=FALSE,
                                    df=8)
# get indices closest to dist_vec_seedlings
seedling_ind <- sapply(dist_vec_seedlings, function(x){which.min(abs(x-sim_seedling_pts$real$predicted$x))})

# number of statistics (8 = mean and sd of seeds and seedlings, e.intercept, x.intercept for seeds and seedlings)
num_stats <- 8 + length(dist_vec_seeds)+length(dist_vec_seeds)

# save "observed" stats
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
obs_stats[1,]<- c(mean(sim_seeds_traps),
                  sd(sim_seeds_traps),
                  sim_seeds_pts$real$predicted$y[seed_ind],
                  sim_seeds_pts$real$x.intercept,
                  sim_seeds_pts$real$e.intercept,
                  mean(sim_seedlings),
                  sd(sim_seedlings),
                  sim_seedling_pts$real$predicted$y[seedling_ind],
                  sim_seedling_pts$real$x.intercept,
                  sim_seedling_pts$real$e.intercept)


saveRDS(obs_stats, paste0("./validation/output/knownparam_stats_splinecorr.RDS"))




## Run ABC
#######################################################################
obs_stats <- readRDS("./validation/output/knownparam_stats_splinecorr.RDS")
sim_params <- readRDS("./simulations/output/params_meta_splinecorr.RDS")
sim_stats <- readRDS("./simulations/output/summary_statistics_meta_splinecorr.RDS")

#remove xint,eint (too many NA values)
obs_stats <- obs_stats %>% dplyr::select(-seeds_xint,-seeds_eint, -seedlings_xint, -seedlings_eint)
sim_stats <- sim_stats %>% dplyr::select(-seeds_xint,-seeds_eint, -seedlings_xint, -seedlings_eint)

## ABC step
abc_results <- abc(target=obs_stats,
    param=sim_params,
    sumstat=sim_stats,
    tol=0.01,
    method = "rejection"
    )

## ABC result prep for plotting
# seed posterior
seed_post <- seed_conversion(abc_results$unadj.values)
stat_seed_post <- post_stats(seed_post)

# seed prior
seed_prior <- seed_conversion(sim_params)
seed_dist <- bind_rows(Posterior=seed_post, Prior=seed_prior, .id="dist")


# establishment posterior
establishment_post <- est_conversion(abc_results$unadj.values)
stat_est_post <- post_stats(establishment_post)

# establishment prior
establishment_prior <-est_conversion(sim_params)
est_dist <- bind_rows(Posterior=establishment_post, Prior=establishment_prior, .id="dist")


# summary statistics
seed_sum_stats <- seed_stats_prep(abc_results$ss)
seedling_sum_stats <- seedling_stats_prep(abc_results$ss)

# observed statistics
seed_obs_stats <- seed_stats_prep(obs_stats)
seedling_obs_stats <- seedling_stats_prep(obs_stats)

# known parameters
seed_known_params <- seed_conversion(known_params)
est_known_params <- est_conversion(known_params)

## Plots
# seed distributions
seed_plot <- (ggplot(seed_dist, aes(value, after_stat(density),fill=dist))
              + geom_histogram(alpha=0.6,position='identity',bins=50)
              #+ geom_density(alpha=0.6,position='identity')
              + geom_vline(data=seed_known_params,aes(xintercept = value,colour='known parameter'))
              + geom_vline(data=stat_seed_post, aes(xintercept=mean, colour='posterior mean'))
              + geom_segment(data=stat_seed_post, inherit.aes=FALSE, aes(x=hdi_lower,
                                                                         xend=hdi_upper,
                                                                         y=0,
                                                                         yend=0,colour='posterior HDI'),
                             size=1)
              + facet_wrap(~key, scales = 'free')
              + scale_colour_manual(name = '',values=c('black','#D93802','red'))
              + scale_fill_manual(name='',values=c('#D95F02','#1B9E77'))
              + ylab("")
              + theme_bw()
              + scale_x_log10()
              + ggtitle("Seed distributions with known parameters")
              + theme(legend.position="bottom", axis.title = element_blank())
)

seed_plot

# establishment distributions
est_plot<- (ggplot(est_dist, aes(value, after_stat(density),fill=dist))
            + geom_histogram(alpha=0.6,position='identity',bins=50)
            #+ geom_density(alpha=0.6,position='identity')
            + geom_vline(data=est_known_params,aes(xintercept = value,colour='known parameter'))
            + geom_vline(data=stat_est_post, aes(xintercept=mean, colour='posterior mean'))
            + geom_segment(data=stat_est_post, inherit.aes=FALSE, aes(x=hdi_lower,
                                                                      xend=hdi_upper,
                                                                      y=0,
                                                                      yend=0,colour='posterior HDI'),
                           size=1)
            + facet_wrap(~key, scales = 'free')
            + scale_colour_manual(name = '',values=c('black','#D93802','red'))
            + scale_fill_manual(name='',values=c('#D95F02','#1B9E77'))
            + ylab("")
            + theme_bw()
            + ggtitle("Establishment distributions with known parameters")
            + theme(legend.position="bottom", axis.title = element_blank())
)
est_plot

ggarrange(seed_plot,est_plot,ncol=1)

