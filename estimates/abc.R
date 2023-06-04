## ABC run with observed data
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
load("./data/pine2.RData")
source("./general/functions.R")
source("./general/abc_local.R")


## Data Prep
######################################################################
## Seeds
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
# create log breaks along spline correlogram
dist_vec_seeds <- exp(seq(log(5),log(50),length.out=6))


## Seedlings
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
# create log breaks along spline correlogram
dist_vec_seedlings <- exp(seq(log(5),log(60),length.out=6))


## ABC
######################################################################
## observed sum stats
obs_stats <- readRDS("./summary_statistics/output/observed_stats.RDS")
## simulations
sim_params <- readRDS("./simulations/output/params_meta_splinecorr.RDS")
sim_stats <- readRDS("./simulations/output/summary_statistics_meta_splinecorr.RDS")

#remove xint,eint
obs_stats <- obs_stats %>% dplyr::select(-seeds_xint,-seeds_eint, -seedlings_xint, -seedlings_eint)
sim_stats <- sim_stats %>% dplyr::select(-seeds_xint,-seeds_eint, -seedlings_xint, -seedlings_eint)

## perform scaling on simulated sum stats
sim_stats <- (sim_stats
              %>% mutate(across(everything(), ~drop(scale.default(.))))
              )
## get scale factors from simulated data
sim_scales <- map_dbl(sim_stats, attr, "scaled:scale")

## scale observed stats
obs_stats <- (obs_stats)/sim_scales

## use local abc function
abc_results <- abc_local(target=obs_stats,
                   param=sim_params,
                   sumstat=sim_stats,
                   tol=0.01,
                   method = "rejection"
)



## Prep data for plotting
######################################################################
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


# remove periods in parameter names
seed_dist <- (seed_dist
              %>% mutate(key=gsub("\\."," ",key,perl=TRUE)))
est_dist <- (est_dist
             %>% mutate(key=gsub("\\."," ",key,perl=TRUE)))
stat_seed_post <- (stat_seed_post
                   %>% mutate(key=gsub("\\."," ",key,perl=TRUE)))
stat_est_post <- (stat_est_post
                  %>% mutate(key=gsub("\\."," ",key,perl=TRUE)))

## Plots
######################################################################
# seed distributions
seed_plot <- (ggplot(seed_dist, aes(value, after_stat(density),fill=dist))
              #+ geom_histogram(alpha=0.6,position='identity',bins=50)
              + geom_density(alpha=0.8,position='identity')
              # + geom_vline(data=seed_known_params,aes(xintercept = value,colour='known parameter'))
              #+ geom_vline(data=stat_seed_post, aes(xintercept=mean, colour='posterior mean'))
              + geom_segment(data=stat_seed_post, inherit.aes=FALSE, aes(x=hdi_lower,
                                                                         xend=hdi_upper,
                                                                         y=0,
                                                                         yend=0,colour='posterior HDI'),
                             size=2,lineend="round",alpha=0.5,guides=FALSE)
              #geom_segment(x=dat.hdi[1],xend=dat.hdi[2],y=0,yend=0,color="blue",size=2,lineend="round")
              + facet_wrap(~key, scales = 'free',nrow=1)
              #+ facet_grid(~key, scales = 'free',labeller=as_labeller(labels))
              + scale_colour_manual(name = '',values=c('#670000'),guide="none")
              + scale_fill_manual(name='',values=c('#670000','darkgrey'))
              + ylab("")
              + theme_bw()
              + ggtitle("Seed dispersal distributions")
              + theme(legend.position="bottom", axis.title = element_blank())
              #+ scale_x_log10()
)

seed_plot

# establishment distributions
est_plot<- (ggplot(est_dist, aes(value, after_stat(density),fill=dist))
            #+ geom_histogram(alpha=0.6,position='identity',bins=50)
            + geom_density(alpha=0.8,position='identity')
            # + geom_vline(data=est_known_params,aes(xintercept = value,colour='known parameter'))
            #+ geom_vline(data=stat_est_post, aes(xintercept=mean, colour='posterior mean'))
            + geom_segment(data=stat_est_post, inherit.aes=FALSE, aes(x=hdi_lower,
                                                                      xend=hdi_upper,
                                                                      y=0,
                                                                      yend=0,colour='posterior HDI'),
                           size=2,lineend="round",alpha=0.5,guides=FALSE)
            + facet_wrap(~key, scales = 'free',nrow=1)
            + scale_colour_manual(name = '',values=c('#670000'),guide="none")
            + scale_fill_manual(name='',values=c('#670000','darkgrey'))
            + ylab("")
            + theme_bw()
            + ggtitle("Plant establishment distributions")
            + theme(legend.position="bottom", axis.title = element_blank())
            #+ scale_x_log10()
)
est_plot
