## Functions used to sample from priors, simulate data, and prep data
#########################################################################

## From: https://github.com/bbolker/bbmisc/blob/master/Rmisc/powexp_prior.R 
## https://cran.r-project.org/web/packages/gnorm/vignettes/gnormUse.html
#' @param lwr lower (loose) bound
#' @param upr upper (loose) bound
#' @param tail_prob prob outside {lwr, upr}
#' @param ctr_prob probability of being in the middle 50\%  of {lwr, upr}
#' @examples
#' get_gnorm()  ## default; gaussian with sd=1
#' get_gnorm(tail_prob=0.01, ctr_prob=0.6)
#' ## set prior approx between 1 and 1000 (log scale)
#' (p <- get_gnorm(lwr=log(1), upr=log(1000), tail_prob=0.01, ctr_prob=0.55))
#' fx <- function(x) do.call("dgnorm", c(list(exp(x)), as.list(p)))
#' curve(fx, from=log(0.5), to=log(1200), n=501)
#' abline(v=c(0,log(1000)), lty=2)
#' brks <- outer(c(1,2,5),c(1,10,100))
#' axis(side=3, at=log(brks), labels=brks)
get_gnorm <- function(lwr=-1, upr=1, tail_prob=2*pnorm(lwr),
                      ctr_prob=abs(diff(pnorm(c(-1,1)*lwr/2)))) {
  require("gnorm")
  ## default tail_prob/ctr_prob assume lwr/upr symmetric around 0 ...
  ## start from Gaussian
  ## desired alpha
  sd <- abs(upr-lwr)/(-2*qnorm(tail_prob/2))
  ## convert to sd (?pgnorm)
  ## conversion factor: sqrt(1/(gamma(3/2)/(gamma(1/2))))
  alpha <- sd*sqrt(2)
  mu <- (lwr+upr)/2 ## symmetric, we don't have to estimate this
  start <- c(alpha=alpha, beta=2)
  tfun <- function(x) {
    ## compute probability within range
    pfun <- function(r) abs(diff(vapply(r,
                                        function(z) do.call("pgnorm",c(list(z, mu=mu), as.list(x))),
                                        FUN.VALUE=numeric(1))))
    tail_obs <- 1-pfun(c(upr,lwr))
    ctr_range <- c((mu+lwr)/2, (mu+upr)/2)
    ctr_obs <- pfun(ctr_range)
    return((tail_prob-tail_obs)^2 + (ctr_prob-ctr_obs)^2)
  }
  return(c(mu=mu,optim(par=start,fn=tfun)$par))
}


if (FALSE) {
  png("dgn.png")
  curve(do.call("dgnorm", c(list(x), as.list(res))), from=-2, to=15, ylab="")
  dev.off()
}

#' Matérn correlation function
#' 
#' @param phi scale parameter
#' @param k shape parameter
#' @param u distance
#' 
corr_matern <- function(phi,k,u){
  ((2^(k-1)*gamma(k))^(-1))*((u/phi)^k)*besselK((u/phi),k)
}

#' log odds function
#' 
logit <- qlogis

#' Seed dispersal priors used in final analysis. 
#' 
#' Draws one sample from each prior, \code{set.seed()} should be used prior to this function
#' call for reproducibility.
#' 
#' @return list containing one parameter set in seed dispersal geostatistical model, in order
#' (Matérn scale, Matérn shape, nugget variance, signal variance, mean, total CV, nugget proportion)
seed_priors <- function(){

  ## seed count mean prior (on log scale)
  mu <- rnorm(n=1,mean=log(10),sd=(log(40)-log(5))/4)
  
  ## Matérn shape prior
  sh <- runif(n=1, min=0.5, max=4)

  ## Matérn reparameterized scale prior
  p_scale<-get_gnorm(lwr=log(3), upr=log(50), tail_prob=0.05, ctr_prob=0.55)
  alpha_log <- rgnorm(1,mu=p_scale[1],alpha=p_scale[2],beta=p_scale[3])
  ## convert to original scale (patch size)
  alpha <- exp(alpha_log)
  ## compute standard Matérn scale parameter (phi)
  sc <- (alpha/(2*sqrt(sh)))
  
  ## nugget proportion prior
  prop_nug_logodds <- rnorm(n=1, mean=logit(0.1), sd=((logit(0.40))-(logit(0.01)))/4)
  ## convert to [0,1] scale
  prop_nug <- exp(prop_nug_logodds)/(1+exp(prop_nug_logodds))
  
  ## total CV prior
  total_cv_log <- rnorm(n=1,mean=log(0.25),sd=(log(0.90)-log(0.1))/4)
  ## convert to original scale
  total_cv <- exp(total_cv_log)
  
  ## solve for nugget and signal variance
  signal_var <- (total_cv)^2*(1-prop_nug)
  nugget_var <- (total_cv)^2*prop_nug
  
  return(as.list(data.frame(sc, sh, nugget_var,signal_var, mu,total_cv,prop_nug)))
}

#' Simulate seedfall
#' 
#' Seed dispersal geostatistical model. \code{set.seed()} should be used prior to this function
#' call for reproducibility.
#' 
#' @param scale Matérn scale (standard parameterization)
#' @param shape Matérn shape
#' @param nugget_var nugget variance on log scale $\tau^2_{log}$
#' @param signal_var signal variance on log scale $\sigma^2_{log}$
#' @param mu mean seed count on log scale
#' @param distance_matrix_list list of distance matrices between each pair of spatial locations (output from \code{stats::dist()})
#' 
#' @return vector of seed counts for each spatial location in distance matrix
simulate_seedfall <- function(scale,shape,nugget_var,signal_var,mu,distance_matrix_list){

  ## compute covariance matrix by plot, at seed locations
  cov_by_plot <- lapply(seq_along(distance_matrix_list),function(i){
    temp_matrix <- apply(as.matrix(distance_matrix_list[[i]]),2,function(x){corr_matern(scale,shape,x)*signal_var})
    diag(temp_matrix) <- signal_var
    temp_matrix
  })
  ## create covariance matrix by combining diagonal blocks
  sigma <-bdiag(cov_by_plot)
  
  ## compute Cholesky of sigma
  U<-chol(sigma%>%as.matrix)
  
  ## draw random normal deviates
  Z <- rnorm(mean=0,sd=1,n=nrow(U))
  N <- rnorm(mean=0,sd=1,n=nrow(U))
  
  ## compute Gaussian random field
  Gaussian_RF <- drop(U%*%Z)
  
  ## seed dispersal geostatistical model
  ## compute vector of Poisson means for each spatial location
  mu_i <- exp(mu + Gaussian_RF + sqrt(nugget_var)*N)

  ## Poisson draws to simulate seeds counts at each spatial location
  sim_seeds<-rpois(n=length(Gaussian_RF),lambda=mu_i)
}

#' Plant establishment priors used in final analysis. 
#' 
#' Draws one sample from each prior, \code{set.seed()} should be used prior to this function
#' call for reproducibility.
#' 
#' @return list containing one parameter set in plant establishment geostatistical model, in order
#' (Matérn scale, Matérn shape, nugget variance, signal variance, mean, total CV, nugget proportion)
establishment_priors <- function(){
  
  ## Matérn shape prior, same as seed prior
  sh <- runif(n=1, min=0.5, max=4)
  
  ## Matérn reparameterized scale prior
  p_scale<-get_gnorm(lwr=log(3), upr=log(50), tail_prob=0.005, ctr_prob=0.55)
  alpha_log <- rgnorm(1,mu=p_scale[1],alpha=p_scale[2],beta=p_scale[3])
  ## convert to original scale (patch size)
  alpha <- exp(alpha_log)
  ## compute standard Matérn scale parameter (phi)
  sc <- (alpha/(2*sqrt(sh)))
  
  ## Mean establishment probability prior (on log odds scale)
  mu <- rnorm(n=1, mean=logit(0.1), sd=((logit(0.5))-(logit(0.01)))/4)
  
  ## total CV prior
  logcv <- rnorm(n=1, mean=log(1), sd=(log(3.5)-log(0.5))/4)
  ## convert to original scale
  total_cv <- exp(logcv)
  
  # nugget proportion, same prior as seeds
  prop_nug_logodds <- rnorm(n=1, mean=logit(0.1), sd=((logit(0.40))-(logit(0.01)))/4)
  ## convert to [0,1] scale
  prop_nug <- exp(prop_nug_logodds)/(1+exp(prop_nug_logodds))
  
  # solve for nugget and signal variance
  signal_var <- (total_cv)^2*(1-prop_nug)
  nugget_var <- (total_cv)^2*prop_nug
  
  return(as.list(data.frame(sc, sh, nugget_var,signal_var, mu,total_cv,prop_nug)))
}

#' Simulate seedling counts
#' 
#' Plant establishment geostatistical model. \code{set.seed()} should be used prior to this function
#' call for reproducibility.
#' 
#' @param scale Matérn scale (standard parameterization)
#' @param shape Matérn shape
#' @param nugget_var nugget variance on log scale $\tau^2_{log}$
#' @param signal_var signal variance on log scale $\sigma^2_{log}$
#' @param mu mean establishment probability on log odds scale
#' @param distance_matrix_list list of distance matrices between each pair of spatial locations (output from \code{stats::dist()})
#' @seed_sims vector of seed counts at each spatial location in \code{distance_matrix_list}
#' 
#' @return vector of seedling counts for each spatial location in distance matrix
simulate_seedlings <- function(scale, shape, mu, signal_var, nugget_var, distance_matrix_list, seed_sims){
  
  ## compute covariance matrix by plot, at seedling locations
  cov_by_plot <- lapply(seq_along(distance_matrix_list),function(i){
    temp_matrix <- apply(as.matrix(distance_matrix_list[[i]]),2,function(x){corr_matern(scale,shape,x)*signal_var})
    diag(temp_matrix) <- signal_var
    temp_matrix
  })
  ## create covariance matrix by combining diagonal blocks
  sigma <-bdiag(cov_by_plot)
  
  ## compute Cholesky of sigma
  U<-chol(sigma%>%as.matrix)
  
  ## draw random normal deviates
  Z <- rnorm(mean=0,sd=1,n=nrow(U))
  N <- rnorm(mean=0,sd=1,n=nrow(U))
  
  ## compute Gaussian random field
  Gaussian_RF <- drop(U%*%Z)

  ## plant establishment geostatistical model
  ## compute vector of log odds establishment probability means for each spatial location
  logit_prob <- mu + Gaussian_RF + sqrt(nugget_var)*N
  ## convert to probability scale
  p_i <- plogis(logit_prob)
  
  ## Binomial draws to simulate seedling counts at each spatial location
  sim_seedlings<-rbinom(n=length(Gaussian_RF), size = seed_sims, prob=p_i)
}



#' Simulations of Seed and Seedlings
#' 
#' Simulating seed and seedling count data using geostatistical models and priors for seed dispersal and plant establishment .
#' The output is expected to be used in Approximate Bayesian Computation.
#' 
#' @param seed_base base random seed number
#' @param i simulation number (iteration)
#' @param num_params total number of parameters from both geostatistical models (seed dispersal and plant establishment)
#' @param num_stats_sp total number of summary statistics to be used in ABC
#' @param dist_vec_seeds vector of distances at which to record autocorrelations on the spline correlogram for the 
#' simulated seed data. Autocorrelations are approximate, \code{ncf::spline.correlog()} evaluates the spline correlogram
#' at 300 points by default. The recorded autocorrelations are the spline correlogram values corresponding to the 
#' closest distances of the 300 points to those in `dist_vec_seeds`
#' @param dist_vec_seedlings vector of distances at which to record autocorrelations on the spline correlogram for the 
#' simulated seedling data. The recorded autocorrelations are approximate, as in `dist_vec_seeds` explanation.
#' @param seeds_u list of distance matrices between each pair of spatial locations for seed traps (output from \code{stats::dist()})
#' @param seed_locations data frame of seed trap locations with fields `x` and `y`
#' @param seedlings_u list of distance matrices between each pair of spatial locations for seed traps (output from \code{stats::dist()})
#' @param seedling_locations data frame of seedling quadrat locations with fields `x` and `y`
#' @param seed_prior_function function that returns seed dispersal model parameters after one draw from seed dispersal priors.
#' Expects output named as in output from \code{seed_priors()}
#' @param est_prior_function function that returns plant establishment model parameters after one draw from plant establishment priors
#' Expects output named as in output from \code{establishment_priors()}
#' 
#' @return data frame of parameters and associated summary statistics from simulated data for plant recruitment model
abc_simulations_splinecorr <- function(seed_base,
                                  i,
                                  num_params,
                                  num_stats_sp,
                                  dist_vec_seeds,
                                  dist_vec_seedlings,
                                  seeds_u, 
                                  seed_locations, 
                                  seedlings_u, 
                                  seedling_locations,
                                  seed_prior_function,
                                  est_prior_function){
  
  ## initialize vectors to store autocorrelation values from spline correlogram
  num_seed_dist <- length(dist_vec_seeds)
  num_seedling_dist <- length(dist_vec_seedlings)
  seed_pts <- vector(length=num_seed_dist)
  seedling_pts <- vector(length=num_seedling_dist)
  ## add field names
  for (j in 1:num_seed_dist){seed_pts[j]<-paste0("seeds_pt_",j)}
  for (k in 1:num_seedling_dist){seedling_pts[k]<-paste0("seedlings_pt_",k)}
  
  ## initialize data frame to save summary statistics of simulated data
  sum_stats_pts <- (matrix(nrow=1,ncol=num_stats_sp)
                    %>% data.frame()
                    %>% setNames(c("seeds_mean_pts","seeds_sd_pts",seed_pts,"seeds_xint","seeds_eint",
                                   "seedlings_mean_pts","seedlings_sd_pts",seedling_pts, "seedlings_xint","seedlings_eint"))
  )
  
  ## initialize data frame to save parameters from priors
  params <- (matrix(nrow=1,ncol=num_params)
             %>% data.frame()
             %>% setNames(c("seed_Matern Scale","seed_Matern Shape","seed_Total CV",
                            "seed_Nugget Proportion","seed_Mean","est_Matern Scale",
                            "est_Matern Shape","est_Total CV","est_Nugget Proportion",
                            "est_Mean"))
  )
  
  ## set seed
  set.seed(seed_base + i)
  ## sample from seed priors
  sp <- seed_prior_function()
  ## save seed parameters
  params[1:5] <- cbind(sp$sc,sp$sh,sp$total_cv, sp$prop_nug,sp$mu)
  
  ## simulate seedfall (at seed trap locations)
  sim_seeds_traps <- tryCatch(simulate_seedfall(sp$sc,sp$sh,sp$nugget_var,sp$signal_var,sp$mu,seeds_u),
                              error = function(e) e)
  ## if simulation was successful (seed counts were generated)
  if (!is.list(sim_seeds_traps)){
    ## compute spline correlogram for simulated seed data
    sim_seeds_pts <- spline.correlog(x=seed_locations$x,
                                     y=seed_locations$y,
                                     z=sim_seeds_traps,
                                     xmax=60,
                                     resamp=0,
                                     save=FALSE,
                                     df=8)
    ## compute closest indices along spline correlogram to dist_vec_seeds
    seed_ind <- sapply(dist_vec_seeds, function(x){which.min(abs(x-sim_seeds_pts$real$predicted$x))})
    
    ## save seed summary statistics
    ## rows are NA if simulation fails
    sum_stats_pts[1:(num_seed_dist+4)] <- c(mean(sim_seeds_traps),
                                            sd(sim_seeds_traps),
                                            sim_seeds_pts$real$predicted$y[seed_ind],
                                            sim_seeds_pts$real$x.intercept,
                                            sim_seeds_pts$real$e.intercept)
  }
  
  ## simulate seedfall (at seedling locations) using same seed prior draws
  sim_seeds_land <- tryCatch(simulate_seedfall(sp$sc,sp$sh,sp$nugget_var,sp$signal_var,sp$mu,seedlings_u),
                             error = function(e) e)
  
  ## if simulation was successful (seed counts were generated)
  if (!is.list(sim_seeds_land)){
    ## sample from establishment priors
    ep <- est_prior_function()
    ## save establishment  parameters
    params[6:10] <- cbind(ep$sc,ep$sh,ep$total_cv,ep$prop_nug,ep$mu)
    
    ## simulate seedlings at seedling locations
    sim_seedlings <- tryCatch(simulate_seedlings(ep$sc,ep$sh,ep$mu,ep$signal_var,ep$nugget_var,
                                                 seedlings_u,
                                                 sim_seeds_land),
                              error = function(e) e)
    
    ## if simulation was successful (seedling counts were generated)
    if (!is.list(sim_seedlings)){
      ## compute spline correlogram for simulated seedling data
      sim_seedling_pts <- spline.correlog(x=seedling_locations$x,
                                          y=seedling_locations$y,
                                          z=sim_seedlings,
                                          xmax=70,
                                          resamp=0,
                                          save=FALSE,
                                          df=8)
      ## compute closest indices along spline correlogram to dist_vec_seedlings
      seedling_ind <- sapply(dist_vec_seedlings, function(x){which.min(abs(x-sim_seedling_pts$real$predicted$x))})
      
      ## save seedling summary statistics
      ## rows are NA if simulation fails
      sum_stats_pts[(num_seed_dist+5):num_stats_sp] <- c(mean(sim_seedlings),
                                                         sd(sim_seedlings),
                                                         sim_seedling_pts$real$predicted$y[seedling_ind],
                                                         sim_seedling_pts$real$x.intercept,
                                                         sim_seedling_pts$real$e.intercept)
    }
    
    
  }
  return(cbind(sum_stats_pts,params))
  
}


#' Simulations of Seed and Seedlings (Old Version)
#' 
#' Simulating seed and seedling count data using geostatistical models and priors for seed dispersal and plant establishment .
#' The output is expected to be used in Approximate Bayesian Computation. Old version using variogram instead of spline correlogram.
#' 
#' @param seed_base base random seed number
#' @param i simulation number (iteration)
#' @param num_params total number of parameters from both geostatistical models (seed dispersal and plant establishment)
#' @param num_stats_vario total number of summary statistics to be used in ABC
#' @param bin_lims_seeds vector of distances bin breaks at which to bin the seed variogram
#' @param bin_lims_seedlings  vector of distances bin breaks at which to bin the seedling variogram
#' @param seeds_u list of distance matrices between each pair of spatial locations for seed traps (output from \code{stats::dist()})
#' @param seed_locations data frame of seed trap locations with fields `x` and `y`
#' @param seedlings_u list of distance matrices between each pair of spatial locations for seed traps (output from \code{stats::dist()})
#' @param seedling_locations data frame of seedling quadrat locations with fields `x` and `y`
#' @param seed_prior_function function that returns seed dispersal model parameters after one draw from seed dispersal priors.
#' Expects output named as in output from \code{seed_priors()}
#' @param est_prior_function function that returns plant establishment model parameters after one draw from plant establishment priors
#' Expects output named as in output from \code{establishment_priors()}
#' 
#' @return data frame of parameters and associated summary statistics from simulated data for plant recruitment model
abc_simulations_vario <- function(seed_base,
                                  i,
                                  num_params,
                                  num_stats_vario,
                                  bin_lims_seeds,
                                  bin_lims_seedlings, 
                                  seeds_u, 
                                  seed_locations, 
                                  seedlings_u, 
                                  seedling_locations,
                                  seed_prior_function,
                                  est_prior_function){

  # adjust bin breaks to number of variogram bins
  num_seed_bins <- length(bin_lims_seeds)-2
  num_seedling_bins <- length(bin_lims_seedlings)-1
  seed_bins <- vector(length=num_seed_bins)
  seedling_bins <- vector(length=num_seedling_bins)
  
  for (j in 1:num_seed_bins){seed_bins[j]<-paste0("seeds_bin_",j)}
  for (k in 1:num_seedling_bins){seedling_bins[k]<-paste0("seedlings_bin_",k)}

  
  sum_stats_bins <- (matrix(nrow=1,ncol=num_stats_vario)
                     %>% data.frame()
                     %>% setNames(c("seeds_mean","seeds_sd",seed_bins,
                                    "seedlings_mean","seedlings_sd",seedling_bins))
  )

  params <- (matrix(nrow=1,ncol=num_params)
             %>% data.frame()
             %>% setNames(c("seed_Matern Scale","seed_Matern Shape","seed_Total CV",
                            "seed_Nugget Proportion","seed_Mean","est_Matern Scale",
                            "est_Matern Shape","est_Total CV","est_Nugget Proportion",
                            "est_Mean"))
  )
  
  ## set seed
  set.seed(seed_base + i)
  ## sample from seed priors
  sp <- seed_prior_function()
  ## save seed parameters
  params[1:5] <- cbind(sp$sc,sp$sh,sp$total_cv, sp$prop_nug,sp$mu)

  # simulate seedfall (at seed trap locations)
  sim_seeds_traps <- tryCatch(simulate_seedfall(sp$sc,sp$sh,sp$nugget_var,sp$signal_var,sp$mu,seeds_u),
                              error = function(e) e)
  if (!is.list(sim_seeds_traps)){
    sim_seeds_bins <- variog(coords = seed_locations,
                             data = sim_seeds_traps,
                             option = "bin",
                             breaks=bin_lims_seeds,
                             messages = FALSE)

    sum_stats_bins[1:(num_seed_bins+2)] <- c(mean(sim_seeds_traps),sd(sim_seeds_traps),sim_seeds_bins$v)

  }
  
  # simulate seedfall (at seedling locations)
  sim_seeds_land <- tryCatch(simulate_seedfall(sp$sc,sp$sh,sp$nugget_var,sp$signal_var,sp$mu,seedlings_u),
                             error = function(e) e)
  
  if (!is.list(sim_seeds_land)){
    ## sample from establishment priors
    ep <- est_prior_function()
    ## save establishment  parameters
    params[6:10] <- cbind(ep$sc,ep$sh,ep$total_cv,ep$prop_nug,ep$mu)
    # simulate seedlings
    sim_seedlings <- tryCatch(simulate_seedlings(ep$sc,ep$sh,ep$mu,ep$signal_var,ep$nugget_var,
                                                 seedlings_u,
                                                 sim_seeds_land),
                              error = function(e) e)
    if (!is.list(sim_seedlings)){
      # compute variogram bins
      sim_seedling_bins <- variog(coords = seedling_locations,
                                  data = sim_seedlings,
                                  option = "bin",
                                  breaks=bin_lims_seedlings,
                                  messages = FALSE)

      
      # rows are NA if simulation fails
      sum_stats_bins[(num_seed_bins+3):num_stats_vario] <- c(mean(sim_seedlings),
                                                             sd(sim_seedlings),
                                                             sim_seedling_bins$v) 
    }
    
    
  }
  return(cbind(sum_stats_bins,params))
  
}


## Data Prep Functions #################

#' Seed Dispersal Parameter Prep
#' 
#' Convert parameters to long format for plotting, and transform some parameter back to original (ecological) scale.
#' @param x matrix of seed dispersal parameters (ex. output from ABC run)
seed_conversion <- function(x){
  (x
   %>% data.frame()
   %>% dplyr::select(starts_with("seed"))
   %>% rename_with(~ gsub("seed_","",.x,perl=TRUE))
   %>% mutate(Mean = exp(Mean))
   %>% mutate(Matern.Scale = 2*Matern.Scale*sqrt(Matern.Shape))
   #%>% mutate(Nugget.Proportion = oob_censor(Nugget.Proportion, range=c(0,1)))
   %>% tidyr::gather())
}

#' Plant Establishment Parameter Prep
#' 
#' Convert parameters to long format for plotting, and transform some parameter back to original (ecological) scale.
#' @param x matrix of plant establishment parameters (ex. output from ABC run)
est_conversion <- function(x){
  (x
   %>% data.frame() 
   %>% dplyr::select(starts_with("est")) 
   %>% rename_with(~ gsub("est_","",.x,perl=TRUE))
   %>% mutate(Mean = plogis(Mean))
   %>% mutate(Matern.Scale = 2*Matern.Scale*sqrt(Matern.Shape))
   %>% tidyr::gather())
}

#' Seed Count Summary Statistic Prep
#' 
#' Convert summary statistics to long format for plotting, and tidy field names.
#' @param x matrix of seed summary statistics (ex. output from ABC run)
seed_stats_prep <- function(x){
  (x
   %>% data.frame() 
   %>% dplyr::select(starts_with("seeds")) 
   %>% rename_with(~ gsub("seeds_","",.x,perl=TRUE))
   %>% tidyr::gather())
}

#' Seedling Count Summary Statistic Prep
#' 
#' Convert summary statistics to long format for plotting, and tidy field names.
#' @param x matrix of seedling summary statistics (ex. output from ABC run)
seedling_stats_prep <- function(x){
  (x
   %>% data.frame() 
   %>% dplyr::select(starts_with("seedlings")) 
   %>% rename_with(~ gsub("seedlings_","",.x,perl=TRUE))
   %>% tidyr::gather())
}

#' Posterior Summary Statistics
#' 
#' Computes mean, median, 2.5% quantiles, 97.5% quantiles, and 95% HDI intervals
#' @param x data frame of posterior parameters in long format, parameter names in `key` field and 
#' parameter value in `value` field. Typically used output from \code{seed_conversion()} or \code{est_conversion()}
post_stats <- function(x){
  (x
   %>% group_by(key) 
   %>% dplyr::summarize(
                        mean = mean(value),
                        median = median(value),
                        q_lower = quantile(value,probs=0.025),
                        q_upper = quantile(value,probs=0.975),
                        hdi_lower = hdi(value,credMass=0.95)[1],
                        hdi_upper = hdi(value,credMass=0.95)[2]
                        ) 
  )
}

