## Validating acceptance criteria using Coverage
######################################################################

## Load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(abc)
library(HDInterval)
library(stringr)
source("./general/abc_local.R")

## Read/prep simulations
######################################################################
full_params <- readRDS("./simulations/output/params_meta_splinecorr.RDS")
full_sum_stats <- readRDS("./simulations/output/summary_statistics_meta_splinecorr.RDS")

# drop xint and eint, and scale by hand (divide by sd)
full_sum_stats <- (full_sum_stats 
                   %>% dplyr::select(-seeds_xint,-seeds_eint, -seedlings_xint, -seedlings_eint)
                   %>% mutate(across(.cols=everything(), ~drop(scale.default(.x))))
)

seed_base <- 111
check_coverage <- list()
cat_plot <- list()

# compute coverage for 500 "observed" data sets
for (i in 1:500){
  # read data
  all_sims <- cbind(full_params,full_sum_stats)
  
  # randomly select one record as observed data
  seed_val <-seed_base + i
  set.seed(seed_val)
  rand_row <- sample(nrow(all_sims %>% drop_na()),1)
  observed_data <-(all_sims %>% drop_na())[rand_row,]
  

  # print iteration number
  if (i%%50==0){
    cat(i,"\n")
  }

  # compute abc for each group of summary statistics
  abc_sims <- abc_local(target=observed_data[1,11:26],
                  param = all_sims[,1:10],
                  sumstat = all_sims[,11:26],
                  tol=0.01,
                  method="rejection")
  
  # ABC posteriors
  accepted_params <- (abc_sims$unadj.values 
                      %>% data.frame 
                      %>% pivot_longer(cols=everything(), names_to="param",values_to="value")
                      %>% group_by(param)
                      %>% dplyr::summarize(
                        hdi_lower = HDInterval::hdi(value,credMass=0.90)[1],
                        hdi_upper = HDInterval::hdi(value,credMass=0.90)[2],
                        mean=mean(value))
  )
  # true parameters
  obs_param <- (observed_data[,1:10] 
                %>% pivot_longer(cols=everything(), names_to="param",values_to="value")
                %>% mutate(across('param',str_replace," ","."))
  )
  cat_plot[[i]] <-(accepted_params
                             %>% left_join(obs_param, by="param") 
                             %>% mutate(mean=mean-value,
                                        hdi_lower=hdi_lower-value,
                                        hdi_upper=hdi_upper-value)
  )
  check_coverage[[i]] <- (accepted_params
                           %>% left_join(obs_param, by="param")
                           %>% rename(true_param=value)
                           %>% ungroup
                           %>% mutate(inside=if_else(true_param >= hdi_lower & true_param <= hdi_upper,1,0))
  )
}
## Save results
saveRDS(check_coverage,"./validation/output/check_coverage.RDS")
saveRDS(cat_plot, "./validation/output/cat_plot.RDS")

## Read results

cat_plot <- readRDS("./validation/output/cat_plot.RDS")
check_coverage <- readRDS("./validation/output/check_coverage.RDS")

## Caterpillar Plot
###########################################################
## HDI credible intervals for each of the 500 runs (in grey)
## Posterior mean - true parameter (in black)
## Ordered in increasing order
results_cat_plot <- (bind_rows(cat_plot) 
                               %>% arrange(mean)
                               %>% group_by(param)
                               %>% mutate(ind=row_number())
                               %>% ungroup())
ggplot(results_cat_plot, aes(x=ind,y=mean))+
  geom_linerange(aes(ymin = hdi_lower, ymax = hdi_upper),colour="grey")+
  geom_hline(yintercept=0)+
  geom_point(alpha=0.3)+
  facet_wrap(~param,scales="free",ncol=5)+
  theme_bw()

## Coverage (90%)
###########################################################
coverage_prob <- (results_cat_plot
                  %>% group_by(param)
                  %>% mutate(n_records=n())
                  %>% ungroup()
                  %>% filter(hdi_lower < 0 & hdi_upper > 0)
                  %>% group_by(param)
                  %>% dplyr::summarize(coverage_prob=n()/n_records)
                  %>% unique()
                  )
results <- (bind_rows(check_coverage)
                      %>% drop_na()
                      %>% group_by(param)
                      %>% dplyr::summarize(percent=sum(inside, na.rm=TRUE)/n())
                      %>% mutate(cov_per=0.9)
)
b_test <- binom.test(0.9*500, 500, p = 0.9)
results <- (results 
            %>% mutate(lower_ci=b_test$conf.int[1],upper_ci=b_test$conf.int[2])
            %>% mutate(param=gsub("est_","Establishment ",param,perl=TRUE),
                       param=gsub("seed_","Seed ",param,perl=TRUE)
                       ,param=gsub("\\."," ",param,perl=TRUE)
                       )

)

ggplot(results, aes(y=param,x=percent))+
  annotate("rect", xmin=b_test$conf.int[1], xmax=b_test$conf.int[2], ymin=-Inf, ymax=Inf, alpha=0.2, fill="blue")+
  #geom_rect(xmin=b_test$conf.int[1],xmax=b_test$conf.int[2],ymin=-Inf,ymax=Inf,alpha=0.3)+
  geom_point()+
  geom_vline(aes(xintercept=cov_per),colour="blue")+
  #geom_ribbon(aes(xmin=lower_ci,xmax=upper_ci,x=percent),colour="blue")+
  #geom_linerange(aes(xmin = lower_ci, xmax = upper_ci))+
  #facet_wrap(~param,scales="free_x",ncol=5)+
  ggtitle("Coverage Probability for Parameters")+
  xlab("%")+
  ylab("")+
  theme_bw()
