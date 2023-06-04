## Validating acceptance criteria using RMSRE
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
source("./general/abc_local.R")


## Read/prep simulations
######################################################################
full_params <- readRDS("./simulations/output/params_meta_splinecorr.RDS")
full_sum_stats <- readRDS("./simulations/output/summary_statistics_meta_splinecorr.RDS")

# drop xint and eint, and scale by hand (divide by sd)
full_sum_stats <- (full_sum_stats 
                   %>% dplyr::select(-seeds_xint,-seeds_eint, -seedlings_xint, -seedlings_eint)
                   %>% mutate(across(.cols=everything(), ~drop(scale(.x,center=FALSE,scale=TRUE))))
)


# join parameters and summary statistics to keep track of parameter
# set groupings
# ignoring the fact that a minimal number of records (< 1%) in sum_stats are NA due to low probability regions
# in the parameter space
sims_full <- cbind(full_params,full_sum_stats)

## Experimental design
######################################################################
# choice of 10 true parameter sets to choose
num_sets = 10
# create batch grouping variable
num_batch <- 100

# split into num_sets parameter sets
sims_full <- (sims_full
              %>% group_by(set = (row_number()-1) %/% (n()/num_sets))
)
# randomly select one record per each parameter set
set.seed(111)
observed_data <- (sims_full 
                  %>% sample_n(1)
)
sims_full <- sims_full %>% anti_join(observed_data)

# randomly select 10 simulations to duplicate to create a total of 1e6
sims_full <- (sims_full 
                 %>% sample_n(1) 
                 %>% union_all(sims_full)
                 %>% arrange(set)
)
# batch sizes
num_batch <- c(1,5,10,50,100,5,25,50,250,500,10,50,100,500,1000)
# acceptance criteria 
ac <- c(rep(0.01,5),rep(0.05,5),rep(0.1,5))
# experimental design
scheme <- (data.frame(num_batch=num_batch,ac=ac)
           %>% mutate(batch_size = 1e5/num_batch,
                      num_abc=ac*batch_size)
)



## Compute RMSRE for each experimental design
######################################################################
rmsre_results <- list()
#std_err_results <- list()
abc_params <- list()
abc_sumstats <- list()
for(i in 1:nrow(scheme)){
  # get scheme values
  num_batch = scheme$num_batch[i]
  ac = scheme$ac[i]
  batch_size = scheme$batch_size[i]
  num_abc = scheme$num_abc[i]
  
  # replicate observed data for the number of batches
  observed_data_reps <- (observed_data
                         %>% slice(rep(1:n(), each = num_batch))  
                         %>% group_by(row_number())
                         
  )
  
  # create list of observed data
  # split observed data into number of batches * parameter sets (each abc run)
  observed_data_reps <- observed_data_reps %>% group_split()
  
  # check on duplicates
  #dups <- updated_sims %>% group_by_all() %>% summarize(n=n()) %>% filter(n>1)
  ## STEPS
  
  
  # for each parameter set
  # split up into num_batch 
  # create batches for each abc run (number of batches times number of parameter sets)
  batches <- (sims_full
              %>% group_by(batch= (row_number()-1) %/% ((n()/num_sets)/num_batch),set)
              %>% group_split()
  )
  
  # split batches into parameters and summary statistics
  sim_params <- batches %>% map(dplyr::select,1:10)
  sim_stats <- batches %>% map(dplyr::select,11:26)
  
  true_params <- observed_data_reps %>% map(dplyr::select,1:10)
  true_stats <- observed_data_reps %>% map(dplyr::select,11:26)
  
  # summarize sets and batch numbers for each batch
  sets_batches <- (batches 
                   %>% map(dplyr::select,27:28)
                   %>% map(unique)
                   #%>% bind_rows()
  )
  
  
  # accepted parameters from abc simulations
  abc_par <- (pmap_dfr(list(true_stats,sim_params,sim_stats,sets_batches),
                        ~ data.frame(abc_local(target=..1,
                                         param=..2,
                                         sumstat = ..3,
                                         tol=ac,
                                         method="rejection")$unadj.values) 
                       %>% mutate(set=..4$set,
                                  batch=..4$batch,
                                  ac=ac,
                                  num_batch=num_batch,
                                  batch_size = batch_size,
                                  num_abc=num_abc)
                        
  ))
  
  # accepted summary stats from abc simulations
  abc_ss <- (pmap_dfr(list(true_stats,sim_params,sim_stats,sets_batches),
                      ~ data.frame(abc_local(target=..1,
                                       param=..2,
                                       sumstat = ..3,
                                       tol=ac,
                                       method="rejection")$ss) 
                      %>% mutate(set=..4$set,
                                 batch=..4$batch,
                                 ac=ac,
                                 num_batch=num_batch,
                                 batch_size = batch_size,
                                 num_abc=num_abc)
                      
  ))
  
  # all_std_err <- (abc_par
  #                 %>% group_by(set, batch)
  #                 #%>% mutate(n=n())
  #                 %>% dplyr::summarize(num_abc=n(),
  #                               across(starts_with(c("seed","est")),sd)
  #                 )
  #                 %>% mutate(across(starts_with(c("seed","est")),~.x/sqrt(num_abc)))
  # )
  
  # split into each abc run (number of batches times the number of parameter sets)
  abc_params_split <- (abc_par  
                       %>% dplyr::select(starts_with(c("seed","est")),set,batch) 
                       %>% group_by(set,batch) 
                       %>% group_split())
  
  # RMSRE computation
  all_rmsre <- lapply(seq_along(abc_params_split),function(k){
    
    rmsre <- setNames(data.frame(matrix(ncol = ncol(true_params[[1]]), nrow = 1)), colnames(true_params[[1]]))
    
    # for each parameter, compute rmsre
    for (j in 1:(ncol(abc_params_split[[k]])-2)){
      # root mean square relative error
      rmsre[,j]<-sqrt((sum(((abc_params_split[[k]][,j]-as.numeric(true_params[[k]][,j]))/(as.numeric(true_params[[k]][,j])))^2))/nrow(abc_params_split[[k]]))
      #sqrt((sum(((h-observed_data$est_Mean[2])/(observed_data$`est_Mean`[2]))^2))/200)
      }
    rmsre
  }) %>% bind_rows() %>% cbind(bind_rows(sets_batches))

  # all_std_err <- lapply(seq_along(all_stats),function(k){
  #   all_stats[[k]][[2]]
  # }) %>% bind_rows() %>% cbind(sets_batches)
  # 
  all_rmsre <- (all_rmsre
                %>% group_by(set)
                %>% dplyr::summarize(across(starts_with(c("seed","est")),mean))
                %>% mutate(num_batch=num_batch,
                           ac=ac,
                           batch_size = batch_size,
                           num_abc = num_abc)
  )
  # all_std_err <- (all_std_err
  #                 %>% group_by(batch)
  #                 %>% dplyr::summarize(across(starts_with(c("seed","est")),mean))
  #                 %>% mutate(num_batch=num_batch,
  #                            ac=ac,
  #                            batch_size = batch_size)
  #                 
  # )

  rmsre_results[[i]] <-all_rmsre
  # std_err_results[[i]]<-all_std_err
  abc_params[[i]]<-abc_par
  abc_sumstats[[i]] <- abc_ss
}
## Save results
rmsre_results <- rmsre_results %>% bind_rows()
#std_err_results <- std_err_results %>% bind_rows()
abc_sumstats <- abc_sumstats %>% bind_rows()
abc_params <- abc_params %>% bind_rows()
saveRDS(rmsre_results,"./validation/output/rmsre_results.RDS")
#saveRDS(std_err_results,"./validation/splinecorr_fullpriors/std_err_results.RDS")


## Read in results
rmsre_results<-readRDS("./validation/output/rmsre_results.RDS")
#std_err_results<-readRDS("./validation/splinecorr_fullpriors/std_err_results.RDS")


# RMSRE averaged across parameter sets and ac combinations
rmsre_ac_set <- (rmsre_results
                 %>% mutate(set=as.factor(set),ac=as.factor(ac),num_abc=as.factor(num_abc))
                 %>% group_by(ac,set)
                 %>% dplyr::summarize(across(starts_with(c("seed","est")),mean))
                 %>% ungroup
                 %>% pivot_longer(cols=starts_with(c("seed","est")),names_to="param",values_to="rmsre")

)
rmsre_ac_set_plot <- ggplot(rmsre_ac_set, aes(x=param,y=rmsre,colour=set,shape=ac))+
  geom_jitter(alpha=0.8)+
  facet_wrap(~param,scales="free",ncol=5)+
  theme_bw()
direct.label(rmsre_ac_set_plot,"last.points")

# RMSRE  averaged across parameter sets and number of accepted abc parameter combinations
rmsre_ac_numabc <- (rmsre_results
                   %>% mutate(set=as.factor(set),ac=as.factor(ac),num_abc=as.factor(num_abc))
                   %>% group_by(ac,num_abc)
                   %>% dplyr::summarize(across(starts_with(c("seed","est")),mean))
                   %>% ungroup
                   %>% pivot_longer(cols=starts_with(c("seed","est")),names_to="param",values_to="rmsre")
                 
)
rmsre_ac_numabc_plot <- ggplot(rmsre_ac_numabc, aes(x=param,y=rmsre,colour=num_abc,shape=ac))+
  geom_jitter(alpha=0.8)+
  facet_wrap(~param,scales="free",ncol=5)+
  theme_bw()
direct.label(rmsre_ac_numabc_plot,"last.points")

# RMSRE averaged across acceptance criterias
rmsre_ac <- (rmsre_results
            %>% mutate(set=as.factor(set),ac=as.factor(ac),num_abc=as.factor(num_abc))
            %>% group_by(ac)
            %>% dplyr::summarize(across(starts_with(c("seed","est")),mean))
            %>% ungroup
            %>% pivot_longer(cols=starts_with(c("seed","est")),names_to="param",values_to="rmsre")
                    
)
rmsre_ac_plot <- ggplot(rmsre_ac, aes(x=param,y=rmsre,colour=ac))+
  geom_jitter(alpha=0.8)+
  facet_wrap(~param,scales="free",ncol=5)+
  theme_bw()
rmsre_ac_plot


