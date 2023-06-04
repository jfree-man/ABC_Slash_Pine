## Visualize SBC results from sbc.R META-Farm run
######################################################################


######################################################################
## Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

## SBC results using default abc::abc scaling (scaling by MAD)
######################################################################
rankstat <- (list.files(path = "./validation/output/sbc/",
                        pattern="*_[0-9]+\\.RDS$",
                        full.names = TRUE)
             %>% map_dfr(readRDS)
)

saveRDS(rankstat,file="./validation/output/sbc/rankstat.RDS")
rankstat <- readRDS("./validation/output/sbc/rankstat.RDS")

## 
data_to_plot <- (rankstat
                 %>% pivot_longer(cols=everything(),
                                  names_to="param",
                                  values_to="value")
                 )

# b_test <- binom.test(x =  round((1/(100+1))*1000), n = 1000, p = (1/(100+1)))
# lower_ci=b_test$conf.int[1]
# upper_ci=b_test$conf.int[2]

# qbinom(0.995, size=1000, p = (1/(100+1)))
# plot(dbinom(0:10, size  = 1000, prob = 1/101))
# qbinom(c(0.005, 0.995), size = 1000, prob = 0.01)


qbinom(c(0.005, 0.995), size = 1000, prob = 1/51)
## N=1000 runs

# N=100000
# L=100
# B=50


# 51 bins, 100000 rank statistics
bw=2
ci_band <- qbinom(c(0.005, 0.995), size = 100000, prob = 1/51)

ggplot(data_to_plot, aes(x=value))+
  geom_histogram(binwidth=bw)+
  geom_density(aes(y=bw * after_stat(count)))+
  annotate("rect", ymin=ci_band[1], ymax=ci_band[2], xmin=-Inf, xmax=Inf, alpha=0.2, fill="#670000")+
  facet_wrap(~param, nrow = 2, ncol=5)+
  theme_bw()


## SBC results using default abc::abc scaling (scaling by MAD), and
## lowered acceptance rate from 1% to 0.5%
######################################################################
rankstat_ac <- (list.files(path = "./validation/output/sbc_ac/",
                        pattern="*_[0-9]+\\.RDS$",
                        full.names = TRUE)
             %>% map_dfr(readRDS)
)

saveRDS(rankstat_ac,file="./validation/output/sbc_ac/rankstat.RDS")
rankstat_ac <- readRDS("./validation/output/sbc_ac/rankstat.RDS")

## 
data_to_plot_ac <- (rankstat_ac
                 %>% pivot_longer(cols=everything(),
                                  names_to="param",
                                  values_to="value")
)


bw=2
ci_band <- qbinom(c(0.005, 0.995), size = 100000, prob = 1/51)

ggplot(data_to_plot_ac, aes(x=value))+
  geom_histogram(binwidth=bw)+
  geom_density(aes(y=bw * after_stat(count)))+
  annotate("rect", ymin=ci_band[1], ymax=ci_band[2], xmin=-Inf, xmax=Inf, alpha=0.2, fill="#670000")+
  facet_wrap(~param, nrow = 2, ncol=5)+
  theme_bw()

## SBC results using no scaling abc_local function
######################################################################
rankstat_ns <- (list.files(path = "./validation/output/sbc_noscale/",
                           pattern="*_[0-9]+\\.RDS$",
                           full.names = TRUE)
                %>% map_dfr(readRDS)
)

saveRDS(rankstat_ns,file="./validation/output/sbc_noscale/rankstat.RDS")
rankstat_ns <- readRDS("./validation/output/sbc_noscale/rankstat.RDS")

## 
data_to_plot_ns <- (rankstat_ns
                    %>% pivot_longer(cols=everything(),
                                     names_to="param",
                                     values_to="value")
)


bw=2
ci_band <- qbinom(c(0.005, 0.995), size = 100000, prob = 1/51)

ggplot(data_to_plot_ns, aes(x=value))+
  geom_histogram(binwidth=bw)+
  geom_density(aes(y=bw * after_stat(count)))+
  annotate("rect", ymin=ci_band[1], ymax=ci_band[2], xmin=-Inf, xmax=Inf, alpha=0.2, fill="#670000")+
  facet_wrap(~param, nrow = 2, ncol=5)+
  theme_bw()


## SBC results using scale.default before abc_local function 
## (subtract mean, divide by sd)
######################################################################
rankstat_sd <- (list.files(path = "./validation/output/sbc_scale.default/",
                           pattern="*_[0-9]+\\.RDS$",
                           full.names = TRUE)
                %>% map_dfr(readRDS)
)

saveRDS(rankstat_sd,file="./validation/output/sbc_scale.default/rankstat.RDS")
rankstat_sd <- readRDS("./validation/output/sbc_scale.default/rankstat.RDS")

## 
data_to_plot_sd <- (rankstat_sd
                    %>% pivot_longer(cols=everything(),
                                     names_to="param",
                                     values_to="value")
)


bw=2
ci_band <- qbinom(c(0.005, 0.995), size = 100000, prob = 1/51)

ggplot(data_to_plot_sd, aes(x=value))+
  geom_histogram(binwidth=bw)+
  geom_density(aes(y=bw * after_stat(count)))+
  annotate("rect", ymin=ci_band[1], ymax=ci_band[2], xmin=-Inf, xmax=Inf, alpha=0.2, fill="#670000")+
  facet_wrap(~param, nrow = 2, ncol=5)+
  theme_bw()



## sbc with ss scaling, no centering (divide by sd only)
## SBC results using scaling no centering before abc_local function 
## (divide by sd only) - centering gets cancelled out in ABC distance
## step, should be identical to previous plot
######################################################################
rankstat_nc <- (list.files(path = "./validation/output/sbc_scaleonly/",
                           pattern="*_[0-9]+\\.RDS$",
                           full.names = TRUE)
                %>% map_dfr(readRDS)
)

saveRDS(rankstat_nc,file="./validation/output/sbc_scaleonly/rankstat.RDS")
rankstat_nc <- readRDS("./validation/output/sbc_scaleonly/rankstat.RDS")

## 
data_to_plot_nc <- (rankstat_nc
                    %>% pivot_longer(cols=everything(),
                                     names_to="param",
                                     values_to="value")
)


bw=2
ci_band <- qbinom(c(0.005, 0.995), size = 100000, prob = 1/51)

ggplot(data_to_plot_nc, aes(x=value))+
  geom_histogram(binwidth=bw)+
  geom_density(aes(y=bw * after_stat(count)))+
  annotate("rect", ymin=ci_band[1], ymax=ci_band[2], xmin=-Inf, xmax=Inf, alpha=0.2, fill="#670000")+
  facet_wrap(~param, nrow = 2, ncol=5)+
  theme_bw()



