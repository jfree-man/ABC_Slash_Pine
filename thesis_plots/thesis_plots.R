## All plots to go in thesis (in order as they appear)
######################################################################

## Load packages, functions
load("./data/pine2.RData")
source("./general/functions.R")
source("./general/abc_local.R")
library(ggplot2)
library(dplyr)
library(ggforce)
library(tikzDevice)
library(latex2exp)
library(purrr)
library(tidyr)
library(HDInterval)
library(ggpubr)
library(RColorBrewer)
library(ncf)
library(RandomFields)
library(TMB)
library(glmmTMB)
library(forcats)
file_location <- "./thesis_plots/output"


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





## 1. Plot of data set ###############################################
## best to keep in one plot object so the count scale is the same
## is there a way to create two plots using the same count scale,
## I recall having dificulty with this before
## add rectangles around plots with labels? (A-H)
## ggforce::facet_zoom doesn't seem to work when using facet_wrap
## be nice to have seeds and seedlings faceted...
## not sure how to label axis, lat long?
######################################################################

all_dat <- ((seeds %>% mutate(type="Seeds"))
            %>% union((seedlings%>% mutate(type="Seedlings")))
            %>% mutate(type=factor(type, levels=c("Seeds","Seedlings")))
)
tikz(file.path(file_location,"dataset.tex"),width = 6, height = 2.5)
ggplot(all_dat, aes(x=x,y=y,size=count,colour=type))+
  geom_point(alpha=0.7)+
  scale_size_area()+
  #scale_size_continuous(range=c(0, 10))+
  scale_colour_manual(name = '',values=c("#7B3F00","#186218"))+
  #name = '',values=c('black','#670000')
  facet_wrap(~type)+
  facet_zoom(xy= plot=="C", zoom.size = 1)+
  #facet_zoom(xlim=c(500,600),ylim=c(200,300))
  xlab("\nx-coordinate (m)")+
  ylab("y-coordinate (m)\n")+
  theme_bw()+
  theme(aspect.ratio = 1,legend.text =element_text(size=10), legend.title=element_text(size=10), legend.spacing.y=unit(0.1,"cm"))+
  #ggtitle("Data Set")+
  labs(size="Count")
dev.off()


## 2. Matérn Correlation function  ##################################
## demonstrate difference in reparameterization
######################################################################

## shape parameters
sh <- c(1.5, 0.5)
## scale parameters
sc <- c(2,10)
## compute reparameterized scale alpha
sc_sh <- (cbind(sc,sh)
          %>% data.frame()
          %>% mutate(alpha=2*sc*sqrt(sh))

)

params <- (sc_sh
           %>% mutate(`(scale, alpha, shape)` = paste0("(",sc,", ",round(alpha,1),", ",sh,")"))
           %>% dplyr::select(`(scale, alpha, shape)`)
           %>% pull()
           
)
## distances
u <- seq(0.1,100,by=0.01)
## compute correlation
corr_u<-apply(sc_sh,1,function(x){corr_matern(x[1],x[2],u)})
dat <- data.frame(cbind(u,corr_u))
names(dat) <- c("u",params)


dat <- (dat
        %>% pivot_longer(cols=-u,names_to="(scale, alpha, shape)", values_to="correlation")
)

sc_sh <- (sc_sh %>% mutate(`(scale, alpha, shape)` = paste0("(",sc,", ",alpha,", ",sh,")"))
          %>% mutate(`(scale, alpha_round, shape)` = paste0("(",sc,", ",round(alpha,1),", ",sh,")"))
)

## interpolate correlation at standard scale parameter (phi)
## and reparamterized scale parameter (alpha)
corr_phi <- list()
corr_alpha <- list()
for (i in 1:(nrow(sc_sh))){
  line <- (dat 
           %>% filter(`(scale, alpha, shape)`==sc_sh$`(scale, alpha_round, shape)`[i])
           %>% unique()
  )
  corr_alpha[[i]] <- approx(x=line$u,
                            y=line$correlation,
                            xout=sc_sh$alpha[i])
  corr_phi[[i]] <- approx(x=line$u,
                          y=line$correlation,
                          xout=sc_sh$sc[i])
}

corr_phi <- bind_rows(corr_phi)
corr_alpha <- bind_rows(corr_alpha)

corr_phi_segments <- (corr_phi
                      %>% mutate(xstart=x,ystart=-Inf,
                                 label=c("$\\phi_1$","$\\phi_2$"))
)
corr_alpha_segments <- (corr_alpha
                        %>% mutate(xstart=x,ystart=-Inf,
                                   label=c("$\\alpha_1$","$\\alpha_2$"))
)

dat$`(scale, alpha, shape)`<-factor(dat$`(scale, alpha, shape)`,levels=rev(sort(unique(dat$`(scale, alpha, shape)`))))
tikz(file.path(file_location,"matern_reparam.tex"),width = 5.5, height = 3)
(ggplot(dat, aes(u,correlation, colour = `(scale, alpha, shape)`))
  + geom_line(lwd=1)
  + theme_bw()
  + geom_point(data=corr_phi, aes(x=x,y=y, shape="Matérn Scale $\\phi$"),inherit.aes=FALSE)
  + geom_point(data=corr_alpha, aes(x=x,y=y, shape="Reparameterized \nMatérn Scale $\\alpha$"),inherit.aes=FALSE)
  + geom_segment(data=corr_phi_segments, aes(x=xstart,y=ystart, xend=x,yend=y),inherit.aes=FALSE, lty="dashed", lineend = "round")
  + geom_segment(data=corr_alpha_segments, aes(x=xstart,y=ystart, xend=x,yend=y),inherit.aes=FALSE, lty="dotted", lineend = "round")
  #+ geom_label(data=corr_phi_segments, aes(x=xstart,y=ystart,label=label),inherit.aes=FALSE)
  + scale_x_continuous(breaks=c(0,corr_phi_segments$xstart[1],corr_alpha_segments$xstart[1],corr_phi_segments$xstart[2],corr_alpha_segments$xstart[2],25,50),
                       labels=c("0","$\\phi_1$","$\\alpha_1$","$\\phi_2$","$\\alpha_2$", "25","50"),limits=c(0,50))
  + labs(# not actually plotting practical range here
    x="distance(m)"
    , colour="$(\\phi, \\alpha, \\kappa)$"
    ,shape = "")
  + guides(shape = guide_legend(order=2), colour = guide_legend(order=1))
  + ylab("correlation\n")
  + scale_color_brewer(palette = "Dark2")
  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)
dev.off()



## 3. GRF and Matern shape #############################################
## Plot of simulated GRFs for various shape parameters(over regular grid)
## 
##
## FIX: have to manually edit tikz output tex file to insert "../thesis_plots/output/"
## so pgfimage can find png raster images from location of tex document (must be a better way)
## In grf_shape.tex, make these changes:
## (line 86) \pgfimage[width= 95.30pt,height= 78.73pt,interpolate=false]{../thesis_plots/output/grf_shape_ras1}};
## (line 164) \pgfimage[width= 95.30pt,height= 78.73pt,interpolate=false]{../thesis_plots/output/grf_shape_ras2}};
## (line 242) \pgfimage[width= 95.30pt,height= 78.73pt,interpolate=false]{../thesis_plots/output/grf_shape_ras3}};
## (line 320) \pgfimage[width= 95.30pt,height= 78.73pt,interpolate=false]{../thesis_plots/output/grf_shape_ras4}};
######################################################################

# shape parameters
shape_1=0.1
shape_2=0.5
shape_3=1.5
shape_4=6

# scale
scale=1
# signal variance
signal_var=0.01
# coordinates
x <- seq(0,5, by=0.1)
y <- seq(0,5, by=0.1)
# Matern function (this should be using the Reparameterized scale so scale=\\alpha)
matern_mod_1 <- RMhandcock(nu=shape_1,scale=scale,var=signal_var)
matern_mod_2 <- RMhandcock(nu=shape_2,scale=scale,var=signal_var)
matern_mod_3 <- RMhandcock(nu=shape_3,scale=scale,var=signal_var)
matern_mod_4 <- RMhandcock(nu=shape_4,scale=scale,var=signal_var)

# Simulate field
set.seed(111)
sim_RF_1 <- RFsimulate(model=matern_mod_1,x=x,y=y)
set.seed(111)
sim_RF_2 <- RFsimulate(model=matern_mod_2,x=x,y=y)
set.seed(111)
sim_RF_3 <- RFsimulate(model=matern_mod_3,x=x,y=y)
set.seed(111)
sim_RF_4 <- RFsimulate(model=matern_mod_4,x=x,y=y)

# combine data
dat_1 <- data.frame(x=x,y=y) %>% expand.grid() %>% mutate(z=c(sim_RF_1$variable1), kappa="$\\kappa=0.1$")
dat_2 <- data.frame(x=x,y=y) %>% expand.grid() %>% mutate(z=c(sim_RF_2$variable1), kappa="$\\kappa=0.5$")
dat_3 <- data.frame(x=x,y=y) %>% expand.grid() %>% mutate(z=c(sim_RF_3$variable1), kappa="$\\kappa=1.5$")
dat_4 <- data.frame(x=x,y=y) %>% expand.grid() %>% mutate(z=c(sim_RF_4$variable1), kappa="$\\kappa=6$")
dat <- union(dat_1,dat_2)%>% union(dat_3) %>% union(dat_4)
#GRFs
tikz(file.path(file_location,"grf_shape.tex"),width = 3, height = 3)
(ggplot(dat,aes(x=x,y=y,fill=z))
  +geom_raster()
  +facet_wrap(~kappa)
  +scale_fill_continuous(type="viridis")
  +scale_x_continuous(expand=c(0,0))
  +scale_y_continuous(expand=c(0,0))
  +xlab("x-coordinate")
  +ylab("y-coordinate")
  +theme_bw()
  +theme(panel.spacing.x=unit(0, "cm"),panel.spacing.y=unit(0, "cm"),
         legend.position = "none",strip.background =element_rect(fill="white"),
         axis.ticks = element_blank(),axis.text=element_blank())
  
)
dev.off()


## 4. Spline Correlograms  ###########################################
## For observed data
######################################################################

## Seeds
seed_locations <- seeds %>% dplyr::select(x,y)
obs_seeds_splinecorr <- spline.correlog(x=seed_locations$x,
                                        y=seed_locations$y,
                                        z=seeds$count,
                                        xmax=60,
                                        resamp=100, # for CI bands
                                        save=TRUE,
                                        df=8)


## Seedlings
seedling_locations <- seedlings %>% dplyr::select(x,y)
obs_seedlings_splinecorr <- spline.correlog(x=seedling_locations$x,
                                            y=seedling_locations$y,
                                            z=seedlings$count,
                                            xmax=70,
                                            resamp=100,
                                            save=TRUE,
                                            df=8)

corr_data <- (tibble(.rows = 600)
              %>% mutate(x=c(obs_seeds_splinecorr$real$predicted$x,obs_seedlings_splinecorr$real$predicted$x),
                         y=c(obs_seeds_splinecorr$real$predicted$y,obs_seedlings_splinecorr$real$predicted$y),
                         y_lower = c(obs_seeds_splinecorr$boot$boot.summary$predicted$y["0.025",],obs_seedlings_splinecorr$boot$boot.summary$predicted$y["0.025",]),
                         y_upper= c(obs_seeds_splinecorr$boot$boot.summary$predicted$y["0.975",],obs_seedlings_splinecorr$boot$boot.summary$predicted$y["0.975",]),
                         type = c(rep("Seeds",300),rep("Seedlings",300)))
              %>% mutate(y = if_else(y < -0.5, -0.5, y), #filter out large confidence band in plot
                         y_lower = if_else(y_lower < -0.5, -0.5, y_lower))
              %>% mutate(y = if_else(y > 2, 2, y), #filter out large confidence band in plot
                         y_upper = if_else(y_upper > 2, 2, y_upper))
              
)
tikz(file.path(file_location,"splinecorr.tex"),width = 5, height = 3)
ggplot(corr_data, aes(x=x,y=y,colour=type))+
  geom_line(lwd=1)+
  geom_ribbon(aes(ymin=y_lower,
                  ymax=y_upper,fill=type),alpha=0.4,colour = NA)+
  geom_hline(yintercept = 0, linetype = "dotted")+
  scale_colour_manual(name = '',values=c("#186218","#7B3F00"))+
  scale_fill_manual(name = '',values=c("#186218","#7B3F00"))+
  theme_bw()+
  ggtitle("Spline Correlogram")+
  xlab("distance(m)")+
  ylab("correlation")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


## 5. Within plot - Distances  #######################################
######################################################################

# assemble data set
plot_u <- union(seeds_u %>% unlist() %>% data.frame() %>% setNames("distance") %>% mutate(type="seed traps"),
                seedlings_u %>% unlist() %>% data.frame() %>% setNames("distance") %>% mutate(type="quadrats"))
plot_u$type=factor(plot_u$type, levels=c("seed traps", "quadrats"))


tikz(file.path(file_location,"interpoint_distances.tex"),width = 6, height = 3)
(ggplot(plot_u, aes(x=distance,fill=type, after_stat(density))) +
    geom_histogram(position="identity", binwidth=1, alpha=0.7) +
    scale_fill_manual(name = '', values=c("#7B3F00", "#186218"))+
    theme_bw() +
    ggtitle("Within Plot Distances") +
    xlab("\ndistance(m)") +
    ylab("density\n")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)
dev.off()


## 6. PPD Mean & CV Parameters  ######################################
######################################################################

## Seeds
set.seed(111)
## draw from seed CV and mean priors
seed_sd <- rnorm(n=1000,mean=log(0.25),sd=(log(0.95)-log(0.05))/4)
# sd(log of seed count) \approx cv of seed count
seed_mu <- rnorm(n=1000, mean = log(10), sd = (log(40)-log(5))/4)
## sample from seed count priors and compute quantiles
seed_count <- list()
q <- list()
for (i in 1:1000){
  seed_count[[i]] <- exp(rnorm(1000,mean=seed_mu[i],sd=exp(seed_sd[i])))
  q[[i]] <- quantile(seed_count[[i]],probs=c(0.1,0.9))
}

sum_q_seed <- (bind_rows(q)
               %>% pivot_longer(cols=everything(),names_to = "quantile",values_to = "value")
               %>% mutate(quantile = gsub("(\\%)","\\\\%",quantile,perl=TRUE)))


## Establishment
mu_i <- list()
sigma_i <- list()
logitp_i <- list()
q <- list()
seed_base <- 111

for (i in 1:1000){
  
  set.seed(seed_base+i)
  
  #sample mu_i
  mu_i[[i]]<-rnorm(n=1, mean=logit(0.1), sd=((logit(0.5))-(logit(0.01)))/4)
  
  #sample sigma_i
  sigma_i[[i]]<-rnorm(n=1, mean=log(1), sd=(log(3.5)-log(0.5))/4)
  
  #sample logit(p_i), compute quantiles
  logitp_i[[i]] <- rnorm(1000,mean=mu_i[[i]],sd=exp(sigma_i[[i]]))
  q[[i]] <- quantile(logitp_i[[i]],probs=c(0.1,0.9))
}
sum_q_est <- (bind_rows(q)
              %>% pivot_longer(cols=everything(),names_to = "quantile",values_to = "value")
              %>% mutate(value = plogis(value))
              %>% mutate(quantile = gsub("(\\%)","\\\\%",quantile,perl=TRUE)))


# quantile plots
seed_quant <- (ggplot(sum_q_seed,aes(x=value,y=after_stat(density),fill=quantile))+
                 geom_density(alpha=0.7)+
                 xlab("\nSeed Count")+
                 ylab("density\n")+
                 scale_fill_brewer(palette = "Dark2")+
                 theme_bw()+
                 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.spacing.x = unit(0.2,"cm"))+
                 guides(fill=guide_legend(title.position="top"))
)

est_quant <- (ggplot(sum_q_est,aes(x=value,y=after_stat(density),fill=quantile))+
                geom_density(alpha=0.7)+
                xlab("\nEstablishment Probability")+
                ylab("")+
                scale_fill_brewer(palette = "Dark2")+
                theme_bw()+
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.spacing.x = unit(0.2,"cm"))+
                guides(fill=guide_legend(title.position="top"))
)

tikz(file.path(file_location,"cv_quant.tex"),width = 6.1, height = 3)
ggarrange(seed_quant, est_quant, nrow=1,common.legend = TRUE)
dev.off()



## 7. AC validation RMSRE Results  ###################################
## RMSRE results for ac grouping only
######################################################################

rmsre_results<-readRDS("./validation/output/rmsre_results.RDS")
# 
rmsre_ac <- (rmsre_results
             %>% mutate(set=as.factor(set),ac=as.factor(ac),num_abc=as.factor(num_abc))
             %>% group_by(ac)
             %>% dplyr::summarize(across(starts_with(c("seed","est")),mean))
             %>% ungroup
             %>% rename_with(~ gsub("seed_","Seed ",.x,perl=TRUE))
             %>% rename_with(~ gsub("est_","Establishment ",.x,perl=TRUE))
             %>% pivot_longer(cols=starts_with(c("Seed","Establishment")),names_to="param",values_to="rmsre")
             %>% mutate(type = if_else(startsWith(param,"Seed"),"Seed","Establishment"))
             %>% mutate(param = gsub("Seed |Establishment ","",param,perl=TRUE))
             %>% mutate(param = if_else(param == "Matern Scale", "Matérn\nReparameterized Scale", param))
             %>% mutate(param = if_else(param == "Matern Shape", "Matérn Shape", param))
             
             
)
rmsre_ac$param <- factor(rmsre_ac$param, levels=c("Matérn\nReparameterized Scale","Matérn Shape", "Mean","Nugget Proportion","Total CV"))
tikz(file.path(file_location,"rmsre.tex"),width = 5.8, height = 2.8)
ggplot(rmsre_ac, aes(y=param,x=rmsre,colour=type,fill=type,shape=ac))+
  geom_point(alpha=0.6)+
  #scale_fill_brewer(palette = "Dark2")+
  #geom_jitter(alpha=0.8)+
  #facet_wrap(~param,scales="free",ncol=5)+
  scale_x_log10()+
  theme_bw()+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("RMSRE")+
  ylab("")+
  scale_colour_manual(name="", labels=c("Establishment"="Plant Establishment",
                                        "Seed"="Seed Dispersal"),values=c("#1B9E77","#D95F02"))+
  scale_fill_manual(name="", labels=c("Establishment"="Plant Establishment",
                                      "Seed"="Seed Dispersal"),values=c("#1B9E77","#D95F02"))+
  scale_shape_manual(name="Acceptance Rate", values=c(21,22,24))+
  theme(legend.position="top",legend.direction = "vertical",legend.box="horizontal",legend.justification = c(0, 1),
        legend.title=element_text(size=9),axis.title.x=element_text(size=9),
        axis.text.y=element_text(size=9),axis.text.x=element_text(size=9),panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  guides(fill="none",colour = guide_legend(override.aes = list(alpha = 1)), shape=guide_legend(order=1))+
  scale_y_discrete(limits=rev)
dev.off()


## 8. Coverage  ######################################################
## 90% coverage results
######################################################################

check_coverage <- readRDS("./validation/output/check_coverage.RDS")

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
                       ,param=gsub("\\."," ",param,perl=TRUE),
                       type = if_else(startsWith(param,"Seed"),"Seed","Establishment"),
                       param = gsub("Seed |Establishment ","",param,perl=TRUE),
                       param = if_else(param == "Matern Scale", "Matérn\nReparameterized Scale", param),
                       param = if_else(param == "Matern Shape", "Matérn Shape", param)
            )

)
results$param <- factor(results$param, levels=c("Matérn\nReparameterized Scale","Matérn Shape", "Mean","Nugget Proportion","Total CV"))
tikz(file.path(file_location,"coverage.tex"),width = 5.8, height = 2.8)
ggplot(results, aes(y=param,x=percent,colour=type))+
  annotate("rect", xmin=b_test$conf.int[1], xmax=b_test$conf.int[2], ymin=-Inf, ymax=Inf, alpha=0.15, fill="#670000")+
  geom_point()+
  scale_colour_manual(name="", labels=c("Establishment"="Plant Establishment  ",
                                        "Seed"="Seed Dispersal"),values=c("#1B9E77","#D95F02"))+
  geom_vline(aes(xintercept=cov_per),colour="#670000")+
  ggtitle("")+
  xlab("Coverage Probability")+
  ylab("")+
  theme_bw()+
  #guides(colour=guide_legend( label.hjust=0.1))+
  theme(axis.text.y=element_text(size=9),axis.text.x=element_text(size=9),axis.title.x=element_text(size=10),
        legend.direction = "vertical",legend.box="horizontal",legend.justification = c(0, 1),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(),legend.position="top")+
  scale_y_discrete(limits=rev)
dev.off()

## 9. SBC uniformity plots ###########################################
## Performed SBC using:
## N = 100000 rank statistics
## b = 51 equi-spaced histogram bins
## L = 100 draws from the posterior (to compute ranks statistics)
## scaled summary statistics manually using scaling by sd
## instead of using the default abc::abc MAD scaling
######################################################################

rankstat_sd <- (list.files(path = "./validation/output/sbc_scaleonly/",
                           pattern="*_[0-9]+\\.RDS$",
                           full.names = TRUE)
                %>% map_dfr(readRDS)
)


sbc_dat <- (rankstat_sd
            %>% rename_with(~ gsub("seed_","Seed ",.x,perl=TRUE))
            %>% rename_with(~ gsub("est_","Establishment ",.x,perl=TRUE))
            %>% pivot_longer(cols=everything(),
                             names_to="param",
                             values_to="value")
            %>% mutate(param=if_else(param=="Seed Matern Scale","Seed Matérn Reparameterized Scale",param),
                                        param=if_else(param=="Establishment Matern Scale","Establishment Matérn Reparameterized Scale",param),
                                        param=if_else(param=="Seed Matern Shape","Seed Matérn Shape",param),
                                        param=if_else(param=="Establishment Matern Shape","Establishment Matérn Shape",param))
)

sbc_dat$param <- factor(x=sbc_dat$param, levels=c("Seed Matérn Reparameterized Scale",
                                                  "Establishment Matérn Reparameterized Scale",
                                                    "Seed Matérn Shape",
                                                  "Establishment Matérn Shape",
                                                    "Seed Mean",
                                                  "Establishment Mean",
                                                    "Seed Nugget Proportion",
                                                  "Establishment Nugget Proportion",
                                                    "Seed Total CV",
                                                    "Establishment Total CV"))
bw=2
ci_band <- qbinom(c(0.005, 0.995), size = 100000, prob = 1/51)

tikz(file.path(file_location,"sbc.tex"),width = 5.9, height = 6.5)
ggplot(sbc_dat, aes(x=value))+
  geom_histogram(binwidth=bw)+
  geom_density(aes(y=bw * after_stat(count)))+
  annotate("rect", ymin=ci_band[1], ymax=ci_band[2], xmin=-Inf, xmax=Inf, alpha=0.2, fill="#670000")+
  facet_wrap(~param, nrow = 5, ncol=2)+
  ylab("")+
  xlab("")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 7.5),strip.background = element_rect(fill="white"),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())
dev.off()






## 10. ABC Parameter Estimates ######################################
## Seemed easiest to plot all parameters separately so some can be 
## log scaled individually
######################################################################

## read in data
obs_stats <- readRDS("./summary_statistics/output/observed_stats.RDS")
sim_params <- readRDS("./simulations/output/params_meta_splinecorr.RDS")
sim_stats <- readRDS("./simulations/output/summary_statistics_meta_splinecorr.RDS")

#remove xint,eint
obs_stats <- obs_stats %>% dplyr::select(-seeds_xint,-seeds_eint, -seedlings_xint, -seedlings_eint)
sim_stats <- sim_stats %>% dplyr::select(-seeds_xint,-seeds_eint, -seedlings_xint, -seedlings_eint)

## perform scaling on simulated sum stats
sim_stats <- (sim_stats
              %>% mutate(across(.cols=everything(), ~drop(scale.default(.x))))
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
# seed posterior
seed_post <- seed_conversion(abc_results$unadj.values)
# Table 4a in thesis - posterior summary statistics
stat_seed_post <- post_stats(seed_post)
# seed prior
seed_prior <- seed_conversion(sim_params)
seed_dist <- bind_rows(Posterior=seed_post, Prior=seed_prior, .id="dist")
# establishment posterior
establishment_post <- est_conversion(abc_results$unadj.values)
# Table 4a in thesis - posterior summary statistics
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
              %>% mutate(key=gsub("\\."," ",key,perl=TRUE))
              %>% mutate(key=if_else(key=="Matern Scale","Matérn Scale",key))
              %>% mutate(key=if_else(key=="Matern Shape","Matérn Shape",key))
              )
est_dist <- (est_dist
             %>% mutate(key=gsub("\\."," ",key,perl=TRUE))
             %>% mutate(key=if_else(key=="Matern Scale","Matérn Scale",key))
             %>% mutate(key=if_else(key=="Matern Shape","Matérn Shape",key))
             )
stat_seed_post <- (stat_seed_post
                   %>% mutate(key=gsub("\\."," ",key,perl=TRUE))
                   %>% mutate(key=if_else(key=="Matern Scale","Matérn Scale",key))
                   %>% mutate(key=if_else(key=="Matern Shape","Matérn Shape",key))
                   )
stat_est_post <- (stat_est_post
                  %>% mutate(key=gsub("\\."," ",key,perl=TRUE))
                  %>% mutate(key=if_else(key=="Matern Scale","Matérn Scale",key))
                  %>% mutate(key=if_else(key=="Matern Shape","Matérn Shape",key))
                  )


## individual plots
seed_plot1 <- (ggplot((seed_dist %>% filter(key %in% c("Matérn Scale"))), aes(value, after_stat(density),fill=dist))
              + geom_density(alpha=0.8,position='identity')
              #+ geom_vline(data=stat_seed_post, aes(xintercept=mean, colour='posterior mean'))
              + geom_segment(data=(stat_seed_post %>% filter(key %in% c("Matérn Scale"))), inherit.aes=FALSE, aes(x=hdi_lower,
                                                                         xend=hdi_upper,
                                                                         y=0,
                                                                         yend=0,colour='posterior HDI'),
                             size=2,lineend="round",alpha=0.5,guides=FALSE)
              #geom_segment(x=dat.hdi[1],xend=dat.hdi[2],y=0,yend=0,color="blue",size=2,lineend="round")
              #+ facet_wrap(~key, scales = 'free',nrow=1)
              + scale_colour_manual(name = '',values=c('#670000'),guide="none")
              + scale_fill_manual(name = '',labels=c("Posterior   ","Prior"),values=c('#670000','darkgrey'))
              + ylab("")
              + theme_bw()
              + ggtitle("Seed dispersal",subtitle="Matérn Reparameterized Scale (m)")
              + theme(legend.position="bottom", axis.title = element_blank(),title =element_text(size=8),legend.text=element_text(size=9),legend.spacing.x = unit(0.4, 'cm'))
              #+ scale_x_log10()
)
seed_plot2 <- (ggplot((seed_dist %>% filter(key %in% c("Matérn Shape"))), aes(value, after_stat(density),fill=dist))
               + geom_density(alpha=0.8,position='identity')
               #+ geom_vline(data=stat_seed_post, aes(xintercept=mean, colour='posterior mean'))
               + geom_segment(data=(stat_seed_post %>% filter(key %in% c("Matérn Shape"))), inherit.aes=FALSE, aes(x=hdi_lower,
                                                                                                                     xend=hdi_upper,
                                                                                                                     y=0,
                                                                                                                     yend=0,colour='posterior HDI'),
                              size=2,lineend="round",alpha=0.5,guides=FALSE)
               #geom_segment(x=dat.hdi[1],xend=dat.hdi[2],y=0,yend=0,color="blue",size=2,lineend="round")
               #+ facet_wrap(~key, scales = 'free',nrow=1)
               + scale_colour_manual(name = '',values=c('#670000'),guide="none")
               + scale_fill_manual(name='',labels=c("Posterior   ","Prior"),values=c('#670000','darkgrey'))
               + ylab("")
               + theme_bw()
               + ggtitle("",subtitle ="Matérn Shape")
               + theme(legend.position="bottom", axis.title = element_blank(),title =element_text(size=8),legend.text=element_text(size=9),legend.spacing.x = unit(0.4, 'cm'))
               #+ scale_x_log10()
)
seed_plot3 <- (ggplot((seed_dist %>% filter(key %in% c("Mean"))), aes(value, after_stat(density),fill=dist))
               + geom_density(alpha=0.8,position='identity')
               #+ geom_vline(data=stat_seed_post, aes(xintercept=mean, colour='posterior mean'))
               + geom_segment(data=(stat_seed_post %>% filter(key %in% c("Mean"))), inherit.aes=FALSE, aes(x=hdi_lower,
                                                                                                             xend=hdi_upper,
                                                                                                             y=0,
                                                                                                             yend=0,colour='posterior HDI'),
                              size=2,lineend="round",alpha=0.5,guides=FALSE)
               #geom_segment(x=dat.hdi[1],xend=dat.hdi[2],y=0,yend=0,color="blue",size=2,lineend="round")
               #+ facet_wrap(~key, scales = 'free',nrow=1)
               + scale_colour_manual(name = '',values=c('#670000'),guide="none")
               + scale_fill_manual(name='',labels=c("Posterior   ","Prior"),values=c('#670000','darkgrey'))
               + ylab("")
               + theme_bw()
               + ggtitle("", subtitle="Mean (seed count)")
               + theme(legend.position="bottom", axis.title = element_blank(),title =element_text(size=8),legend.text=element_text(size=9),legend.spacing.x = unit(0.4, 'cm'))
               + scale_x_log10()
)

seed_plot4 <- (ggplot((seed_dist %>% filter(key %in% c("Nugget Proportion"))), aes(value, after_stat(density),fill=dist))
               + geom_density(alpha=0.8,position='identity')
               #+ geom_vline(data=stat_seed_post, aes(xintercept=mean, colour='posterior mean'))
               + geom_segment(data=(stat_seed_post %>% filter(key %in% c("Nugget Proportion"))), inherit.aes=FALSE, aes(x=hdi_lower,
                                                                                                               xend=hdi_upper,
                                                                                                               y=0,
                                                                                                               yend=0,colour='posterior HDI'),
                              size=2,lineend="round",alpha=0.5,guides=FALSE)
               #geom_segment(x=dat.hdi[1],xend=dat.hdi[2],y=0,yend=0,color="blue",size=2,lineend="round")
               #+ facet_wrap(~key, scales = 'free',nrow=1)
               + scale_colour_manual(name = '',values=c('#670000'),guide="none")
               + scale_fill_manual(name='',labels=c("Posterior   ","Prior"),values=c('#670000','darkgrey'))
               + ylab("")
               + theme_bw()
               + ggtitle("", subtitle="Nugget Proportion")
               + theme(legend.position="bottom", axis.title = element_blank(),title =element_text(size=8),legend.text=element_text(size=9),legend.spacing.x = unit(0.4, 'cm'))
               + scale_x_log10()
)

seed_plot5 <- (ggplot((seed_dist %>% filter(key %in% c("Total CV"))), aes(value, after_stat(density),fill=dist))
               + geom_density(alpha=0.8,position='identity')
               #+ geom_vline(data=stat_seed_post, aes(xintercept=mean, colour='posterior mean'))
               + geom_segment(data=(stat_seed_post %>% filter(key %in% c("Total CV"))), inherit.aes=FALSE, aes(x=hdi_lower,
                                                                                                               xend=hdi_upper,
                                                                                                               y=0,
                                                                                                               yend=0,colour='posterior HDI'),
                              size=2,lineend="round",alpha=0.5,guides=FALSE)
               #geom_segment(x=dat.hdi[1],xend=dat.hdi[2],y=0,yend=0,color="blue",size=2,lineend="round")
               #+ facet_wrap(~key, scales = 'free',nrow=1)
               + scale_colour_manual(name = '',values=c('#670000'),guide="none")
               + scale_fill_manual(name='',labels=c("Posterior   ","Prior"),values=c('#670000','darkgrey'))
               + ylab("")
               + theme_bw()
               + ggtitle("", subtitle="Total CV")
               + theme(legend.position="bottom", axis.title = element_blank(),title =element_text(size=8),legend.text=element_text(size=9),legend.spacing.x = unit(0.4, 'cm'))
               + scale_x_log10()
)

est_plot1<- (ggplot((est_dist %>% filter(key %in% c("Matérn Scale"))), aes(value, after_stat(density),fill=dist))
            + geom_density(alpha=0.8,position='identity')
            #+ geom_vline(data=stat_est_post, aes(xintercept=mean, colour='posterior mean'))
            + geom_segment(data=(stat_est_post %>% filter(key %in% c("Matérn Scale"))), inherit.aes=FALSE, aes(x=hdi_lower,
                                                                      xend=hdi_upper,
                                                                      y=0,
                                                                      yend=0,colour='posterior HDI'),
                           size=2,lineend="round",alpha=0.5,guides=FALSE)
            #+ facet_wrap(~key, scales = 'free',nrow=1)
            + scale_colour_manual(name = '',values=c('#670000'),guide="none")
            + scale_fill_manual(name='',labels=c("Posterior   ","Prior"),values=c('#670000','darkgrey'))
            + ylab("")
            + theme_bw()
            + ggtitle("Plant establishment",subtitle="Matérn Reparameterized Scale (m)")
            + theme(legend.position="bottom", axis.title = element_blank(), title =element_text(size=8),legend.text=element_text(size=9),legend.spacing.x = unit(0.4, 'cm'))
            #+ scale_x_log10()
)

est_plot2<- (ggplot((est_dist %>% filter(key %in% c("Matérn Shape"))), aes(value, after_stat(density),fill=dist))
             + geom_density(alpha=0.8,position='identity')
             #+ geom_vline(data=stat_est_post, aes(xintercept=mean, colour='posterior mean'))
             + geom_segment(data=(stat_est_post %>% filter(key %in% c("Matérn Shape"))), inherit.aes=FALSE, aes(x=hdi_lower,
                                                                                                                xend=hdi_upper,
                                                                                                                y=0,
                                                                                                                yend=0,colour='posterior HDI'),
                            size=2,lineend="round",alpha=0.5,guides=FALSE)
             #+ facet_wrap(~key, scales = 'free',nrow=1)
             + scale_colour_manual(name = '',values=c('#670000'),guide="none")
             + scale_fill_manual(name='',labels=c("Posterior   ","Prior"),values=c('#670000','darkgrey'))
             + ylab("")
             + theme_bw()
             + ggtitle("",subtitle ="Matérn Shape")
             + theme(legend.position="bottom", axis.title = element_blank(),title =element_text(size=8),legend.text=element_text(size=9),legend.spacing.x = unit(0.4, 'cm'))
             #+ scale_x_log10()
)
est_plot3<- (ggplot((est_dist %>% filter(key %in% c("Mean"))), aes(value, after_stat(density),fill=dist))
             + geom_density(alpha=0.8,position='identity')
             #+ geom_vline(data=stat_est_post, aes(xintercept=mean, colour='posterior mean'))
             + geom_segment(data=(stat_est_post %>% filter(key %in% c("Mean"))), inherit.aes=FALSE, aes(x=hdi_lower,
                                                                                                                xend=hdi_upper,
                                                                                                                y=0,
                                                                                                                yend=0,colour='posterior HDI'),
                            size=2,lineend="round",alpha=0.5,guides=FALSE)
             #+ facet_wrap(~key, scales = 'free',nrow=1)
             + scale_colour_manual(name = '',values=c('#670000'),guide="none")
             + scale_fill_manual(name='',labels=c("Posterior   ","Prior"),values=c('#670000','darkgrey'))
             + ylab("")
             + theme_bw()
             + ggtitle("", subtitle="Mean (probability)")
             + theme(legend.position="bottom", axis.title = element_blank(),title =element_text(size=8),legend.text=element_text(size=9),legend.spacing.x = unit(0.4, 'cm'))
             + scale_x_log10()
)
est_plot4<- (ggplot((est_dist %>% filter(key %in% c("Nugget Proportion"))), aes(value, after_stat(density),fill=dist))
             + geom_density(alpha=0.8,position='identity')
             #+ geom_vline(data=stat_est_post, aes(xintercept=mean, colour='posterior mean'))
             + geom_segment(data=(stat_est_post %>% filter(key %in% c("Nugget Proportion"))), inherit.aes=FALSE, aes(x=hdi_lower,
                                                                                                                xend=hdi_upper,
                                                                                                                y=0,
                                                                                                                yend=0,colour='posterior HDI'),
                            size=2,lineend="round",alpha=0.5,guides=FALSE)
             #+ facet_wrap(~key, scales = 'free',nrow=1)
             + scale_colour_manual(name = '',values=c('#670000'),guide="none")
             + scale_fill_manual(name='',labels=c("Posterior   ","Prior"),values=c('#670000','darkgrey'))
             + ylab("")
             + theme_bw()
             + ggtitle("", subtitle="Nugget Proportion")
             + theme(legend.position="bottom", axis.title = element_blank(),title =element_text(size=8),legend.text=element_text(size=9),legend.spacing.x = unit(0.4, 'cm'))
             + scale_x_log10()
)
est_plot5<- (ggplot((est_dist %>% filter(key %in% c("Total CV"))), aes(value, after_stat(density),fill=dist))
             + geom_density(alpha=0.8,position='identity')
             #+ geom_vline(data=stat_est_post, aes(xintercept=mean, colour='posterior mean'))
             + geom_segment(data=(stat_est_post %>% filter(key %in% c("Total CV"))), inherit.aes=FALSE, aes(x=hdi_lower,
                                                                                                                xend=hdi_upper,
                                                                                                                y=0,
                                                                                                                yend=0,colour='posterior HDI'),
                            size=2,lineend="round",alpha=0.5,guides=FALSE)
             #+ facet_wrap(~key, scales = 'free',nrow=1)
             + scale_colour_manual(name = '',values=c('#670000'),guide="none")
             + scale_fill_manual(name='',labels=c("Posterior   ","Prior"),values=c('#670000','darkgrey'))
             + ylab("")
             + theme_bw()
             + ggtitle("", subtitle="Total CV")
             + theme(legend.position="bottom", axis.title = element_blank(),title =element_text(size=8),legend.text=element_text(size=9),legend.spacing.x = unit(0.4, 'cm'))
             + scale_x_log10()
)


tikz(file.path(file_location,"abc_parameters.tex"),width = 5.5, height = 8)
ggarrange(seed_plot1, est_plot1, ggplot() + theme_void(), ggplot() + theme_void(),
          seed_plot2, est_plot2, ggplot() + theme_void(), ggplot() + theme_void(),
          seed_plot3, est_plot3, ggplot() + theme_void(), ggplot() + theme_void(),
          seed_plot4, est_plot4, ggplot() + theme_void(), ggplot() + theme_void(),
          seed_plot5, est_plot5, ncol=2, nrow=9, heights=c(1,-0.01,1,-0.01,1,-0.01,1,-0.01,1),common.legend=TRUE)
dev.off()


## 11. Frequentist Profile Likelihood #############################################
######################################################################

# make coordinates factors
seed_data <- seeds
seed_data$pos <- numFactor(seed_data$x, seed_data$y)
seed_data$plot_group <- factor(seed_data$plot)
seed_data$obs <- factor(1:nrow(seed_data))

fit.mat <- glmmTMB(count ~ mat(pos + 0 | plot_group),family=poisson(link="log"),data=seed_data)
# add nugget (random effects)
fit.mat.nugget <- update(fit.mat,.~.+(1|obs))

# Table 5 in thesis - MLE estimates
# in order, mean, sd, range, smoothness, and nugget
exp(fit.mat.nugget$sdr$par.fixed)
# mean
exp(2.05309641 + c(-1,1)*1.96*0.2895127)
# sd
exp(-0.01288684 + c(-1,1)*1.96*0.1963664)
# range
exp(5.87611822 + c(-1,1)*1.96*1.8595513)
# smoothness
exp(-1.53081872 + c(-1,1)*1.96*0.4780515)
# nugget
exp(-7.91728480 + c(-1,1)*1.96*2788.3747517)#huge standard error

# profile scale
tt_scale <- tmbprofile(fit.mat$obj, 3)

tt_nugget <- tmbprofile(fit.mat.nugget$obj, 5
                        # ,parm.range=seq(1e-4,1,by=0.01)# doesn't seem to matter what range I specify
)
# remove extreme observations for plotting
tt_nugget <- tt_nugget[-c(1:6),]


tikz(file.path(file_location,"profile_likelihood.tex"),width = 3, height = 3)
ggplot(tt_scale, aes(x=exp(theta),y=value))+
  geom_line()+
  #geom_hline(yintercept=min(tt_scale$value)+1.92,colour="red")+
  theme_bw()+
  xlab("\nMatérn Scale (m)")+
  ylab("Profile Negative Log-likelihood\n")+
  scale_x_log10()
dev.off()



tikz(file.path(file_location,"profile_likelihood_nugget.tex"),width = 3, height = 3)
ggplot(tt_nugget, aes(x=exp(theta),y=value))+
  geom_line()+
  theme_bw()+
  xlab("\nNugget Variance")+
  ylab("Profile Negative Log-likelihood\n")+
  scale_x_log10()
dev.off()

## Not Included
######################################################################

## 0. Spatial trends by coordinates ##################################
## Not sure this really tells me anything
## Plot J for seeds - higher counts at higher x and y coords
## all other seed plots ?
## Plot K for seedlings - decreasing counts for larger y coords
## all other seedling plots ?
## No obvious trends - excluding from thesis
######################################################################

seed_data <- all_dat %>% filter(type=="Seeds")

ggplot(seed_data, aes(x=x,y=count)) +
  geom_point() + 
  geom_smooth() +
  facet_wrap(~plot, scales = "free")+
  xlab("x coord")

ggplot(seed_data, aes(x=y,y=count)) +
  geom_point() + 
  geom_smooth() +
  facet_wrap(~plot, scales = "free")+
  xlab("y coord")


seedling_data <- all_dat %>% filter(type=="Seedlings")

ggplot(seedling_data, aes(x=x,y=count)) +
  geom_point() + 
  geom_smooth() +
  facet_wrap(~plot, scales = "free")+
  xlab("x coord")

ggplot(seedling_data, aes(x=y,y=count)) +
  geom_point() + 
  geom_smooth() +
  facet_wrap(~plot, scales = "free")+
  xlab("y coord")


### Trying ABC plots using facets but scaling individually
# ## Plots
# # seed distributions
# # likely issues with oob stuff that I can't figure out right now
# library(ggh4x)
# scales_seeds <- list(
#   NULL,
#   NULL,
#   scale_x_log10(),
#   scale_x_log10(),
#   scale_x_log10()
# )
# library(scales)
# seed_dist <- seed_dist %>% filter(oob_discard)
# seed_plot <- (ggplot(seed_dist, aes(value, after_stat(density),fill=dist))
#               + geom_density(alpha=0.8,position='identity')
#               #+ geom_vline(data=stat_seed_post, aes(xintercept=mean, colour='posterior mean'))
#               + geom_segment(data=stat_seed_post, inherit.aes=FALSE, aes(x=hdi_lower,
#                                                                          xend=hdi_upper,
#                                                                          y=0,
#                                                                          yend=0,colour='posterior HDI'),
#                              size=2,lineend="round",alpha=0.5,guides=FALSE)
#               #geom_segment(x=dat.hdi[1],xend=dat.hdi[2],y=0,yend=0,color="blue",size=2,lineend="round")
#               + facet_wrap(~key, scales = 'free',ncol=1)
#               
#               + scale_colour_manual(name = '',values=c('#670000'),guide="none")
#               + scale_fill_manual(name='',values=c('#670000','darkgrey'))
#               + ylab("")
#               + theme_bw()
#               + ggtitle("Seed dispersal distributions")
#               + theme(legend.position="bottom", axis.title = element_blank())
#               #+ scale_x_log10()
#               + facetted_pos_scales(x=scales_seeds,y=NULL)
# )
# 
# seed_plot
#seed_plot + facetted_pos_scales(x=scales)
# 
# 
# # establishment distributions
# est_plot<- (ggplot(est_dist, aes(value, after_stat(density),fill=dist))
#             + geom_density(alpha=0.8,position='identity')
#             #+ geom_vline(data=stat_est_post, aes(xintercept=mean, colour='posterior mean'))
#             + geom_segment(data=stat_est_post, inherit.aes=FALSE, aes(x=hdi_lower,
#                                                                       xend=hdi_upper,
#                                                                       y=0,
#                                                                       yend=0,colour='posterior HDI'),
#                            size=2,lineend="round",alpha=0.5,guides=FALSE)
#             + facet_wrap(~key, scales = 'free',ncol=1)
#             + scale_colour_manual(name = '',values=c('#670000'),guide="none")
#             + scale_fill_manual(name='',values=c('#670000','darkgrey'))
#             + ylab("")
#             + theme_bw()
#             + ggtitle("Plant establishment distributions")
#             + theme(legend.position="bottom", axis.title = element_blank())
#             #+ scale_x_log10()
# )


