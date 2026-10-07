# R-code for Multiview Distance Regularized S-map (MDR S-map)
rm(list = ls());gc()
# Set path
OPath <- 'C:\\data\\Cooperation\\Prodigio\\Code\\Demo_RScript_20250529\\'
setwd(OPath)

library(rEDM) # Package for EDM analyses (e.g., simplex projection and CCM) were conducted in rEDM package ver. 1.2.3
library(doParallel)
library(parallel)
library(foreach)
library(Kendall)
library(MASS)
library(dplyr)
library(glmnet)

library(ggplot2)
library(cowplot)
library(ggpattern)
source('./Code/Function_MDR_inferM.R')
seed <- 49563
set.seed(seed)
SaveFile <- F#T # (T/F) for saving files

#############################################################
#### Data Preparation 
# Read dataset
otufing=read.csv('./Data/Perturbed_CH4crnotu_verc_Species_D47.csv', header=T,stringsAsFactors = F)
otufing[,grep('OTU_',colnames(otufing))]=otufing[,grep('OTU_',colnames(otufing))]*repmat(apply(otufing[,grep('OTU_',colnames(otufing))],1,sum,na.rm=T),1,ncol(otufing[,grep('OTU_',colnames(otufing))]))^-1
ch4 <- as.numeric(otufing[,1])

(da.name <- 'Demo_CH4_otucrn_Perturbed_D47')
do <- otufing

ndo <- nrow(do)
nin <- ndo # library sample size

# Excluding rare species with no clear temporal dynamics 
dto <- apply(do,1,sum) # total abundance
dpo <- do*repmat(dto^-1,1,ncol(do))
(nsp <- ncol(do))                            # number of selected species

# Target for the estimation X -> CH4
infm <- matrix(1,nsp,nsp)
infm[-1,1] <- 0
infm[-1,-1] <- 0

# Excluding the inference from missing values and zero-abundance/concentration
cabs <- 1*(do==0);
cabs <- as.matrix(cabs[-nrow(cabs),]*cabs[-1,])
cabs[is.na(cabs)] <- 0

cna <- 1*(is.na(do));
cna <- as.matrix(cna[-nrow(cna),]+cna[-1,])
cna[cna>1] <- 1


do.mean <- apply(do[1:nin,],2,mean,na.rm=T)  # mean abundance in in-sample
do.sd <- apply(do[1:nin,],2,sd,na.rm=T)      # SD of abundance in in-sample
d <- do[1:(nin-1),]                          # In-sample dataset at time t
d_tp1 <- do[2:(nin),]                        # In-sample dataset at time t+1
ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t
ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1

# Out-sample
d.test <- dt_tp1 <- dst_tp1 <- ds.test <- NULL

# Compiled data at time t 
ds.all <- rbind(ds,ds.test)

#############################################################
# Find the optimal embedding dimension & nonlinearity parameter for each variable 
# based on univariate simplex projection and S-map, respectively

# Univariate simplex projection
Emax <- 10
cri <- 'rmse' # model selection 
Ed <- NULL
forecast_skill_simplex <- NULL
for(i in 1:ncol(ds)){
  spx.i <- simplex(ds[,i],E=2:Emax)
  Ed <- c(Ed,spx.i[which.min(spx.i[,cri])[1],'E'])
  forecast_skill_simplex <- c(forecast_skill_simplex,spx.i[which.min(spx.i[,cri])[1],'rho'])
}
Ed # The optimal embedding dimension for each variable
forecast_skill_simplex # Forecast skills for each variable based on simplex projection

######################################################################
# Find causal variables by CCM analysis for multiview embedding
# Warning: It is time consuming for calculating the causation for each node
# CCM causality test for all node pairs 
ccm.out <- ccm.fast.inferM(ds,inferM=infm, Lmin=5,Lint=10, silience=T)
ccm.sig <- ccm.out[['ccm.sig']]
ccm.rho <- ccm.out[['ccm.rho']]

######################################################################
# Perform multiview embedding analysis for each node
# Warning: It is time consuming for running multview embedding for each nodes
do.multiview <- F
if(do.multiview){
  esele_lag <- esim.lag.demo2(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax,inferM=infm, seed=seed)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(SaveFile){write.csv(esele_lag,paste('./Output/eseleLag_',da.name,'_nin',nin,'_demo.csv',sep=''),row.names=F)}
}

esele <- read.csv(paste('./Output/eseleLag_',da.name,'_nin',nin,'_demo.csv',sep=''),header=T)

####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo2(ds, ds.all, esele, inferM=infm)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]


if(T){
  (da.name <- paste(da.name,'NA',sep='_'))
  cna[apply(cna==1, 1, any), ] <- 1
}

######## Leave-one-out cross-validation for finding the optimal parameters for MDR S-map analysis
### Warning: The cross-validation is the most time-consuming step in MDR S-map requiring massive computations and .  
### Thus, we recommend dividing job into smaller parts (sub.da>1)  or used parallel computation (parall=T, ncore>=1)
### The example showing below divided the parameter space into five parts and ran independently (sub.da=5).
do.MDR.CV <- F 
### The parameter cv.unit determines the precision of selected parameters and strongly influences computation time.
### In our cases, we used cv.unit=0.025 to obtain more precise estimations
### This parameter may be adjusted to 0.05 or even 0.1, depending on how sensitive the results to parameter precision. 
cv.unit <- 0.025 
alpha.so <- seq(0, 1, cv.unit);            # Sequence of alpha
sub.da <- 1                                # Divide the computation job into five parts 
afsp <- eqsplit(1:length(alpha.so),sub.da) # Divide the parameter space based on alpha parameter
alf <- 1 # R-code for  Multiview Distance Regularized S-map (MDR S-map)

# Cross-validation of MDR analysis    
if(do.MDR.CV){ #slow process; 3204 seconds in my computer
  alpha.s <- alpha.so[afsp[alf,1]:afsp[alf,2]] # Subset parameter pace
  cv.ind <- cv.MDR.infM3(ds, ds_tp1, dmatrix.list=dmatrix.train.mvx,rm_zero=T,absence_m=cabs,na_m=cna, inferM=infm,
                         parall=F, ncore=1, keep_intra=T,alpha.seq=alpha.s)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(SaveFile){write.csv(cv.ind,paste('./Output/',da.name,'_nin',nin,'_cvunit',cv.unit,'_alph',alpha.s[1]*100,'_cvout_Nmvx_Rallx.csv',sep=''),row.names=F)}
}
# Repeat the CV under different parameter subsets by changing alf=2,3..,sub.da                                  
# Run CV in the first parameter subset 

cv.ind <- read.csv(paste('./Output/',da.name,'_nin',nin,'_cvunit',cv.unit,'_alph',alpha.so[afsp[alf,1]]*100,
                         '_cvout_Nmvx_Rallx.csv',sep=''),header=T) 
# Select the optimal parameter set with the minimal MSE
paracv.demo <- secv.demo(cv.ind)

############################################################
# Fitting MDR S-map based on the parameters selected by CV
ptype <- 'enet'                           # enet:elastic-net regularization

# Fitting the MDR S-map
smap.demo <- MDRsmap.infM3(paracv=paracv.demo,ptype=ptype,keep_intra=T,out.sample=F,
                           ds,ds_tp1,ds.test,dst_tp1,inferM=infm,
                           dmatrix.list=dmatrix.train.mvx,rm_zero=T,absence_m=cabs,na_m=cna,
                           dmatrix.test.list=dmatrix.test.mvx)

# Save forecast skills
nr.out <- smap.demo[['nr.out']];
# To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
if(SaveFile){
  write.csv(nr.out,paste('./Output/',da.name,'_nin',nin,'_cvunit',cv.unit,'_',ptype,'_nrout_Nmvx_Rallx_demo_NEW.csv',sep=''),row.names=F)
  # Save interaction Jacobian matrices at all time points
  write.csv(smap.demo[['jcof']],paste('./Output/',da.name,'_nin',nin,'_cvunit',cv.unit,'_',ptype,'_jcof_Nmvx_Rallx_demo_NEW.csv',sep=''),row.names=F)
}


nr.out
pred <- apply(filter(smap.demo[['jcof']],Insample==1)[,-c(1:3)]*cbind(rep(1),ds),1,sum)

tipp <- 38
tipp2 <- 44

# Summary statistics of MDR S-map for reconstrucuted dynamic interactome
nr.out.a <- read.csv(paste('./Output/Perturbed_otucrn_CH4_D47_Species_NA_cvunit0.025_enet_nrout_Nmvx_Rallx_demo.csv',sep=''),header=T)
domisp <- as.character(read.csv(paste('./Output/domisp_list.csv',sep=''),header=T)[,1])

nr.out.a$domisp <- factor(1*c(colnames(ds)[-1]%in%domisp),
                          levels = c(0, 1),
                          labels = c("Non-Cascade", "Cascade"))


# ============================================================
# Fig. S11 - Part A: Violin plot
# ============================================================

x0 <- 1
w  <- 0.4

# Compute density and restrict to (0, 1)
dens <- density(nr.out.a$rho, from = 0, to = 1, n = 512, bw = 0.02)
d <- dens$y / max(dens$y) * w

violin_df <- data.frame(
  x = c(x0 + d, rev(x0 - d)),
  y = c(dens$x, rev(dens$x))
)

# Single-row data frame for the CH4 point
ch4_df <- data.frame(x = x0, y = 0.95, label = "CH4")

set.seed(123)

pA <- ggplot(nr.out.a, aes(x = x0, y = rho)) +
  # Violin body (single closed polygon)
  geom_polygon(data = violin_df,
               aes(x = x, y = y),
               fill = "lightgrey", color = "black", alpha = 0.6) +
  # Boxplot (narrow, centered)
  geom_boxplot(width = 0.1, fill = "white", color = "black",
               outlier.shape = NA, alpha = 0.8) +
  # Jittered points colored by domisp
  geom_jitter(aes(color = domisp),
              width = 0.35, height = 0, size = 1.0, alpha = 0.7) +
  # CH4 point
  geom_point(data = ch4_df,
             aes(x = x, y = y, color = label),
             size = 3, shape = 18) +
  scale_color_manual(
    name = NULL,
    values = c("Non-Cascade" = "black",
               "Cascade"     = "blue",
               "CH4"         = "red")
  ) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0.3, 1.7)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(x = "", y = expression(rho),
       title = expression(paste("Prediction skill (", rho, ")",
                                " distribution"))) +
  theme_bw() +
  theme(
    plot.title       = element_text(hjust = 0.5, size = 11),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    legend.position      = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.direction     = "vertical",
    legend.text          = element_text(size = 8)
  )

# ============================================================
# Part B: Time series plot with hatched shading
# ============================================================

tipp  <- 38
tipp2 <- 44

rho_val <- round(cor(ds_tp1[, 1], pred, use = "complete.obs"), 3)

ts_df <- data.frame(
  Day  = 1:46,
  Obs  = ds_tp1[, 1],
  Pred = as.numeric(pred)
)

pB <- ggplot(ts_df, aes(x = Day)) +
  # Hatched rectangle for the perturbation window
  geom_rect_pattern(
    aes(xmin = tipp, xmax = tipp2, ymin = -5, ymax = 5),
    pattern         = "stripe",
    pattern_angle   = 45,
    pattern_density = 0.15,
    pattern_spacing = 0.02,
    pattern_color   = "orange",
    pattern_fill    = NA,
    fill            = NA,
    color           = "white"
  ) +
  # Orange vertical line at the start of the window
  geom_vline(xintercept = tipp, color = "orange") +
  # Observation line
  geom_line(aes(y = Obs, color = "Observation"), linewidth = 0.6) +
  geom_point(aes(y = Obs, color = "Observation"), size = 2.0) +
  # Prediction line
  geom_line(aes(y = Pred, color = "Prediction"), linewidth = 0.6) +
  geom_point(aes(y = Pred, color = "Prediction"), size = 2.0) +
  scale_color_manual(
    name = NULL,
    values = c("Observation" = "black", "Prediction" = "blue")
  ) +
  scale_x_continuous(breaks = seq(0, 46, 10)) +
  coord_cartesian(ylim = c(-3, 2)) +
  labs(x = "Day",
       y = expression(paste("Standardized CH"[4], " production"))) +
  annotate("text", x = 1, y = -2.1,
           label = paste0("Prediction skill = ", rho_val),
           hjust = 0, size = 3.5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black"),
    legend.position      = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.direction     = "vertical",
    legend.text          = element_text(size = 8)
  )

# ============================================================
# Combine with cowplot and add A / B labels
# ============================================================

g_combined <- plot_grid(
  pA, pB,
  ncol           = 2,
  rel_widths     = c(1, 1.6),
  align          = "h",
  axis           = "tb",
  labels         = c("A", "B"),
  label_size     = 16,
  label_fontface = "bold"
)

win.graph(10, 5)
g_combined

# Optional: save to file
# ggsave("FigS11_plot.pdf", g_combined, width = 14, height = 5, dpi = 600)

