# R-code for Multiview Distance Regularized S-map (MDR S-map)
rm(list = ls());gc()
# Set path
OPath <- 'D:/data/Cooperation/Prodigio/Demo_RScript'
setwd(OPath)

library(rEDM) # Package for EDM analyses (e.g., simplex projection and CCM) were conducted in rEDM package ver. 1.2.3
library(doParallel)
library(parallel)
library(foreach)
library(Kendall)
library(MASS)
library(dplyr)
library(glmnet)
source('./Code/Function_MDR_inferM.R')
seed <- 49563
set.seed(seed)
SaveFile <- T #T/F for saving files

#############################################################
#### Data Preparation 
# Read dataset
out.sample <- F # T/F for out-of-sample forecast
if(out.sample){nout <- 2}else{nout <- 0}  # number of out-of-sample

otufing=read.csv('./Data/Perturbed_otucrn_ch4_verc_Species.csv', header=T,stringsAsFactors = F)
otufing[,grep('OTU_',colnames(otufing))]=otufing[,grep('OTU_',colnames(otufing))]*repmat(apply(otufing[,grep('OTU_',colnames(otufing))],1,sum,na.rm=T),1,ncol(otufing[,grep('OTU_',colnames(otufing))]))^-1
ch4 <- as.numeric(otufing[,1])

(da.name <- 'Demo_CH4_otucrn_Perturbed')
do <- otufing

ndo <- nrow(do)
nin <- ndo-nout # library sample size

# Excluding rare species with no clear temporal dynamics 
dto <- apply(do,1,sum) # total abundance
dpo <- do*repmat(dto^-1,1,ncol(do))
(nsp <- ncol(do))                            # number of selected species

# Target for the estimation X -> CH4
infm=matrix(1,nsp,nsp)
infm[-1,1]=0
infm[-1,-1]=0

# Excluding the inference from missing values and zero-abundance/concentration
cabs=1*(do==0);
cabs=as.matrix(cabs[-nrow(cabs),]*cabs[-1,])
cabs[is.na(cabs)]=0

cna=1*(is.na(do));
cna=as.matrix(cna[-nrow(cna),]+cna[-1,])
cna[cna>1]=1


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
do.CCM <- T
ccm.out <- ccm.fast.inferM(ds,inferM=infm, Lmin=5,Lint=10, silience=T)
ccm.sig <- ccm.out[['ccm.sig']]
ccm.rho <- ccm.out[['ccm.rho']]

######################################################################
# Perform multiview embedding analysis for each node
# Warning: It is time consuming for running multview embedding for each nodes
do.multiview <- T
if(do.multiview){
  esele_lag <- esim.lag.demo2(ds,ccm.rho,ccm.sig,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax,inferM=infm)
  # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
  if(SaveFile){write.csv(esele_lag,paste('./Output/eseleLag_',da.name,'_nin',nin,'_demo.csv',sep=''),row.names=F)}
}

esele <- read.csv(paste('./Output/eseleLag_',da.name,'_nin',nin,'_demo.csv',sep=''),header=T)

####################################################
## The computation of multiview distance
dmatrix.mv <- mvdist.demo2(ds, ds.all, esele, inferM=infm)
dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]

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
smap.demo <- MDRsmap.infM3(paracv=paracv.demo,ptype=ptype,keep_intra=T,out.sample=out.sample,
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
