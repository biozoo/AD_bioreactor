############################################################
### DEV analysis
rm(list = ls());gc()
# Set path
OPath <- 'C:/data/Cooperation/Prodigio/Code/Demo_RScript_20250529'
setwd(OPath)
source('./Code/Function_DEV_V2.R')

library(rEDM) # EDM analyses (e.g., simplex projection and CCM) were conducted in rEDM package ver. 1.2.3
library(MASS)
library(dplyr)

seed <- 49563
set.seed(seed)
SaveFile <- F #T# (T/F) for saving files
interpolation <- T
tipp <- 38+1
tipp2 <- 44+1
ind <- (1:tipp2)
################################################################################
# DEV analysis for metabolites
################################################################################
#### Data Preparation
crnotu.p <- read.csv('./Data/OLR_Perturbed_crnotu_verc_CH4_Species_D47.csv',header=T,stringsAsFactors = F)
vfa_rep <- read.csv('.//Data//Reactor_data_envi.csv',header=T,stringsAsFactors = F)[,c('HAc_O2','HPro_O2','i.HBu_O2','HBu_O2','iHVal_O2','HVal_O2','HCa_O2')]

crn.p <- crnotu.p[,-c(grep('OTU_',colnames(crnotu.p)))]

if(interpolation){
  crn.p.na <- apply(crn.p,2,EDM_NA_fill)
  crn.p <- crn.p.na
}

(da.name1 <- 'OLR_R3_crn_verc_D47')
ds <- crn.p[ind,]
nin <- nrow(ds) # library sample size


### Determine the optimal DEV parameters using cross-validation ##
if(T){
  zeroint=apply(ds,2,zeropad)# find the length of zero-padding region in the time series
  dev_p  <-  NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    wd.i <- 15
    if(zeroint[i]>wd.i){wd.i=zeroint[i]+2}
    if(i<=7){
      param.opt.a <- dev_cvls_v2(ds[,i], E = seq(1,10,by=1), window_size = seq(wd.i, 40, by=2), replicate_data= vfa_rep[,i] )
    }else{
      param.opt.a <- dev_cvls_v2(ds[,i], E = seq(1,10,by=1), window_size = seq(wd.i, 40, by=2) )
    }
    (para.optt <- param.opt.a[[1]])
    dev_p <- rbind(dev_p, c(variable=colnames(ds)[i], para.optt))
    rm(list = c('param.opt.a','para.optt'));gc()
    time.used <- proc.time() - t.begin
    cat("effect", i, "in sample completed:", time.used[3],"sec\n")
  }
  if(SaveFile) write.csv(dev_p, paste('./Output/', da.name1,'_DEV_Optparameter.csv',sep=''),row.names=F)
}

##################################################
### Calculate DEV based on the optimal parameters##
if(T){
  dev_crn_p <- read.csv(paste('./Output/', da.name1,'_DEV_Optparameter.csv',sep=''),header=T)
  dev_p <- dev_crn_p
  DEVresult <-NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    para.optt <- as.numeric(dev_p[i,-1])
    names(para.optt) <- colnames(dev_p)[-1]
    ### Performing DEV algorithm
    if(i<=7){
      devts <- DEV_perform_ls2(ds[,i],param.opt=para.optt, replicate_data= vfa_rep[,i], pred_rep1_only=T)
    }else{
      devts <- DEV_perform_ls2(ds[,i],param.opt=para.optt)
    }
    
    DEVresult <- rbind(DEVresult, data.frame(variable=colnames(ds)[i], devts))
    rm(list = c('param.opt.a', 'devts'));gc()
    time.used <- proc.time() - t.begin
    cat("effect", i, "in sample completed:", time.used[3],"sec\n")
  }
  if(SaveFile) write.csv(DEVresult,paste('./Output/', da.name1,'_DEV_TimeSeries.csv',sep=''),row.names=F)
}

dev.crn=read.csv(paste('./Output/', da.name1,'_DEV_TimeSeries.csv',sep=''),header=T)


################################################################################
# DEV analysis for microbial OTUs
################################################################################
#### Data Preparation 

otu.p <- crnotu.p[,grep('OTU_',colnames(crnotu.p))]
(da.name2 <- 'OLR_R3_otu_verc_D47')

### DEV analysis for microbial OTU
ds <- otu.p[ind,]
nin <- nrow(ds) # library sample size


### Determine the optimal DEV parameters using cross-validation ##
if(T){
  zeroint=apply(ds,2,zeropad)# find the length of zero-padding region in the time series
  
  dev_p  <-  NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    wd.i <- 15
    if(zeroint[i]>wd.i){wd.i=zeroint[i]+2}
    param.opt.a <- dev_cvls_v2(ds[,i], E = seq(1,10,by=1), window_size = seq(wd.i, 40, by=2) )
    (para.optt <- param.opt.a[[1]])
    dev_p <- rbind(dev_p, c(variable=colnames(ds)[i], para.optt))
    rm(list = c('param.opt.a','para.optt'));gc()
    time.used <- proc.time() - t.begin
    cat("effect", i, "in sample completed:", time.used[3],"sec\n")
  }
  if(SaveFile) write.csv(dev_p, paste('./Output/', da.name2,'_DEV_Optparameter.csv',sep=''),row.names=F)
}

##################################################
### Calculate DEV based on the optimal parameters##
if(T){
  dev_otu_p <- read.csv(paste('./Output/', da.name2,'_DEV_Optparameter.csv',sep=''),header=T)
  dev_p <- dev_otu_p
  DEVresult <-NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    para.optt <- as.numeric(dev_p[i,-1])
    names(para.optt) <- colnames(dev_p)[-1]
    ### Performing DEV algorithm
    devts <- DEV_perform_ls2(ds[,i],param.opt=para.optt)
    DEVresult <- rbind(DEVresult, data.frame(variable=colnames(ds)[i], devts))
    rm(list = c('devts'));gc()
    time.used <- proc.time() - t.begin
    cat("effect", i, "in sample completed:", time.used[3],"sec\n")
  }
  if(SaveFile) write.csv(DEVresult,paste('./Output/', da.name2,'_DEV_TimeSeries.csv',sep=''),row.names=F)
}

dev.crn <- read.csv(paste('./Output/', da.name1,'_DEV_TimeSeries.csv',sep=''),header=T)
dev.otu <- read.csv(paste('./Output/', da.name2,'_DEV_TimeSeries.csv',sep=''),header=T)

dev.crnotu <- rbind(dev.crn, dev.otu)
if(SaveFile) write.csv(dev.crnotu, './Output/DEV_D47_TS_crnotu.csv',row.names=F)


################################################################################
# DEV analysis for network properties
################################################################################
#### Data Preparation 
netprop.o <- read.csv('./Output/demo_network_prop_ccm_perturbed_D47.csv',header=T)
netprop.p <- filter(netprop.o, Reactor=='Perturbed')

(da.name3 <- 'OLR_R3_network_D47')

do <- matrix(NA, 46, NROW(netprop.p))
for(i in 1:NROW(netprop.p)) do[,i] <- as.numeric(netprop.p[i, -c(1:4)])

nin <- nrow(do) # library sample size

if(interpolation){
  do.na <- apply(do,2,EDM_NA_fill)
  do <- do.na
}
colnames(do) <- netprop.p[,'Property']

ds <- do[ind, ] #

### Determine the optimal DEV parameters using cross-validation ##
if(T){
  zeroint=apply(do,2,zeropad)# find the length of zero-padding region in the time series
  dev_p  <-  NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    wd.i <- 15
    if(zeroint[i]>wd.i){wd.i=zeroint[i]+2}
    param.opt.a <- dev_cvls_v2(ds[,i], E = seq(1,10,by=1), window_size = seq(wd.i, 40, by=2) )
    (para.optt <- param.opt.a[[1]])
    dev_p <- rbind(dev_p, c(variable=colnames(ds)[i], para.optt))
    rm(list = c('param.opt.a','para.optt'));gc()
    time.used <- proc.time() - t.begin
    cat("effect", i, "in sample completed:", time.used[3],"sec\n")
  }
  if(SaveFile) write.csv(dev_p, paste('./Output/', da.name3,'_DEV_Optparameter.csv',sep=''),row.names=F)
}

##################################################
### Calculate DEV based on the chosen parameters##
if(T){
  dev_net_p <- read.csv(paste('./Output/', da.name3,'_DEV_Optparameter.csv',sep=''),header=T)
  dev_p <- dev_net_p
  DEVresult <-NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    para.optt <- as.numeric(dev_p[i,-1])
    names(para.optt) <- colnames(dev_p)[-1]
    ### Performing DEV algorithm
    devts <- DEV_perform_ls2(ds[,i],param.opt=para.optt)
    DEVresult <- rbind(DEVresult, data.frame(variable=colnames(ds)[i], devts))
    rm(list = c('param.opt.a', 'devts'));gc()
    time.used <- proc.time() - t.begin
    cat("effect", i, "in sample completed:", time.used[3],"sec\n")
  }
  if(SaveFile) write.csv(DEVresult,paste('./Output/', da.name3,'_DEV_TimeSeries.csv',sep=''),row.names=F)
}

dev.net <- read.csv(paste('./Output/', da.name3,'_DEV_TimeSeries.csv',sep=''),header=T)
