############################################################
### DEV analysis
rm(list = ls());gc()
# Set path
OPath <- 'D:/data/Cooperation/Prodigio/Demo_RScript'
setwd(OPath)
source('./Code/Function_DEV.R')

library(rEDM) # EDM analyses (e.g., simplex projection and CCM) were conducted in rEDM package ver. 1.2.3
library(MASS)
library(dplyr)

seed <- 49563
set.seed(seed)
SaveFile <- F #T/F for saving files
interpolation <- T

################################################################################
# DEV analysis for metabolites
################################################################################
#### Data Preparation
crnotu.p <- read.csv('./Data/Demo_Perturbed_otucrn_abundance.csv',header=T,stringsAsFactors = F)

crn.p <- crnotu.p[,-grep('OTU_',colnames(crnotu.p))]

if(interpolation){
  crn.p.na <- apply(crn.p,2,EDM_NA_fill)
  crn.p <- crn.p.na
}

(da.name1 <- 'OLR_R3_crn_verc_ls_sfill')
do <- crn.p
nin <- nrow(do) # library sample size

# Standardization
ds <- apply(do,2,scale) #

# Identify the length of zero-padding region in the time series
zeroint=apply(do,2,zeropad)

if(F){
  dev_p  <-  NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    wd.i <- 15
    if(zeroint[i]>wd.i){wd.i=zeroint[i]+2}
    param.opt.a <- dev_cvls(ds[,i], window_size = seq(wd.i, 30, by=1))
    (para.optt <- param.opt.a[[1]])
    dev_p <- rbind(dev_p, c(variable=colnames(ds)[i], para.optt))
    rm(list = c('param.opt.a','para.optt'));gc()
    time.used <- proc.time() - t.begin
    cat("effect", i, "in sample completed:", time.used[3],"sec\n")
  }
  if(SaveFile) write.csv(dev_p, paste('./Output/', da.name1,'_DEV_Optparameter.csv',sep=''),row.names=F)
}

##################################################
### Note that S-Map calculations for large window sizes and/or time series need a lot of computational time
### start algorithm
if(F){
  dev_crn_p <- read.csv(paste('./Output/', da.name1,'_DEV_Optparameter.csv',sep=''),header=T)
  dev_p <- dev_crn_p
  DEVresult <-NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    para.optt <- as.matrix(dev_p[i,-1])
    ### Performing DEV algorithm
    devts <- DEV_perform_ls(ds[,i],param.opt=para.optt)
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
(da.name2 <- 'OLR_R3_otu_verc_ls_s')

### DEV analysis for microbial OTU
do <- otu.p
nin <- nrow(do) # library sample size

# Standardization
ds <- apply(do,2,scale) #

# Identify the length of zero-padding region in the time series
zeroint=apply(do,2,zeropad)

if(F){
  dev_p  <-  NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    wd.i <- 15
    if(zeroint[i]>wd.i){wd.i=zeroint[i]+2}
    param.opt.a <- dev_cvls(ds[,i], window_size = seq(wd.i, 30, by=1))
    (para.optt <- param.opt.a[[1]])
    dev_p <- rbind(dev_p, c(variable=colnames(ds)[i], para.optt))
    rm(list = c('param.opt.a','para.optt'));gc()
    time.used <- proc.time() - t.begin
    cat("effect", i, "in sample completed:", time.used[3],"sec\n")
  }
  if(SaveFile) write.csv(dev_p, paste('./Output/', da.name2,'_DEV_Optparameter.csv',sep=''),row.names=F)
}

##################################################
### Note that S-Map calculations for large window sizes and/or time series need a lot of computational time
### start algorithm
if(F){
  dev_otu_p <- read.csv(paste('./Output/', da.name2,'_DEV_Optparameter.csv',sep=''),header=T)
  dev_p <- dev_otu_p
  DEVresult <-NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    para.optt <- as.matrix(dev_p[i,-1])
    ### Performing DEV algorithm
    devts <- DEV_perform_ls(ds[,i],param.opt=para.optt)
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
if(SaveFile) write.csv(dev.crnotu, './Output/DEV_ls_sfill_TS_crnotu.csv',row.names=F)

################################################################################
# DEV analysis for network properties
################################################################################
#### Data Preparation 
netprop.o <- read.csv('./Output/demo_network_prop_perturbed_new.csv',header=T)
netprop.p <- filter(netprop.o, Reactor=='Perturbed')

(da.name3 <- 'OLR_R3_network_ls_sfill_new')

do <- matrix(NA, 45, NROW(netprop.p))
for(i in 1:NROW(netprop.p)) do[,i] <- as.numeric(netprop.p[i, -c(1:4)])

nin <- nrow(do) # library sample size

if(interpolation){
  do.na <- apply(do,2,EDM_NA_fill)
  do <- do.na
}
colnames(do) <- netprop.p[,'Property']

# Standardization
ds <- apply(do,2,scale) #

# Identify the length of zero-padding region in the time series
zeroint=apply(do,2,zeropad)

if(F){
  dev_p  <-  NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    wd.i <- 15
    if(zeroint[i]>wd.i){wd.i=zeroint[i]+2}
    param.opt.a <- dev_cvls(ds[,i], window_size = seq(wd.i, 30, by=1))
    (para.optt <- param.opt.a[[1]])
    dev_p <- rbind(dev_p, c(variable=colnames(ds)[i], para.optt))
    rm(list = c('param.opt.a','para.optt'));gc()
    time.used <- proc.time() - t.begin
    cat("effect", i, "in sample completed:", time.used[3],"sec\n")
  }
  if(SaveFile) write.csv(dev_p, paste('./Output/', da.name3,'_DEV_Optparameter.csv',sep=''),row.names=F)
}

##################################################
### Note that S-Map calculations for large window sizes and/or time series need a lot of computational time
### start algorithm
if(F){
  dev_net_p <- read.csv(paste('./Output/', da.name3,'_DEV_Optparameter.csv',sep=''),header=T)
  dev_p <- dev_net_p
  DEVresult <-NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time()
    ##Performing cross-validation for selecting the optimal parameter set
    para.optt <- as.matrix(dev_p[i,-1])
    ### Performing DEV algorithm
    devts <- DEV_perform_ls(ds[,i],param.opt=para.optt)
    DEVresult <- rbind(DEVresult, data.frame(variable=colnames(ds)[i], devts))
    rm(list = c('param.opt.a', 'devts'));gc()
    time.used <- proc.time() - t.begin
    cat("effect", i, "in sample completed:", time.used[3],"sec\n")
  }
  if(SaveFile) write.csv(DEVresult,paste('./Output/', da.name3,'_DEV_TimeSeries.csv',sep=''),row.names=F)
}

dev.net=read.csv(paste('./Output/', da.name3,'_DEV_TimeSeries.csv',sep=''),header=T)
