# R code for performing DEV early warning analysis in Fig. 1
rm(list = ls());gc()
# Set path
OPath <- 'D:/data/Cooperation/Prodigio/Demo_RScript'
setwd(OPath)

source('./Code/Function_DEV.R')

library(rEDM) # Some EDM analyses (e.g., simplex projection and CCM) were conducted in rEDM package ver. 1.2.3
library(dplyr)
library(MASS)
library(mgcv)
library(gratia)
seed <- 49563
set.seed(seed)
SaveFile <- T #T/F for saving files
interpolation <- T

tipp <- 37
tipp2 <- 43
nd <- 46
#############################################################
#### Data Preparation 
(da.name <- 'DEV_CH4_ls')
reactor_ch4_ews <- read.csv( './Data/reactor_ch4_ews.csv', header=T)

do <- reactor_ch4_ews[,c('CH4_control', 'CH4_perturbed', 'OTU_Richness_control', 'OTU_Richness_perturbed')]
nin <- nrow(do) # library sample size

if(interpolation) do[,'OTU_Richness_control'] <- EDM_NA_fill(reactor_ch4_ews[,c('OTU_Richness_control')])

# Standardization
ds <- apply(do,2,scale) #

# Identify the length of zero-paddings region in the time series
zeroint=apply(do,2,zeropad)
##################################################
### DEV analysis
### Note that S-Map calculations for large window sizes and/or time series need a lot of computational time
### start algorithm
DEVresult = para.optR = NULL
for(i in 1:ncol(ds)){
  t.begin <- proc.time()
  
  ##Performing cross-validation for selecting the optimal parameter set
  wd.i=15
  if(zeroint[i]>wd.i){wd.i=zeroint[i]+1}
  param.opt.a <- dev_cvls(ds[,i], window_size = seq(wd.i, 30, by=1))
  para.optt=as.matrix(param.opt.a[[1]])
  
  ### Performing DEV algorithm
  devts=DEV_perform_ls(ds[,i],param.opt=para.optt)
  DEVresult = rbind(DEVresult, data.frame(variable=colnames(ds)[i], devts))
  para.optR = rbind(para.optR, c(variable=colnames(ds)[i], para.optt))
  rm(list = c('param.opt.a', 'devts'));gc()
  write.csv(DEVresult,paste('./Output/',da.name,'_DEV_TimeSeries.csv',sep=''),row.names=F)
  write.csv(para.optR,paste('./Output/',da.name,'_DEV_Optparameter.csv',sep=''),row.names=F)
  time.used <- proc.time() - t.begin
  cat("effect", i, "in sample completed:", time.used[3],"sec\n")
}


################################################################################
## (Fig. 1) Early warning signal analysis for CH4 production 
################################################################################
# Plot the proxy of organic loading
yrg1 <- range(c(reactor_ch4_ews[,c('TCOD_control','TCOD_perturbed')]),na.rm=T)
win.graph(160,50)
par(mfrow=c(1,2),mar=c(4,4,1,4))
plot(TCOD_control~Day, data=reactor_ch4_ews, type='o',ylab='TCOD (mg/L)',xlab='Day',ylim=yrg1,
     main='Control reactor',xlim=c(1,(nd+1)),lwd=2,pch=1, cex=1.5, cex.lab=1.5, cex.axis=1.5)
plot(TCOD_perturbed~Day, data=reactor_ch4_ews, type='o',ylab='TCOD (mg/L)',xlab='Day',ylim=yrg1,
     main='Perturbed reactor',xlim=c(1,(nd+1)),lwd=2,pch=1, cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(h=1,col='grey')
polygon(x=c(tipp,tipp2,tipp2,tipp),y=rep(yrg1*c(0.5,1.5),each=2), col="orange",density=15, border="white")

# Plot the CH4 production and DEV
win.graph(160,70)
par(mfrow=c(1,2),mar=c(4,4,1,4))
yrg2 <- range(c(do[,c('CH4_control','CH4_perturbed')]),na.rm=T)

plot(do[,'CH4_control']~c(1:nd),type='o',ylab='CH4 (%)',xlab='Day',ylim=yrg2,
     main='Control reactor',xlim=c(1,nd),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
par(new=T)
plot(Abs~I(t_max),filter(DEVresult,variable=='CH4_control'), type='o',col='blue',xlim=c(1,nd),ylim=c(0.4,8),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
mtext('',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='', col.axis = "blue", cex.axis=1.5)
abline(h=1,col='grey')

plot(do[,'CH4_perturbed']~c(1:nd), type='o',ylab='',xlab='Day',ylim=yrg2,
     main='Perturbed reactor',xlim=c(1,nd),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
par(new=T)
plot(Abs~I(t_max),filter(DEVresult,variable=='CH4_perturbed'), type='o',col='blue',xlim=c(1,nd),ylim=c(0.4,1.17),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
mtext('DEV',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='DEV', col.axis = "blue", cex.axis=1.5)
abline(h=1,col='grey')
polygon(x=c(tipp,tipp2,tipp2,tipp),y=c(-10,-10,10,10), col="orange",density=15, border="white")

################################################################################
## Microbial OTU Richness (Fig. S1)
win.graph(160,70)
par(mfrow=c(1,2),mar=c(4,4,1,4))
yrg2 <- range(c(do[,c('OTU_Richness_control','OTU_Richness_perturbed')]),na.rm=T)
yrgd <- range(c(filter(DEVresult,variable%in%c('OTU_Richness_control','OTU_Richness_perturbed'))[,'Abs']),na.rm=T)

plot(do[,'OTU_Richness_control']~c(1:nd),type='o',ylab='Microbial OTU Richness',xlab='Day',ylim=yrg2,
     main='Control reactor',xlim=c(1,nd),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
par(new=T)
plot(Abs~I(t_max),filter(DEVresult,variable=='OTU_Richness_control'), type='o',
     ylim=yrgd, col='blue',xlim=c(1,nd),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
mtext('',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='', col.axis = "blue", cex.axis=1.5)
abline(h=1,col='grey')

plot(do[,'OTU_Richness_perturbed']~c(1:nd), type='o',ylab='',xlab='Day',ylim=yrg2,
     main='Perturbed reactor',xlim=c(1,nd),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
par(new=T)
plot(Abs~I(t_max),filter(DEVresult,variable=='OTU_Richness_perturbed'), type='o',
     ylim=yrgd, col='blue',xlim=c(1,nd),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
mtext('DEV',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='DEV', col.axis = "blue", cex.axis=1.5)
abline(h=1,col='grey')
polygon(x=c(tipp,tipp2,tipp2,tipp),y=c(-10,-10,10,10), col="orange",density=15, border="white")



## Automatic detection of bifurcation event  
udev3 <- unique(DEVresult[,'variable'])
dev.mx <- matrix(NA,nd,length(udev3))
for(i in 1:length(udev3)){
  da.t <- filter(DEVresult, variable==udev3[i])
  dev.mx[(1:nd)%in%da.t[,'t_max'],i] <- da.t[,'Abs']
}
colnames(dev.mx) <- udev3

windows()
crip.ch4 <- NULL
ind.tip <- tipp;#nd;#
for(i in 1:NCOL(dev.mx)){
  crip.i=NULL
  da.t2 = data.frame(time=1:ind.tip,value=dev.mx[1:ind.tip,i])
  crip.i <- find.bif2(da.t = da.t2, ind.tip, ploty=T)
  if(NROW(crip.i)>0) crip.ch4=rbind(crip.ch4, data.frame(variable=rep(i),crip.i))
}  

(crip.ch4 <- data.frame(ID=colnames(dev.mx)[crip.ch4[,1]],crip.ch4))


