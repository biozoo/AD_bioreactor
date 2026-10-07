# R code for performing DEV early warning analysis in Fig. 1
rm(list = ls());gc()
# Set path
OPath <- 'C:/data/Cooperation/Prodigio/Code/Demo_RScript_20250529'
setwd(OPath)

source('./Code/Function_DEV_V2.R')

library(rEDM) # Some EDM analyses (e.g., simplex projection and CCM) were conducted in rEDM package ver. 1.2.3
library(dplyr)
library(MASS)
library(mgcv)
library(gratia)
seed <- 49563
set.seed(seed)
SaveFile <- F #T/F for saving files
interpolation <- T

tipp <- 38+1
tipp2 <- 44+1
nd <- 46+1
#############################################################
#### Data Preparation 
(da.name <- 'DEV_CH4_indR_compiled')
reactor_ch4 <- read.csv( './Data/Reactor_data_envi.csv', header=T)
richness <- read.csv( './Data/OTU_richness.csv', header=T)
vname <- c('CH4',"COD_rm")
reaname <- c('C1', 'O1', 'C2', 'O2')
vname2 <- c(paste(rep(vname,4),rep(reaname,each=2),sep='_'))

do <- reactor_ch4[,vname2]#apply(reactor_ch4[,vname2],2,scale)#  

nin <- nrow(do) # library sample size

if(interpolation) {
  for(i in 1:NCOL(do))  do[,i] <- EDM_NA_fill(do[,i])
}

# Standardization
ind <- 1:tipp2#nd#
do1 <- do[ind, 1:length(reaname)]; do2 <- do[ind, (1+length(reaname)):(2*length(reaname))]
colnames(do1) <- colnames(do2) <- unlist(strsplit(colnames(do)[1:4],'1')) 

# Identify the length of zero-paddings region in the time series
zeroint <- apply(rbind(do1,rep(NA,length(vname)),do2),2,zeropad)
##################################################
### DEV analysis
### Note that S-Map calculations for large window sizes and/or time series need a lot of computational time
### start algorithm
DEVresult = para.optR = NULL
for(i in 1:ncol(do1)){
  t.begin <- proc.time()
  
  ##Performing cross-validation for selecting the optimal parameter set
  wd.i <- 15
  if(zeroint[i]>=(wd.i-1)){wd.i=zeroint[i]+2}
  param.opt.a <- dev_cvls_v2(do1[,i], E = seq(1,10,by=1),
                             window_size = seq(wd.i, 40, by=2), replicate_data= do2[,i])
  para.optt <- param.opt.a[[1]]
  param.opt.a2 <- dev_cvls_v2(do2[,i], E = seq(1,10,by=1),
                             window_size = seq(wd.i, 40, by=2), replicate_data= do1[,i])
  para.optt2 <- param.opt.a2[[1]]
  ### Performing DEV algorithm
  devts <- DEV_perform_ls2(do1[,i],param.opt=para.optt, replicate_data= do2[,i], pred_rep1_only=T)
  devts2 <- DEV_perform_ls2(do2[,i],param.opt=para.optt2, replicate_data= do1[,i], pred_rep1_only=T)
  DEVresult <- rbind(DEVresult, data.frame(variable=paste(colnames(do1)[i],'1',sep=''), devts),
                     data.frame(variable=paste(colnames(do2)[i],'2',sep=''), devts2))
  para.optR <- rbind(para.optR, c(variable=paste(colnames(do1)[i],'1',sep=''), para.optt),
                     c(variable=paste(colnames(do2)[i],'2',sep=''), para.optt))
  
  rm(list = c('param.opt.a', 'devts')); gc()
  time.used <- proc.time() - t.begin
  cat("effect", i, "in sample completed:", time.used[3],"sec\n")
}

do.r <- richness #apply(richness,2,scale)
if(interpolation) {
  for(i in 1:NCOL(do.r))  do.r[,i] <- EDM_NA_fill(do.r[,i])
}

DEVresult.r = para.optR.r = NULL
for(i in 1:ncol(do.r)){
  t.begin <- proc.time()
  
  ##Performing cross-validation for selecting the optimal parameter set
  wd.i <- 15
  if(zeroint[i]>=(wd.i-1)){wd.i=zeroint[i]+2}
  param.opt.a <- dev_cvls_v2(do.r[ind,i], E = seq(1,10,by=1),
                             window_size = seq(wd.i, 40, by=2))
  para.optt <- param.opt.a[[1]]
  
  ### Performing DEV algorithm
  devts <- DEV_perform_ls2(do.r[ind,i], param.opt = para.optt)
  DEVresult.r <- rbind(DEVresult.r, data.frame(variable=colnames(do.r)[i], devts))
  para.optR.r <- rbind(para.optR.r, c(variable=colnames(do.r)[i], param.opt.a[[1]]))
  
  rm(list = c('param.opt.a', 'devts')); gc()
  time.used <- proc.time() - t.begin
  cat("effect", i, "in sample completed:", time.used[3],"sec\n")
}

DEVresult <- rbind(DEVresult, DEVresult.r)
para.optR <- rbind(para.optR, para.optR.r)

colnames(para.optR)[-1] <- colnames(para.optt)
para.optR <- data.frame(variable=unlist(para.optR[,1]), apply(para.optR[,-1],2,as.numeric))
DEVresult[,'t_max'] <-DEVresult[,'t_max'] -1
if(SaveFile){
  write.csv(DEVresult,paste('./Output/',da.name,'_DEV_TimeSeries.csv',sep=''), row.names=F)
  write.csv(para.optR,paste('./Output/',da.name,'_DEV_Optparameter.csv',sep=''), row.names=F)
}  

# Plot the proxy of organic loading
yrg1 <- range(c(reactor_ch4[,c('CODt_C1','CODt_C2','CODt_O1','CODt_O1')]),na.rm=T)
win.graph(160,50)
par(mfrow=c(1,2),mar=c(4,4,1,4))
plot(CODt_C2~Day, data=reactor_ch4, type='o',ylab='TCOD (mg/L)',xlab='Day',ylim=yrg1,col='grey',xlim=c(0,(nd-1)),
     main='Control reactor',lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
lines(CODt_C1~Day, data=reactor_ch4, type='o',lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)

plot(CODt_O2~Day, data=reactor_ch4, type='o',ylab='',xlab='Day',ylim=yrg1,col='grey',xlim=c(0,(nd-1)),
     main='Perturbed reactor',lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
lines(CODt_O1~Day, data=reactor_ch4, type='o',lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(h=1,col='grey')
polygon(x=c(tipp,tipp2,tipp2,tipp)-1,y=rep(yrg1*c(0.5,1.5),each=2), col="orange",density=15, border="white")
abline(v=tipp-1, col="orange")

# Plot the CH4 production and DEV
###############################
##### Methane production (ml/d)
###############################
win.graph(200,90)
par(mfrow=c(1,2),mar=c(4,4,1,4))
# Control reactors
name.t <- paste("CH4",reaname[c(1,3)],sep='_')
yrg2 <- range(c(do[,paste("CH4",reaname,sep='_')]),na.rm=T)

plot(do[,name.t[2]]~c(0:(nd-1)),type='o',ylab='CH4 production (ml/d)',xlab='Day',ylim=yrg2, col='grey',
     main='',xlim=c(0,(nd-1)),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
lines(do[,name.t[1]]~c(0:(nd-1)),lwd=2,pch=16, cex=1.5,type='o',ylim=yrg2)
par(new=T)
dev.t <- filter(DEVresult,variable==name.t[2])
plot(Abs~I(t_max),dev.t, type='o',col='lightblue',xlim=c(0,(nd-1)),ylim=range(c(filter(DEVresult,variable%in%name.t)[,'Abs'],1),na.rm=T),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
mtext('',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='', col.axis = "blue", cex.axis=1.5)
dev.t <- filter(DEVresult,variable==name.t[1])
lines(Abs~I(t_max),dev.t, type='o',col='blue',xlim=c(0,(nd-1)),lwd=2,pch=16, cex=1.5)
abline(h=1,col='grey')

# Perturbed reactors
name.t <- paste("CH4",reaname[c(2,4)],sep='_')

plot(do[,name.t[2]]~c(0:(nd-1)),type='o',ylab='CH4 production (ml/d)',xlab='Day',ylim=yrg2, col='grey',
     main='',xlim=c(0,(nd-1)),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
lines(do[,name.t[1]]~c(0:(nd-1)),lwd=2,pch=16, cex=1.5,type='o',ylim=yrg2)

par(new=T)
dev.t <- filter(DEVresult,variable==name.t[2])
plot(Abs~I(t_max),dev.t, type='o',col='lightblue',xlim=c(0,(nd-1)),ylim=c(0.9, 1.1),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
#range(c(filter(DEVresult,variable%in%name.t)[,'Abs'],1),na.rm=T)
mtext('',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='', col.axis = "blue", cex.axis=1.5)

dev.t <- filter(DEVresult,variable==name.t[1])
lines(Abs~I(t_max),dev.t, type='o',col='blue',xlim=c(0,(nd-1)),lwd=2,pch=16, cex=1.5)
polygon(x=c(tipp,tipp2,tipp2,tipp)-1,y=c(-10,-10,10,10), col="orange",density=15, border="white")
abline(v=tipp-1, col="orange")

abline(h=1,col='grey')

# Plot the COD removals and DEV individual reactors
win.graph(200,90)
par(mfrow=c(1,2),mar=c(4,4,1,4))
name.t <- paste("COD_rm",reaname[c(1,3)],sep='_')
yrg2 <- range(c(do[,paste("COD_rm",reaname,sep='_')]),na.rm=T)

plot(do[,name.t[2]]~c(0:(nd-1)),type='o',ylab='COD removal efficiency',xlab='Day',ylim=yrg2, col='grey',
     main='',xlim=c(0,(nd-1)),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
lines(do[,name.t[1]]~c(0:(nd-1)),lwd=2,pch=16, cex=1.5,type='o',ylim=yrg2)

name.t <- paste("COD_rm",reaname[c(2,4)],sep='_')

plot(do[,name.t[2]]~c((1:nd)-1),type='o',ylab='COD removal efficiency',xlab='Day',ylim=yrg2, col='grey',
     main='',xlim=c(0,(nd-1)),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
lines(do[,name.t[1]]~c(0:(nd-1)),lwd=2,pch=16, cex=1.5,type='o',ylim=yrg2)

par(new=T)
dev.t <- filter(DEVresult,variable==name.t[2])
plot(Abs~I(t_max),dev.t, type='o',col='lightblue',xlim=c(0,(nd-1)),ylim=range(c(filter(DEVresult,variable%in%name.t)[,'Abs'],1),na.rm=T),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
mtext('',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='', col.axis = "blue", cex.axis=1.5)
dev.t <- filter(DEVresult,variable==name.t[1])
lines(Abs~I(t_max),dev.t, type='o',col='blue',xlim=c(0,(nd-1)),lwd=2,pch=16, cex=1.5)

abline(h=1,col='grey')
polygon(x=c(tipp,tipp2,tipp2,tipp)-1,y=c(-10,-10,10,10), col="orange",density=15, border="white")
abline(v=tipp-1, col="orange")

#####################
##### OTU richness
#####################
# Control reactors
win.graph(200,90)
par(mfrow=c(1,2),mar=c(4,4,1,4))
yrg2 <- range(c(do.r),na.rm=T)
plot(do.r[,1]~c(0:(nd-1)),type='o',ylab='OTU richness',xlab='Day',ylim=yrg2, 
     main='',xlim=c(0,(nd-1)),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
par(new=T)
dev.t <- filter(DEVresult,variable==colnames(do.r)[1])
plot(Abs~I(t_max),dev.t, type='o',col='blue',xlim=c(0,(nd-1)),ylim=range(c(filter(DEVresult,variable%in%colnames(do.r)[1])[,'Abs'],1),na.rm=T),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
mtext('',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='', col.axis = "blue", cex.axis=1.5)
abline(h=1,col='grey')

# Perturbed reactors
plot(do.r[,2]~c((1:nd)-1),type='o',ylab='OTU richness',xlab='Day',ylim=yrg2, 
     main='',xlim=c(0,(nd-1)),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
par(new=T)
dev.t <- filter(DEVresult,variable==colnames(do.r)[2])
plot(Abs~I(t_max),dev.t, type='o',col='blue',xlim=c(0,(nd-1)),ylim=range(c(filter(DEVresult,variable%in%colnames(do.r)[2])[,'Abs'],1),na.rm=T),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
mtext('',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='', col.axis = "blue", cex.axis=1.5)

abline(h=1,col='grey')
polygon(x=c(tipp,tipp2,tipp2,tipp)-1,y=c(-10,-10,10,10), col="orange",density=15, border="white")
abline(v=tipp-1, col="orange")

## Automatic detection of bifurcation event  
udev3 <- unique(DEVresult[,'variable'])
dev.mx <- matrix(NA,nd,length(udev3))
for(i in 1:length(udev3)){
  da.t <- filter(DEVresult, variable==udev3[i])
  dev.mx[(0:(nd-1))%in%da.t[,'t_max'],i] <- da.t[,'Abs']
}
colnames(dev.mx) <- udev3

dev.mx <- dev.mx[, apply(!is.na(dev.mx),2,sum)>10]

crip.ch4 <- NULL
ind.tip <- tipp2;#
for(i in 1:NCOL(dev.mx)){
  crip.i <- NULL
  da.t2 <- data.frame(time=1:ind.tip,value=dev.mx[1:ind.tip,i])
  crip.i <- find.bif2(da.t = da.t2, ind.tip, ploty=F)
  if(NROW(crip.i)>0) crip.ch4=rbind(crip.ch4, data.frame(variable=rep(i),crip.i))
}  

(crip.ch4 <- data.frame(ID=colnames(dev.mx)[crip.ch4[,1]],crip.ch4))


################################################################################
## Tipping point of methane production is defined as the first statistically
## significant transition after the bifurcation warning time, based on the
## smoothed trajectory and its first-derivative 95% confidence interval.
## The operational tipping point is defined as the first sampling day after
## the estimated transition time.
ch4 <- do[,'CH4_O1']
da.t <- data.frame(time=1:NROW(do), value=ch4)
g.t1 <- gam(value~s(time),data=da.t)
newd <- seq(0,46,0.01)

## Estimate the smoothed methane trajectory and its first derivative
yhat <- predict.gam(g.t1, newdata=data.frame(time=newd))
fd <- fderiv(g.t1, newdata = data.frame(time=newd))

## Calculate 95% confidence intervals for the first derivative
ci <- as.matrix(confint(fd, type = "confidence")[,c('lower','est','upper')])
yhatd1 <- ci[,'est']
smo <- cbind(time=newd, yhat, ci)

## Detect statistically significant transitions by identifying sign changes
## in the first derivative whose confidence interval excludes zero
pk.i <- data.frame(signch2(yci=smo[,c('lower','est','upper')],x=unlist(c(newd))))

## Identify the first statistically significant transition occurring after the identified bifurcation event
## Here, the transition time is defined as the point at which the first derivative 
## becomes significantly negative following a bifurcation event identified by DEV analysis 
## (see Early warning signal analysis based on the dynamical eigenvalue (DEV) analysis)

(tipping <- filter(pk.i, sign!=0 & (ts_i>=filter(crip.ch4,ID=='CH4_O1')[,'BFwarn']))[1,'ts_i'])

## Define the operational tipping point as the first sampling day after the estimated transition time
(tipping.integer <- ceiling(tipping))

