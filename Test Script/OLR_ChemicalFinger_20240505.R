# R code for performing MDR S-map anaysis
rm(list = ls());gc()
# Set path
Root_Path <- 'D:'
Sub_Path <- '\\data\\Cooperation\\Prodigio\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)
source(paste(Root_Path ,'\\data\\Feitsui\\other.r',sep=''))
source(paste(Root_Path ,Sub_Path,'\\Code\\DEV_function_20230601.R',sep=''))

library(vegan)
library(rEDM) # Some EDM analyses (e.g., simplex projection and CCM) were conducted in rEDM package ver. 1.2.3
library(doParallel)
library(parallel)
library(foreach)
library(Kendall)
library(MASS)
library(dplyr)
library(glmnet)
library(zyp)
library(Kendall)
library(igraph)
library(NADA2)
library(diffusionMap)
library(earlywarnings)
library(Rtsne)
library(mdatools)
library(bigmemory)

source(paste(Root_Path,'\\data\\Cooperation\\Pear\\Code\\Demo_MDR_function_inferM.R',sep=''))
source(paste(Root_Path,'\\data\\Cooperation\\Pear\\Code\\Demo_MDR_function.R',sep=''))
seed <- 49563
set.seed(seed)
SaveFile <- T #T/F for saving files
interp <- T 

mag_scale=function(x,mgf=0.5){
  if(prod(x)>=0){
    if(x[2]>0){y=x;y[1]=y[1]*(1.0-mgf);y[2]*(1.0+mgf)}else{
      y=x;y[1]=y[1]*(1.0+mgf);y[2]*(1.0-mgf)}
  }else{y=x*(1.0+mgf)}
  return(y)
}

envi12 <- read.csv(paste(Root_Path,Sub_Path,'\\Data\\R1_chemical fingerprints_final_XK_20231007.csv',sep=''),header=T,stringsAsFactors = F)
envi34 <- read.csv(paste(Root_Path,Sub_Path,'\\Data\\R3R4_chemical fingerprints_final_XK_20231007.csv',sep=''),header=T,stringsAsFactors = F)
#otu12=read.csv('./Data/Mgnus_compo12_shock_tau1_23.csv',header=T)
#otu34=read.csv('./Data/Mgnus_compo34_shock_tau1_23.csv',header=T)
otu12=read.csv('./Data/OLR1_otu_verc_ch4p_Species_tau132.csv',header=T)
otu34=read.csv('./Data/OLR3_otu_verc_ch4p_Species_tau132.csv',header=T)

otu12.2=read.csv('./Data/OLR1_otu_ver2_Species_tau132.csv',header=T)
otu34.2=read.csv('./Data/OLR3_otu_ver2_Species_tau132.csv',header=T)

otu12.g=read.csv('./Data/Mgnus_compo12_Genus_tau1_23.csv',header=T)
otu34.g=read.csv('./Data/Mgnus_compo34_Genus_tau1_23.csv',header=T)
otu12.f=read.csv('./Data/Mgnus_compo12_Family_tau1_23.csv',header=T)
otu34.f=read.csv('./Data/Mgnus_compo34_Family_tau1_23.csv',header=T)

dor12 <- read.csv('.//Data//Table.EDM.Long.Shocks_NMBU_R1R2.csv',header=T,stringsAsFactors = F,sep=';',dec=',', fileEncoding="latin1")
dor34 <- read.csv('.//Data//Table.EDM.Long.Shocks_NMBU_R3R4.csv',header=T,stringsAsFactors = F,sep=';',dec=',', fileEncoding="latin1")
intervals=cbind(diff(dor12[,'Days']),diff(dor34[,'Days']))[-1,]
colnames(intervals)=c('R12','R34')

dat12=dsplit(dor12[,'Data.of.sampling'],nco=3,sepp='\\.',c('Day','Month','Year'))[,3:1]
dat34=dsplit(dor34[,'Data.of.sampling'],nco=3,sepp='\\.',c('Day','Month','Year'))[,3:1]
Date=as.Date(pastt(dat34,sepp='-'))
Days=as.numeric(Date-Date[1]+1)

shockd12=filter(dor12,Status%in%paste('Shock day', 1:46))
shockd34=filter(dor34,Status%in%paste('Shock day', 1:46))

Days=Days[dor34[,'Status']%in%paste('Shock day', 1:46)]
Date=Date[dor34[,'Status']%in%paste('Shock day', 1:46)]

envivb=c("T..ºC.","pH","redox","OLR..gCOD.Ld.","TS..mg.L.","ash..mg.L.",
         "VS..mg.L.","VS....","CO2....","CH4....","H2....","TCOD..mg.L.",      
         "SCOD..mg.L.","N.NH4...mg.L.","N.NH3..mg.L.","Qbiogas..mL.day.")

envivb34=c("Hac","Hpro","isoHBu","Hbu","isoHVal","Hval","Hca")

info12=shockd12[,1:5]
compo12=shockd12[,!colnames(shockd12)%in%c(envivb,colnames(info12))]  

info34=shockd34[,1:5]
compo34=shockd34[,!colnames(shockd34)%in%c(envivb,envivb34,colnames(info34))]  


olrall=read.csv(paste(Root_Path,Sub_Path,'\\Data\\OLR_shock_all_data.csv',sep=''),header=T,stringsAsFactors = F)

envi34 <- filter(envi34,reactor=='R3')
envi34 <- envi34[grep('Shock day',envi34[,'status']),]
#envi34[,-c(1:11)][is.na(envi34[,-c(1:11)])]=0
#envi34[apply(envi34[,-c(1:11)]==0,1,all),-c(1:11)]=NA

envi12 <- filter(envi12,reactor=='R1')
envi12 <- envi12[grep('Shock day',envi12[,'status']),]
envi12t=matrix(NA,NROW(envi34),NCOL(envi34));envi12t=data.frame(envi12t);colnames(envi12t)=colnames(envi34)
envi12t[as.Date(envi34[,'date'],'%m/%d/%Y')%in%as.Date(envi12[,'date'],'%m/%d/%Y'),]=envi12
envi12t[,c('date','status')]=envi34[,c('date','status')];
envi12t[,'reactor']=envi12t[1,'reactor'];
envi12=envi12t

#envi12[24,c("TS..mg.L.",'ash..mg.L.',"VS..mg.L.")]=envi12[24,c("TS..mg.L.",'ash..mg.L.',"VS..mg.L.")]/10


olr34=filter(olrall,Reactor=='R3')
olr34=olr34[as.Date(olr34[,'Date'],'%d/%m/%Y')%in%as.Date(envi34[,'date'],'%m/%d/%Y'),]

olr12=olr34;olr12[,-1]=NA
olr12t=filter(olrall,Reactor=='R1')
olr12[as.Date(olr12[,'Date'],'%d/%m/%Y')%in%as.Date(olr12t[,'Date'],'%d/%m/%Y'),-1]=olr12t[as.Date(olr12t[,'Date'],'%d/%m/%Y')%in%as.Date(olr12[,'Date'],'%d/%m/%Y'),-1]
olr12[,'Reactor']='R1';olr12[,'Type']='Control1'     


envi12[,c('CH4_percentage','CH4_mL_day')]=cbind(round(shockd12[,c('CH4....')],2), olr12[,c('mL.CH4.d')])
envi34[,c('CH4_percentage','CH4_mL_day')]=cbind(round(shockd34[,c('CH4....')],2), olr34[,c('mL.CH4.d')])
olr12[,c('CH4....')]=envi12[,c('CH4_percentage')]
olr34[,c('CH4....')]=envi34[,c('CH4_percentage')]
#envi34=rbind(envi34[1:4,],rep(NA),envi34[-c(1:4),])
#cbind(ch4p_olr=olr34[,'CH4....'],ch4_olr=olr34[,'mL.CH4.d'],ch4_finger=envi34[,'mL.CH4.d'])
nsp=NCOL(envi34)
nd=NROW(envi34)


###################
## Plot environmental variables
envivb2=c('pH','ORP','T..ºC.','NH4...mgN.L.','Qbiogas..mL.day.', 'mL.CH4.d')
tind=1:46
tipp=37
tipp2=43
win.graph(90,60);par(mfcol=c(2,3),mar=c(4,4,1,1))
for(i in 1:length(envivb2)){
  
  y.t=cbind(olr34,shockd34)[,envivb2[i]]
  y.t1=cbind(olr12,shockd12)[,envivb2[i]]
  plot(y.t~tind,type='o',xlab='Day',ylab=envivb2[i],pch=16,
       cex=2,lwd=2,ylim=range(c(y.t,y.t1),na.rm=T))
  lines(y.t1~tind,type='o',col='grey',pch=16,cex=2,lwd=2)
  polygon(x=c(tipp,tipp2,tipp2,tipp),y=rep(mag_scale(range(c(y.t,y.t1),na.rm=T)),each=2), col="orange",density=15, border="white")
  if(i==2)legend('topleft',c('Control','Experiment'),bty='n',lty=c(1,1),col=c('black','grey'),pch=c(16,1),cex=2)
}


if(F){
  win.graph(40,90);
  par(mfrow=c(3,1),mar=c(4,4,1,1))
  #for(i in 1:NCOL(stable.t)){
  for(i in c(1,2,4)){
    x.t=c((1:(tipp2-1))+0.5)
    y.t=stable.t[,i]
    y.t1=stable.t1[,i]
    if(i==1){x.t=x.t[-c(1:8)];y.t=y.t[-c(1:8)]}
    plot(y.t~x.t,ylab=colnames(stable.t)[i],type='o',xlab='Time',pch=16,cex=2,lwd=2,ylim=range(c(y.t,y.t1),na.rm=T))
    if(i!=1){
      lines(y.t1~x.t,type='o',pch=1,cex=2,lwd=2,col='grey')
    } 
    abline(h=1,lty=2)
    polygon(x=c(tipp,tipp2,tipp2,tipp),y=rep(mag_scale(range(c(y.t,y.t1),na.rm=T)),each=2), col="orange",density=15, border="white")
    legend('topleft',c('Control','Experiment'),bty='n',lty=c(1,1),col=c('black','grey'),pch=c(16,1))
  }
}




finger12=envi12[,-c(1:9)]
finger34=envi34[,-c(1:9)]

frich12=apply(finger12!=0,1,sum)
frich34=apply(finger34!=0,1,sum)
orich12=apply(compo12!=0,1,sum)
orich34=apply(compo34!=0,1,sum)
tipp=37;tipp2=43; tind.o=1:nd
win.graph(80,100);par(mfcol=c(2,1),mar=c(4,4,1,1))
yrg.t=range(c(orich12,orich34),na.rm=T)
plot(orich12~tind.o,type='o',ylab='OTU richness',xlab='Day',ylim=yrg.t,lwd=2)
lines(orich34~tind.o,type='o',col='blue',lwd=2)
polygon(x=c(tipp,tipp2,tipp2,tipp),y=rep(mag_scale(yrg.t),each=2), col="orange",density=15, border="white")

yrg.t=range(c(frich12,frich34),na.rm=T)
plot(frich12~tind.o,type='o',ylab='Chemical fingerprint richness',xlab='Day',ylim=yrg.t,lwd=2)
lines(frich34~tind.o,type='o',col='blue',lwd=2)
polygon(x=c(tipp,tipp2,tipp2,tipp),y=rep(mag_scale(yrg.t),each=2), col="orange",density=15, border="white")

richo=cbind(frich12,frich34,orich12,orich34, envi12[,c('CH4_percentage','CH4_mL_day')], envi34[,c('CH4_percentage','CH4_mL_day')])
colnames(richo)=c('frich12','frich34','orich12','orich34','ch4p12','ch4m12','ch4p34','ch4m34')
write.csv(richo,'olr_richness.csv',row.names = F)

ind.t34=which(apply(finger34!=0,2,sum,na.rm=T)>=((nd-sum(is.na(finger34[,1])))*(2/3)))
ind.t12=which(apply(finger12!=0,2,sum,na.rm=T)>=((nd-sum(is.na(finger12[,1])))*(2/3)))

ind.t=intersect(ind.t34,ind.t12)#ind.t34#

finger12=finger12[,ind.t]
finger34=finger34[,ind.t]

write.csv(finger12,'./Data/OLR1_crn_ver2_Species_tau132.csv',row.names=F)
write.csv(finger34,'./Data/OLR3_crn_ver2_Species_tau132.csv',row.names=F)



ind.t34=which(apply(finger34!=0,2,sum,na.rm=T)>=((nd-sum(is.na(finger34[,1])))*(2/3)))
ind.t12=which(apply(finger12!=0,2,sum,na.rm=T)>=((nd-sum(is.na(finger12[,1])))*(2/3)))

ada12=apply(cbind(envi12[,c('CH4_percentage','CH4_mL_day')],finger12),2,scale)
ada34=apply(cbind(envi34[,c('CH4_percentage','CH4_mL_day')],finger34),2,scale)
ada34.vfa=apply(cbind(envi34[,c('CH4_percentage','CH4_mL_day')],olr34[,c('HAc','HPro','i.HBu','HBu','iHVal','HVal','HCa')]),2,scale)

inferM=matrix(0, nsp+2, nsp+2);inferM[1:2,-c(1:2)]=1
inferM.vfa=matrix(0, NCOL(ada34.vfa), NCOL(ada34.vfa));inferM.vfa[1:2,-c(1:2)]=1
if(F){
  ccm.r34=ccm.fast.inferM(ds=ada34,inferM,Lmin=5,Lint=10)
  ccm.r34.vfa=ccm.fast.inferM(ds=ada34.vfa,inferM.vfa,Lmin=5,Lint=10)
  ccm.r12=ccm.fast.inferM(ds=ada12,inferM,Lmin=5,Lint=10)
  write.csv(ccm.r12[[1]],'ccm_finger_ver2_t12_w.csv',row.names=F)
  write.csv(ccm.r12[[2]],'ccm_finger_ver2_t12_p.csv',row.names=F)
  write.csv(ccm.r34[[1]],'ccm_finger_ver2_t34_w.csv',row.names=F)
  write.csv(ccm.r34[[2]],'ccm_finger_ver2_t34_p.csv',row.names=F)
  write.csv(ccm.r34.vfa[[1]],'ccm_VFA_ver2_t34_w.csv',row.names=F)
  write.csv(ccm.r34.vfa[[2]],'ccm_VFA_ver2_t34_p.csv',row.names=F)
  sum(ccm.r34[[2]][1,])
}

ccm.r34w=read.csv('ccm_finger_ver2_t34_w.csv',header=T)
ccm.r34s=read.csv('ccm_finger_ver2_t34_p.csv',header=T)
(spind.p=which(ccm.r34s[1,]==1)-2)
#(spind.ml=which(ccm.r34s[2,]==1)-2)

ccm.r12w=read.csv('ccm_finger_ver2_t12_w.csv',header=T)
ccm.r12s=read.csv('ccm_finger_ver2_t12_p.csv',header=T)
(spind.p1=which(ccm.r12s[1,]==1)-2)
#(spind.ml1=which(ccm.r12s[2,]==1)-2)

(spind=spind.p)#intersect(spind.p1,spind.p))

write.csv(finger12[,spind],'./Data/OLR1_finger_verc_ch4p_Species_tau132.csv',row.names=F)
write.csv(finger34[,spind],'./Data/OLR3_finger_verc_ch4p_Species_tau132.csv',row.names=F)

crn12=finger12[,spind]
crn34=cbind(olr34[,c('HAc', 'HPro', 'i.HBu', 'HBu', 'iHVal', 'HVal', 'HCa')],finger34[,spind])

write.csv(crn12,'./Data/OLR1_crn_verc_ch4p_Species_tau132.csv',row.names=F)
write.csv(crn34,'./Data/OLR3_crn_verc_ch4p_Species_tau132.csv',row.names=F)

if(F){
  crn12.v2=read.csv('./Data/OLR1_crn_verc_ch4p_Species_tau132.csv',header=T,fill=T)
  crn34.v2=read.csv('./Data/OLR3_crn_verc_ch4p_Species_tau132.csv',header=T,fill=T)
  crn12.v1=read.csv('./Data/OLR1_crn_tau132.csv',header=T,fill=T)
  crn34.v1=read.csv('./Data/OLR3_crn_tau132.csv',header=T,fill=T)
  
  intersect(colnames(crn34.v1),colnames(crn34.v2))
  dim(crn34.v2)
  dim(crn34.v1)
}


otufing12.2=cbind(ch4p=olr12[,c('CH4....')],finger12,otu12.2)
otufing34.2=cbind(olr34[,c('CH4....','HAc', 'HPro', 'i.HBu', 'HBu', 'iHVal', 'HVal', 'HCa')],finger34,otu34.2)

write.csv(otufing12.2,'./Data/OLR1_otufinger_ver2_ch4p_Species_tau132.csv',row.names=F)
write.csv(otufing34.2,'./Data/OLR3_otufinger_ver2_ch4p_Species_tau132.csv',row.names=F)

ada34.ind=cbind(olr34[,c('CH4....','mL.CH4.d', 'HAc', 'HPro', 'i.HBu', 'HBu', 'iHVal', 'HVal', 'HCa')],finger34[,spind.p])
ada12.ind=cbind(olr12[,c( 'CH4....', 'mL.CH4.d')],finger12[,intersect(spind.p,ind.t12)])
otufing12=cbind(ch4p=olr12[,c('CH4....')],finger12[,intersect(spind.p,ind.t12)],otu12)
otufing34=cbind(olr34[,c('CH4....','HAc', 'HPro', 'i.HBu', 'HBu', 'iHVal', 'HVal', 'HCa')],finger34[,spind.p],otu34)

write.csv(otufing12,'./Data/OLR1_otufinger_verc_ch4p_Species_tau132.csv',row.names=F)
write.csv(otufing34,'./Data/OLR3_otufinger_verc_ch4p_Species_tau132.csv',row.names=F)

write.csv(otufing12[,-1],'./Data/OLR1_otucrn_verc_Species_tau132.csv',row.names=F)
write.csv(otufing34[,-1],'./Data/OLR3_otucrn_verc_Species_tau132.csv',row.names=F)

write.csv(otufing12.2[,-1],'./Data/OLR1_otucrn_ver2_Species_tau132.csv',row.names=F)
write.csv(otufing34.2[,-1],'./Data/OLR3_otucrn_ver2_Species_tau132.csv',row.names=F)


############################################################
### DEV analysis
source(paste(Root_Path ,'\\data\\Feitsui\\other.r',sep=''))
source(paste(Root_Path ,Sub_Path,'\\Code\\DEV_function_20230601.R',sep=''))

#############################################################
#### Data Preparation 
# Read dataset
out.sample <- F # T/F for out-of-sample forecast
if(out.sample){nout <- 2}else{nout <- 0}  # number of out-of-sample

#(da.name <- 'OLR_R1_tau1p_finger')
#do <- ada12.ind
(da.name <- 'OLR_R3_crn_verc')
do <- crn34

#(da.name <- 'OLR_R1_otu_verc')
#do <- otu12


ndo <- nrow(do)
nin <- ndo-nout # library sample size

# Excluding rare species with no clear temporal dynamics 
dto <- apply(do,1,sum) # total abundance
dpo <- do*repmat(dto^-1,1,ncol(do))

# In-sample
do.mean <- apply(do[1:nin,],2,mean,na.rm=T)  # mean abundance in in-sample
do.sd <- apply(do[1:nin,],2,sd,na.rm=T)      # SD of abundance in in-sample
ds <- (do-repmat(do.mean,nrow(do),1))*repmat(do.sd,nrow(do),1)^-1 # Normalized in-sample dataset at time t

cons=apply(is.na(ds),2,all)
ds=ds[,!cons]
do=do[,!cons]

zeropad=function(x){
  zint=which(x==0)
  z1=zL=zint[1]
  zc=0
  zpd=1
  zlist=list()
  if(length(zint)>1){
    for(i in 2:length(zint)){
      if((zint[i]-z1)<=1){
        zL=c(zL,zint[i])
      }else{
        if(length(zL)>1){zc=zc+1; zlist[[zc]]=zL}
        zL=zint[i]
      }
      z1=zint[i]
    }
    zlist[[zc+1]]=zL
    
  }else{zlist=NULL}
  return(zlist)
}

(zeroint=apply(do,2,function(x){max(unlist(lapply(zeropad(x),length)))}))
#for(i in 1:ncol(do)) zeropad(do[,i])
##################################################
### DEV

### Note that S-Map calculations for large window sizes and/or time series need a lot of computational time
### start algorithm
DEVresult = para.optR = NULL
for(i in 1:ncol(ds)){
  #for(i in 1:3){
  t.begin <- proc.time()
  ##Performing cross-validation for selecting the optimal parameter set
  #if(i!=3){
  wd.i=15
  if(zeroint[i]>wd.i){wd.i=zeroint[i]+1}
  param.opt.a <- dev_cv2(ds[,i], window_size = seq(wd.i, 30, by=1))
  
  para.optt=param.opt.a[[1]]
  ### Performing DEV algorithm
  devts=DEV_perform(ds[,i],param.opt=para.optt)
  
  DEVresult = rbind(DEVresult, data.frame(variable=colnames(ds)[i], devts))
  para.optR = rbind(para.optR, c(variable=colnames(ds)[i], para.optt))
  rm(list = c('param.opt.a', 'devts'));gc()
  time.used <- proc.time() - t.begin
  cat("effect", i, "in sample completed:", time.used[3],"sec\n")
}

write.csv(DEVresult,paste(da.name,'_DEV_TimeSeries.csv',sep=''),row.names=F)
write.csv(para.optR,paste(da.name,'_DEV_Optparameter.csv',sep=''),row.names=F)

