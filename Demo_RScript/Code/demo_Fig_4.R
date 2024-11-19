# R code for plotting Fig. 4 
rm(list = ls());gc()
# Set path
OPath <- 'D:\\data\\Cooperation\\Prodigio\\'
setwd(OPath)

library(factoextra)
library(igraph)
library(dplyr)
library(mgcv)
library(gratia)
library(EnvStats)

seed <- 49563
set.seed(seed)
SaveFile <- T #T/F for saving files

tipp=37
tipp2=43
nd=46
tind.o=1:46


#####################################################################################################################################
# Data loading
clasdo <- read.csv('./Demo_RScript/Data/demo_otu_taxa.csv',header=T)
rownames(clasdo)=clasdo[,'OTU']
crnotu1 <- read.csv('./Demo_RScript/Data/Demo_OLR1_otucrn_verc_Species_tau132.csv',header=T,stringsAsFactors = F)
crnotu3 <- read.csv('./Demo_RScript/Data/Demo_OLR3_otucrn_verc_Species_tau132.csv',header=T,stringsAsFactors = F)
reactor_ch4_ews <- read.csv( './Demo_RScript/Data/reactor_ch4_ews.csv', header=T)

# Calculate percentage data for microbial OTU 
otu1 <- crnotu1[,grep('OTU_',colnames(crnotu1))]
otu3 <- crnotu3[,grep('OTU_',colnames(crnotu3))]

# Extract metabolomic data
crn1 <- crnotu1[,-grep('OTU_',colnames(crnotu1))]
crn3 <- crnotu3[,-grep('OTU_',colnames(crnotu3))]
rownames(otu3)=rownames(crn3)=1:nd
crnotu3s=apply(crnotu3,2,scale)

nd=NROW(otu3)
nsp=NCOL(otu3)
ncrn=NCOL(crn1)

nco1=ncrn+nsp
nco3=ncrn+nsp+7

ind.co=data.frame(Type=c(rep('VFA',7),rep('Fingerprint',ncrn), clasdo[,'Functional_guild']),
                    Name=colnames(crnotu3),Abb=c(rep(1,7),rep(2,ncrn),
                                                 match(clasdo[,'Functional_guild'],c("Methanogen","Acetogen","Other_taxa"))+2))
rownames(ind.co)=ind.co[,'Name']

newd=seq(0,tipp,0.01)
################################################################################
## Function for detecting significant sign changes in the trend of a time series
signch2=function(oy,y,x,yci=NULL){
  n=length(x)
  ay=abs(y)
  may=mean(ay,na.rm=T)
  ind.t=which(apply(embed(y,2),1,prod)<0)+1
  ys=sign((sign(yci[,'lower']*yci[,'upper'])==1)*yci[,'lower'])
  ind.ts0=which(sign(y)[-1]-sign(y)[-length(y)]!=0)+1
  ind.ts=which(ys[-1]-ys[-length(ys)]!=0)
  if(length(ind.ts)==0){
    yss=data.frame(ts_i=1,ts_f=max(x),sign=0,extreme=0,significance=0,psign_ch=0,sign_tp=0)
    if(length(ind.ts0)!=0){yss['psign_ch']=1}}else{
      yss=cbind(c(1,ind.ts+1),c(ind.ts, length(ys)),c(ys[ind.ts],ys[1+ind.ts[length(ind.ts)]]),rep(0),rep(0),rep(1),rep(0))
      for(i in 1:NROW(yss)) {
        sign.tp=unique(sign(y[yss[i,1]:yss[i,2]]))
        if(length(sign.tp)==1){yss[i,7]=sign.tp}
        if(i>1 & i<NROW(yss)){
          if(any(c(yss[i,1]:yss[i,2])%in%ind.t)){
            yss[i,4]=max(ind.t[ind.t%in%c(yss[i,1]:yss[i,2])])
            yss[i,5]=sign(yss[i+1,3])*((yss[i-1,3]*yss[i+1,3])==-1)
          }
        }
      }
      
      yss[,c(1:2)]=cbind(x[yss[,1]],x[yss[,2]])
      yss[yss[,4]!=0,4]=x[yss[yss[,4]!=0,4]]
      
      colnames(yss)=c('ts_i','ts_f','sign','extreme','significance','psign_ch','sign_tp')
    }
  return(yss)
}

periodf=function(x){
  ust=unique(x)
  period.a=NULL
  for(ui in 1:length(ust)){
    y=which(x==ust[ui])
    if(length(y)!=1){
      emb2=embed(y,2)
      period.ind=emb2[,1]-emb2[,2]
      indi=which(period.ind!=1)
      period.i=matrix(c(emb2[1,2], max(emb2)),nrow=1,ncol=2)
      if(length(indi)>0){
        for(ki in 1:length(indi)){
          period.i[ki,2]=emb2[indi[ki],2]
          period.i=rbind(period.i,c(emb2[indi[ki],1], NA))
        }
        period.i[NROW(period.i),2]=max(emb2)
      }
      colnames(period.i)=c('initial','final')
      period.a=rbind(period.a, data.frame(state=rep(ust[ui]),period.i)) 
    }
  }
  return(period.a)
}

find.bif=function(da.t, sig=1, above.pre=T){
  ne2=sum(!is.na(da.t[,'value']))
  g.t <- gam(value~s(time, k = trunc(ne2/2)),data=da.t)
  xmin=da.t[which(!is.na(da.t[,'value']))[1],'time']
  newd.i=seq(xmin,tipp2,0.01)
  
  fd <- fderiv(g.t, newdata = data.frame(time=newd.i))
  # Calculate 95% confidence intervals for derivatives
  cid <- as.matrix(confint(fd, type = "confidence")[,c('lower','est','upper')])
  ci <- confint(g.t,parm='s(time)',data=data.frame(time=newd.i))[,c('lower','est','upper')]+g.t$coefficients["(Intercept)"]
  yhat=predict.gam(g.t,newdata=data.frame(time=newd.i))
  yhatd1=cid[,'est']
  
  crxp=as.numeric(which((sign(yhat-1)[-length(yhat)]-sign((yhat-1))[-1])<0))
  exmp=which((sign(yhatd1)[-length(yhatd1)]-sign((yhatd1))[-1])>0)
  warnn.u=filter(periodf((ci[,'lower']>1)*1),state==1)
  warnn.l=filter(periodf((ci[,'upper']<1)*1),state==1)
  
  if(length(crxp)!=0){
    ind.t=NULL
    for(cxi in 1:length(crxp)){
      ind1=any(warnn.u[,'initial']>crxp[cxi])
      ind2=any(warnn.l[,'final']<crxp[cxi])
      ind3=!any(yhat[newd.i<newd.i[crxp[cxi]]]>1)
      if(sig==1){ind.t=c(ind.t,ind1)}else if(sig==2){ind.t=c(ind.t,ind1&ind2)}
      if(above.pre){ind.t[cxi]=ind.t[cxi]&ind3}
    }  
    crxp=as.matrix(crxp)[ind.t,]#[ind1,]#
  }
  
  if(length(crxp)!=0){
    warp=warnn.u[,'initial']
    warp.i=NULL
    for(cx.i in 1:length(crxp)){
      warp.t=sort(warp[warp>crxp[cx.i]])[1]
      if(length(warp.t)==0){warp.t=NA}
      if(cx.i>1){if(!(is.na(warp.i[cx.i-1])|is.na(warp.t))){if(warp.t==warp.i[cx.i-1]){warp.i[cx.i-1]=NA}}}
      warp.i=c(warp.i,warp.t)
    }
    
    exmp.i=NULL
    for(wx.i in 1:length(warp.i)){
      if(!is.na(warp.i[wx.i])){
        exmp.t=sort(exmp[exmp>warp.i[wx.i]])[1]
        if(length(exmp.t)==0){exmp.t=max(yhat[(1:length(yhat))>warp.i[wx.i]])}
        exmp.i=c(exmp.i,exmp.t)
      }else{exmp.i=c(exmp.i,NA)}
    }
    
    crip.i=data.frame(BFcross=newd.i[as.numeric(crxp)], BFwarn=newd.i[as.numeric(warp.i)],BFextrem=newd.i[as.numeric(exmp.i)])
    crip.i=crip.i[(!is.na(crip.i[,'BFwarn'])),]
    if(NROW(crip.i)>0) {crip.i=data.frame(order=1:NROW(crip.i),crip.i)}
    
  }else{crip.i=NULL}
  return(crip.i)
}

mag_scale=function(x,mgf=0.5){
  if(prod(x)>=0){
    if(x[2]>0){y=x;y[1]=y[1]*(1.0-mgf);y[2]*(1.0+mgf)}else{
      y=x;y[1]=y[1]*(1.0+mgf);y[2]*(1.0-mgf)}
  }else{y=x*(1.0+mgf)}
  return(y)
}

warn_turn <- function(x,y,tmp){
  da.t1=data.frame(time=x, value=y)
  g.t1 <- gam(value~s(time),data=da.t1)
  yhat=predict.gam(g.t1,newdata=data.frame(time=newd))
  
  turn.i=rep(0,7)
  names(turn.i)=c('sign_f','turn','warn_time','change','change_time','extinction','psign_ch')#,'extinction_time')
  turn.ind=which(tmp[,'significance']!=0)
  turn.i['sign_f']=as.numeric(tmp[NROW(tmp), c('sign')])
  sign_tp.i=unique(tmp[NROW(tmp):(NROW(tmp)-1),'sign_tp'])
  if(length(sign_tp.i)==1){turn.i[c(1)]=sign_tp.i}
  turn.i['psign_ch']=any(tmp[, c('psign_ch')]==1)*1
  if(length(turn.ind)>0){
    turn.i[c('turn','warn_time')]=as.numeric(tmp[max(turn.ind), c('significance','extreme')])
  }else if(any(tmp[,'sign']!=0)){
    turn.i[c('change','change_time')]=as.numeric(tmp[min(which(tmp[,'sign']!=0)),c('sign','ts_i')])
  }
  if(da.t1[nrow(da.t1),'value']==0){
    turn.i['extinction']=1;ett=da.t1[which(da.t1[,'value']==0),'time'];
    ett=ett[ett>tmp[NROW(tmp),'ts_i']]
  }

  k.t <- kendallTrendTest(da.t1[,'value'])
  kti=c(k.t$estimate[c('tau')], ifelse(k.t$p.value<0.05,1,0));
  names(kti)=c('tau', 'sig_tau')
  
  vtmp=c(turn.i,kti)
  tp.i=NULL
  if(vtmp['turn']>0){tp.i='B1'}else if(vtmp['turn']<0){tp.i='B2'}else if(vtmp['tau']>0& vtmp['sig_tau']==1 & vtmp['change_time']==0){
    tp.i='A1';vtmp['change_time']=1}else if(vtmp['tau']<0& vtmp['sig_tau']==1& vtmp['change_time']==0){tp.i='A2';vtmp['change_time']=1
    }else if(vtmp['change']>0){if(vtmp['sign_f']==1 & vtmp['change_time']>1){tp.i='C1'}else{tp.i='D1'}
    }else if(vtmp['change']<0){if(vtmp['sign_f']==-1 & vtmp['change_time']>1){tp.i='C2'}else{tp.i='D2'}}else{if(vtmp['psign_ch']==0){tp.i='E1'}else{tp.i='E2'}}
  vtmp['warn_time'] <- vtmp['warn_time']+vtmp['change_time'] 
  
  return(list(tp.i, vtmp))
}

##################################################################################
## Time pattern analysis clusters
mem.all <- read.csv('demo_cluster_all_verc.csv',header=T)
cluID=apply(mem.all[,c('type','cluster')],1,paste,collapse = '-')
da.tmp=NULL
for(i in 1:NROW(mem.all)){
  da.t=data.frame(time=1:tipp2, value=as.numeric(mem.all[i,-c(1:2)]))
  da.t1=filter(da.t,time<=tipp);da.t1[,'value']=scale(da.t1[,'value'])
  g.t1 <- gam(value~s(time),data=da.t1)
  
  yhat=predict.gam(g.t1,newdata=data.frame(time=newd))
  yhatd1 <- derivatives(g.t1, type = "central",data=data.frame(time=newd), order=1)$derivative
  fd <- fderiv(g.t1, newdata = data.frame(time=newd))
  # Calculate 95% confidence intervals for derivatives
  ci <- as.matrix(confint(fd, type = "confidence")[,c('lower','est','upper')])
  pk.i=signch2(oy=yhat,y=yhatd1,yci=ci,x=unlist(c(newd)))
  da.tmp=rbind(da.tmp, data.frame(variable=rep(i,NROW(pk.i)),ID=cluID[i],pk.i))
}

clu.tmp=tpp=NULL
for(i in 1:NROW(mem.all)){
  turn.i <- warn_turn(x=1:tipp, y=as.numeric(mem.all[i,-c(1:2)][1:tipp]), tmp=filter(da.tmp,variable==i))
  tpp=c(tpp,turn.i[[1]])
  clu.tmp=rbind(clu.tmp, turn.i[[2]])
}
clu.tmp=data.frame(ID=cluID,Type=mem.all[,1],Ttype= tpp, clu.tmp)
rownames(clu.tmp)= cluID

ind.w=which(clu.tmp[,'Ttype']%in%c('B1','B2','C1','C2'))


myColor <- c(colorRampPalette(c("lightblue","blue"))(sum(mem.all[,'type']=='OTU')),
  colorRampPalette(c("pink",'red' ))(sum(mem.all[,'type']=='CRN')),
  colorRampPalette(c("lightyellow2","orange" ))(sum(mem.all[,'type']=='NET')))
windows(100,70)
plot(as.numeric(mem.all[1,-c(1:2)][1:tipp2])~c(1:tipp2),type='n',ylim=c(-2.05,2.05),xlim=c(1,tipp+1),ylab='Smoothed cluster average',xlab='Time')
for(i in 1:NROW(mem.all)){
  mem.i <- mem.all[i,]
  da.t=data.frame(time=1:tipp, value=scale(as.numeric(mem.i[-c(1:2)][1:tipp])))
  g.t1 <- gam(value~s(time),data=da.t)
  yhat <- predict.gam(g.t1,newdata=data.frame(time=newd))
  lines(yhat~newd,col=myColor[i],lwd=2) 
}

polygon(x=c(tipp,tipp2,tipp2,tipp),y=rep(mag_scale(c(-1.5,1.5)),each=2), col="orange",density=15, border="white")
legend('bottomright',c('Microbial OTUs','Metabolites','Network properties'),col=c('blue','red','orange'),lty=1,lwd=2,bty='n')

pch.t=clu.tmp[ind.w,'sign_f']
points((rep(2.05, length(ind.w))*-pch.t)~clu.tmp[ind.w,'warn_time'], 
       pch=ifelse(pch.t>0,24,25), bg= myColor[ind.w], col=myColor[ind.w],cex=1.5)
text(clu.tmp[ind.w,'warn_time'],rep(1.8, length(ind.w))*-pch.t, 
     clu.tmp[ind.w,'ID'], col=myColor[ind.w],cex=0.8)


################################################################################
### CCM results
ccm.otu=read.csv('./Demo_RScript/Data/demo_ccm_OTU.csv', header=T)
ccm.crn=read.csv('./Demo_RScript/Data/demo_ccm_metabolite.csv', header=T)
ccm.net=read.csv('./Demo_RScript/Output/demo_ccm_Network.csv', header=T)

################################################################################
## Computed the smoothed estimate of variable trends in microbiome & metabolome
newd=seq(0,tipp,0.01)
nnd=length(newd)
if(F){
  smo.crnotu3=matrix(NA, nco3*nnd ,5)
  for(i in 1:nco3){
    cat("\r",round(i/nco3*100,1),'% of the network nodes')
    da.t=data.frame(time=tind.o, value=crnotu3[,i])
    da.t1=filter(da.t,time<=tipp)
    g.t1 <- gam(value~s(time),data=da.t1)
    yhat=predict.gam(g.t1,newdata=data.frame(time=newd))
    fd <- fderiv(g.t1, newdata = data.frame(time=newd))
    # Calculate 95% confidence intervals for derivatives
    ci <- as.matrix(confint(fd, type = "confidence")[,c('lower','est','upper')])
    yhatd1 <- ci[,'est']
    smo.crnotu3[(nnd*(i-1)+1):(nnd*i),]=cbind(time=newd,yhat, ci)
  }
  smo.crnotu3=data.frame(ID=rep(colnames(crnotu3), each=nnd),smo.crnotu3)
  colnames(smo.crnotu3) <- c('ID','time','yhat','lower','est','upper')
  write.csv(smo.crnotu3, './Demo_RScript/Output/smo_otucrn3_37.csv',row.names = F)
}

smo.crnotu3=read.csv('./Demo_RScript/Output/smo_otucrn3_37.csv', header = T)

################################################################################
## Bifurcation analysis for microbial OTUs
dev_crnotu <- read.csv('./Demo_RScript/Data/DEV_TS_crnotu.csv',header = T)
devf3=filter(dev_crnotu,variable%in%colnames(otu3))

udev3=unique(devf3[,'variable'])
dev.mx=matrix(NA,nd,length(udev3))
for(i in 1:length(udev3)){
  da.t=filter(devf3, variable==udev3[i])
  dev.mx[(1:nd)%in%da.t[,'t_max'],i]=da.t[,'Abs']
}
colnames(dev.mx)=udev3

crip=NULL
j=1
for(i in 1:NCOL(dev.mx)){
  crip.i=NULL
  da.t2 = data.frame(time=1:tipp2,value=dev.mx[1:tipp2,i])
  crip.i <- find.bif(da.t = da.t2, sig=1)
  if(NROW(crip.i)>0) crip=rbind(crip, data.frame(variable=rep(i),crip.i))
}  

crip.otu=data.frame(ID=colnames(dev.mx)[crip[,1]],crip)

################################################################################
### Bifurcation analysis for metabolites
devf3=filter(dev_crnotu,variable%in%colnames(crn3))

udev3=unique(devf3[,'variable'])
dev.mx=matrix(NA,nd,length(udev3))
for(i in 1:length(udev3)){
  da.t=filter(devf3, variable==udev3[i])
  dev.mx[(1:nd)%in%da.t[,'t_max'],i]=da.t[,'Abs']
}
colnames(dev.mx)=udev3

crip=NULL
j=1
for(i in 1:NCOL(dev.mx)){
  crip.i=NULL
  da.t2 = data.frame(time=1:tipp2,value=dev.mx[1:tipp2,i])
  crip.i <- find.bif(da.t = da.t2, sig=1)
  if(NROW(crip.i)>0) crip=rbind(crip, data.frame(variable=rep(i),crip.i))
}  

crip.crn=data.frame(ID=colnames(dev.mx)[crip[,1]],crip)

################################################################################
### Bifurcation analysis for network properties
devf3=read.csv('./Demo_RScript/Output/DEV_demo_Network_Perturbed.csv',header=T)
udev3=unique(devf3[,'variable'])

dev.mx=matrix(NA,nd,length(udev3))
for(i in 1:length(udev3)){
  da.t=filter(devf3, variable==udev3[i])
  dev.mx[(1:nd)%in%da.t[,'t_max'],i]=da.t[,'Abs']
}
colnames(dev.mx)=udev3

crip=NULL
j=1
for(i in 1:NCOL(dev.mx)){
  crip.i=NULL
  da.t2 = data.frame(time=1:tipp2,value=dev.mx[1:tipp2,i])
  crip.i <- find.bif(da.t = da.t2, sig=1)
  if(NROW(crip.i)>0) crip=rbind(crip, data.frame(variable=rep(i),crip.i))
}  
crip.net=data.frame(ID=udev3[crip[,'variable']],crip)

################################################################################
# Temporal pattern analysis for microbiome & metabolome
################################################################################
co.tmp=tpp=NULL
for(i in 1:nco3){
  cat("\r",round(i/nco3*100,1),'% of the network nodes')
  smo.i=filter(smo.crnotu3,ID==colnames(crnotu3)[i])
  yhat=smo.i[,'yhat']
  yhatd1=smo.i[,'est']
  pk.i=signch2(oy=yhat,y=yhatd1,yci=smo.i[,c('lower','est','upper')],x=unlist(c(newd)))
  da.tmp <- data.frame(variable=rep(i,NROW(pk.i)),ID=colnames(crnotu3)[i],pk.i)
  turn.i <- warn_turn(x=1:tipp, y=crnotu3[1:tipp,i], tmp=da.tmp)
  tpp=c(tpp,turn.i[[1]])
  co.tmp=rbind(co.tmp, turn.i[[2]])
}
co.tmp=data.frame(ID=colnames(crnotu3),Type=ind.co[,'Type'],Ttype= tpp, co.tmp)
rownames(co.tmp)= colnames(crnotu3)

write.csv(co.tmp, './Demo_RScript/Output/Localized_Transition_otucrn3.csv',row.names = F)

###############################################################################
# Temporal pattern analysis for network properties
################################################################################
netprop.o=read.csv('./Demo_RScript/Data/demo_network_prop_perturbed.csv',header=T)
netprop.p=filter(netprop.o, Reactor=='Perturbed')
netprop=t(netprop.p[,-c(1:4)]);colnames(netprop)=netprop.p[,'Property'];rownames(netprop)=NULL

net.tmp=tpp=NULL
for(i in 1:NCOL(netprop)){
  cat("\r",round(i/NCOL(netprop)*100,1),'% of the network nodes')
  da.t1=data.frame(time=1:tipp, value=netprop[1:tipp,i])
  g.t1 <- gam(value~s(time),data=da.t1)
  yhat=predict.gam(g.t1,newdata=data.frame(time=newd))
  yhatd1 <- derivatives(g.t1, type = "central",data=data.frame(time=newd), order=1)$derivative
  fd <- fderiv(g.t1, newdata = data.frame(time=newd))
  # Calculate 95% confidence intervals for derivatives
  ci <- as.matrix(confint(fd, type = "confidence")[,c('lower','est','upper')])
  pk.i=signch2(oy=yhat,y=yhatd1,yci=ci,x=unlist(c(newd)))
  da.tmp <- data.frame(variable=rep(i,NROW(pk.i)),ID=netprop.p[,'Property'][i],pk.i)
  turn.i <- warn_turn(x=1:tipp, y=da.t1[,'value'], tmp=da.tmp)
  tpp=c(tpp, turn.i[[1]])
  net.tmp=rbind(net.tmp, turn.i[[2]])
}

net.tmp=data.frame(ID=netprop.p[,'Property'],Type=netprop.p[,'Property_type'],Ttype= tpp, net.tmp)
rownames(net.tmp)= netprop.p[,'Property']

############################################################################################
# Detect localized transitions in microbiome, metabolome, and dynamical interactome
# microbial OTUs
crip.otu2 <- filter(crip.otu,order==1, BFwarn <=tipp)
otu.mvs <- co.tmp[grep('OTU_',co.tmp[,'ID']),]
fin.mvs <- data.frame(otu.mvs,BFcross=NA,BFwarn=NA,BFextrem=NA, CCM=ccm.otu[,'Perturbed_OTU'])
fin.mvs[fin.mvs[,'ID']%in%crip.otu2[,'ID'],c('BFcross', 'BFwarn', 'BFextrem')]=crip.otu2[,c('BFcross', 'BFwarn', 'BFextrem')]
fin.mvs[,'Type'] <- factor(clasdo[,'Functional_guild'], level=c('Methanogen','Acetogen', 'Other_taxa'))
fin.mvs <- fin.mvs%>%arrange(Type, desc(warn_time))%>%filter(Ttype%in%c('B1','B2','C1','C2'))
otu.g <- ind.co[fin.mvs[,'ID'],]

win.graph(50,90)
par(mar=c(4,6,1,1))
ytick <- nrow(fin.mvs):1
fg.pos1 <- c(min(ytick[fin.mvs[,'Type']=='Methanogen']),min(ytick[fin.mvs[,'Type']=='Acetogen']))
fg.pos2 <- c(max(ytick[fin.mvs[,'Type']=='Methanogen']),max(ytick[fin.mvs[,'Type']=='Acetogen']),max(ytick[fin.mvs[,'Type']=='Other_taxa']))
plot(y=ytick,x=fin.mvs[,'warn_time'],xlim=c(1,tipp),
     col=ifelse(fin.mvs[,'Ttype']%in%c('B1','C1'),4,2), ,pch=16,cex= (fin.mvs[,'CCM'])*2,ylab='',xlab='Time', yaxt="n")
abline(h=fg.pos1-0.8,col='grey')
text(rep(8,3), fg.pos2-c(-1.5,1.5,1.5),c('Methanogen','Acetogen', 'Other_taxa'))

points(y=ytick,x=fin.mvs[,'BFwarn'],col='blue',pch=17,cex= 0.75)
axis(side=2, at=nrow(fin.mvs):1, labels = FALSE)
text(par("usr")[1], ytick,  cex=0.4,labels = clasdo[fin.mvs[,'ID'], 'finest_taxa_name'], srt = 0, pos = 2, xpd = TRUE)

warn.otu <- fin.mvs %>% filter(!is.na(BFwarn)) %>% dplyr::select(ID, Type, Ttype, warn_time, BFwarn, CCM)

# metabolites
crip.crn2=filter(crip.crn,order==1, BFwarn <=tipp)
crn.mvs=co.tmp[-grep('OTU_',co.tmp[,'ID']),]
fin.mvs=data.frame(crn.mvs,BFcross=NA,BFwarn=NA,BFextrem=NA, CCM=ccm.crn[,'Perturbed_finger'])
fin.mvs[fin.mvs[,'ID']%in%crip.crn2[,'ID'],c('BFcross', 'BFwarn', 'BFextrem')]=crip.crn2[,c('BFcross', 'BFwarn', 'BFextrem')]
fin.mvs=fin.mvs%>%arrange(desc(Type), desc(warn_time))%>%filter(Ttype%in%c('B1','B2','C1','C2'))#,'D1','D2'))
crn.g=ind.co[fin.mvs[,'ID'],]
warn.crn=fin.mvs[!is.na(fin.mvs[,'BFwarn']),]

win.graph(50,90)
par(mar=c(4,6,1,1))
ytick=nrow(fin.mvs):1
fg.pos=min(ytick[fin.mvs[,'Type']=='VFA'])
plot(y=ytick,x=fin.mvs[,'warn_time'],xlim=c(1,tipp),
     col=ifelse(fin.mvs[,'Ttype']%in%c('B1','C1'),4,2),pch=16,cex= (fin.mvs[,'CCM'])*2,ylab='',xlab='Time', yaxt="n")
abline(h=fg.pos-0.8,col='grey')
text(rep(7,2), c(max(ytick)+2,max(ytick)-8),c('VFA','Fingerprint'))

points(y=ytick,x=fin.mvs[,'BFwarn'],col='blue',pch=17,cex= 0.75)
axis(side=2, at=nrow(fin.mvs):1, labels = FALSE)
text(par("usr")[1], ytick,  cex=0.6, 
     labels = fin.mvs[,'ID'], srt = 0, pos = 2, xpd = TRUE)

warn.crn <- fin.mvs %>% filter(!is.na(BFwarn)) %>% dplyr::select(ID, Type, Ttype, warn_time, BFwarn, CCM)

write.csv(rbind(warn.crn,warn.otu), './Demo_RScript/Output/warn_time_CRNOTU.csv',row.names = F)


# Network proprties
crip.net2=filter(crip.net,order==1, BFwarn <=tipp)
fin.mvs=data.frame(net.tmp,BFcross=NA,BFwarn=NA,BFextrem=NA,CCM=as.numeric(ccm.net))
fin.mvs[fin.mvs[,'ID']%in%crip.net2[,'ID'],c('BFcross', 'BFwarn', 'BFextrem')]=crip.net2[,c('BFcross', 'BFwarn', 'BFextrem')]
fin.mvs[,'Type']=factor(fin.mvs[,'Type'], level=c('Centrality','FG_Net', 'Feedback', 'IS_type','Topology'))
fin.mvs=fin.mvs%>%arrange(Type, desc(warn_time))%>%filter(Ttype%in%c('B1','B2','C1','C2'))

win.graph(50,90)
par(mar=c(4,8.5,1,1))
ytick=nrow(fin.mvs):1
fg.pos1=c(min(ytick[fin.mvs[,'Type']=="Centrality"]),min(ytick[fin.mvs[,'Type']=="FG_Net"]),min(ytick[fin.mvs[,'Type']=="IS_type"]))
fg.pos2=c(max(ytick[fin.mvs[,'Type']=="Centrality"]),max(ytick[fin.mvs[,'Type']=="FG_Net"]),max(ytick[fin.mvs[,'Type']=="IS_type"]),max(ytick[fin.mvs[,'Type']=="Topology"]))
plot(y=ytick,x=fin.mvs[,'warn_time'],xlim=c(1,tipp),
     col=ifelse(fin.mvs[,'Ttype']%in%c('B1','C1'),4,2), pch=16, cex= fin.mvs[,'CCM']*2, ylab='',xlab='Time', yaxt="n")
abline(h=fg.pos1-0.8,col='grey')
text(rep(8,4), fg.pos2-c(0,0.8,0.8,0.8),c('Centrality','FG_Net', 'IS_type','Topology'))
legend('bottomright',c('Turning point', 'Bifurcation point'),bty='n',col=c('black','blue'),pch=c(1,17))

points(y=ytick,x=fin.mvs[,'BFwarn'],col='blue',pch=17,cex= 1)

axis(side=2, at=nrow(fin.mvs):1, labels = FALSE)
text(par("usr")[1], ytick,  cex=0.75, labels = fin.mvs[,'ID'], srt = 0, pos = 2, xpd = TRUE)

