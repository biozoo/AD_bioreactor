# R code for plotting Fig. 4 
rm(list = ls());gc()
# Set path
OPath <- 'D:/data/Cooperation/Prodigio/Code/Demo_RScript_20250529'
setwd(OPath)

library(factoextra)
library(igraph)
library(dplyr)
library(mgcv)
library(MASS)
library(gratia)
library(EnvStats)

seed <- 49563
set.seed(seed)
SaveFile <- F #T/F for saving files

tipp <- 38+1
tipp2 <- 44+1
nd <- 47
tind.o <- 1:47

newd <- seq(1,tipp,0.01)
nnd <- length(newd)
source('./Code/Function_DEV_V2.R')

#####################################################################################################################################
# Data loading
clasdo <- read.csv('./Data/demo_otu_taxa_D47.csv',header=T)
rownames(clasdo)=clasdo[,'OTU']
crnotu <- read.csv('./Data/OLR_Perturbed_crnotu_verc_CH4_Species_D47.csv',header=T,stringsAsFactors = F)

# Calculate percentage data for microbial OTU 
otu3 <- crnotu[,grep('OTU_',colnames(crnotu))]

# Extract metabolomic data
crn3 <- crnotu[,-grep('OTU_',colnames(crnotu))]
rownames(otu3) <- rownames(crn3) <- 1:nd

nd <- NROW(otu3)
nsp <- NCOL(otu3)
ncrn <- NCOL(crn3)

#nco1=ncrn+nsp
nco3 <- ncrn+nsp

ind.co <- data.frame(Type=c(rep('VFA',7),rep('Fingerprint',ncrn-7), clasdo[,'Functional_guild']),
                    Name=colnames(crnotu),Abb=c(rep(1,7),rep(2,ncrn-7),
                                                 match(clasdo[,'Functional_guild'],c("Methanogen","Acetogen","Other_taxa"))+2))
rownames(ind.co) <- ind.co[,'Name']



################################################################################
## Function for detecting significant sign changes in the trend of a time series
signch2 <- function(x,yci=NULL){
  y <- yci[,'est']
  n <- length(x)
  ind.t <- which(apply(embed(y,2),1,prod)<0)+1
  ys <- sign((sign(yci[,'lower']*yci[,'upper'])==1)*yci[,'lower'])
  ind.ts0 <- which(sign(y)[-1]-sign(y)[-length(y)]!=0)+1
  ind.ts <- which(ys[-1]-ys[-length(ys)]!=0)
  if(length(ind.ts)==0){
    yss <- data.frame(ts_i=1,ts_f=max(x),sign=0,extreme=0,extreme_sign=0,sign_change=0)
    if(length(ind.ts0)!=0){yss['sign_change'] <- 1}
    }else{
      yss <- cbind(ts_i=c(1,ind.ts+1),
                   ts_f=c(ind.ts,length(ys)),
                   sign=c(ys[ind.ts],ys[1+ind.ts[length(ind.ts)]]),
                   extreme=rep(0),extreme_sign=rep(0),sign_change=rep(1))
      for(i in 1:NROW(yss)) {
        if(i>1 & i<NROW(yss)){
          if(any(c(yss[i,'ts_i']:yss[i,'ts_f'])%in%ind.t)){
            if((yss[i-1,'sign']*yss[i+1,'sign'])<0){
              yss[i,'extreme'] <- x[max(ind.t[ind.t%in%c(yss[i,'ts_i']:yss[i,'ts_f'])])]
              yss[i,'extreme_sign'] <- sign(yss[i+1,'sign'])*((yss[i-1,'sign']*yss[i+1,'sign'])==-1)
            }
          }
        }
      }
      yss[,c('ts_i','ts_f')] <- cbind(x[yss[,'ts_i']], x[yss[,'ts_f']])
    }
  return(yss)
}

mag_scale <- function(x,mgf=0.5){
  if(prod(x)>=0){
    if(x[2]>0){y <- x;y[1] <- y[1]*(1.0-mgf);y[2]*(1.0+mgf)}else{
      y <- x;y[1] <- y[1]*(1.0+mgf);y[2]*(1.0-mgf)}
  }else{y <- x*(1.0+mgf)}
  return(y)
}

warn_turn2=function(y, yci, pk){
  k.t <- kendallTrendTest(y)
  kti <- c(k.t$estimate[c('tau')], ifelse(k.t$p.value<0.05,1,0));
  names(kti) <- c('tau', 'sig_tau')
  
  turn.i <- c(turn=0,warn_time=0)
  tp.i <- NULL
  np <- NROW(pk)
  sign.f <- pk[np,'sign']
  
  if(np==1){
    turn.i['turn'] <- 0
    if(kti['tau']>0){tp.i <- 'A1'}else{tp.i <- 'A2'}
  }else{
    if(T){
      if(sign.f==0){
        ind.f <- which((yci[,'x']>=pk[np,'ts_i'])&(yci[,'x']<=pk[np,'ts_f']))
        if(all(sign(yci[ind.f,'est'])==pk[np-1,'sign'])){sign.f <- pk[np,'sign'] <- pk[np-1,'sign']}
      }}
    if(sign.f==0){
      if(pk[np-1, 'sign']==1){tp.i <- 'C1'}else{tp.i <- 'C2'}
    }else{
      ind.op <- which(pk[, 'sign']==-sign.f)
      ind.sp <- which(pk[, 'sign']==sign.f)
      if(length(ind.op)!=0){ # if any sign change compared to the trend of the last period
        turn.i['warn_time'] <- max(pk[(1:np)>ind.op, 'extreme']) # find the last extreme value
        ifelse(sign.f==1, tp.i <- 'D1', tp.i <- 'D2')
      }else{ # monotonic trends
        if(pk[1, 'sign']==0){ifelse(sign.f==1, tp.i <- 'B1',tp.i <- 'B2')}
        if(pk[1, 'sign']!=0){ifelse(sign.f==1, tp.i <- 'A1',tp.i <- 'A2')}
        if(tp.i %in% c('B1','B2')){turn.i['warn_time'] <- pk[min(which(pk[, 'sign']==sign.f)),'ts_i']}
      }
    }
  }  
  if(tp.i == 'A1' & kti['sig_tau']==0){tp.i <- 'E1'}
  if(tp.i == 'A2' & kti['sig_tau']==0){tp.i <- 'E2'}
  if(tp.i%in%c('B1','D1')){turn.i['turn'] <- 1}
  if(tp.i%in%c('B2','D2')){turn.i['turn'] <- -1}
  return(list(tp.i, c(turn.i, kti)))
}

##################################################################################
## Time pattern analysis clusters
mem.all <- read.csv('./Output/demo_cluster_all_verc_D47.csv',header=T)
cluID <- apply(mem.all[,c('type','cluster')],1,paste,collapse = '-')
cluID <- gsub('NET-', 'N',gsub('CRN-', 'C',gsub('OTU-', 'M',cluID)))

clu.tmp <- tpp <- NULL
for(i in 1:NROW(mem.all)){
  da.t <- data.frame(time=1:tipp2, value=as.numeric(mem.all[i,-c(1:2)]))
  da.t1 <- filter(da.t,time<=tipp);da.t1[,'value']=scale(da.t1[,'value'])
  g.t1 <- gam(value~s(time),data=da.t1)
  
  yhat <- predict.gam(g.t1,newdata=data.frame(time=newd))
  yhatd1 <- derivatives(g.t1, type = "central",data=data.frame(time=newd), order=1)$derivative
  fd <- fderiv(g.t1, newdata = data.frame(time=newd))
  # Calculate 95% confidence intervals for derivatives
  ci <- as.matrix(confint(fd, type = "confidence")[,c('lower','est','upper')])
  pk.i <- signch2(yci=ci, x=unlist(c(newd)))
  turn.i <- warn_turn2(y=as.numeric(mem.all[i,-c(1:2)][1:tipp]), yci=data.frame(x=newd,ci), pk=pk.i)
  tpp <- c(tpp,turn.i[[1]])
  clu.tmp <- rbind(clu.tmp, turn.i[[2]])
}

clu.tmp <- data.frame(ID=cluID,Type=mem.all[,1],Ttype= tpp, clu.tmp)
rownames(clu.tmp) <- cluID


myColor <- c(colorRampPalette(c("lightblue","blue"))(sum(mem.all[,'type']=='OTU')),
  colorRampPalette(c("pink",'red' ))(sum(mem.all[,'type']=='CRN')),
  colorRampPalette(c("lightyellow2","orange" ))(sum(mem.all[,'type']=='NET')))

clu.tmp <- data.frame(clu.tmp, myColor)
clu.turn <- filter(clu.tmp,Ttype%in%c('B1','B2','D1','D2'))

windows(100,70)
plot(as.numeric(mem.all[1,-c(1:2)][1:tipp2])~c((1:tipp2)-1),
     type='n',ylim=c(-2.05,2.05),xlim=c(1,tipp2)-1,ylab='Smoothed cluster average',xlab='Time')
for(i in 1:NROW(mem.all)){
  mem.i <- mem.all[i,]
  da.t <- data.frame(time=(1:tipp), value=scale(as.numeric(mem.i[-c(1:2)][1:tipp])))
  g.t1 <- gam(value~s(time),data=da.t)
  yhat <- predict.gam(g.t1,newdata=data.frame(time=newd))
  lines(yhat~I(newd-1), col=myColor[i], lwd=2) 
}

polygon(x=c(tipp,tipp2,tipp2,tipp)-1,y=rep(mag_scale(c(-1.5,1.5)),each=2), col="orange",density=15, border="white")
legend('bottomright',c('Microbial OTUs','Metabolites','Network properties'),col=c('blue','red','orange'),lty=1,lwd=2,bty='n')
abline(v=tipp-1,col='orange')

pch.t <- clu.turn[,'turn']
points((rep(2.05, NROW(clu.turn))*-pch.t)~c(clu.turn[,'warn_time']-1), 
       pch=ifelse(pch.t>0,24,25), bg= clu.turn[,'myColor'], col=clu.turn[,'myColor'],cex=1.5)
text((clu.turn[,'warn_time']-1),rep(1.8, NROW(clu.turn))*-pch.t, 
     clu.turn[,'ID'], col=clu.turn[,'myColor'],cex=0.8)


################################################################################
### loading CCM results
ccm.otu <- read.csv('./Output/Demo_ccm_OTU_CH4_Species_D47.csv', header=T)
ccm.crn <- read.csv('./Output/Demo_ccm_metabolite_CH4_D47.csv', header=T)

ccm.otu.c <- ccm.otu[,'CCM_Control']
ccm.crn.c <- ccm.crn[,'CCM_Control']

ccm.otu <- ccm.otu[,'CCM_Perturbed']
ccm.crn <- ccm.crn[,'CCM_Perturbed']

ccm.crnotu <- c(ccm.crn,ccm.otu)
names(ccm.crnotu) <- colnames(crnotu)

ccm.net <- read.csv('./Output/Demo_ccm_Net_CH4_D47.csv', header=T)

################################################################################
## Computed the smoothed estimate of variable trends in microbiome & metabolome
if(T){
  smo.crnotu3 <- matrix(NA, nco3*nnd ,5)
  for(i in 1:nco3){
    cat("\r",round(i/nco3*100,1),'% of the network nodes')
    da.t <- data.frame(time=tind.o, value=crnotu[,i])
    da.t1 <- filter(da.t,time<=tipp)
    g.t1 <- gam(value~s(time),data=da.t1)
    yhat <- predict.gam(g.t1, newdata=data.frame(time=newd))
    fd <- fderiv(g.t1, newdata = data.frame(time=newd))
    # Calculate 95% confidence intervals for derivatives
    ci <- as.matrix(confint(fd, type = "confidence")[,c('lower','est','upper')])
    yhatd1 <- ci[,'est']
    smo.crnotu3[(nnd*(i-1)+1):(nnd*i),] <- cbind(time=newd,yhat, ci)
  }
  smo.crnotu3 <- data.frame(ID=rep(colnames(crnotu), each=nnd),smo.crnotu3)
  colnames(smo.crnotu3) <- c('ID','time','yhat','lower','est','upper')
  write.csv(filter(smo.crnotu3, time%in%(1:nd)), './Output/smo_otucrn3_D47.csv',row.names = F)
}

################################################################################
## Bifurcation analysis for microbial OTUs
dev_crnotu <- read.csv('./Output/DEV_D47_TS_crnotu.csv',header = T)
ind.tip <- tipp2#nd#

devf3=filter(dev_crnotu,variable%in%colnames(otu3))

udev3 <- unique(devf3[,'variable'])
dev.mx <- matrix(NA,nd,length(udev3))
for(i in 1:length(udev3)){
  da.t <- filter(devf3, variable==udev3[i])
  dev.mx[(1:nd)%in%da.t[,'t_max'],i] <- da.t[,'Abs']
}
colnames(dev.mx) <- udev3

crip <- NULL
j <- 1
for(i in 1:NCOL(dev.mx)){
  crip.i <- NULL
  da.t2 <- data.frame(time=1:ind.tip ,value=dev.mx[1:ind.tip,i])
  if(sum(!is.na(da.t2[,'value']))>2){
    crip.i <- find.bif2(da.t = da.t2, ind.tip)
  }
  if(NROW(crip.i)>0) crip=rbind(crip, data.frame(variable=rep(i),crip.i))
}  

crip.otu <- data.frame(ID=colnames(dev.mx)[crip[,1]],crip)

################################################################################
### Bifurcation analysis for metabolites
devf3 <- filter(dev_crnotu,variable%in%colnames(crn3))

udev3 <- unique(devf3[,'variable'])
dev.mx <- matrix(NA,nd,length(udev3))
for(i in 1:length(udev3)){
  da.t <- filter(devf3, variable==udev3[i])
  dev.mx[(1:nd)%in%da.t[,'t_max'],i]=da.t[,'Abs']
}
colnames(dev.mx) <- udev3

crip <- NULL
j <- 1
for(i in 1:NCOL(dev.mx)){
  crip.i <- NULL
  da.t2 <- data.frame(time=1:ind.tip ,value=dev.mx[1:ind.tip ,i])
  if(sum(!is.na(da.t2[,'value']))>2){
    crip.i <- find.bif2(da.t = da.t2, ind.tip)
  }
  if(NROW(crip.i)>0) crip <- rbind(crip, data.frame(variable=rep(i),crip.i))
}  

crip.crn <- data.frame(ID=colnames(dev.mx)[crip[,1]],crip)

################################################################################
### Bifurcation analysis for network properties
devf3 <- read.csv('./Output/OLR_R3_network_D47_DEV_TimeSeries.csv',header=T)
udev3 <- unique(devf3[,'variable'])

dev.mx <- matrix(NA,nd,length(udev3))
for(i in 1:length(udev3)){
  da.t <- filter(devf3, variable==udev3[i])
  dev.mx[(1:nd)%in%da.t[,'t_max'],i] <- da.t[,'Abs']
}
colnames(dev.mx) <- udev3

crip <- NULL
j <- 1
for(i in 1:NCOL(dev.mx)){
  crip.i <- NULL
  da.t2 <- data.frame(time=1:ind.tip ,value=dev.mx[1:ind.tip ,i])
  if(sum(!is.na(da.t2[,'value']))>2){
    crip.i <- find.bif2(da.t = da.t2, ind.tip)
  }
  if(NROW(crip.i)>0) crip=rbind(crip, data.frame(variable=rep(i),crip.i))
}  
crip.net <- data.frame(ID=udev3[crip[,'variable']],crip)

################################################################################
# Temporal pattern analysis for microbiome & metabolome
################################################################################
co.tmp <- tpp <- NULL
ntmp <- ntmp2 <- NULL
for(i in 1:nco3){
  #cat("\r",round(i/nco3*100,1),'% of the network nodes')
  smo.i <- filter(smo.crnotu3,ID==colnames(crnotu)[i])
  yhat <- smo.i[,'yhat']
  yhatd1 <- smo.i[,'est']
  pk.i=signch2(yci=smo.i[,c('lower','est','upper')],x=unlist(c(newd)))
  ntmp <- c(ntmp, NROW(pk.i))
  if(NROW(pk.i)==6){ntmp2 <- c(ntmp2, paste(pk.i[,'sign'],collapse = '-'))}
  da.tmp <- data.frame(variable=rep(i,NROW(pk.i)),ID=colnames(crnotu)[i],pk.i)
  #turn.i <- warn_turn(x=1:tipp, y=crnotu[1:tipp,i], tmp=da.tmp)
  turn.i <- warn_turn2(y=crnotu[1:tipp,i], yci=data.frame(x=newd,ci), pk=da.tmp)
  
  tpp <- c(tpp,turn.i[[1]])
  co.tmp <- rbind(co.tmp, turn.i[[2]])
}
co.tmp <- data.frame(ID=colnames(crnotu),Type=ind.co[,'Type'],Ttype= tpp, co.tmp)
rownames(co.tmp)= colnames(crnotu)

write.csv(co.tmp, './Output/Localized_Transition_crnotu_D47.csv',row.names = F)

###############################################################################
# Temporal pattern analysis for network properties
################################################################################
netprop.o=read.csv('./Output/demo_network_prop_ccm_perturbed_D47.csv',header=T)
netprop.p=filter(netprop.o, Reactor=='Perturbed')
netprop=t(netprop.p[,-c(1:4)]);colnames(netprop)=netprop.p[,'Property'];rownames(netprop)=NULL

net.tmp <- tpp <- NULL
for(i in 1:NCOL(netprop)){
  cat("\r",round(i/NCOL(netprop)*100,1),'% of the network nodes')
  da.t1=data.frame(time=1:tipp, value=netprop[1:tipp,i])
  g.t1 <- gam(value~s(time),data=da.t1)
  yhat=predict.gam(g.t1,newdata=data.frame(time=newd))
  yhatd1 <- derivatives(g.t1, type = "central",data=data.frame(time=newd), order=1)$derivative
  fd <- fderiv(g.t1, newdata = data.frame(time=newd))
  # Calculate 95% confidence intervals for derivatives
  ci <- as.matrix(confint(fd, type = "confidence")[,c('lower','est','upper')])
  pk.i=signch2(yci=ci,x=unlist(c(newd)))
  da.tmp <- data.frame(variable=rep(i,NROW(pk.i)),ID=netprop.p[,'Property'][i],pk.i)
  #turn.i <- warn_turn(x=1:tipp, y=da.t1[,'value'], tmp=da.tmp)
  turn.i <- warn_turn2(y=da.t1[,'value'], yci=data.frame(x=newd,ci), pk=da.tmp)
  
  tpp=c(tpp, turn.i[[1]])
  net.tmp=rbind(net.tmp, turn.i[[2]])
}

net.tmp=data.frame(ID=netprop.p[,'Property'],Type=netprop.p[,'Property_type'],Ttype= tpp, net.tmp)
rownames(net.tmp)= netprop.p[,'Property']

############################################################################################
# Detect localized transitions in microbiome, metabolome, and dynamical interactome
# microbial OTUs
crip.otu2 <- filter(crip.otu ,BFcross <= ind.tip)
otu.mvs <- co.tmp[grep('OTU_',co.tmp[,'ID']),]
fin.mvs <- data.frame(otu.mvs,BFcross=NA,BFwarn=NA,BFmax=NA, CCM=ccm.otu)
fin.mvs[fin.mvs[,'ID']%in%crip.otu2[,'ID'],c('BFcross', 'BFwarn', 'BFmax')]=crip.otu2[,c('BFcross', 'BFwarn', 'BFmax')]
fin.mvs[,'Type'] <- factor(clasdo[,'Functional_guild'], level=c('Methanogen','Acetogen', 'Other_taxa'))
fin.mvs <- fin.mvs%>%arrange(Type, desc(warn_time))%>%filter(Ttype%in%c('B1','B2','D1','D2'))
otu.g <- ind.co[fin.mvs[,'ID'],]

win.graph(150,90)
par(mfcol=c(1,3),mar=c(4,6,1,1))
ytick <- nrow(fin.mvs):1
fg.pos1 <- c(min(ytick[fin.mvs[,'Type']=='Methanogen']),min(ytick[fin.mvs[,'Type']=='Acetogen']))
fg.pos2 <- c(max(ytick[fin.mvs[,'Type']=='Methanogen']),max(ytick[fin.mvs[,'Type']=='Acetogen']),max(ytick[fin.mvs[,'Type']=='Other_taxa']))
plot(y=ytick,x=fin.mvs[,'warn_time']-1,xlim=c(1,tipp2)-1,
     col=ifelse(fin.mvs[,'Ttype']%in%c('B1','D1'),'#33c4ff','#fa8379'), ,pch=16,cex= (fin.mvs[,'CCM'])*3,ylab='',xlab='Time', yaxt="n")
abline(h=fg.pos1-0.8,col='grey',lty=2)
text(rep(8,3), fg.pos2-c(-1.5,1.5,1.5),c('Methanogen','Acetogen', 'Other taxa'),cex=2)

points(y=ytick,x=fin.mvs[,'BFcross']-1,col='blue',pch=17,cex= 1.5)
axis(side=2, at=nrow(fin.mvs):1, labels = FALSE)
text(par("usr")[1], ytick,  cex=0.4,labels = clasdo[fin.mvs[,'ID'], 'finest_taxa_name'], srt = 0, pos = 2, xpd = TRUE)

warn.otu <- fin.mvs %>% filter(!is.na(BFcross)) %>% dplyr::select(ID, Type, Ttype, warn_time, BFcross, CCM)

detail_otu <- data.frame(clasdo[fin.mvs[,'ID'], ],fin.mvs[,c('warn_time','BFcross','CCM')]);colnames(detail_otu) <- c(colnames(clasdo), 'Transition','Bifurcation','CCM')
detail_otu

# metabolites
crip.crn2 <- filter(crip.crn , BFcross <= ind.tip)
crn.mvs <- co.tmp[-grep('OTU_',co.tmp[,'ID']),]
fin.mvs <- data.frame(crn.mvs,BFcross=NA,BFwarn=NA,BFmax=NA, CCM=ccm.crn)
fin.mvs[fin.mvs[,'ID']%in%crip.crn2[,'ID'],c('BFcross', 'BFwarn', 'BFmax')]=crip.crn2[,c('BFcross', 'BFwarn', 'BFmax')]
fin.mvs <- fin.mvs%>%arrange(desc(Type), desc(warn_time))%>%filter(Ttype%in%c('B1','B2','D1','D2'))#,'D1','D2'))
crn.g <- ind.co[fin.mvs[,'ID'],]
warn.crn <- fin.mvs[!is.na(fin.mvs[,'BFcross']),]

#win.graph(50,90)
#par(mar=c(4,6,1,1))
ytick <- nrow(fin.mvs):1
fg.pos <- min(ytick[fin.mvs[,'Type']=='VFA'])
plot(y=ytick,x=fin.mvs[,'warn_time']-1,xlim=c(1,tipp2)-1,
     col=ifelse(fin.mvs[,'Ttype']%in%c('B1','D1'),'#33c4ff','#fa8379'),pch=16,cex= (fin.mvs[,'CCM'])*3,ylab='',xlab='Time', yaxt="n")
abline(h=fg.pos-0.8,col='grey',lty=2)
text(rep(7,2), c(max(ytick)+2,max(ytick)-8),c('VFA','Chemical \nFingerprint'),cex=2)

points(y=ytick,x=fin.mvs[,'BFcross']-1,col='blue',pch=17,cex= 1.5)
axis(side=2, at=nrow(fin.mvs):1, labels = FALSE)
text(par("usr")[1], ytick,  cex=0.6, 
     labels = fin.mvs[,'ID'], srt = 0, pos = 2, xpd = TRUE)

warn.crn <- fin.mvs %>% filter(!is.na(BFcross)) %>% dplyr::select(ID, Type, Ttype, warn_time, BFcross, CCM)

warn.crnotu <- rbind(warn.crn,warn.otu)

write.csv(warn.crnotu, './Output/warn_time_CRNOTU_D47.csv',row.names = F)

detail_crn <- data.frame(fin.mvs[,c('ID','Type')],fin.mvs[,c('warn_time','BFcross','CCM')]);colnames(detail_crn) <- c(c('Chemical_formula','Type','Transition','Bifurcation','CCM'))
detail_crn

# Network properties
crip.net2 <- filter(crip.net , BFcross <= ind.tip)
fin.mvs <- data.frame(net.tmp,BFcross=NA,BFwarn=NA,BFmax=NA,CCM=as.numeric(c(ccm.net[,'CCM_Perturbed'])))
fin.mvs[fin.mvs[,'ID']%in%crip.net2[,'ID'],c('BFcross', 'BFwarn', 'BFmax')]=crip.net2[,c('BFcross', 'BFwarn', 'BFmax')]
fin.mvs[,'Type'] <- factor(fin.mvs[,'Type'], level=c('Centrality','Guild_Int', 'Feedback', 'Int_type','Topology'))
fin.mvs <- fin.mvs%>%arrange(Type, desc(warn_time))%>%filter(Ttype%in%c('B1','B2','D1','D2'))

#win.graph(50,90)
#par(mar=c(4,8.5,1,1))
ytick=nrow(fin.mvs):1
fg.pos1=c(min(ytick[fin.mvs[,'Type']=="Centrality"]),min(ytick[fin.mvs[,'Type']=="Guild_Int"]), min(ytick[fin.mvs[,'Type']=="Feedback"]), min(ytick[fin.mvs[,'Type']=="Int_type"]))
fg.pos2=c(max(ytick[fin.mvs[,'Type']=="Centrality"]),max(ytick[fin.mvs[,'Type']=="Guild_Int"]), max(ytick[fin.mvs[,'Type']=="Feedback"]), max(ytick[fin.mvs[,'Type']=="Int_type"]),max(ytick[fin.mvs[,'Type']=="Topology"]))
plot(y=ytick,x=fin.mvs[,'warn_time']-1,xlim=c(1,tipp2)-1,
     col=ifelse(fin.mvs[,'Ttype']%in%c('B1','D1'),'#33c4ff','#fa8379'), 
     pch=16, cex= fin.mvs[,'CCM']*3, ylab='',xlab='Time', yaxt="n")
abline(h=fg.pos1-0.8,col='grey',lty=2)
text(rep(8,4), fg.pos2-c(0,0.8,0.8,0.8),c("Centrality", 'Guild Interaction', "Feedback", "Interaction type", 'Topology'),cex=2)
legend('bottomright',c('Turning point', 'Bifurcation point'),bty='n',col=c('black','blue'),pch=c(1,17),cex=1.5)

points(y=ytick,x=fin.mvs[,'BFcross']-1,col='blue',pch=17,cex= 1.5)

axis(side=2, at=nrow(fin.mvs):1, labels = FALSE)
text(par("usr")[1], ytick,  cex=0.75, labels = fin.mvs[,'ID'], srt = 0, pos = 2, xpd = TRUE)

detail_net <- data.frame(fin.mvs[,c('ID','Type')],fin.mvs[,c('warn_time','BFcross','CCM')]);colnames(detail_net) <- c(c('Property','Type','Transition','Bifurcation','CCM'))
detail_net
