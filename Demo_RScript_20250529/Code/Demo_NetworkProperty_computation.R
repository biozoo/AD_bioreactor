# R code for calculating Network property
rm(list = ls());gc()
# Set path
OPath <- 'D:/data/Cooperation/Prodigio/Code/Demo_RScript_20250529'
setwd(OPath)

library(igraph)
library(dplyr)
library(mgcv)
library(EnvStats)

source('./Code/Function_network.R')
seed <- 49563
set.seed(seed)
SaveFile <- T#F #/F for saving files

#############################################################
#### Data Preparation
tipp <- 38
tipp2 <- 44
nd <- 46
# Data loading
clasdo <- read.csv('./Data/demo_otu_taxa_D47.csv',header=T)
rownames(clasdo) <- clasdo[,'OTU']
crnotu1 <- read.csv('./Data/OLR_Control_crnotu_verc_CH4_Species_D47.csv',header=T,stringsAsFactors = F)
crnotu3 <- read.csv('./Data/OLR_Perturbed_crnotu_verc_CH4_Species_D47.csv',header=T,stringsAsFactors = F)

# Calculate percentage data for microbial OTU 
otu1 <- crnotu1[,grep('OTU_',colnames(crnotu1))]
otu3 <- crnotu3[,grep('OTU_',colnames(crnotu3))]

# Extract metabolomic data
crn1 <- crnotu1[,-grep('OTU_',colnames(crnotu1))]
crn3 <- crnotu3[,-grep('OTU_',colnames(crnotu3))]

rownames(otu3) <- rownames(crn3) <- 0:46

nsp <- NCOL(otu3)
ncrn <- NCOL(crn1)

nco1 <- ncrn+nsp
nco3 <- ncrn+nsp+7

crnotu3s <- apply(crnotu3,2,scale)

ind.co <- data.frame(Type=c(rep('VFA',7),rep('Fingerprint',ncrn), clasdo[,'Functional_guild']),
                  Name=colnames(crnotu3),Abb=c(rep(1,7),rep(2,ncrn),
                                               match(clasdo[,'Functional_guild'],c("Methanogen","Acetogen","Other_taxa"))+2))

ind.co <- mutate(ind.co, Abb1= 1*(Abb<=2)+2*(Abb>2))

rownames(ind.co) <- ind.co[,'Name']

Netpro <- c('Guild_Int','Topology','Centrality','Int_type','Feedback')

ind1.co <- ind.co[-c(1:7),'Abb']
ind3.co <- ind.co[,'Abb']

ind11.co <- ind.co[-c(1:7),'Abb1']
ind31.co <- ind.co[,'Abb1']

################################################################################################
### Network analysis
daname.c <- 'Control_otucrn_CH4_D47_Species_NA'
daname.p <- 'Perturbed_otucrn_CH4_D47_Species_NA'  
co1.ins <- read.csv(paste('./Output/',daname.c,'_cvunit0.025_enet_jcof_Nmvx_Rallx_demo.csv',sep=''),stringsAsFactors = F)
co3.ins <- read.csv(paste('./Output/',daname.p,'_cvunit0.025_enet_jcof_Nmvx_Rallx_demo.csv',sep=''),stringsAsFactors = F)

tind <- sort(unique(co3.ins[,'time']))

# Interaction array 
co.ary1 <- coa.ary1 <- array(NA,dim=c(nco1,nco1,length(tind)))
co.ary3 <- coa.ary3 <- array(NA,dim=c(nco3,nco3,length(tind)))

for(i in tind){
  ints.i11 <- as.matrix(filter(co1.ins,time==i)[,-c(1:4)])
  ints.i31 <- as.matrix(filter(co3.ins,time==i)[,-c(1:4)])
  co.ary1[,,i] <- ints.i11
  co.ary3[,,i] <- ints.i31
  }

# Time integrated
co.ins3 <- apply(co.ary3,1:2,mean,na.rm=T)
co.ins1 <- apply(co.ary1,1:2,mean,na.rm=T)

### Computation of network properties
###############################################################################
### Block interaction matrix between functional components

if(T){
  co.bary1 <- blockj.a(ism.a=co.ary1, member=ind.co[-c(1:7),'Abb'])
  co.bary3 <- blockj.a(ism.a=co.ary3, member=ind.co[,'Abb'])
  co1.bary1 <- blockj.a(ism.a=co.ary1, member=ind11.co)
  co1.bary3 <- blockj.a(ism.a=co.ary3, member=ind31.co)
  
  bj.co11 <- filter(blockins(co1.bary1),sign==0)
  bj.co31 <- filter(blockins(co1.bary3),sign==0)
  bj.co1 <- filter(blockins(co.bary1),sign==0);bj.co1[,c('cause','effect')] <- bj.co1[,c('cause','effect')]+1
  bj.co3 <- filter(blockins(co.bary3),sign==0)
  
  bj.da <- rbind(data.frame(treat='Control',level=1,bj.co11),
              data.frame(treat='Perturbed',level=1,bj.co31),
              data.frame(treat='Control',level=2,bj.co1),
              data.frame(treat='Perturbed',level=2,bj.co3)
  )
  
  Gnet.c <- Gnet.p <- g.name <- NULL
  for(cau.i in 1:5){
    for(eff.j in 1:5){
      Gnet.p <- cbind(Gnet.p, as.numeric(filter(bj.co3, cause==cau.i, effect==eff.j)[,'value']))
      if(cau.i==1|eff.j==1){
        Gnet.c <- cbind(Gnet.c, rep(NA,nd))
      }else{Gnet.c <- cbind(Gnet.c, as.numeric(filter(bj.co1, cause==cau.i, effect==eff.j)[,'value']))}
      g.name <- c(g.name, paste('i', eff.j, '_j', cau.i,sep=''))
    }
  }
  Gnet <- data.frame(Reactor= rep(c('Perturbed', 'Control'),each=nd),
             Day=rep(1:nd,2),rbind(Gnet.p, Gnet.c))
  colnames(Gnet)[-c(1:2)] <- g.name
  
  if(SaveFile){
    write.csv(bj.da,'./Output/Demo_bjda_co_verc_Species.csv',row.names=F)
    write.csv(Gnet,'./Output/demo_GuildNet.csv',row.names=F)
  } 
}

################################################################################
## Topological analysis
if(T){
  topo1.co <- block.prop(ism=co.ary1,member=ind1.co, exclude=which(ind11.co==1))
  topo11.co <- block.prop(ism=co.ary1,member=ind11.co, exclude=which(ind11.co==1),centrality=F)
  topo3.co <- block.prop(ism=co.ary3,member=ind3.co, exclude=which(ind31.co==1))
  topo31.co <- block.prop(ism=co.ary3,member=ind31.co, exclude=which(ind31.co==1),centrality=F)
  
  pname <- cbind(
    c(rep('Topology',3),rep('Centrality',5),rep('Int_type',9)),
    c("sdstrength.all","modularity","deig",
      "strength.in","strength.out",'strength.all',"hub" ,"auth",
      "Cooperation","Competition","Exploitated","Exploitator",
      "Commensalism1","Commensalism0","Ammensalism-1","Ammensalism0","Interaction_Diversity"
    ))
  
  topp11 <- data.frame(treat = 'Control', level=1, datransf(topo.o=topo11.co,ind.t=1:2))
  topp31 <- data.frame(treat = 'Perturbed', level=1, datransf(topo.o=topo31.co,ind.t=1:2))
  
  topp1 <- data.frame(treat = 'Control', level=2, datransf(topo.o=topo1.co,ind.t=1:2))
  topp3 <- data.frame(treat = 'Perturbed', level=2, datransf(topo.o=topo3.co,ind.t=1:2))
  
  topoda <- rbind(topp11,  topp31,  topp1,  topp3)
  topoda <- data.frame(topoda, type=pname[match(topoda[,'variable'],pname[,2]),1])
  topoda <- filter(topoda, variable!='strength.all') #Exclude redundant centrality indices
  
  if(SaveFile) write.csv(topoda,'./Output/Demo_topoda_co_verc_Species.csv',row.names=F)
}

################################################################################
# Feedback analysis 
################################################################################
# Searching for Feedbacks  
lp <- 3
if(T){
  co.fdm3 <- FindCycles.lapply(co.ins3, L=lp)
  if(SaveFile) write.csv(co.fdm3,paste('./Output/Demo_co3.fdm3_verc_L3_Species.csv',sep=''),row.names=F)
  
  co.fdm1 <- FindCycles.lapply(co.ins1, L=lp)
  if(SaveFile) write.csv(co.fdm1,paste('./Output/Demo_co1.fdm3_verc_L3_Species.csv',sep=''),row.names=F)
}

# Computation of feedback weight  
# Perturbed reactor
if(T){
  fdm <- read.csv('./Output/Demo_co3.fdm3_verc_L3_Species.csv',header=T)[,-4]
  member <- ind3.co
  umem <- sort(unique(member))
  fdm2 <- NULL 
  for(i in umem){
    mem.i <- which(member==umem[i])
    fdm2 <- cbind(fdm2, apply(fdm,1,function(x){return(any(x%in%mem.i)*1)}))
  }
  
  # All feedback loops
  fdw3.co3.all <- fdwf(ism=co.ary3, fdm=fdm, L=lp, memm=F)
  
  # Feedback within chemical subnetwork
  indx <- ((fdm2[,3]==0)+(fdm2[,4]==0)+(fdm2[,5]==0))==3
  fdw3.co3.c <- fdwf(ism=co.ary3, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # Feedback within microbial subnetwork
  indx <- ((fdm2[,1]==0)+(fdm2[,2]==0))==2
  fdw3.co3.o <- fdwf(ism=co.ary3, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves biological-microbial interaction
  indx <- ((fdm2[,1]+fdm2[,2])>0)+((fdm2[,3]+fdm2[,4]+fdm2[,5])>0)
  fdw3.co3.co <- fdwf(ism=co.ary3, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves methanogen
  indx <- fdm2[,3]
  fdw3.co3.m <- fdwf(ism=co.ary3, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves acetogen
  indx <- fdm2[,4]
  fdw3.co3.a <- fdwf(ism=co.ary3, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves methanogen & acetogen
  indx <- (fdm2[,3]+fdm2[,4])==2
  fdw3.co3.am <- fdwf(ism=co.ary3, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves methanogen & acetogen & other OTU
  indx <- (fdm2[,3]+fdm2[,4]+fdm2[,4])==3
  fdw3.co3.amo <- fdwf(ism=co.ary3, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves methanogen & acetogen & chemicals
  indx <- (((fdm2[,1]+fdm2[,2])>0)+fdm2[,3]+fdm2[,4])==3
  fdw3.co3.cam <- fdwf(ism=co.ary3, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves methanogen & chemicals
  indx <- (((fdm2[,1]+fdm2[,2])>0)+fdm2[,3])==2
  fdw3.co3.cm <- fdwf(ism=co.ary3, fdm=fdm, L=lp, memm=T,ind.t=indx)
  rm(fdm);rm(fdm2);gc()
  
  # Control reactor
  fdm <- read.csv('./Output/Demo_co1.fdm3_verc_L3_Species.csv',header=T)[,-4]
  member <- ind1.co
  umem <- sort(unique(member))
  fdm2 <- NULL 
  for(i in umem){
    mem.i <- which(member==umem[i])
    fdm2 <-  cbind(fdm2, apply(fdm,1,function(x){return(any(x%in%mem.i)*1)}))
  }
  fdm2 <- cbind(rep(0),fdm2)
  
  # all feedback
  fdw3.co1.all <- fdwf(ism=co.ary1, fdm=fdm, L=lp, memm=F)
  
  # Feedback within chemical subnetwork
  indx <- ((fdm2[,3]==0)+(fdm2[,4]==0)+(fdm2[,5]==0))==3
  fdw3.co1.c <- fdwf(ism=co.ary1, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # Feedback within microbial subnetwork
  indx <- ((fdm2[,1]==0)+(fdm2[,2]==0))==2
  fdw3.co1.o <- fdwf(ism=co.ary1, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves biological-microbial interaction
  indx <- ((fdm2[,1]+fdm2[,2])>0)+((fdm2[,3]+fdm2[,4]+fdm2[,5])>0)
  fdw3.co1.co <- fdwf(ism=co.ary1, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves methanogen
  indx <- fdm2[,3]
  fdw3.co1.m <- fdwf(ism=co.ary1, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves acetogen
  indx <- fdm2[,4]
  fdw3.co1.a <- fdwf(ism=co.ary1, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves methanogen & acetogen
  indx <- (fdm2[,3]+fdm2[,4])==2
  fdw3.co1.am <- fdwf(ism=co.ary1, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves methanogen & acetogen & other OTU
  indx <- (fdm2[,3]+fdm2[,4]+fdm2[,4])==3
  fdw3.co1.amo <- fdwf(ism=co.ary1, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves methanogen & acetogen & chemicals
  indx <- (((fdm2[,1]+fdm2[,2])>0)+fdm2[,3]+fdm2[,4])==3
  fdw3.co1.cam <- fdwf(ism=co.ary1, fdm=fdm, L=lp, memm=T,ind.t=indx)
  
  # feedback involves methanogen & chemicals
  indx <- (((fdm2[,1]+fdm2[,2])>0)+fdm2[,3])==2
  fdw3.co1.cm <- fdwf(ism=co.ary1, fdm=fdm, L=lp, memm=T,ind.t=indx)
  rm(fdm);rm(fdm2);gc()
  
  fdda <- rbind(
    data.frame(treat='Perturbed', type='Feedback', level=rep(c(rep(1,4),rep(2,5)),each=nd), block=rep(c(0,12,11,22,3,4,34,1234,123),each=nd),L=3, sign=0, time=rep(tind,9), 
               value=c(fdw3.co3.all,fdw3.co3.co,fdw3.co3.c,fdw3.co3.o,fdw3.co3.m,fdw3.co3.a,fdw3.co3.am,fdw3.co3.cam,fdw3.co3.cm)),
    data.frame(treat='Control', type='Feedback', level=rep(c(rep(1,4),rep(2,5)),each=nd), block=rep(c(0,12,11,22,3,4,34,1234,123),each=nd),L=3, sign=0, time=rep(tind,9), 
               value=c(fdw3.co1.all,fdw3.co1.co,fdw3.co1.c,fdw3.co1.o,fdw3.co1.m,fdw3.co1.a,fdw3.co1.am,fdw3.co1.cam,fdw3.co1.cm))
  )  
  if(SaveFile) write.csv(fdda,'./Output/Demo_fdda_co_verc_Species.csv',row.names=F)
}

#################################################################################################
bjda <- read.csv('./Output/Demo_bjda_co_verc_Species.csv',header=T)
topoda <- read.csv('./Output/Demo_topoda_co_verc_Species.csv',header=T)
fdda <- read.csv('./Output/Demo_fdda_co_verc_Species.csv',header=T)

bjda.nm <- apply(bjda[,c('level','type','sign','cause','effect')], 1, paste, collapse='-')
bjda.unm <- unique(bjda.nm)

topo.nm <- apply(topoda[,c('level','type','variable','member')], 1, paste, collapse='-')
topo.unm <- unique(topo.nm)

fdda.nm <- NULL
for(i in 1:NROW(fdda)) {
  fdda.nm <- c(fdda.nm, paste(fdda[i,c('level','type','sign','block')], collapse='-'))}
fdda.unm <- unique(fdda.nm)

netda <- rbind(
  data.frame(bjda[,c('treat','level','type')],name=bjda.nm,bjda[,c('time','value')]),
  data.frame(topoda[,c('treat','level','type')],name=topo.nm,topoda[,c('time','value')]),
  data.frame(fdda[,c('treat','level','type')],name=fdda.nm,fdda[,c('time','value')])
)

if(SaveFile) write.csv(netda,'./Output/Demo_netda_co_verc_Species_D47.csv',row.names=F)

###################################################################################################
## Comparison between the control and perturbed reactor
netda <- read.csv('./Output/Demo_netda_co_verc_Species_D47.csv',header=T)

netda.c <- filter(netda,treat=='Control')
netda34 <- filter(netda,treat=='Perturbed')

nam3 <- unique(netda34[,'name'])

eff.nam <- NULL
for(i in 1:length(nam3)){
  da.t12 <- filter(netda.c,name==nam3[i])
  da.t34 <- filter(netda34,name==nam3[i])
  eff.nam[i] <- (!all(is.na(da.t34[,'value'])))*1+(NROW(da.t12)>0 & !all(is.na(da.t34[,'value'])))*2
}

nam3 <- nam3[eff.nam>0]
eff.nam <- eff.nam[eff.nam>0]

netda.p <- filter(netda34,name%in%nam3)

nd <- max(netda.p[,'time'])
gamda <- gamda3 <- NULL
t.begin <- proc.time()
for(i in 1:length(nam3)){
  da.t12 <- filter(netda.c,name==nam3[i])
  da.t34 <- filter(netda.p,name==nam3[i])
  da.t3 <- filter(da.t34,time<=tipp)
  g.t <- gam(value~s(time),data=da.t3)
  k.t <- kendallTrendTest(da.t3[,'value'])
  gamda3 <- rbind(gamda3, c(unlist(c(rsq=summary(g.t)$r.sq, devex=summary(g.t)$dev.expl, as.data.frame(summary(g.t)$s.table), k.t$estimate, pk=k.t$p.value) ) ) )
  
  if(eff.nam[i]==3){
    da.t0 <- da.t0r <- da.t34;da.t0[,'treat'] <- 0;da.t0r[,'treat'] <- 0
    da.t0[,'value'] <- scale(cinter(da.t34[,'value']))-scale(cinter(da.t12[,'value']))
    da.t0r[,'value'] <- scale(da.t34[,'value'])-scale(da.t12[,'value'])
    da.t1 <- filter(da.t0,time<=tipp)
    if(!all(is.na(da.t1[,'value']))){
      g.t1 <- gam(value~s(time),data=da.t1)
      gamda <- rbind(gamda, c(unlist(c(rsq=summary(g.t1)$r.sq, devex=summary(g.t1)$dev.expl, as.data.frame(summary(g.t1)$s.table)))))
    }else{gamda <- rbind(gamda,rep(NA));}
  }else{
    gamda <- rbind(gamda,rep(NA));
  }
  time.used <- proc.time() - t.begin
  cat("\r",paste('property', i, round(i/length(nam3)*100,2),"%", round(time.used[3],2),"sec"))
}

gamda <- data.frame(name=nam3,gamda)%>%mutate(significance=(p.value<0.05)*1)
gamda3 <- data.frame(name=nam3,gamda3)
gamda[is.na(gamda[,'significance']),'significance'] <- 0

ind.t <- which(((gamda[,'p.value']<0.05)|eff.nam==1) & (gamda3[,'p.value']<0.05))
uvs <- nam3[ind.t]

netda.pv <- filter(netda.p,name%in%uvs)
netda.cv <- filter(netda.c,name%in%uvs)

netype <- c('Chem','Bio', 'Chem_Bio')
fg <- list(c('Chem', 'Bio'),
           c('VFA','Fingerprint','Methanogen','Acetogen','Other_taxa'))

Property_type <- Property <- Network_type <- c()
value.c <- value.p <- NULL
for(i in 1:length(uvs)){
  net.pi <- filter(netda.pv, name==uvs[i])
  net.ci <- filter(netda.cv, name==uvs[i])
  value.p <- rbind(value.p,net.pi[,'value'])
  if(NROW(net.ci)==0) {value.c <- rbind(value.c,rep(NA,nd))}else{value.c <- rbind(value.c,net.ci[,'value'])}
  
  uvs.i <- unlist(strsplit(uvs[i],'-'))
  Property_type <- c(Property_type, uvs.i[2])
  fg.i <- fg[[as.numeric(uvs.i[1])]]
  ni <- length(uvs.i)
  subnet.i <- uvs.i[ni]
  if(Property_type[i]=='Guild_Int'){
    subnet.i <- c(uvs.i[ni-1], subnet.i)
    if(subnet.i[1]==subnet.i[2]){
      subnet.i <- subnet.i[1]
      Property_type[i] <- 'Topology'
      uvs.i[3] <- 'Mean_Int_strength'}
    }else if(Property_type[i]=='Feedback'){
      subnet.j <- c()
      for(j in 1:nchar(subnet.i)){subnet.j[j] <- substring(subnet.i,j,j)}
      subnet.i <- subnet.j
  }
  subnet.i <- as.numeric(subnet.i)
  involve.i <- unique(fg.i[subnet.i])
  
  chem.i <- any(involve.i%in%c('Chem','VFA', 'Fingerprint'))
  bio.i <- any(involve.i%in%c('Bio','Methanogen','Acetogen','Other_taxa'))
  if(any(subnet.i==0)){chem.i <- bio.i <- T}
  Network_type[i] <- netype[chem.i*1 + bio.i*2]
  if(Property_type[i]=='Guild_Int'){prop.i <- paste(involve.i,collapse=' -> ')}
  if(Property_type[i]%in%c('Topology','Centrality','Int_type')){prop.i <- paste(c(involve.i,uvs.i[3]),collapse='_')}
  if(Property_type[i]%in%'Feedback'){
    if(any(subnet.i!=0)){prop.i <- paste('Feedback', paste(involve.i, collapse='-'))}else{prop.i <- paste('Feedback', 'Whole')}}
  Property[i] <- prop.i
}

value.a <- rbind(value.p, value.c)
colnames(value.a) <- paste('Day',1:nd, sep='_')
netprop <- data.frame(Reactor=rep(c('Perturbed', 'Control'),each=length(uvs)), Property=rep(Property,2), 
           Property_type=rep(Property_type,2), Network_type=rep(Network_type,2), value.a)

if(SaveFile) write.csv(netprop,'./Output/demo_network_prop_perturbed_D47.csv',row.names=F)

#######################
## Causal analysis
source('./Code/Function_MDR_inferM.R')
library(rEDM)
library(Kendall)
netprop.o <- read.csv('./Output/demo_network_prop_perturbed_D47.csv',header=T)
netprop.p <- filter(netprop.o, Reactor=='Perturbed')
prop.name <- netprop.p[,'Property']
netprop.c <- t(filter(netprop.o, Reactor=='Control')[,-c(1:4)])
netprop.p <- t(netprop.p[,-c(1:4)]) 
colnames(netprop.p) <- colnames(netprop.c) <- c(prop.name)
rownames(netprop.p) <- rownames(netprop.c) <- NULL
nsp <- NCOL(netprop.p)

envi <- read.csv('.//Data//Reactor_data_envi.csv',header=T,stringsAsFactors = F)

inferM <- matrix(0, nsp+1, nsp+1);inferM[1,-1] <- 1
ccm.c <- ccm.fast.inferM(ds=apply(cbind(CH4=envi[-nd,'CH4_C1'], netprop.c),2,scale),inferM,Lmin=5,Lint=10)
ccm.p <- ccm.fast.inferM(ds=apply(cbind(CH4=envi[-nd,'CH4_O1'], netprop.p),2,scale),inferM,Lmin=5,Lint=10)

(spind <- which(ccm.p[[2]][1,]==1)-1)
ccm.pw <- ccm.p[[1]][1,-1][spind]
ccm.cw <- ccm.c[[1]][1,-1][spind]

ccm.net <- data.frame(Variable=prop.name[spind],
                      CCM_Control=ccm.cw, 
                      CCM_Perturbed=ccm.pw)

if(SaveFile) {
  write.csv(ccm.net,'./Output/Demo_ccm_Net_CH4_D47.csv',row.names=F)
  write.csv(filter(netprop, Property%in%prop.name[spind]),'./Output/demo_network_prop_ccm_perturbed_D47.csv',row.names=F)
  }
