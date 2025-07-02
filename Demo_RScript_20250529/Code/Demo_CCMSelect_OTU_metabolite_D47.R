# R code for performing MDR S-map analysis
rm(list = ls());gc()
# Set path
Root_Path <- 'D:'
Sub_Path <- '\\data\\Cooperation\\Prodigio\\Code\\Demo_RScript_20250529\\'
OPath <- paste(Root_Path,Sub_Path,sep='')
setwd(OPath)
source('./Code/Function_MDR_inferM.R')

library(rEDM) # Some EDM analyses (e.g., simplex projection and CCM) were conducted in rEDM package ver. 1.2.3
library(doParallel)
library(parallel)
library(foreach)
library(Kendall)
library(zyp)
library(dplyr)

seed <- 49563
set.seed(seed)
SaveFile <- F #T#/F for saving files
type.ch4 <- 'CH4'

################################################################################
### Select OTU
clasotu <- read.csv('.//Data//AR_taxonomy.csv',header=T,stringsAsFactors = F, fileEncoding="latin1")
rownames(clasotu) <- clasotu[,'OTU']
envi <- read.csv('.//Data//Reactor_data_envi.csv',header=T,stringsAsFactors = F)

chlorop <- filter(clasotu, Order=='o__Chloroplast')[,'OTU']

compo.c <- read.csv('./Data/OLR_Control_otu_D47_Species.csv',header=T)
compo.p <- read.csv('./Data/OLR_Perturbed_otu_D47_Species.csv',header=T)

compo.cinfo <- compo.c[,1:5]
compo.c <- compo.c[,-c(1:5)]  
compo.pinfo <- compo.p[,1:5]
compo.p <- compo.p[,-c(1:5)]  

nd <- NROW(compo.p)

(TCount <- rowSums(compo.p))

pcri <- 2/3#
doind.c <- which(apply(compo.c!=0,2,sum,na.rm=T)>=(nd*pcri)& !names(compo.c)%in%chlorop)
doind.p <- which(apply(compo.p!=0,2,sum,na.rm=T)>=(nd*pcri)& !names(compo.c)%in%chlorop)

(doind <- intersect(doind.c,doind.p))#
nsp <- length(doind)
otu.c <- compo.c[,doind]
otu.p <- compo.p[,doind]

otu.cs <- apply(data.frame(envi[,c('CH4_C1')],otu.c),2,scale)
otu.ps <- apply(data.frame(envi[,c('CH4_O1')],otu.p),2,scale)

inferM <- matrix(0, nsp+1, nsp+1);inferM[1,-1] <- 1
ccm.c <- ccm.fast.inferM(ds=otu.cs,inferM,Lmin=5,Lint=10)
ccm.p <- ccm.fast.inferM(ds=otu.ps,inferM,Lmin=5,Lint=10)

(spind <- which(ccm.p[[2]][1,]==1)-1)
ccm.pw <- ccm.p[[1]][1,-1][spind]
ccm.cw <- ccm.c[[1]][1,-1][spind]

ccm.otu <- data.frame(Variable=colnames(otu.p)[spind],CCM_Control=ccm.cw, CCM_Perturbed=ccm.pw)

if(SaveFile){
  write.csv(ccm.otu,paste('./Output/Demo_ccm_OTU',type.ch4,'Species_D47.csv',sep='_'),row.names=F)
  write.csv(otu.c[,spind],paste('./Data/OLR_Control_otu_verc',type.ch4,'Species_D47.csv',sep='_'),row.names=F)
  write.csv(otu.p[,spind],paste('./Data/OLR_Perturbed_otu_verc',type.ch4,'Species_D47.csv',sep='_'),row.names=F)
}

meth.ls=list(
  Phylum=c('Euryarchaeota','Halobacterota','Thermoplasmatota'),
  Class=c('Methanobacteria','Thermococci','Methanomicrobia','Methanosarcinia','Thermoplasmata'),
  Order=c('Methanosarciniales', 'Methanomassiliicoccales','Methanobacteriales','Methanofastidiosales','Methanomicrobiales'),
  Family=c('Methanobacteriaceae','Methanothermobacteriaceae','Methanofastidiosaceae','Methanomicrobiaceae','Methanocorpusculaceae',
           'Methanoregulaceae','Methanospirillaceae','Methanosarcinaceae','Methanosaetaceae','Methanomassiliicoccaceae'),
  Genus=c('Methanobacterium','Methanobrevibacter','Methanothermobacter','Ca_Methanofastidiosum','Methanoculleus',
          'Methanocorpusculum','Methanolinea','Methanospirillum','Methanosarcina','Methanomethylovorans','Methanothrix','Methanomassiliicoccus')
)

aceto.ls=list(
  Class=c('Syntrophia','Syntrophorhabdia','Syntrophomonadia'),
  Order=c('Syntrophales', 'Syntrophorhabdales','Syntrophomonadales'),	
  Family=c('Syntrophaceae','Smithellaceae','Syntrophomonadaceae','Syntrophorhabdaceae','Syntrophotaleaceae','Pelotomaculaceae'),
  Genus=c('Syntrophus', 'Smithella', 'Syntrophobacter', 'Syntrophorhabdus', 'Syntrophotalea', 'Syntrophomonas', 'midas_g_1312', 'Pelotomaculum', 'Anaerovorax')
)

clasdo <- NULL
for(i in 1:NROW(ccm.otu)){clasdo <- rbind(clasdo, filter(clasotu,OTU%in%ccm.otu[i,'Variable']))}
otu.nam.s <- unlist(lapply(clasdo[,'Genus'],function(x){strsplit(x,'g__')[[1]][2]}))
otu.nam.f <- unlist(lapply(clasdo[,'Family'],function(x){strsplit(x,'f__')[[1]][2]}))

func_otu <- c('Methanogen','Acetogen','Other_taxa')
ind.otu <- rep(3,NROW(ccm.otu))-
  2*otu.nam.f%in%meth.ls[['Family']]-
  1*otu.nam.s%in%aceto.ls[['Genus']]


clasdo <- data.frame(clasdo, Finest_taxon=apply(clasdo, 1, function(x){
  ind <- which((!grepl('midas',x[1:6]))*(!x[1:6]=='')==1)
  if(length(ind)>0){y <- x[max(ind)]}else{y <- x[8]}
  return(y)}))
clasdo <- data.frame(Functional_guild=func_otu[ind.otu], clasdo)

richness.c <- apply(compo.c[,!colnames(compo.c)%in%chlorop]!=0,1,sum)
richness.p <- apply(compo.p[,!colnames(compo.p)%in%chlorop]!=0,1,sum)

if(SaveFile){
  write.csv(data.frame(OTU_rich_C1=richness.c, OTU_rich_O1=richness.p),'./Data/OTU_richness.csv',row.names=F)
  write.csv(clasdo,'./Data/demo_otu_taxa_D47.csv',row.names=F)
}

################################################################################
### select metabolite
VFA <- envi[,c('HAc_O1', 'HPro_O1', 'i.HBu_O1', 'HBu_O1', 'iHVal_O1', 'HVal_O1', 'HCa_O1')]
finger.co <- read.csv('./Data/OLR_Control_finger_D47.csv',header=T)
finger.po <- read.csv('./Data/OLR_Perturbed_finger_D47.csv',header=T)


if(F){
  # Function for analyzing chemical formula 
  chemfor <- function(nam.o, element= c('C','H','N','O','S')){
    nam.i <- nam.o
    for(i in 1:length(element)) nam.i <- sub(element[i], paste('-',element[i],'-',sep=''), nam.i)
    
    namx <- matrix(0, ncol=5,nrow=length(nam.i))
    for(i in 1:length(element)){
      ind.i <- grep(element[i],nam.i)
      nub <- unlist(lapply(strsplit(nam.i[ind.i],paste('-',element[i],'-',sep='')),function(x) {unlist(strsplit(x[2],'-'))[1]}))
      nub[nub==""] <- 1
      nub[is.na(nub)] <- 1
      nub <- as.numeric(nub)
      namx[ind.i,i] <- nub
    }
    
    colnames(namx) <- element
    return(namx)
  }
  
  elenum <- chemfor(colnames(finger.co),element= c('C','H','N','O','S'))
  wa_sc.c <- apply(finger.co*repmat(elenum[,'S'],NROW(finger.co),1),1,sum)/apply(finger.co*repmat(elenum[,'C'],NROW(finger.co),1),1,sum)
  wa_sc.p <- apply(finger.po*repmat(elenum[,'S'],NROW(finger.po),1),1,sum)/apply(finger.po*repmat(elenum[,'C'],NROW(finger.po),1),1,sum)
  
  win.graph(100,80)
  plot(wa_sc.c~c(0:(nd-1)),type='o',pch=16,ylab='wa(S/C)',col='grey',xlab='Day',lwd=2,cex=2,ylim=c(0,0.09))
  lines(wa_sc.p~c(0:(nd-1)),type='o',pch=16,lwd=2,cex=2)
  abline(v=38,col='orange')
}


nd <- NROW(finger.po)
pcri <- 2/3
ind.tp <- which(apply(finger.po!=0,2,sum,na.rm=T)>=((nd-sum(is.na(finger.po[,1])))*pcri))
ind.tc <- which(apply(finger.co!=0,2,sum,na.rm=T)>=((nd-sum(is.na(finger.co[,1])))*pcri))

ind.t <- intersect(ind.tp,ind.tc)#

finger.c <- finger.co[,ind.t]
finger.p <- finger.po[,ind.t]

nsp <- NCOL(finger.p)

ada.c <- apply(cbind(envi[,c('CH4_C1')],finger.c),2,scale)
ada.p <- apply(cbind(envi[,c('CH4_O1')],finger.p),2,scale)
ada.pvfa <- apply(cbind(envi[,c('CH4_O1')],VFA),2,scale)

inferM <- matrix(0, nsp+1, nsp+1);inferM[1,-1] <- 1
inferM.vfa <- matrix(0, NCOL(ada.pvfa), NCOL(ada.pvfa));inferM.vfa[1,-1] <- 1

ccm.rp <- ccm.fast.inferM(ds=ada.p,inferM,Lmin=5,Lint=10)
ccm.vfa <- ccm.fast.inferM(ds=ada.pvfa,inferM.vfa,Lmin=5,Lint=10)
ccm.rc <- ccm.fast.inferM(ds=ada.c,inferM,Lmin=5,Lint=10)

(spind.tp <- which(ccm.rp[[2]][1,]==1)-1)
(spind.tc <- which(ccm.rc[[2]][1,]==1)-1)
(spind.vfa <- which(ccm.vfa[[2]][1,]==1)-1)

ccm.rpw <- c(ccm.vfa[[1]][1,-1][spind.vfa],ccm.rp[[1]][1,-1][spind.tp])
ccm.rcw <- ccm.rc[[1]][1,-1][spind.tp]


(spind <- spind.tp)#intersect(spind.p1,spind.p))

ccm.crn <- data.frame(Variable=c(colnames(VFA),colnames(finger.p)[spind]),
                      CCM_Control=c(rep(NA,length(ccm.rpw)-length(ccm.rcw)),ccm.rcw), 
                      CCM_Perturbed=ccm.rpw)
if(SaveFile){
  write.csv(ccm.crn,paste('./Output/Demo_ccm_metabolite',type.ch4,'D47.csv',sep='_'),row.names=F)
  write.csv(finger.c[,spind],paste('./Data/OLR_Control_crn_verc',type.ch4,'D47.csv',sep='_'),row.names=F)
  write.csv(cbind(VFA[,spind.vfa],finger.p[,spind]),paste('./Data/OLR_Perturbed_crn_verc',type.ch4,'D47.csv',sep='_'),row.names=F)
}

# Combines microbiome and metabolome data 

crn.c <- read.csv(paste('./Data/OLR_Control_crn_verc',type.ch4,'D47.csv',sep='_'), header=T, stringsAsFactors = F)
crn.p <- read.csv(paste('./Data/OLR_Perturbed_crn_verc',type.ch4,'D47.csv',sep='_'), header=T, stringsAsFactors = F)

otu.c <- read.csv(paste('./Data/OLR_Control_otu_verc',type.ch4,'Species_D47.csv',sep='_'), header=T, stringsAsFactors = F)
otu.p <- read.csv(paste('./Data/OLR_Perturbed_otu_verc',type.ch4,'Species_D47.csv',sep='_'), header=T, stringsAsFactors = F)


data.p <- data.frame(CH4=envi[,c('CH4_O1')],crn.p, otu.p)


if(SaveFile){
  write.csv(data.frame(crn.c, otu.c),paste('./Data/OLR_Control_crnotu_verc',type.ch4,'Species_D47.csv',sep='_'),row.names=F)
  write.csv(data.frame(crn.p, otu.p),paste('./Data/OLR_Perturbed_crnotu_verc',type.ch4,'Species_D47.csv',sep='_'),row.names=F)
  write.csv(data.p,'./Data/Perturbed_CH4crnotu_verc_Species_D47.csv',row.names=F)
}
