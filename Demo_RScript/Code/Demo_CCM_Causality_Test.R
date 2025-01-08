# R code for performing MDR S-map anaysis
rm(list = ls());gc()
# Set path
OPath <- 'D:/data/Cooperation/Prodigio/Demo_RScript'
setwd(OPath)

library(dplyr)
library(rEDM)
library(Kendall)

source('./Code/Function_MDR_inferM.R')

seed <- 49563
set.seed(seed)
SaveFile <- F #T/F for saving files
run.ccm <- F

reactor_ch4_ews <- read.csv( './Data/reactor_ch4_ews.csv', header=T)

# CCM causal analysis for microbial OTU

otu.da <- read.csv('./Data/OLR_otu_verc_ch4p_Species_tau132.csv',header=T)
otu.p <- filter(otu.da, Reactor=='Perturbed')[,-c(1:2)]
otu.c <- filter(otu.da, Reactor=='Control')[,-c(1:2)]

nsp <- NCOL(otu.p)
inferM <- matrix(0, nsp+1, nsp+1);
inferM[1,-1] <- 1

if(run.ccm){
  ccm.otu.c <- ccm.fast.inferM(ds=apply(cbind(reactor_ch4_ews[,'CH4_control'],otu.c),2,scale),inferM,Lmin=5,Lint=10)
  ccm.otu.p <- ccm.fast.inferM(ds=apply(cbind(reactor_ch4_ews[,'CH4_perturbed'],otu.p),2,scale),inferM,Lmin=5,Lint=10)
  ccm.otu <- data.frame(ID=colnames(otu.p), ccm_OTU_control=ccm.otu.c[[1]][1,-1],ccm_OTU_perturbed=ccm.otu.p[[1]][1,-1])
}

# CCM causal analysis for metabolite
crn.da <- read.csv('./Data/OLR_crn_verc_ch4p_Species_tau132.csv',header=T)
crn.p <- filter(crn.da, Reactor=='Perturbed')[,-c(1:2)]
crn.c <- filter(crn.da, Reactor=='Control')[,-c(1:2)][,-c(1:7)]

# Compiled metabolomic & microbiom data
crnotu.c <- cbind(crn.c, otu.c)
crnotu.p <- cbind(crn.p, otu.p)
crnotu.p.ch4 <- cbind(CH4=reactor_ch4_ews[,'CH4_perturbed'], crn.p, otu.p)
crnotu.c.ch4 <- cbind(CH4=reactor_ch4_ews[,'CH4_control'], crn.c, otu.c)
crnotu.c[,grep('OTU_',colnames(crnotu.c))] <- crnotu.c[,grep('OTU_',colnames(crnotu.c))]*repmat(apply(crnotu.c[,grep('OTU_',colnames(crnotu.c))],1,sum),1,ncol(crnotu.c[,grep('OTU_',colnames(crnotu.c))]))^-1
crnotu.p[,grep('OTU_',colnames(crnotu.p))] <- crnotu.p[,grep('OTU_',colnames(crnotu.p))]*repmat(apply(crnotu.p[,grep('OTU_',colnames(crnotu.p))],1,sum),1,ncol(crnotu.p[,grep('OTU_',colnames(crnotu.p))]))^-1

nsp <- NCOL(crn.p)
inferM <- matrix(0, nsp+1, nsp+1);
inferM[1,-1] <- 1
inferM.c <- inferM[,-c(2:8)]

if(run.ccm){
  ccm.crn.c <- ccm.fast.inferM(ds=apply(cbind(reactor_ch4_ews[,'CH4_control'],crn.c),2,scale),inferM.c,Lmin=5,Lint=10)
  ccm.crn.p <- ccm.fast.inferM(ds=apply(cbind(reactor_ch4_ews[,'CH4_perturbed'],crn.p),2,scale),inferM,Lmin=5,Lint=10)
  ccm.crn <- data.frame(ID=colnames(crn.p), ccm_metabolite_control=c(rep(NA,7),ccm.crn.c[[1]][1,-1]),ccm_metabolite_perturbed=ccm.crn.p[[1]][1,-1])
}

# CCM causal analysis for network property
net.o <- filter(read.csv('./Output/demo_network_prop_perturbed_new.csv',header=T),Reactor=='Perturbed')
net.p <- t(net.o[,-c(1:4)])
colnames(net.p) <- net.o[,'Property']
  
nsp <- NCOL(net.p)
inferM <- matrix(0, nsp+1, nsp+1);
inferM[1,-1] <- 1

if(run.ccm){ccm.net.p <- ccm.fast.inferM(ds=apply(cbind(reactor_ch4_ews[-1,'CH4_perturbed'],net.p),2,scale),inferM,Lmin=5,Lint=10)}

if(SaveFile){
  write.csv(ccm.otu,'./Output/demo_ccm_otu.csv',row.names = F)
  write.csv(ccm.crn,'./Output/demo_ccm_metabolite.csv',row.names = F)
  write.csv(ccm.net.p[[1]][1,-1],'./Output/demo_ccm_Network.csv',row.names = F)
  
  write.csv(crnotu.c, './Data/Demo_Control_otucrn_abundance.csv', row.names = F)
  write.csv(crnotu.p, './Data/Demo_Perturbed_otucrn_abundance.csv', row.names = F)
  write.csv(crnotu.p.ch4, './Data/Perturbed_otucrn_ch4_verc_Species.csv', row.names = F)
  write.csv(crnotu.c.ch4, './Data/Control_otucrn_ch4_verc_Species.csv', row.names = F)
}

