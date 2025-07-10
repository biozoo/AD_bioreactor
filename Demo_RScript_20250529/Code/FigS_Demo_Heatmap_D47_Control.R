# R code for plotting Fig. 1-3
rm(list = ls());gc()
# Set path
OPath <- 'D:/data/Cooperation/Prodigio/Code/Demo_RScript_20250529'
setwd(OPath)

library(pheatmap)
library(factoextra)
library(igraph)
library(EnvStats)
library(dplyr)

seed <- 49563
set.seed(seed)
SaveFile <- T #T/F for saving files

tipp <- 38+1
tipp2 <- 44+1
nd <- 47
tind.o <- 0:(nd-1)


#####################################################################################################################################
#####################################################################################################################################
# Data loading
clasdo <- read.csv('./Data/demo_otu_taxa_D47.csv',header=T)
rownames(clasdo) <- clasdo[,'OTU']
crnotu1 <- read.csv('./Data/OLR_Control_crnotu_verc_CH4_Species_D47.csv',header=T,stringsAsFactors = F)
crnotu3 <- read.csv('./Data/OLR_Perturbed_crnotu_verc_CH4_Species_D47.csv',header=T,stringsAsFactors = F)
envi <- read.csv('.//Data//Reactor_data_envi.csv',header=T,stringsAsFactors = F)

# Calculate percentage data for microbial OTU 
otu1 <- crnotu1[,grep('OTU_',colnames(crnotu1))]
otu3 <- crnotu3[,grep('OTU_',colnames(crnotu3))]

# Extract metabolomic data
crn1 <- crnotu1[,-grep('OTU_',colnames(crnotu1))]
crn3 <- crnotu3[,-grep('OTU_',colnames(crnotu3))]

rownames(otu3) <- rownames(crn3) <- 1:nd

otu.nam <- colnames(otu3)
nsp <- length(otu.nam)

nals <- list(c('Chem','Bio'),
          c('VFA','Fingerprint','Methanogen','Acetogen','Other taxa')
)

######################################################################
### CCM results
ccm.otu <- read.csv('./Output/Demo_ccm_OTU_CH4_Species_D47.csv', header=T)
ccm.crn <- read.csv('./Output/Demo_ccm_metabolite_CH4_D47.csv', header=T)

# function for plot time series by clusters 
plot.den <- function(da, out.o=NULL, memb=NULL, cut.k, mem.out=F, plot_p=T, yrg=NULL){
  if(is.null(memb)){
    memb <- data.frame(variable=colnames(da[,out.o$tree_row[["order"]]]),member=rep(0))
    for(i in 1:cut.k) {
      memb.i <- cutree(out.o$tree_row, k=cut.k)
      ind.i <- (memb.i[names(memb.i)==memb[memb[,'member']==0,'variable'][1]])
      memb[memb[,'variable']%in%names(which(memb.i==ind.i)),'member'] <- i
    }
  }
  if(is.null(yrg)){da.yrg <- range(da,na.rm=T)}else{da.yrg <- yrg}
  nd.t <- NROW(da)
  if(plot_p){
    #win.graph(cut.k*10*2,cut.k*10)
    #par(mfrow=c(ifelse(cut.k<=4,1,2),min(cut.k,4)),mar=c(4,3,1,1))
    win.graph(15,cut.k*10)
    par(mfrow=c(cut.k,1),mar=c(3,3,1,1))
  }
  clumean <- NULL
  for(i in 1:cut.k){
    da.i <- as.matrix(da[,memb[memb[,'member']==i,'variable']])
    if(!all(is.na(c(da.i)))){
      if(plot_p){
        plot(c(da.i)~rep(1:nd.t,NCOL(da.i)),type='n', xaxt = ifelse(i==cut.k,"s","n"),
             xlab=ifelse(i==cut.k,'Time',''),ylim=da.yrg,
             ylab='values',main=paste('cluster',i))
        for(j in 1:NCOL(da.i)){
          lines(da.i[,j]~c(1:nd.t),type='l',col='grey80')
        }
        lines(apply(da.i,1,mean,na.rm=T)~c(1:nd.t), type='l', lwd=2)
        polygon(x=c(tipp,tipp2,tipp2,tipp),y=rep(mag_scale(c(-10,10)),each=2), col="orange",density=15, border="white")
      }
      clumean <- rbind(clumean,apply(da.i,1,mean,na.rm=T))
    }
  }
  if(mem.out){return(list(memb,clumean))}
}

mag_scale <- function(x,mgf=0.5){
  if(prod(x)>=0){
    if(x[2]>0){y <- x;y[1] <- y[1]*(1.0-mgf);y[2]*(1.0+mgf)}else{
      y <- x;y[1] <- y[1]*(1.0+mgf);y[2]*(1.0-mgf)}
  }else{y <- x*(1.0+mgf)}
  return(y)
}

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


###############################################################################
## Fig S time series of environmental parameters
win.graph(200,100)
par(mfcol=c(2,3),mar=c(4,4,1,4))
# Control reactors
name.a <- c('pH','ORP','Temp','NH4','Qbiogas','CH4p')
reac <- c('C1','C2','O1','O2')
xlab.a <- c('pH','ORP','Temperature (oC)','NH4 (mg N/L)','Biogas production (ml/Day)','Methane production (%)')

for(i in 1:length(name.a)){
  (name.t <- paste(name.a[i],reac,sep='_'))
  yrg2 <- range(c(envi[,name.t]),na.rm=T)
  
  plot(envi[,name.t[2]]~c((1:nd)-1),type='o',ylab=xlab.a[i],xlab='Day',ylim=yrg2, col='grey',
       main='',lwd=2,pch=17, cex=2, cex.lab=1.5, cex.axis=1.5)
  lines(envi[,name.t[1]]~c((1:nd)-1),type='o',col='black',
        main='',lwd=2,pch=16, cex=2, cex.lab=1.5, cex.axis=1.5)
  lines(envi[,name.t[4]]~c((1:nd)-1),type='o',col='lightgreen',
        main='',lwd=2,pch=17, cex=2, cex.lab=1.5, cex.axis=1.5)
  lines(envi[,name.t[3]]~c((1:nd)-1),type='o',col='darkgreen',
        main='',lwd=2,pch=16, cex=2, cex.lab=1.5, cex.axis=1.5)
  
  abline(v=tipp-1,col='orange',lwd=2)
  polygon(x=c(tipp,tipp2,tipp2,tipp)-1,y=rep(c(0.1,2)*yrg2,each=2), col="orange",density=15, border="white")
  if(i==3) { legend('topleft',c('Control C1','Control C2','Perturbed O1','Perturbed O2'),cex=2,
                    pch=c(16,17,16,17),col=c('black','grey','darkgreen','lightgreen'),lty=1,lwd=2,bty='n')}
}


graphics.off()

################################################################################
## Fig. 2A Heatmap for microbial OTUs
################################################################################
annotation_row <- data.frame(
  Function = factor(clasdo[,'Functional_guild']),
  CCM=as.numeric(ccm.otu[,'CCM_Control'])
)
rownames(annotation_row) <- colnames(otu3)

ann_colors  <-  list(
  Function = c(Methanogen = "#E7298A", Acetogen= "#7570B3", Other_taxa= "grey80"),
  CCM=c("white", "darkgreen")
)

ind <- (1:tipp2)
da.to <- otu3[ind,]
da.t <- apply(da.to,2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
colnames(da.t) <- colnames(da.to) <- rownames(annotation_row)

# Assign cluster number  
cut.k <- 7
out <- pheatmap(t(da.t),cutree_rows = cut.k, cluster_cols = F, clustering_method = "complete",
                border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors, silent=T)
mem.otu <- plot.den(da=da.t, out.o=out, cut.k=cut.k, mem.out=T, plot_p=F)

memb.i <- cutree(out$tree_row, k=cut.k)
ord <- out$tree_row$order
memb.o <- memb.i[ord]

# Reordered clusters according to peak timings
reo <- order(apply(mem.otu[[2]][,1:tipp],1,which.max))
memb.reo <- NULL
for(i in 1:cut.k){
  mem.oi <- filter(mem.otu[[1]],member==reo[i])
  mem.oi[,'member'] <- i
  memb.reo <- rbind(memb.reo, mem.oi)
}
gap.o <- c();for(i in 1:cut.k){gap.o[i] <- max(which(memb.reo[,'member']==i))}

annotation_row.o <- annotation_row
annotation_row <- annotation_row[memb.reo[,'variable'],]
memb.o2 <- data.frame(cluster=memb.reo['member'], clasdo[paste('OTU_',memb.reo[,'variable'],sep=''),], annotation_row)

# Plot reordered cluster in control reactor
win.graph(120,140)
da.t <- apply(otu1[ind,],2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
da.t3r <- da.t[,memb.reo[,'variable']]
outr <- pheatmap(t(da.t3r),cutree_rows = cut.k, cluster_cols = F, cluster_row = F, gaps_row=sort(gap.o),clustering_method =  "complete",
                 border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors)
mem.otu <- plot.den(da=da.t, memb=memb.reo, cut.k=cut.k, mem.out=T)

### Plot cluster based on the grouping of control reactor
annotation_row <- data.frame(
  Function = factor(clasdo[,'Functional_guild']),
  CCM=as.numeric(ccm.otu[,'CCM_Control'])
)
rownames(annotation_row) <- colnames(otu3)

ann_colors  <-  list(
  Function = c(Methanogen = "#E7298A", Acetogen= "#7570B3", Other_taxa= "grey80"),
  CCM=c("white", "darkgreen")
)

ind <- 1:tipp2
da.to <- otu1[ind,]
da.t <- apply(da.to,2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
colnames(da.t) <- colnames(da.to) <- rownames(annotation_row)

# Assign cluster number  
cut.k <- 7
out <- pheatmap(t(da.t),cutree_rows = cut.k, cluster_cols = F, clustering_method = "complete",
                border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors, silent=T)
mem.otu <- plot.den(da=da.t, out.o=out, cut.k=cut.k, mem.out=T, plot_p=F)

memb.i <- cutree(out$tree_row, k=cut.k)
ord <- out$tree_row$order
memb.o <- memb.i[ord]

# Reordered clusters according to peak timings
reo <- order(apply(mem.otu[[2]][,1:tipp],1,which.max))
memb.reo <- NULL
for(i in 1:cut.k){
  mem.oi <- filter(mem.otu[[1]],member==reo[i])
  mem.oi[,'member'] <- i
  memb.reo <- rbind(memb.reo, mem.oi)
}
gap.o <- c();for(i in 1:cut.k){gap.o[i] <- max(which(memb.reo[,'member']==i))}

annotation_row.o <- annotation_row
annotation_row <- annotation_row[memb.reo[,'variable'],]
memb.o2 <- data.frame(cluster=memb.reo['member'], clasdo[paste('OTU_',memb.reo[,'variable'],sep=''),], annotation_row)

# Plot reordered cluster 
win.graph(120,140)
da.t <- apply(otu1[ind,],2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
da.t3r <- da.t[,memb.reo[,'variable']]
outr <- pheatmap(t(da.t3r),cutree_rows = cut.k, cluster_cols = F, cluster_row = F, gaps_row=sort(gap.o),clustering_method =  "complete",
                 border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors)
mem.otu <- plot.den(da=da.t, memb=memb.reo, cut.k=cut.k, mem.out=T)

graphics.off()

################################################################################
## (Fig. 2B) Heatmap for chemical fingerprint
################################################################################
#Analysis of chemical formula 
elenum <- chemfor(colnames(crn3[,-c(1:7)]),element= c('C','H','N','O','S'))
#C number: HAc: acetic acid 2; HPro: propionic acid 3; HBu: butyric acid 4; iHVal: iso-valeric acid 5; HVal: valeric acid 5; HCa: caproic acid 6;
fin.mvs <- data.frame(Cn=c(c(2, 3, 4, 4, 5, 5, 6),elenum[,'C']),NCr=c(rep(0,7),elenum[,'N']/elenum[,'C']),SCr=c(rep(0,7),elenum[,'S']/elenum[,'C']))

res.crn <- colnames(crn3)

ele.t <- rep('CH',nrow(elenum));ele.t[elenum[,'N']!=0] <- 'N';ele.t[elenum[,'S']!=0] <- 'S';
annotation_row <- data.frame(
  C_number=fin.mvs[,'Cn'],
  NC_ratio=fin.mvs[,'NCr'],
  SC_ratio=fin.mvs[,'SCr'],
  CCM=as.numeric(ccm.crn[,'CCM_Control'])
)
rownames(annotation_row) <- colnames(crn3)

ann_colors  <-  list(
  C_number = c("white", "black"),
  NC_ratio = c("white", "orange"),
  SC_ratio = c("white", "cyan"),
  CCM=c("white", "darkgreen")
)

# Assign cluster number  
cut.k <- 6
ind <- 1:(tipp2)

da.t <- apply(crn3[ind,],2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
out <- pheatmap(t(da.t),cutree_rows = cut.k, cluster_cols = F, clustering_method =  "complete",
                border_color = "grey100", fontsize_row=8, annotation_row=annotation_row, annotation_colors  = ann_colors, silent=T)
mem.crn <- plot.den(da=da.t, out.o=out, cut.k=cut.k,mem.out=T, plot_p = F)
rownames(da.t) <- (1:tipp2)-1

# Reordered clusters according to peak timings
reo <- order(apply(mem.crn[[2]][,1:tipp],1,which.max))
memb.reo <- NULL
for(i in 1:cut.k){
  mem.oi <- filter(mem.crn[[1]],member==reo[i])
  mem.oi[,'member'] <- i
  memb.reo <- rbind(memb.reo, mem.oi)
}
gap.o <- c();for(i in 1:cut.k){gap.o[i] <- max(which(memb.reo[,'member']==i))}

annotation_row.o <- annotation_row
annotation_row <- annotation_row[memb.reo[,'variable'],]

win.graph(150,140)
da.t <- apply(cbind(matrix(NA,nd,7),crn1)[ind,],2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
colnames(da.t)[1:7] <- colnames(crn3)[1:7]
da.t3r <- da.t[,memb.reo[,'variable']]
outr <- pheatmap(t(da.t3r),cutree_rows = cut.k, cluster_cols = F, cluster_row = F, gaps_row=sort(gap.o),clustering_method =  "complete",
              border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors )
mem.crn <- plot.den(da=da.t, memb=memb.reo, cut.k=cut.k, mem.out=T)

### Plot cluster based on the grouping of control reactor
fin.mvs <- fin.mvs[-c(1:7),]
annotation_row <- data.frame(
  C_number=fin.mvs[,'Cn'],
  NC_ratio=fin.mvs[,'NCr'],
  SC_ratio=fin.mvs[,'SCr'],
  CCM=as.numeric(ccm.crn[-c(1:7),'CCM_Control'])
)
rownames(annotation_row) <- colnames(crn1)

ann_colors  <-  list(
  C_number = c("white", "black"),
  NC_ratio = c("white", "orange"),
  SC_ratio = c("white", "cyan"),
  CCM=c("white", "darkgreen")
)

# Assign cluster number  
cut.k <- 6
ind <- (1:tipp2)
da.t <- apply(crn1[ind,],2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
out <- pheatmap(t(da.t),cutree_rows = cut.k, cluster_cols = F, clustering_method =  "complete",
                border_color = "grey100", fontsize_row=8, silent=T)
mem.crn <- plot.den(da=da.t, out.o=out, cut.k=cut.k,mem.out=T, plot_p = F)

# Reordered clusters according to peak timings
reo <- order(apply(mem.crn[[2]][,1:tipp],1,which.max))
memb.reo <- NULL
for(i in 1:cut.k){
  mem.oi <- filter(mem.crn[[1]],member==reo[i])
  mem.oi[,'member'] <- i
  memb.reo <- rbind(memb.reo, mem.oi)
}

memb.reo <- filter(memb.reo,variable%in%colnames(crn1))

gap.o <- c();for(i in 1:cut.k){gap.o[i] <- max(which(memb.reo[,'member']==i))}

annotation_row.o <- annotation_row
annotation_row <- annotation_row[memb.reo[,'variable'],]

win.graph(150,140)
da.t3r <- da.t[,memb.reo[,'variable']]
outr <- pheatmap(t(da.t3r),cutree_rows = cut.k, cluster_cols = F, cluster_row = F, gaps_row=sort(gap.o),clustering_method =  "complete",
                 border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors )
mem.crn <- plot.den(da=da.t, memb=memb.reo, cut.k=cut.k, mem.out=T)

graphics.off()

##################################################################################
## Fig. 3A Heatmap for network properties
##################################################################################
netprop.o <- read.csv('./Output/demo_network_prop_ccm_perturbed_D47.csv',header=T)
netprop.p <- filter(netprop.o, Reactor=='Perturbed')
netprop.ci <- filter(netprop.o, Reactor=='Control')
netprop <- t(netprop.p[,-c(1:4)]);colnames(netprop)=netprop.p[,'Property'];rownames(netprop)=NULL
netprop.c <- t(netprop.ci[,-c(1:4)]);colnames(netprop.c) <- netprop.ci[,'Property']; rownames(netprop.c) <- NULL

annotation_row <- data.frame(
  Property = factor(netprop.p[,'Property_type']),
  Type = factor(netprop.p[,'Network_type'])
)
rownames(annotation_row) <- netprop.p[,'Property']

ann_colors  <-  list(
  Property = c(Guild_Int= 'cyan2', Topology = "#7570B3", Centrality= "pink", Int_type= "#66A61E", Feedback = "brown"),
  Type = c(Chem="blue4",Bio='yellow2',Chem_Bio='darkgreen')
)

cut.k <- 6
ind <- 1:(tipp2-1)
da.t <- apply(netprop[ind,],2,scale)
rownames(da.t) <- (1:(tipp2-1))-1

da.tmax <- da.t
da.tmax[da.t>=2.0]=2.0;da.tmax[da.t<=-2.0]=-2.0;

out <- pheatmap(t(da.tmax),cutree_rows = cut.k, cluster_cols = F,clustering_method =  "complete", 
             border_color = "grey100", fontsize_row=8, annotation_row=annotation_row, annotation_colors  = ann_colors , silent=T)
mem.net <- plot.den(da=da.t, out.o=out, cut.k=cut.k, mem.out=T, yrg=c(-3,3), plot_p = F)

# Reordered clusters according to peak timings
reo <- order(apply(mem.net[[2]][,1:tipp],1,which.max))
memb.reo <- NULL
for(i in 1:cut.k){
  mem.oi <- filter(mem.net[[1]],member==reo[i])
  mem.oi[,'member'] <- i
  memb.reo <- rbind(memb.reo, mem.oi)
}
gap.o <- c();for(i in 1:cut.k){gap.o[i] <- max(which(memb.reo[,'member']==i))}

annotation_row.o <- annotation_row
annotation_row <- annotation_row[memb.reo[,'variable'],]

win.graph(120,140)
da.t <- apply(netprop.c[ind,],2,scale)
da.tmax <- da.t
da.tmax[da.t>=2.0] <- 2.0; da.tmax[da.t<=-2.0] <- -2.0;

da.t3r <- da.tmax[,memb.reo[,'variable']]
outr <- pheatmap(t(da.t3r),cutree_rows = cut.k, cluster_cols = F, cluster_row = F, gaps_row=sort(gap.o),clustering_method =  "complete",
                 border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors )
mem.net <- plot.den(da=da.tmax, memb=memb.reo, cut.k=cut.k, mem.out=T)

### Plot cluster based on the grouping of control reactor
annotation_row <- data.frame(
  Property = factor(netprop.ci[, 'Property_type']),
  Type = factor(netprop.ci[, 'Network_type'])
)
rownames(annotation_row) <- netprop.ci[,'Property']

ann_colors  <-  list(
  Property = c(Guild_Int= 'cyan2', Topology = "#7570B3", Centrality= "pink", Int_type= "#66A61E", Feedback = "brown"),
  Type = c(Chem="blue4",Bio='yellow2',Chem_Bio='darkgreen')
)

da.t <- apply(netprop.c[ind,!apply(is.na(netprop.c),2,all)],2,scale)
da.tmax <- da.t
da.tmax[da.t>=2.0] <- 2.0; da.tmax[da.t<=-2.0] <- -2.0;

out <- pheatmap(t(da.tmax),cutree_rows = cut.k, cluster_cols = F,clustering_method =  "complete", 
                border_color = "grey100", fontsize_row=8, annotation_row=annotation_row, annotation_colors  = ann_colors , silent=T)
mem.net <- plot.den(da=da.t, out.o=out, cut.k=cut.k, mem.out=T, yrg=c(-3,3), plot_p = F)

# Reordered clusters according to peak timings
reo <- order(apply(mem.net[[2]][,(1:tipp)],1,which.max))
memb.reo <- NULL
for(i in 1:cut.k){
  mem.oi <- filter(mem.net[[1]],member==reo[i])
  mem.oi[,'member'] <- i
  memb.reo <- rbind(memb.reo, mem.oi)
}
gap.o <- c();for(i in 1:cut.k){gap.o[i] <- max(which(memb.reo[,'member']==i))}

annotation_row.o <- annotation_row
annotation_row <- annotation_row[memb.reo[,'variable'],]

win.graph(120,140)
da.t3r <- da.tmax[,memb.reo[,'variable']]
outr <- pheatmap(t(da.t3r),cutree_rows = cut.k, cluster_cols = F, cluster_row = F, gaps_row=sort(gap.o),clustering_method =  "complete",
                 border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors )
mem.net <- plot.den(da=da.tmax, memb=memb.reo, cut.k=cut.k, mem.out=T)

