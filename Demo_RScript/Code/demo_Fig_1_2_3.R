# R code for plotting Fig. 1-3
rm(list = ls());gc()
# Set path
OPath <- 'D:\\data\\Cooperation\\Prodigio\\'
setwd(OPath)

library(pheatmap)
library(factoextra)
library(igraph)
library(EnvStats)

seed <- 49563
set.seed(seed)
SaveFile <- T #T/F for saving files

tipp=37
tipp2=43
tind.o=1:46
nd=46

#####################################################################################################################################
#####################################################################################################################################
# Data loading
clasdo <- read.csv('./Demo_RScript/Data/demo_otu_taxa.csv',header=T)
crnotu1 <- read.csv('./Demo_RScript/Data/OLR1_otucrn_verc_Species_tau132.csv',header=T,stringsAsFactors = F)
crnotu3 <- read.csv('./Demo_RScript/Data/OLR3_otucrn_verc_Species_tau132.csv',header=T,stringsAsFactors = F)

# Calculate percentage data for microbial OTU 
otu1.s <- crnotu1[,grep('OTU_',colnames(crnotu1))]
otu3.s <- crnotu3[,grep('OTU_',colnames(crnotu3))]
otu1.s <- otu1.s*repmat(apply(otu1.s,1,sum),1,ncol(otu1.s))^-1
otu3.s <- otu3.s*repmat(apply(otu3.s,1,sum),1,ncol(otu3.s))^-1

# Extract metabolomic data
crn1.o <- crnotu1[,-grep('OTU_',colnames(crnotu1))]
crn3.o <- crnotu3[,-grep('OTU_',colnames(crnotu3))]

rownames(otu3.s)=rownames(crn3.o)=1:nd

otu.nam <- colnames(otu3.s)
nsp <- length(otu.nam)

nals=list(c('Chem','Bio'),
          c('VFA','Fingerprint','Methanogen','Acetogen','Other taxa')
)

################################################################################
#### Interactomic Data
# Chemical-Microbial interaction network
co1.ins=read.csv('./Demo_RScript/Data/OLR1_otucrn_verc_Species_nin46_cvunit0.025_enet_jcof_Nmvx_Rallx_demo_NEW.csv',stringsAsFactors = F)
co3.ins=read.csv('./Demo_RScript/Data/OLR3_otucrn_verc_Species_nin46_cvunit0.025_enet_jcof_Nmvx_Rallx_demo_NEW.csv',stringsAsFactors = F)

######################################################################
### CCM results
ccm.otu=read.csv('./Demo_RScript/Data/demo_ccm_OTU.csv', header=T)
ccm.crn=read.csv('./Demo_RScript/Data/demo_ccm_metabolite.csv', header=T)

# function for plot time series by clusters 
plot.den=function(da, out.o=NULL, memb=NULL, cut.k, mem.out=F, plot_p=T){
  if(is.null(memb)){
    memb=data.frame(variable=colnames(da[,out.o$tree_row[["order"]]]),member=rep(0))
    for(i in 1:cut.k) {
      memb.i=cutree(out.o$tree_row, k=cut.k)
      ind.i=(memb.i[names(memb.i)==memb[memb[,'member']==0,'variable'][1]])
      memb[memb[,'variable']%in%names(which(memb.i==ind.i)),'member']=i
    }
  }
  
  nd.t=NROW(da)
  if(plot_p){
    win.graph(cut.k*10*2,cut.k*10)
    par(mfrow=c(ifelse(cut.k<=4,1,2),min(cut.k,4)),mar=c(4,3,1,1))
  }
  clumean=NULL
  for(i in 1:cut.k){
    da.i=as.matrix(da[,memb[memb[,'member']==i,'variable']])
    if(!all(is.na(c(da.i)))){
      if(plot_p){
        plot(c(da.i)~rep(1:nd.t,NCOL(da.i)),type='n', xlab=ifelse(i==cut.k,'Time',''),ylab='values',main=paste('cluster',i))
        for(j in 1:NCOL(da.i)){
          lines(da.i[,j]~c(1:nd.t),type='l',col='grey80')
        }
        lines(apply(da.i,1,mean,na.rm=T)~c(1:nd.t), type='l', lwd=2)
        polygon(x=c(tipp,tipp2,tipp2,tipp),y=rep(mag_scale(c(-10,10)),each=2), col="orange",density=15, border="white")
      }
      clumean=rbind(clumean,apply(da.i,1,mean,na.rm=T))
    }
  }
  if(mem.out){return(list(memb,clumean))}
}

mag_scale=function(x,mgf=0.5){
  if(prod(x)>=0){
    if(x[2]>0){y=x;y[1]=y[1]*(1.0-mgf);y[2]*(1.0+mgf)}else{
      y=x;y[1]=y[1]*(1.0+mgf);y[2]*(1.0-mgf)}
  }else{y=x*(1.0+mgf)}
  return(y)
}

# Function for analyzing chemical formula 
chemfor=function(nam.o,element= c('C','H','N','O','S')){
  nam.i=nam.o
  for(i in 1:length(element)) nam.i=sub(element[i], paste('-',element[i],'-',sep=''), nam.i)
  
  namx=matrix(0, ncol=5,nrow=length(nam.i))
  for(i in 1:length(element)){
    ind.i=grep(element[i],nam.i)
    nub=unlist(lapply(strsplit(nam.i[ind.i],paste('-',element[i],'-',sep='')),function(x) {unlist(strsplit(x[2],'-'))[1]}))
    nub[nub==""]=1
    nub[is.na(nub)]=1
    nub=as.numeric(nub)
    namx[ind.i,i]=nub
  }
  
  colnames(namx)=element
  return(namx)
}

################################################################################
## (Fig. 1) Early warning signal analysis for CH4 production 
################################################################################
reactor_ch4_ews <- read.csv( './Demo_RScript/Data/reactor_ch4_ews.csv', header=T)

yrg1=range(c(reactor_ch4_ews[,c('TCOD_control','TCOD_perturbed')]),na.rm=T)
win.graph(160,50)
par(mfrow=c(1,2),mar=c(4,4,1,4))
plot(TCOD_control~Day, data=reactor_ch4_ews, type='o',ylab='TCOD (mg/L)',xlab='Day',ylim=yrg1,
     main='Control reactor',xlim=c(1,(nd+1)),lwd=2,pch=1, cex=1.5, cex.lab=1.5, cex.axis=1.5)
plot(TCOD_perturbed~Day, data=reactor_ch4_ews, type='o',ylab='TCOD (mg/L)',xlab='Day',ylim=yrg1,
     main='Perturbed reactor',xlim=c(1,(nd+1)),lwd=2,pch=1, cex=1.5, cex.lab=1.5, cex.axis=1.5)
abline(h=1,col='grey')
polygon(x=c(tipp,tipp2,tipp2,tipp),y=rep(yrg1*c(0.5,1.5),each=2), col="orange",density=15, border="white")

win.graph(160,70)
par(mfrow=c(1,2),mar=c(4,4,1,4))
yrg2=range(c(reactor_ch4_ews[,c('CH4_control','CH4_perturbed')]),na.rm=T)

plot(CH4_control~Day, data=reactor_ch4_ews,type='o',ylab='CH4 (%)',xlab='Day',ylim=yrg2,
     main='Control reactor',xlim=c(1,nd),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
par(new=T)
plot(DEV_CH4_control~I(Day+0.5), data=reactor_ch4_ews, type='o',col='blue',xlim=c(1,nd),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
mtext('',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='', col.axis = "blue", cex.axis=1.5)
#dev_sum2(gas12)
abline(h=1,col='grey')

plot(CH4_perturbed~Day, data=reactor_ch4_ews,type='o',ylab='',xlab='Day',ylim=yrg2,
     main='Perturbed reactor',xlim=c(1,nd),lwd=2,pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.5)
par(new=T)
plot(DEV_CH4_perturbed~I(Day+0.5), data=reactor_ch4_ews, type='o',col='blue',xlim=c(1,nd),ylab='',xlab='',axes=F,lwd=2,pch=16, cex=1.5)
mtext('DEV',side=4, col='blue',line = 3, cex=1.5)
axis(4,col='blue',ylab='DEV', col.axis = "blue", cex.axis=1.5)
#dev_sum2(gas12)
abline(h=1,col='grey')
polygon(x=c(tipp,tipp2,tipp2,tipp),y=c(-10,-10,10,10), col="orange",density=15, border="white")

graphics.off()

################################################################################
## Fig. 2A Heatmap for microbial OTUs
################################################################################
annotation_row=data.frame(
  Function = factor(clasdo[,'Functional_guild']),
  CCM=as.numeric(ccm.otu[,'Perturbed_OTU'])
)
rownames(annotation_row)=colnames(otu3.s)
rownames(annotation_row)=gsub('OTU_','',rownames(annotation_row))

ann_colors = list(
  Function = c(Methanogen = "#E7298A", Acetogen= "#7570B3", Other_taxa= "grey80"),
  CCM=c("white", "darkgreen")
)

ind=1:tipp2
da.to=otu3.s[ind,]
da.t=apply(da.to,2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
colnames(da.t)=colnames(da.to)=rownames(annotation_row)

# Assign cluster number  
cut.k=7
out=pheatmap(t(da.t),cutree_rows = cut.k, cluster_cols = F, 
             border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors, silent=T)
mem.otu=plot.den(da=da.t, out.o=out, cut.k=cut.k, mem.out=T, plot_p=F)

memb.i=cutree(out$tree_row, k=cut.k)
ord=out$tree_row$order
memb.o=memb.i[ord]

# Reordered clusters according to peak timings
reo=c(7,2,1,3:6)
memb.reo=NULL
for(i in 1:cut.k){
  mem.oi=filter(mem.otu[[1]],member==reo[i])
  mem.oi[,'member']=i
  memb.reo=rbind(memb.reo, mem.oi)
}
gap.o=c();for(i in 1:cut.k){gap.o[i]=max(which(memb.reo[,'member']==i))}

annotation_row.o=annotation_row
annotation_row=annotation_row[memb.reo[,'variable'],]
memb.o2=data.frame(cluster=memb.reo['member'], clasdo[paste('OTU_',memb.reo[,'variable'],sep=''),], annotation_row)

# Plot reordered cluster 
win.graph(120,140)
da.t3r=da.t[,memb.reo[,'variable']]
outr=pheatmap(t(da.t3r),cutree_rows = cut.k, cluster_cols = F, cluster_row = F, gaps_row=sort(gap.o),
              border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors)
mem.otu=plot.den(da=da.t, memb=memb.reo, cut.k=cut.k, mem.out=T)

# Check the assigned cluster number by elbow method
windows()
fviz_nbclust(t(da.t), 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 12          # max number of clusters to consider
) + 
  labs(title="Elbow Method for HC") +
  geom_vline(xintercept = cut.k, linetype = 2)

graphics.off()

################################################################################
## (Fig. 2B) Heatmap for chemical fingerprint
################################################################################

#Analysis of chemical formula 
elenum=chemfor(colnames(crn3.o[,-c(1:7)]),element= c('C','H','N','O','S'))
#C number: HAc: acetic acid 2; HPro: propionic acid 3; HBu: butyric acid 4; iHVal: iso-valeric acid 5; HVal: valeric acid 5; HCa: caproic acid 6;
fin.mvs=data.frame(Cn=c(c(2, 3, 4, 4, 5, 5, 6),elenum[,'C']),NCr=c(rep(0,7),elenum[,'N']/elenum[,'C']),SCr=c(rep(0,7),elenum[,'S']/elenum[,'C']))

res.crn=colnames(crn3.o)

ele.t=rep('CH',nrow(elenum));ele.t[elenum[,'N']!=0]='N';ele.t[elenum[,'S']!=0]='S';
annotation_row=data.frame(
  C_number=fin.mvs[,'Cn'],
  NC_ratio=fin.mvs[,'NCr'],
  SC_ratio=fin.mvs[,'SCr'],
  CCM=as.numeric(ccm.crn[,'Perturbed_finger'])
)
rownames(annotation_row)=colnames(crn3.o)

ann_colors = list(
  C_number = c("white", "black"),
  NC_ratio = c("white", "orange"),
  SC_ratio = c("white", "cyan"),
  CCM=c("white", "darkgreen")
)

# Assign cluster number  
cut.k=6
ind=1:tipp2
da.t=apply(crn3.o[ind,],2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
out=pheatmap(t(da.t),cutree_rows = cut.k, cluster_cols = F, 
             border_color = "grey100", fontsize_row=8, annotation_row=annotation_row, annotation_colors  = ann_colors, silent=T)
mem.crn=plot.den(da=da.t, out.o=out, cut.k=cut.k,mem.out=T, plot_p = F)

# Reordered clusters according to peak timings
reo=c(1:2,6,3:5)
memb.reo=NULL
for(i in 1:cut.k){
  mem.oi=filter(mem.crn[[1]],member==reo[i])
  mem.oi[,'member']=i
  memb.reo=rbind(memb.reo, mem.oi)
}
gap.o=c();for(i in 1:cut.k){gap.o[i]=max(which(memb.reo[,'member']==i))}

annotation_row.o=annotation_row
annotation_row=annotation_row[memb.reo[,'variable'],]

win.graph(150,140)
da.t3r=da.t[,memb.reo[,'variable']]
outr=pheatmap(t(da.t3r),cutree_rows = cut.k, cluster_cols = F, cluster_row = F, gaps_row=sort(gap.o),
              border_color = "grey100", fontsize_row=5, annotation_row=annotation_row, annotation_colors  = ann_colors )
mem.crn=plot.den(da=da.t, memb=memb.reo, cut.k=cut.k, mem.out=T)

# Check the assigned cluster number by elbow method
windows()
fviz_nbclust(t(da.t), 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 12          # max number of clusters to consider
) + 
  labs(title="Elbow Method for HC") +
  geom_vline(xintercept = cut.k, linetype = 2) 

graphics.off()

################################################################################
## Fig. 3A Heatmap for network properties
################################################################################
netprop.o=read.csv('./Demo_RScript/Data/demo_network_prop_perturbed.csv',header=T)
netprop.p=filter(netprop.o, Reactor=='Perturbed')
netprop=t(netprop.p[,-c(1:4)]);colnames(netprop)=netprop.p[,'Property'];rownames(netprop)=NULL

annotation_row=data.frame(
  Property = factor(netprop.p[,'Property_type']),
  Type = factor(netprop.p[,'Network_type'])
)
rownames(annotation_row)=netprop.p[,'Property']

ann_colors = list(
  Property = c(FG_Net= 'cyan2', Topology = "#7570B3", Centrality= "pink", IS_type= "#66A61E", Feedback = "brown"),
  Type = c(Chem="blue4",Bio='yellow2',Chem_Bio='darkgreen')
)

win.graph(150,140)
cut.k=7
ind=1:tipp2-1#
da.t=apply(netprop[ind,],2,scale)
da.tmax=da.t
da.tmax[da.t>=2.5]=2.5;da.tmax[da.t<=-2.5]=-2.5;

out=pheatmap(t(da.tmax),cutree_rows = cut.k, cluster_cols = F, 
             border_color = "grey100", fontsize_row=8, annotation_row=annotation_row, annotation_colors  = ann_colors )
mem.net=plot.den(da=da.t, out.o=out, cut.k=cut.k,mem.out=T)

windows()
fviz_nbclust(t(da.t), 
             FUNcluster = hcut,  # hierarchical clustering
             method = "wss",     # total within sum of square
             k.max = 12          # max number of clusters to consider
) + 
  labs(title="Elbow Method for HC") +
  geom_vline(xintercept = cut.k, linetype = 2)      # 在 X=3的地方 

mem.all=rbind(data.frame(type='OTU',cluster=1:NROW(mem.otu[[2]]),mem.otu[[2]]),
              data.frame(type='CRN',cluster=1:NROW(mem.crn[[2]]),mem.crn[[2]]),
              data.frame(type='NET',cluster=1:NROW(mem.net[[2]]),cbind(mem.net[[2]]),X43=rep(NA)))

if(SaveFile) write.csv(mem.all,'demo_cluster_all_verc.csv',row.names=F)

################################################################################
## Fig. 3B Time-varying interaction networks among Functional Guilds
################################################################################
# about 10-day averaged interaction strengths
periodn=cbind(c(1,11,21,31),c(10,20,30,37))
net.tmax=45
indx=c('VFA','Fingerprint','Methanogen','Acetogen','Other taxa')

fgnet=read.csv('./Demo_RScript/Data/demo_fgnet.csv',header=T)
fgnet.c=filter(fgnet,Reactor=='Control')
fgnet.p=filter(fgnet,Reactor=='Perturbed')

block.c=block.p=array(0, c(5,5,NROW(periodn)))
for(i in 1:NROW(periodn)){
  block.c[,,i]=matrix(as.numeric(apply(fgnet.c[,-c(1:2)][periodn[i,1]:periodn[i,2],],2,mean,na.rm=T)),5,5)
  block.p[,,i]=matrix(as.numeric(apply(fgnet.p[,-c(1:2)][periodn[i,1]:periodn[i,2],],2,mean,na.rm=T)),5,5)  
}

blocknet=list(block.c,block.p)
vcol=c('#DBB100', '#66A61E','#E7298A', "#7570B3", "grey80")#rainbow(length(indx))

win.graph(80,40)
par(mfrow=c(2,4),mar=c(1,1,1,1))
for(k in 1:length(blocknet)){
  for(i in 1:NROW(periodn)){
    g.a=t(blocknet[[k]][,,i])
    g.a[is.na(g.a)]=0
    diag(g.a)=0
    g <- graph_from_adjacency_matrix(g.a, mode = "directed", weighted = T)
    ga.p=ga.n=g.a
    ga.p[g.a<=0]=0;ga.n[g.a>=0]=0;
    ga.i=(ga.p>0)-(ga.n<0)
    ga.m=c(g.a[ga.i!=0])
    ga.i=c(ga.i[ga.i!=0])
    
    yrg.p=range(ga.p,na.rm=T)
    yrg.n=range(ga.n,na.rm=T)
    
    V(g)$color=vcol
    V(g)$label <- indx
    V(g)$label.cex <- 1
    
    E(g)$width=abs(E(g)$weight)*10^3#
    E(g)$color[E(g)$weight>0] <- "blue"
    E(g)$color[E(g)$weight<0] <- "red"
    E(g)$arrow.size=0.5
    
    #Color scaling function
    c_scale.p <- colorRamp(c('white','blue'))
    c_scale.n <- colorRamp(c('white','red'))
    
    coords <- layout_in_circle(g, order = 1:length(indx))
    
    plot(g, layout = coords, main=paste('Day', 1+(i-1)*10,'-', ifelse(i<4,i*10,tipp)), vertex.size=25, edge.curved = 0.2)
    if(i==1& k==1){legend('bottomright', indx,pch=rep(16,length(indx)), col=vcol, bty='n',cex=1)
      legend('bottomleft', c('Positive interaction','Negtive interaction'),
             lty=c(1,1), col=c('blue','red'), bty='n',cex=1)
    }
  }# end of i
}# end of k

################################################################################
### Fig. 3C-E Network instability indices & strength of feedback loops
################################################################################
netprop.co=filter(netprop.o, Reactor=='Control')
netprop.c=t(netprop.co[,-c(1:4)]);colnames(netprop.c)=netprop.co[,'Property'];rownames(netprop.c)=NULL

instab_property=c('VFA_D_eigenvalue','Methano_D_eigenvalue','Feedback  Bio L=3')
stable.p=netprop[1:(tipp2-1),instab_property]
stable.c=netprop.c[1:(tipp2-1),instab_property]

win.graph(90,30);
par(mfrow=c(1,3),mar=c(4,4,1,1))
for(i in 1:length(instab_property)){
  x.t=c((1:(tipp2-1))+0.5)
  y.p=stable.p[,i]
  y.c=stable.c[,i]
  if(i==1){x.t=x.t[-c(1:8)];y.p=y.p[-c(1:8)]} # DEV of VFA subnetwork was zero during Day 1-8 because of very low VFAs concentration. 
  plot(y.p~x.t,ylab=colnames(stable.p)[i],type='o',xlab='Time',pch=16,cex=2,lwd=2,ylim=range(c(y.p,y.c),na.rm=T))
  if(i!=1){
    lines(y.c~x.t,type='o',pch=1,cex=2,lwd=2,col='grey')
  } 
  abline(h=1,lty=2)
  polygon(x=c(tipp,tipp2,tipp2,tipp),y=rep(mag_scale(range(c(y.p,y.c),na.rm=T)),each=2), col="orange",density=15, border="white")
  legend('topleft',c('Perturbed','Control'),bty='n',lty=c(1,1),col=c('black','grey'),pch=c(16,1))
}


