# R code for plotting Fig. 5
rm(list = ls());gc()
# Set path
OPath <- 'D:\\data\\Cooperation\\Prodigio\\'
setwd(OPath)

library(factoextra)
library(igraph)
library(mgcv)
library(gratia)
library(dplyr)

seed <- 49563
set.seed(seed)
SaveFile <- T #T/F for saving files

tipp <- 37
tipp2 <- 43
nd <- 46
tind.o <- 1:nd

################################################################################
## Functions
## R equivalent of repmat (matlab)
repmat = function(X,m=1,n=1){
  if(!is.null(dim(X))){X=as.matrix(X)}else{X=as.numeric(X)}
  if(is.matrix(X)){
    mx = dim(X)[1]
    nx = dim(X)[2]
    y=matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }else if(m>1&n==1) {
    y=matrix(rep(X,each=m),m,length(X))
    
  }else if(m==1&n>1){
    y=matrix(rep(X,n),length(X),n)
  }
  return(y)
}

mytriangle <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1 / 200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(
    x = coords[, 1], y = coords[, 2], bg = vertex.color,
    stars = cbind(vertex.size, vertex.size, vertex.size),
    add = TRUE, inches = FALSE
  )
}
# clips as a circle
add_shape("triangle",
          clip = shapes("circle")$clip,
          plot = mytriangle
)

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

rownames(otu3) <- rownames(crn3) <- 1:nd

nsp <- NCOL(otu3)
ncrn <- NCOL(crn1)

nco1 <- ncrn+nsp
nco3 <- ncrn+nsp+7

crnotu3s <- apply(crnotu3,2,scale)

ind.co=data.frame(Type=c(rep('VFA',7),rep('Fingerprint',ncrn), clasdo[,'Functional_guild']),
                  Name=colnames(crnotu3),Abb=c(rep(1,7),rep(2,ncrn),
                                               match(clasdo[,'Functional_guild'],c("Methanogen","Acetogen","Other_taxa"))+2))
rownames(ind.co)=ind.co[,'Name']

# Loading interactomic data
mco3.ins <- read.csv('./Demo_RScript/Data/CH4_crnotu_verc_OLR3_Species_nin46_cvunit0.025_enet_jcof_Nmvx_Rallx_demo_NEW.csv',stringsAsFactors = F)
co3.ins <- read.csv('./Demo_RScript/Data/OLR3_otucrn_verc_Species_nin46_cvunit0.025_enet_jcof_Nmvx_Rallx_demo_NEW.csv',stringsAsFactors = F)

co.ary3 <- array(NA,dim=c(nco3,nco3,nd-1))
for(i in 1:(nd-1)){
  ints.i31=as.matrix(filter(co3.ins,time==i)[,-c(1:4)])
  co.ary3[,,i]=ints.i31
}

################################################################################
### CCM results
ccm.otu <- read.csv('./Demo_RScript/Data/demo_ccm_OTU.csv', header=T)
ccm.crn <- read.csv('./Demo_RScript/Data/demo_ccm_metabolite.csv', header=T)
ccm.crnotu <- c(ccm.crn[,'Perturbed_finger'],ccm.otu[,'Perturbed_OTU'])
names(ccm.crnotu) <- colnames(crnotu3)
################################################################################
# Contributions for each species to the CH4 critical transition at the tipping point
da.i=data.frame(value=reactor_ch4_ews[,'CH4_perturbed'],time=1:nd)
g.t <- gam(value~s(time),data=da.i)
newd.i=seq(0,nd,0.01)
yhat=predict.gam(g.t,newdata=data.frame(time=newd.i))
yhatd1 <- derivatives(g.t, type = "central",data=data.frame(time=newd.i), order=1)$derivative

# Contribution of each node on the state changes in CH4 based on the dynamical interactome and time series of state variables 
impc=mco3.ins[,-c(1:5)]*crnotu3s[-nd,]
cotr=rbind(rep(NA),impc[-1,]-impc[-(nd-1),])
# Sign of the local trend (t, t+1) defined as the estimated first derivative evaluated at t
sdch4s=sign(yhatd1[newd.i%in%(seq(1,nd-1,1))])
# Examine whether the sign of local trend equals to the sign of contribution  
drcotrs=sign(cotr)*sdch4s
# Corrected contributions depending on whether the sign of contribution facilitate(+) or buffer (-) the observed trend in CH4
cotrds=drcotrs*abs(cotr)

# Cumulative contribution
ind=(tipp):(tipp2-1)#(tipp+1)#
cotrdr0=cotrds[ind,];
ind.r=which(as.numeric(cotrdr0[1,])>0)
cotrdr0=cotrdr0[,ind.r];rownames(cotrdr0)=paste('D',ind,sep='-')

# Examine the contributions for all the state changes among network nodes. 
# This is a similar process repeating the computation presented for CH4 

# smoothed trend estimate for all the variables 
smoin.crnotu3=filter(read.csv('smo_otucrn3_37.csv', header = T), time%in%(1:nd))

crnotu3.smo=matrix(smoin.crnotu3[,'yhat'],nrow=tipp)
crnotu3.smod=matrix(smoin.crnotu3[,'est'],nrow=tipp)

cotrd.ary=cotrds.ary=array(NA, dim=c(NCOL(crnotu3),NCOL(crnotu3),tipp))
for(i in 1:NCOL(crnotu3)){
  impc.i=t(co.ary3[i,,])*crnotu3s[-nd,]
  cotr.i=rbind(rep(NA),impc.i[-1,]-impc.i[-(nd-1),])
  cotr.i=cotr.i[1:tipp,]
  
  sdxis=sign(crnotu3.smod[,i])
  drcotrs.i=sign(cotr.i)*sdxis
  cotrds.i=drcotrs.i*abs(cotr.i)
  cotrds.ary[i,,]=t(cotrds.i)
}

mv.tmp <- read.csv('./Demo_RScript/Output/Localized_Transition_otucrn3.csv',header = T)
warn.crnotu <- read.csv('./Demo_RScript/Output/warn_time_CRNOTU.csv',header = T)
rownames(mv.tmp)=mv.tmp[,'ID']

domino=NULL
for(i in 1:NCOL(crnotu3)){
  mv.tmp.t0=filter(mv.tmp, ID==colnames(crnotu3)[i])
  cht.t=as.numeric(mv.tmp.t0[c('warn_time')])
  tarid=as.character(mv.tmp.t0['ID'])
  
  if(mv.tmp.t0['Ttype']%in%c('C1','C2','B1','B2')){
    ind.r=as.numeric(ceiling(cht.t)):as.numeric(1+ceiling(cht.t)); if(length(ind.r)==1){ind.r=ind.r:(ind.r+1)}
    indi=which(apply(cotrds.ary[i,,ind.r]>0,1,all))
    mv.tmp.t2=mv.tmp[indi,]
    mv.tmp.t2=mv.tmp.t2%>%mutate(contribution=apply(matrix(cotrds.ary[i,indi,ind.r],ncol=2),1,sum))%>%
      filter(Ttype%in%c('B1','B2','C1','C2'))%>%arrange(desc(warn_time))%>%filter((warn_time)<trunc(cht.t))
    tarid2=mv.tmp.t2[,'ID']
    contri=mv.tmp.t2[,'contribution']
  }else{tarid2=NA;contri=NA}
  if(length(tarid2)==0){tarid2=NA;contri=NA}
  domino=rbind(domino, data.frame(Variable=rep(i),ID=rep(tarid),Ttype=rep(mv.tmp.t0['Ttype']),Time=rep(cht.t),Domino=tarid2,Contribution=contri))
}

domino=filter(domino, ID!=Domino)
domino=filter(domino, ID%in%warn.crnotu[,1] & Domino%in%warn.crnotu[,1])
dim(domino)

ch4.contr=mv.tmp%>%filter(ID%in%colnames(cotrdr0),ID%in%warn.crnotu[,1])%>%arrange(desc(warn_time))
domino=rbind(data.frame(Variable=0,ID='CH4',Ttype='',Time=tipp,Domino=ch4.contr[,'ID'],
                        Contribution=as.numeric(cotrdr0[1,ch4.contr[,'ID']])),domino)

################################################################################
## Reconstruct various paths involved in the cascading mechanism
dp.ls <- list();dp.ls[[1]]=cbind(ind=rep('CH4'), ID=ch4.contr[,'ID'])
for(i in 2:30){
  dpath.i <- dp.ls[[i-1]][,2]
  ud.i <- unique(dpath.i)
  ud.i <- ud.i[!is.na(ud.i)]
  dj2 <- NULL
  for(j in 1:length(ud.i)){
    dj <- filter(domino,ID==ud.i[j])[,'Domino']
    if(length(dj)!=0 & all(!is.na(dj))){
      dj2=rbind(dj2, cbind(ind=rep(ud.i[j],each=length(dj)),ID=rep(dj,length(ud.i[j]))))
    }else{
      dj2=rbind(dj2, cbind(ind=rep(ud.i[j],each=1),ID=rep(NA,1)))
    }
  }
  dp.ls[[i]]=dj2;print(paste(i,NROW(dj2)))
  if(all(is.na(dj2[,2]))){break}else{rm(dj2);gc()}
}

# All individual cascading paths 
dpi1 <- dp.ls[[1]]
for(i in 2:(length(dp.ls))){
  dpi=as.matrix(dp.ls[[i]])
  ud.t=unique(dpi[,1])
  path.i=NULL
  for(j in 1:length(ud.t)){
    dc.j=dpi[dpi[,1]==ud.t[j],2]
    dc.j1=matrix(dpi1[which(dpi1[,NCOL(dpi1)]==ud.t[j]),],ncol=i)
    path.i=rbind(path.i, cbind(repmat(matrix(dc.j1,ncol=i),length(dc.j),1),as.character(rep(dc.j,each=NROW(dc.j1)))))
  }  
  colnames(path.i)=NULL
  na.ind <- which(is.na(dpi1[,NCOL(dpi1)]))
  
  if(length(na.ind)>0) {
    dpi1[na.ind,NCOL(dpi1)]=NA
    path.na=cbind(dpi1[na.ind,],rep(NA,length(na.ind)))
    colnames(path.na)=colnames(path.i)
    path.i=rbind(path.i, path.na)
  }
  path.i=path.i
  dpi1=path.i
}
dpi1 <- dpi1[,-NCOL(dpi1)]

# Network nodes involved in the cascading process
domisp=intersect(colnames(crnotu3),unique(c(dpi1[,-1])))

# Timing for localized transition during the cascading process
dptime=matrix(NA,NROW(dpi1),NCOL(dpi1))
dptime[,1]=tipp
timm=NULL
for(i in 1:length(domisp)){
  mv.i=filter(mv.tmp,ID==domisp[i])
  timm=c(timm, as.numeric(mv.i['warn_time']))
  dptime[which(dpi1==domisp[i])]=timm[i]
}
names(timm)=domisp
(timm)

# Extract each link during the cascading process
domi.link=NULL
domi.netw=matrix(0,length(domisp)+1,length(domisp)+1);
colnames(domi.netw) <- rownames(domi.netw) <- c('CH4', domisp)
for(i in 2:NCOL(dpi1)){
  dpi.i=dpi1[,c(i,i-1)]
  dpi.i=dpi.i[apply(!is.na(dpi.i),1,all),]
  if(!is.matrix(dpi.i)){dpi.i=matrix(dpi.i,1,2)}
  ub.i=apply(dpi.i,1,function(x){y=paste(x,collapse = '->');return(y)})
  domi.link=c(domi.link, ub.i)
  for(j in 1:NROW(dpi.i)){
    domino.j=filter(domino, Domino==dpi.i[j,1],ID==dpi.i[j,2])
    domi.netw[dpi.i[j,2], dpi.i[j,1]]=as.numeric(domino.j['Contribution'])
  }
}
(domi.link=unique(domi.link))

################################################################################
# Plot cascading network presented in Fig. 5A
g=graph.adjacency(t(domi.netw), mode="directed",weighted =T)

# Node shape (triangle & circle labeling the increasing and decreasing trends before the tipping point
mvp.bi=rep(0,length(domisp));mvp.bi[grep('1',mv.tmp[domisp,'Ttype'])]=1;mvp.bi[grep('2',mv.tmp[domisp,'Ttype'])]=2
shape.bi <- c('triangle','circle')[c(2, mvp.bi)]
# Network node layout based on x:timing for localized transition & y:dissimilarity in cascading processes among nodes
dpath <- matrix(0, length(domisp),length(domisp))
colnames(dpath) <- rownames(dpath) <- domisp
for(i in 1:length(domisp)){
  sp.i=domisp[i]
  for(j in 1:length(domisp)){
    sp.j=domisp[j]
    if(i>j){
      ij.union=sum(apply(1*(dpi1==sp.i | dpi1==sp.j),1,sum, na.rm=T)>=1)
      ij.inter=sum(apply(1*(dpi1==sp.i | dpi1==sp.j),1,sum, na.rm=T)==2)
      dpath[i,j]=1-(ij.inter/ij.union)
    }
  }
}
dpath[upper.tri(dpath)]=dpath[lower.tri(dpath)]
pcop=ape::pcoa(dpath, correction="none", rn=NULL)
pathd1=(pcop$vectors[,1]-pcop$vectors[names(which.min(abs(pcop$vectors[,1]-mean(pcop$vectors[,1])))),1])

mylay <- cbind(c(tipp,timm), c(-15,(pathd1*20 - mvp.bi*10)))

ucol=c(CH4='#FF0000', VFA='#DBB100', Fingerprint='#66A61E',
       Methanogen='#E7298A', Acetogen="#7570B3", Other_taxa="grey80")
vcol=c(ucol[1], ucol[mv.tmp[domisp,'Type']])

domisp2=c('CH4',domisp);
domisp2[grep('OTU_',domisp2)]=clasdo[domisp2[grep('OTU_',domisp2)],'finest_taxa_name']

V(g)$color=vcol
V(g)$label <- domisp2
V(g)$label.cex <- 0.7
V(g)$size <- c(1, ccm.crnotu[domisp])*8 # Node size is scaled according to the strength of causal influences on CH4
V(g)$shape <- shape.bi

E(g)$arrow.size=0.3
# Edge color presented the quantitative contribution of one node to the localized transition of the other nodes
colfunc <- colorRampPalette(c("grey90", "black"))
col.code=colfunc(102)[-1]
w.per=round(ecdf(E(g)$weight)(E(g)$weight)*100)+1
E(g)$color=col.code[w.per]

windows(160,160);
plot(g, layout=mylay, vertex.label.family="serif",cex.lab=5,edge.label.family="Palatino",
      vertex.label.cex=0.7,edge.label.cex=2, xlab = "Transition time", ylab='')
legend('topright', names(ucol), pch=rep(16,length(ucol)), col=ucol, bty='n',cex=1)
axis(side=1,at=round(2*(seq(0,40,10)-min(mylay[,1]))/(max(mylay[,1])-min(mylay[,1]))-1,2), labels =seq(0,40,10))

tipp.poly=round(2*(c(0,tipp,tipp2)-min(mylay[,1]))/(max(mylay[,1])-min(mylay[,1]))-1,2)
polygon(x=tipp.poly[c(2:3,3:2)],y=rep(c(-1.1,-0.9),each=2), col="orange",density=15, border="white")


################################################################################
# Fig. S9. Include reinforced drivers showing monotonic increasing/decreasing trends
domino.mn=NULL
for(i in 1:length(domisp)){
  mv.tmp.t0=filter(mv.tmp, ID==domisp[i])
  cht.t=as.numeric(mv.tmp.t0[c('warn_time')])
  tarid=as.character(mv.tmp.t0['ID'])
  
  ind.r=as.numeric(ceiling(cht.t)):as.numeric(1+ceiling(cht.t)); if(length(ind.r)==1){ind.r=ind.r:(ind.r+1)}
  indi=which(apply(cotrds.ary[i,,ind.r]>0,1,all))
  mv.tmp.t2=mv.tmp[indi,]%>%mutate(contribution=apply(matrix(cotrds.ary[i,indi,ind.r],ncol=2),1,sum))%>%filter(Ttype%in%c('A1','A2'))
  tarid2=mv.tmp.t2[,'ID']
  contri=mv.tmp.t2[,'contribution']
  
  if(length(tarid2)==0){}else{
  domino.mn=rbind(domino.mn, data.frame(Variable=rep(i),ID=rep(tarid),Ttype=rep(mv.tmp.t0['Ttype']),
                                        Time=rep(cht.t),Domino=tarid2,Contribution=contri))
  }
}

ch4.mn=filter(mv.tmp,ID%in%colnames(cotrdr0),Ttype%in%c('A1','A2'))

domino.mn=rbind(data.frame(Variable=0,ID='CH4',Ttype=ch4.mn[,'Ttype'],Time=tipp,Domino=ch4.mn[,'ID'],
                           Contribution=as.numeric(cotrdr0[1,ch4.mn[,'ID']])),domino.mn)
domisp.mn=intersect(colnames(crnotu3),unique(unlist(c(domino.mn[,c('Domino')])))) 

domi.netw.mn=matrix(0,length(domisp)+length(domisp.mn)+1,length(domisp)+length(domisp.mn)+1)
colnames(domi.netw.mn)=rownames(domi.netw.mn)=c('CH4',domisp,domisp.mn)
domi.link.mn=apply(domino.mn[,c('Domino','ID')],1,paste,collapse = '->')
domi.link2=c(domi.link,domi.link.mn)
for(i in 1:length(domi.link2)){
  link.i=unlist(strsplit(domi.link2[i],'->'))
  domino.i=filter(rbind(domino,domino.mn), Domino==link.i[1],ID==link.i[2])
  domi.netw.mn[link.i[2],link.i[1]]=as.numeric(domino.i['Contribution'])
}

g2=graph.adjacency(t(domi.netw.mn), mode="directed",weighted =T)

mylay2 <- cbind(c(mylay[,1],rep(-2,length(domisp.mn))),
                c(mylay[,2], seq(min(mylay[,2]),max(mylay[,2]),length.out=length(domisp.mn))[rank(-ind.co[c(domisp.mn),'Abb']+seq(0,0.1, length.out=length(domisp.mn)))])
)

mn.bi=rep(0,length(domisp.mn));mn.bi[grep('1',mv.tmp[domisp.mn,'Ttype'])]=1;mn.bi[grep('2',mv.tmp[domisp.mn,'Ttype'])]=2

domisp3=c('CH4',domisp,domisp.mn)
domisp3[grep('OTU_',domisp3)]=clasdo[domisp3[grep('OTU_',domisp3)],'finest_taxa_name']

vcol2=ucol[c(1,ind.co[c(domisp,domisp.mn),'Abb']+1)]
shape.bi2 <- c('triangle','circle')[c(2, mvp.bi, mn.bi)]

colfunc2 <- colorRampPalette(c("#FFF2F4",'pink', "red"))
col.code2=colfunc2(102)[-1]

V(g2)$color=vcol2
V(g2)$label <- domisp3
V(g2)$label.cex <- 0.7
V(g2)$size <- c(1, ccm.crnotu[c(domisp,domisp.mn)])*8
V(g2)$shape <- shape.bi2

E(g2)$arrow.size=0.3
E(g2)$width=c(rep(1.2, length(E(g)$weight)),rep(0.8, length(E(g2)$weight)-length(E(g)$weight)))
w.per=c(round(ecdf(E(g2)$weight)(E(g2)$weight)*100))
E(g2)$color=c(col.code[w.per[1:length(E(g)$weight)]+1],col.code2[w.per[-c(1:length(E(g)$weight))]+1])

windows(160,160);
plot(g2, layout=mylay2, vertex.label.family="serif",cex.lab=5,
     edge.label.family="Palatino", vertex.label.cex=0.7,edge.label.cex=2, xlab = "Warning time", ylab='')
legend('topright', names(ucol), pch=rep(16,length(ucol)), col=ucol, bty='n',cex=1)
axis(side=1,at=round(2*(seq(0,40,10)-min(mylay2[,1]))/(max(mylay2[,1])-min(mylay2[,1]))-1,2), labels =seq(0,40,10))

tipp.poly=round(2*(c(0,tipp,tipp2)-min(mylay2[,1]))/(max(mylay2[,1])-min(mylay2[,1]))-1,2)
polygon(x=tipp.poly[c(2:3,3:2)],y=rep(c(-1.1,-0.9),each=2), col="orange",density=15, border="white")


graphics.off()
################################################################################
# Fig. 5B Standardized levels of proteins involved in CH4 production
small_func <- read.csv('./Demo_RScript/Data/Demo_AD_small_summedLFQ_scaled.csv',header=T)

annotation_row <- data.frame(
  Pathway=small_func[,'AD_Pathway']
)

win.graph(100,60)
cut.k=4
da.t <- apply(t(small_func[,-c(1)]),2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
colnames(da.t)=small_func[,'AD_Pathway']
out <- pheatmap(t(da.t),cluster_cols = F, cluster_rows = F, angle_row =180, cellwidth = 100, cellheight = 100,
                border_color = "white", fontsize_row=15, fontsize_col=15, labels_col=paste('Day ',c(-10,3,11,19,27,35,43),sep=''))
