# R code for plotting Fig. 5
rm(list = ls());gc()
# Set path
OPath <- 'D:/data/Cooperation/Prodigio/Demo_RScript'
setwd(OPath)

library(factoextra)
library(igraph)
library(mgcv)
library(gratia)
library(dplyr)
library(pheatmap)
library(fields)

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


#########################################topological comparison
wtopo_node=function(X){
  n=NROW(X)
  if(any(is.na(X))){
    topov=matrix(NA,NCOL(X),5)
    colnames(topov)=c('strength.in','strength.out','strength.all','hub','auth')
  }else{
    diag(X)=0
    x.d.uw=graph.adjacency(t(X), mode="directed",weighted =T)
    topov=cbind(strength.in=graph.strength(x.d.uw, mode = "in"),
                strength.out=graph.strength(x.d.uw, mode = "out"),
                strength.all=graph.strength(x.d.uw, mode = "all"),
                hub=hub.score(x.d.uw)$vector,
                auth=authority_score(x.d.uw)$vector
    )
  }
  return(topov)
}


#####################################################################################################################################
# Data loading
clasdo <- read.csv('./Data/demo_otu_taxa.csv',header=T)
rownames(clasdo)=clasdo[,'OTU']
crnotu1 <- read.csv('./Data/Demo_Control_otucrn_abundance.csv',header=T,stringsAsFactors = F)
crnotu3 <- read.csv('./Data/Demo_Perturbed_otucrn_abundance.csv',header=T,stringsAsFactors = F)
reactor_ch4_ews <- read.csv( './Data/reactor_ch4_ews.csv', header=T)

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
mco3.ins <- read.csv("./Output/Demo_CH4_otucrn_Perturbed_nin46_cvunit0.025_enet_jcof_Nmvx_Rallx_demo_NEW.csv",stringsAsFactors = F)[,-c(1:5)]

mco3.ins2 <- read.csv("D:/data/Cooperation/Prodigio/Demo_RScript_old/Output/Demo_CH4_otucrn_Perturbed_nin46_cvunit0.025_enet_jcof_Nmvx_Rallx_demo_NEW.csv",stringsAsFactors = F)[,-c(1:5)]

co1.ins <- read.csv('./Output/Control_otucrn_verc_Species_nin46_cvunit0.025_enet_jcof_Nmvx_Rallx_demo.csv',stringsAsFactors = F)
co3.ins <- read.csv('./Output/Perturbed_otucrn_verc_Species_nin46_cvunit0.025_enet_jcof_Nmvx_Rallx_demo.csv',stringsAsFactors = F)

co.ary3 <- array(NA,dim=c(nco3,nco3,nd-1))
for(i in 1:(nd-1)){
  ints.i31=as.matrix(filter(co3.ins,time==i)[,-c(1:4)])
  co.ary3[,,i]=ints.i31
}

################################################################################
### CCM results
ccm.otu <- read.csv('./Output/demo_ccm_otu.csv', header=T)
ccm.crn <- read.csv('./Output/demo_ccm_metabolite.csv', header=T)

ccm.otu.c <- ccm.otu[,'ccm_OTU_control']
ccm.crn.c <- ccm.crn[,'ccm_metabolite_control']

ccm.otu <- ccm.otu[,'ccm_OTU_perturbed']
ccm.crn <- ccm.crn[,'ccm_metabolite_perturbed']

ccm.crnotu <- c(ccm.crn,ccm.otu)

names(ccm.crnotu) <- colnames(crnotu3)
################################################################################
# Contributions for each species to the CH4 critical transition at the tipping point
da.i=data.frame(value=reactor_ch4_ews[,'CH4_perturbed'],time=1:nd)
g.t <- gam(value~s(time),data=da.i)
newd.i=seq(0,nd,0.01)
yhat=predict.gam(g.t,newdata=data.frame(time=newd.i))
yhatd1 <- derivatives(g.t, type = "central",data=data.frame(time=newd.i), order=1)$derivative

# Contribution of each node on the state changes in CH4 based on the dynamical interactome and time series of state variables 
impc=mco3.ins*crnotu3s[-nd,]
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
smo.crnotu3 <- filter(read.csv('./Output/smo_otucrn3_s.csv', header = T), time%in%(1:nd))
#smoin.crnotu3=filter(read.csv('smo_otucrn3_37.csv', header = T), time%in%(1:nd))

crnotu3.smo=matrix(smo.crnotu3[,'yhat'],nrow=tipp)
crnotu3.smod=matrix(smo.crnotu3[,'est'],nrow=tipp)

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

mv.tmp <- read.csv('./Output/Localized_Transition_crnotu_s.csv',header = T)
warn.crnotu <- read.csv('./Output/warn_time_CRNOTU_ls_sfill.csv',header = T)

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
g <- graph.adjacency(t(domi.netw), mode="directed",weighted =T)

# Node shape (triangle & circle labeling the increasing and decreasing trends before the tipping point
mvp.bi <- rep(0,length(domisp)); mvp.bi[grep('1',mv.tmp[domisp,'Ttype'])] <- 1; mvp.bi[grep('2',mv.tmp[domisp,'Ttype'])] <- -1 #2
mvp.bi <- c(-1, mvp.bi) 
shape.bi <- c('triangle','circle')[ifelse(mvp.bi==1,1,2)]

path.num <- rep(0,1+length(domisp))
for(i in 1:length(domisp)){
  path.i <- all_simple_paths(g,domisp[i] ,'CH4', mode="out")
  path.node <- table(unlist(path.i))
  path.num.i <- rep(0,1+length(domisp))
  path.num.i[as.numeric(names(path.node))] <- as.numeric(path.node)
  path.num <- path.num + path.num.i
}

cen.i <- max(path.num[-1])-path.num[-1]
cen.i <- rank(cen.i, ties.method ='first')*c(mvp.bi[-1])

mylay <- cbind(c(tipp,timm), c(0, cen.i))

ucol <- c(CH4='#FF0000', VFA='#DBB100', Fingerprint='#66A61E',
       Methanogen='#E7298A', Acetogen="#7570B3", Other_taxa="grey80")
vcol <- c(ucol[1], ucol[mv.tmp[domisp,'Type']])

domisp2 <- c('CH4',domisp);
domisp2[grep('OTU_',domisp2)]=clasdo[domisp2[grep('OTU_',domisp2)],'finest_taxa_name']

V(g)$color <- vcol
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
#postscript(paste('D:/data/Cooperation/Prodigio/Manuscript/Figure/Domino effect_verc_ls.eps',sep=''), fonts=c("serif", "Palatino"))
plot(g, layout=mylay, vertex.label.family="serif",cex.lab=5,edge.label.family="Palatino",
      vertex.label.cex=0.7,edge.label.cex=2, xlab = "Transition time", ylab='')
legend('topright', names(ucol), pch=rep(16,length(ucol)), col=ucol, bty='n',cex=1)
axis(side=1,at=round(2*(seq(0,40,10)-min(mylay[,1]))/(max(mylay[,1])-min(mylay[,1]))-1,2), labels =seq(0,40,10))

tipp.poly=round(2*(c(0,tipp,tipp2)-min(mylay[,1]))/(max(mylay[,1])-min(mylay[,1]))-1,2)
polygon(x=tipp.poly[c(2:3,3:2)],y=rep(c(-1.1,-0.9),each=2), col="orange",density=15, border="white")

#Display the graph and the legend.
ordered <- order(col.code)
color <- vector(length = length(ordered),mode="double")
for (i in 1:length(ordered)) color[ ordered[i] ] <- col.code [ i ]

image.plot(legend.only=T, zlim=range(E(g)$weight), col=col.code, 
           legend.shrink= 0.25, legend.cex=0.75, smallplot=c(0.025,0.05,0.35,0.5))
#dev.off()


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

mn.bi=rep(0,length(domisp.mn));mn.bi[grep('1',mv.tmp[domisp.mn,'Ttype'])]=1;mn.bi[grep('2',mv.tmp[domisp.mn,'Ttype'])]=-1

path.num.mn <- rep(0,length(domisp.mn))
for(i in 1:length(domisp.mn)){
  path.i <- all_simple_paths(g2, domisp.mn[i], 'CH4', mode="out")
  path.node <- table(unlist(path.i))
  path.node <- path.node[as.numeric(names(path.node))>(length(domisp)+1)]
  path.num.i <- rep(0,length(domisp.mn))
  path.num.i[as.numeric(names(path.node))-(length(domisp)+1)] <- as.numeric(path.node)
  path.num.mn <- path.num.mn + path.num.i
}

cen.i <- max(path.num.mn)-path.num.mn
cen.i <- rank(cen.i)*c(mn.bi)


lay22 <- seq(min(mylay[,2]),max(mylay[,2]),length.out=length(domisp.mn))[rank(cen.i, ties.method ='first')]
mylay2 <- cbind(c(mylay[,1],rep(-2,length(domisp.mn))),c(mylay[,2], lay22))


domisp3=c('CH4',domisp,domisp.mn)
domisp3[grep('OTU_',domisp3)]=clasdo[domisp3[grep('OTU_',domisp3)],'finest_taxa_name']

vcol2=ucol[c(1,ind.co[c(domisp,domisp.mn),'Abb']+1)]

shape.bi2 <- c('triangle','circle')[ifelse(c(mvp.bi,mn.bi)==1,1,2)]

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
#postscript(paste('./Manuscript/Figure/Domino effect_verc_ls_mn.eps',sep=''), fonts=c("serif", "Palatino"))
plot(g2, layout=mylay2, vertex.label.family="serif",cex.lab=5,
     edge.label.family="Palatino", vertex.label.cex=0.7,edge.label.cex=2, xlab = "Warning time", ylab='')
legend('topright', names(ucol), pch=rep(16,length(ucol)), col=ucol, bty='n',cex=1)
axis(side=1,at=round(2*(seq(0,40,10)-min(mylay2[,1]))/(max(mylay2[,1])-min(mylay2[,1]))-1,2), labels =seq(0,40,10))

tipp.poly=round(2*(c(0,tipp,tipp2)-min(mylay2[,1]))/(max(mylay2[,1])-min(mylay2[,1]))-1,2)
polygon(x=tipp.poly[c(2:3,3:2)],y=rep(c(-1.1,-0.9),each=2), col="orange",density=15, border="white")

#Display the graph and the legend.
ordered <- order(col.code)
color <- vector(length = length(ordered),mode="double")
for (i in 1:length(ordered)) color[ ordered[i] ] <- col.code [ i ]

image.plot(legend.only=T, zlim=range(E(g2)$weight[1:length(domisp2)]), col=col.code, 
           legend.shrink= 0.25, legend.cex=0.75, smallplot=c(0.025,0.05,0.55,0.7))

ordered2 <- order(col.code2)
color2 <- vector(length = length(ordered2),mode="double")
for (i in 1:length(ordered2)) color2[ ordered2[i] ] <- col.code2 [ i ]

image.plot(legend.only=T, zlim=range(E(g2)$weight[-c(1:length(domisp2))]), col=col.code2, 
           legend.shrink= 0.25, legend.cex=0.75, smallplot=c(0.025,0.05,0.35,0.5))

#dev.off()

graphics.off()
################################################################################
# Fig. 5B Standardized levels of proteins involved in CH4 production
small_func <- read.csv('./Data/Demo_AD_small_summedLFQ_scaled.csv',header=T)

annotation_row <- data.frame(
  Pathway=small_func[,'AD_Pathway']
)

win.graph(100,60)
cut.k=4
da.t <- apply(t(small_func[,-c(1)]),2,function(x){y=(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))})
colnames(da.t)=small_func[,'AD_Pathway']
out <- pheatmap(t(da.t),cluster_cols = F, cluster_rows = F, angle_row =180, cellwidth = 100, cellheight = 100,
                border_color = "white", fontsize_row=15, fontsize_col=15, labels_col=paste('Day ',c(-10,3,11,19,27,35,43),sep=''))



###################################
### properties of domino nodes
domino.ccm.o=data.frame(type=c(1*colnames(otu3)%in%domisp+2*colnames(otu3)%in%domisp.mn+3*!colnames(otu3)%in%c(domisp,domisp.mn)),ccm=ccm.otu)
domino.ccm.c=data.frame(type=c(1*colnames(crn3)%in%domisp+2*colnames(crn3)%in%domisp.mn+3*!colnames(crn3)%in%c(domisp,domisp.mn)),ccm=ccm.crn)

domino.ccm.o12=data.frame(type=c(1*colnames(otu1)%in%domisp+2*colnames(otu1)%in%domisp.mn+3*!colnames(otu1)%in%c(domisp,domisp.mn)),ccm=ccm.otu.c)
domino.ccm.c12=data.frame(type=c(1*colnames(crn1)%in%domisp+2*colnames(crn1)%in%domisp.mn+3*!colnames(crn1)%in%c(domisp,domisp.mn)),ccm=ccm.crn.c[-c(1:7)])

win.graph(70,150)
par(mfrow=c(5,2),mar=c(1,4,1,1))

# Causal strength on CH4 for those nodes involved in the cascading mechanisms 
# Causal strength of microbial OTUs on CH4 in the control reactor
boxplot(ccm~type, domino.ccm.o12, lwd = 2, xlab='Type',ylab='CCM skill',ylim=c(0,1.05))
stripchart(ccm~type, vertical = TRUE, data = domino.ccm.o12, cex=0.5,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
kruskal.test(ccm~type,domino.ccm.o12)
dunn.test::dunn.test(x=list(filter(domino.ccm.o12, type==1)[,'ccm'],filter(domino.ccm.o12, type==2)[,'ccm'],filter(domino.ccm.o12, type==3)[,'ccm']))

# Causal strength of microbial OTUs on CH4 in the perturbed reactor
boxplot(ccm~type, domino.ccm.o, lwd = 2, xlab='Type',ylab='CCM skill',ylim=c(0,1.05))
stripchart(ccm~type, vertical = TRUE, data = domino.ccm.o, cex=0.5,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
kruskal.test(ccm~type,domino.ccm.o)
dunn.test::dunn.test(x=list(filter(domino.ccm.o, type==1)[,'ccm'],filter(domino.ccm.o, type==2)[,'ccm'],filter(domino.ccm.o, type==3)[,'ccm']))

# Causal strength of metabolite on CH4 in the control reactor
boxplot(ccm~type, domino.ccm.c12, lwd = 2, xlab='Type',ylab='CCM skill',ylim=c(0,1.05))
stripchart(ccm~type, vertical = TRUE, data = domino.ccm.c12, cex=0.5,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
kruskal.test(ccm~type,domino.ccm.c12)

# Causal strength of metabolite on CH4 in the perturbed reactor
boxplot(ccm~type, domino.ccm.c, lwd = 2, xlab='Type',ylab='CCM skill',ylim=c(0,1.05))
stripchart(ccm~type, vertical = TRUE, data = domino.ccm.c, cex=0.5,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
kruskal.test(ccm~type,domino.ccm.c)

# Centrality analysis for the nodes involved in the cascading mechanism 
da.p=abs(co3.ins)
dai.p=da.p[,c(1:4)]
da.p=da.p[,-c(1:4)]

da.c=abs(co1.ins)
dai.c=da.c[,c(1:4)]
da.c=da.c[,-c(1:4)]

ind.abb=ind.co[,'Abb']

centv.p <- centv.c <- NULL
for(i in 1:tipp){
  co.p=as.matrix(da.p[dai.p[,'time']==i,])
  topov.pi=wtopo_node(co.p)
  centv.p=rbind(centv.p,data.frame(time=i, variable=1:NCOL(co.p), topov.pi))
  
  co.c=as.matrix(da.c[dai.c[,'time']==i,])
  topov.ci=wtopo_node(co.c)
  centv.c=rbind(centv.c,data.frame(time=i, variable=1:NCOL(co.c), topov.ci))
}

centsp.p <- data.frame(domino=1*(ind.co[,'Name']%in%domisp)+2*(ind.co[,'Name']%in%domisp.mn),aggregate(centv.p ,list(centv.p[,'variable']),median,na.rm=T)[,-c(1:2)])
centsp.p[centsp.p [,'domino']==0,'domino'] <- 3
aggregate(centsp.p,list(centsp.p[,'domino']),median)

centsp.c <- data.frame(domino=1*(ind.co[-c(1:7),'Name']%in%domisp)+2*(ind.co[-c(1:7),'Name']%in%domisp.mn),aggregate(centv.c ,list(centv.c[,'variable']),median,na.rm=T)[,-c(1:2)])
centsp.c[centsp.c [,'domino']==0,'domino'] <- 3
aggregate(centsp.c,list(centsp.c[,'domino']),median)

da.box.c <- centsp.c[ind.abb>2,]
da.box.p <- centsp.p[ind.abb>2,]
# Node centrality of those involved in the cascading process happened in the perturbed reactor
# Averaged total interaction strength (inward+outward) in the control reactor
vb <- 'strength.all'
formula.t <- as.formula(paste(vb,'domino',sep='~'))

boxplot(formula.t,da.box.c, lwd = 2, xlab='Type',ylab=vb)#,ylim=c(0.2,1.05))
stripchart(formula.t, vertical = TRUE, data = da.box.c, cex=0.5,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
kruskal.test(formula.t,da.box.c)
dunn.test::dunn.test(da.box.c[,vb], as.factor(da.box.c[,'domino']))
# Magnitude of interaction strength (inward+outward) in the perturbed reactor
boxplot(formula.t,da.box.p, lwd = 2, xlab='Type',ylab=vb)#,ylim=c(0.2,1.05))
stripchart(formula.t, vertical = TRUE, data = da.box.p, cex=0.5,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
kruskal.test(formula.t,da.box.p)
dunn.test::dunn.test(da.box.p[,vb], as.factor(da.box.p[,'domino']))

# Hub index in the control reactor
vb <- 'hub'
formula.t <- as.formula(paste(vb,'domino',sep='~'))
boxplot(formula.t,da.box.c, lwd = 2, xlab='Type',ylab=vb)#,ylim=c(0.2,1.05))
stripchart(formula.t,, vertical = TRUE, data = da.box.c, cex=0.5,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
kruskal.test(formula.t,da.box.c)
dunn.test::dunn.test(da.box.c[,vb], as.factor(da.box.c[,'domino']))
# Hub index in the perturbed reactor
boxplot(formula.t,da.box.p, lwd = 2, xlab='Type',ylab=vb)#,ylim=c(0.2,1.05))
stripchart(formula.t, vertical = TRUE, data = da.box.p, cex=0.5,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
kruskal.test(formula.t,da.box.p)
dunn.test::dunn.test(da.box.p[,vb], as.factor(da.box.p[,'domino']))

# Authority index in the control reactor
vb <- 'auth'
formula.t <- as.formula(paste(vb,'domino',sep='~'))
boxplot(formula.t,da.box.c, lwd = 2, xlab='Type',ylab=vb)#,ylim=c(0.2,1.05))
stripchart(formula.t,, vertical = TRUE, data = da.box.c, cex=0.5,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
kruskal.test(formula.t,da.box.c)
dunn.test::dunn.test(da.box.c[,vb], as.factor(da.box.c[,'domino']))
# Authority in the perturbed reactor
boxplot(formula.t,da.box.p, lwd = 2, xlab='Type',ylab=vb)#,ylim=c(0.2,1.05))
stripchart(formula.t, vertical = TRUE, data = da.box.p, cex=0.5,
           method = "jitter", add = TRUE, pch = 20, col = 'blue')
kruskal.test(formula.t,da.box.p)
dunn.test::dunn.test(da.box.p[,vb], as.factor(da.box.p[,'domino']))
