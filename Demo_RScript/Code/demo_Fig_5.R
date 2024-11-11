# R code for plotting Fig. 5
rm(list = ls());gc()
# Set path
OPath <- 'D:\\data\\Cooperation\\Prodigio\\'
setwd(OPath)

library(factoextra)
library(igraph)

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
crnotu1 <- read.csv('./Demo_RScript/Data/OLR1_otucrn_verc_Species_tau132.csv',header=T,stringsAsFactors = F)
crnotu3 <- read.csv('./Demo_RScript/Data/OLR3_otucrn_verc_Species_tau132.csv',header=T,stringsAsFactors = F)
reactor_ch4_ews <- read.csv( './Demo_RScript/Data/reactor_ch4_ews.csv', header=T)

# Calculate percentage data for microbial OTU 
otu1.s <- crnotu1[,grep('OTU_',colnames(crnotu1))]
otu3.s <- crnotu3[,grep('OTU_',colnames(crnotu3))]
otu1.s <- otu1.s*repmat(apply(otu1.s,1,sum),1,ncol(otu1.s))^-1
otu3.s <- otu3.s*repmat(apply(otu3.s,1,sum),1,ncol(otu3.s))^-1

# Extract metabolomic data
crn1.o <- crnotu1[,-grep('OTU_',colnames(crnotu1))]
crn3.o <- crnotu3[,-grep('OTU_',colnames(crnotu3))]

rownames(otu3.s)=rownames(crn3.o)=1:nd

nd=NROW(otu3.s)
nsp=NCOL(otu3.s)
ncrn=NCOL(crn1.o)

nco1=ncrn+nsp
nco3=ncrn+nsp+7

crnotu3=cbind(crn3.o,otu3.s)
crnotu1=cbind(crn1.o,otu1.s)
crnotu3s=apply(crnotu3,2,scale)


# Loading interactomic data
mco3.ins=read.csv('./Demo_RScript/Data/CH4_crnotu_verc_OLR3_Species_nin46_cvunit0.025_enet_jcof_Nmvx_Rallx_demo_NEW.csv',stringsAsFactors = F)
co3.ins=read.csv('./Demo_RScript/Data/OLR3_otucrn_verc_Species_nin46_cvunit0.025_enet_jcof_Nmvx_Rallx_demo_NEW.csv',stringsAsFactors = F)

co.ary3=coa.ary3=array(NA,dim=c(nco3,nco3,nd))
for(i in 1:(nd-1)){
  ints.i31=as.matrix(filter(co3.ins,time==i)[,-c(1:4)])
  co.ary3[,,i]=ints.i31
}


################################################################################
# Contributions for each species to the CH4 critical transition at the tipping point
da.i=data.frame(value=reactor_ch4_ews[,'CH4_perturbed'],time=1:nd)
g.t <- gam(value~s(time),data=da.i)
newd.i=seq(0,nd,0.01)
yhat=predict.gam(g.t,newdata=data.frame(time=newd.i))
yhatd1 <- derivatives(g.t, type = "central",data=data.frame(time=newd.i), order=1)$derivative

# Contribution of each node on the state changes in CH4 
impc=mco3.ins[,-c(1:5)]*crnotu3s[-nd,]
cotr=rbind(rep(NA),impc[-1,]-impc[-(nd-1),])
# Sign of the local trend (t, t+1) defined as the estimated first derivative evaluated at t
sdch4s=sign(yhatd1[newd.i%in%(seq(1,nd-1,1))])
# Examine whether the sign of local trend equals to the sign of contribution  
drcotrs=sign(cotr)*sdch4s
# Corrected contributions depending on whether the sign of contribution facilitate(+) or buffer (-) the observed trend in CH4
cotrds=drcotrs*abs(cotr)


# Examine the contribution
smoin.crnotu3=filter(read.csv('smo_otucrn3_37.csv', header = T), time%in%(1:nd))

crnotu3.smo=matrix(smoin.crnotu3[,'yhat'],nrow=tipp)
crnotu3.smod=matrix(smoin.crnotu3[,'est'],nrow=tipp)

cotrd.ary=cotrds.ary=array(NA, dim=c(tipp,NCOL(crnotu3),NCOL(crnotu3)))
for(i in 1:NCOL(crnotu3)){
  impc.i=t(co.ary3[i,,])*crnotu3s[-nd,]
  cotr.i=rbind(rep(NA),impc.i[-1,]-impc.i[-(nd-1),])
  cotr.i=cotr.i[1:tipp,]
  
  sdxis=sign(crnotu3.smod[,i])
  drcotrs.i=sign(cotr.i)*sdxis
  cotrds.i=drcotrs.i*abs(cotr.i)
  cotrds.ary[,,i]=cotrds.i
}

domino=NULL
for(i in 1:NCOL(crnotu3)){
  mv.tmp.t0=filter(mv.tmp, ID==colnames(crnotu3)[i])
  cht.t=sum(mv.tmp.t0[c('turn_time','change_time')])
  tarid=as.character(mv.tmp.t0['ID'])
  
  if(mv.tmp.t0['Ttype']%in%c('C1','C2','B1','B2') & !is.na(cht.t)){
    if(step.c==1){
      ind.r=as.numeric(ceiling(cht.t))
      indi=which(cotrds.ary[ind.r,,i]>0)
    }else{
      ind.r=as.numeric(ceiling(cht.t)):as.numeric(step.c-1+ceiling(cht.t)); if(length(ind.r)==(step.c-1)){ind.r=ind.r:(ind.r+1)}
      indi=which(apply(cotrds.ary[ind.r,,i]>0,2,all))
    }
    
    mv.tmp.t2=mv.tmp[indi,]
    mv.tmp.t2=mv.tmp.t2%>%mutate(contribution=apply(matrix(cotrds.ary[ind.r,indi,i],nrow=length(ind.r)),2,sum))%>%
      filter(Ttype%in%c('B1','B2','C1','C2'))%>%arrange(desc(turn_time+change_time))%>%filter((turn_time+change_time)<trunc(cht.t))
    tarid2=mv.tmp.t2[,'ID']
    contri=mv.tmp.t2[,'contribution']
  }else{tarid2=NA;contri=NA}
  if(length(tarid2)==0){tarid2=NA;contri=NA}
  domino=rbind(domino, data.frame(Variable=rep(i),ID=rep(tarid),Ttype=rep(mv.tmp.t0['Ttype']),Time=rep(cht.t),Domino=tarid2,Contribution=contri))
}
domino=filter(domino, ID!=Domino)
domino=filter(domino, ID%in%warn.otucrn[,1] & Domino%in%warn.otucrn[,1])

ch4.curv=filter(mv.tmp,ID%in%colnames(cotrdr0),Ttype%in%c('B1','B2','C1','C2'))
domino=rbind(data.frame(Variable=0,ID='CH4',Ttype='C2',Time=tipp,Domino=ch4.curv[,'ID'],Contribution=as.numeric(cotrdr0[1,ch4.curv[,'ID']])),domino)

mv.tmp.t0=mv.tmp.t
mv.tmp.t0=mv.tmp.t0%>%filter(ID!='CH4')%>%filter(Ttype%in%c('B1','B2','C1','C2'))%>%filter(turn_time+change_time>=2)%>%arrange(desc(turn_time+change_time))
chtime=apply(mv.tmp.t0[,c('turn_time','change_time')],1,sum)

dpath=cbind(rep('CH4'), mv.tmp.t0[,'ID'])
dpath.k0=dpath
dp.ls=list()
dp.ls[[1]]=data.frame(ind=1,ID=dpath[,2])
for(i in 2:30){
  dpath.kt=dp.ls[[i-1]][,2]#dpath.k0[,k]
  ud.t=unique(dpath.kt)
  if(any(is.na(ud.t))){ud.t <- ud.t[!is.na(ud.t)]}
  d.t2=NULL#dpath.k1=NULL
  for(k in 1:length(ud.t)){
    d.t=filter(domino,ID==ud.t[k])[,'Domino']
    indp=which(dpath.kt==ud.t[k])
    if( length(d.t)!=0 & all(!is.na(d.t))){
      #dpath.k1=rbind(dpath.k1,cbind(repmat(matrix(dpath.k0[i,],nrow=1),length(d.t),1),d.t))
      d.t2=rbind(d.t2, cbind(ind=rep(ud.t[k],each=length(d.t)),ID=rep(d.t,length(ud.t[k]))))
    }else{
      d.t2=rbind(d.t2, cbind(ind=rep(ud.t[k],each=1),ID=rep(NA,1)))
    }
  }
  print(paste(i,NROW(d.t2)))
  dp.ls[[i]]=as.matrix(d.t2)
  if(all(is.na(d.t2[,2]))){break}else{
    rm(d.t2);gc();
  }
}

dp.ls0 <- dp.ls
dpi1 <- as.matrix(cbind(rep('CH4'), dp.ls[[1]][,2]))

for(i in 2:length(dp.ls)){
  dpi=as.matrix(dp.ls[[i]])
  ud.t=unique(dpi[,1])
  dpath.i=NULL
  for(j in 1:length(ud.t)){
    dc.j=dpi[dpi[,1]==ud.t[j],2]
    dc.j1=matrix(dpi1[which(dpi1[,NCOL(dpi1)]==ud.t[j]),],ncol=i)
    dpath.i=rbind(dpath.i, cbind(repmat(matrix(dc.j1,ncol=i),length(dc.j),1),as.character(rep(dc.j,each=NROW(dc.j1)))))
  }  
  colnames(dpath.i)=NULL
  na.ind <- which(is.na(dpi1[,NCOL(dpi1)]))
  
  if(length(na.ind)>0) {
    dpi1[na.ind,NCOL(dpi1)]=NA
    dpath.na=cbind(dpi1[na.ind,],rep(NA,length(na.ind)))
    colnames(dpath.na)=colnames(dpath.i)
    dpath.i=rbind(as.matrix(dpath.i), as.matrix(dpath.na))
  }
  
  dpath.i=as.matrix(dpath.i)
  dpi1=dpath.i
  rm(dpath.i);gc();
  
  print(i)
}
dpi1=dpi1[,-NCOL(dpi1)]


step.n=NCOL(dpi1)
dptime=matrix(NA,NROW(dpi1),NCOL(dpi1))
dptime[,1]=tipp
spc=c(colnames(crn3.o),colnames(otu3)[spind6])
domisp=intersect(colnames(crnotu3),unique(c(dpi1[,-1])))
sort(table(c(dpi1[,-1])),decreasing = T)
timm=NULL
for(i in 1:length(domisp)){
  mv.i=filter(mv.tmp,ID==domisp[i])
  timm=c(timm, sum(mv.i[c('turn_time','change_time')]))
  dptime[which(dpi1==domisp[i])]=timm[i]
  
}
names(timm)=domisp
domisp2=c('CH4',domisp)
dpi2=dpi1;dpi2[which(is.na(dpi1))]=""
dpath=matrix(0, length(domisp),length(domisp))
for(i in 1:length(domisp)){
  sp.i=domisp[i]
  for(j in 1:length(domisp)){
    sp.j=domisp[j]
    if(i>j){
      ij.union=sum(apply(1*(dpi2==sp.i | dpi2==sp.j),1,sum)>=1)
      ij.inter=sum(apply(1*(dpi2==sp.i | dpi2==sp.j),1,sum)==2)
      dpath[i,j]=1-(ij.inter/ij.union)
    }
  }
}
colnames(dpath)=rownames(dpath)=domisp
dpath[upper.tri(dpath)]=dpath[lower.tri(dpath)]
pcop=ape::pcoa(dpath, correction="none", rn=NULL)

pathd1=(pcop$vectors[,1]-pcop$vectors[names(which.min(abs(pcop$vectors[,1]-mean(pcop$vectors[,1])))),1])

dpi.sta=data.frame(summ=apply(dptime,1,sum,na.rm=T),maxx=apply(dptime[,-1],1,max,na.rm=T),minn=apply(dptime,1,min,na.rm=T),lengthh=apply(dptime,1,function(x){sum(!is.na(x))}))

dosp.c=matrix(0,length(domisp), step.n-1)
for(i in 2:(step.n)){
  tb.i=table(dpi1[,i])
  dosp.c[domisp%in%names(tb.i),i-1]=as.numeric(tb.i)
}
rownames(dosp.c)=domisp
avpath=rowSums((dosp.c/rowSums(dosp.c))*repmat(c(1:(NCOL(dpi1)-1)),length(domisp),1))
domisp.s=apply(dosp.c,2,function(x){domisp[order(x,decreasing = T)[1:3]]})
tio.t=unique(c(domisp.s))


domi.link=NULL
for(i in 2:NCOL(dpi1)){
  dpi.i=dpi1[,c(i,i-1)]
  dpi.i=dpi.i[apply(!is.na(dpi.i),1,all),]
  if(!is.matrix(dpi.i)){dpi.i=matrix(dpi.i,1,2)}
  ub.i=apply(dpi.i,1,function(x){y=paste(x,collapse = '->');return(y)})
  domi.link=c(domi.link, ub.i)
}

colfunc <- colorRampPalette(c("grey90", "black"))
col.code=colfunc(102)[-1]
colfunc2 <- colorRampPalette(c("#FFF2F4",'pink', "red"))
col.code2=colfunc2(102)[-1]#c(colorRampPalette(c('white'))(1),colfunc2(100))


domi.link=unique(domi.link)
domi.mvp=filter(mv.tmp,ID%in%domisp)
mvp.bi=rep(0,length(domisp));mvp.bi[grep('1',domi.mvp[,'Ttype'])]=1;mvp.bi[grep('2',domi.mvp[,'Ttype'])]=2
domi.taxa= data.frame(filter(taxa.otu,OTU%in%domisp),ind3.co2[c(domisp[grep('OTU_',domisp)]),])
domi.net=domi.netw=matrix(0,length(domisp)+1,length(domisp)+1)
colnames(domi.net)=rownames(domi.net)=colnames(domi.netw)=rownames(domi.netw)=c('CH4',domisp)
for(i in 1:length(domi.link)){
  link.i=unlist(strsplit(domi.link[i],'->'))
  domino.i=filter(domino, Domino==link.i[1],ID==link.i[2])
  domi.net[link.i[2],link.i[1]]=1
  domi.netw[link.i[2],link.i[1]]=as.numeric(domino.i['Contribution'])
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

shape=c('triangle','circle')

g=graph.adjacency(t(domi.netw), mode="directed",weighted =T)

mylay <- layout.auto(g)
mf1=0.5;mf2=20
pathd2=(pathd1- mvp.bi*mf1)*mf2
mylay[,1]=c(tipp,timm);mylay[,2]=c(-mf1*mf2*1.5,pathd2);

indx=c('CH4','VFA','Fingerprint','Methanogen','Acetogen','Other_OTU')
ucol=c('#FF0000', '#DBB100', '#66A61E','#E7298A', "#7570B3", "grey80")#rainbow(length(indx))
names(ind3.co)=colnames(crnotu3)
vcol=ucol[c(1,ind3.co[domisp]+1)]

domisp2=gsub('OTU_','',domisp)
domisp3=domisp
domisp3[grep('OTU_',domisp)]=domi.taxa[,'fine.taxa']

V(g)$color=vcol
V(g)$label <- c('CH4',domisp3)
V(g)$label.cex <- 0.7
V(g)$size <- c(1,ccm.crnotu3[domisp])*8# 4
V(g)$shape <- shape[c(2,mvp.bi)]


#E(g)$weight=abs(E(g)$weight)
E(g)$arrow.size=0.3

w.per=round(ecdf(E(g)$weight)(E(g)$weight)*100)
E(g)$color=col.code[w.per+1]

windows(160,160);
#postscript(paste('./Manuscript/Figure/Domino effect_verc_',step.ch4,step.c,'.eps',sep=''), fonts=c("serif", "Palatino"))
plot(g, layout=mylay, vertex.label.family="serif",cex.lab=5, #vertex.shape= shape[c(2,mvp.bi)],
     edge.label.family="Palatino", vertex.label.cex=0.7,edge.label.cex=2, xlab = "Transition time", ylab='')
legend('bottomright', indx, pch=rep(16,length(indx)), col=ucol, bty='n',cex=1)

axis(side=1,at=round(2*(seq(0,40,10)-min(mylay[,1]))/(max(mylay[,1])-min(mylay[,1]))-1,2), labels =seq(0,40,10))
box()
#dev.off()

domino.ls=data.frame(as.matrix(ind3.co[domisp]),timm);colnames(domino.ls)=c('Type','Time')
domino.ls[,1]=indx[-1][domino.ls[,1]]
domino.ls=data.frame(domino.ls,fine.taxa=c(rep(NA,  min(grep('OTU',rownames(domino.ls)))-1), domi.taxa[,'fine.taxa']))
arrange(domino.ls, Type, Time)


# Monotonic driver
domino.mn=NULL
for(i in 1:length(domisp)){
  mv.tmp.t0=filter(mv.tmp, ID==domisp[i])
  cht.t=sum(mv.tmp.t0[c('turn_time','change_time')])
  tarid=as.character(mv.tmp.t0['ID'])
  if(step.c==1){
    ind.r=as.numeric(ceiling(cht.t))
    indi=which(cotrds.ary[ind.r,,i]>0)
  }else{
    ind.r=as.numeric(ceiling(cht.t)):as.numeric(step.c-1+ceiling(cht.t)); if(length(ind.r)==(step.c-1)){ind.r=ind.r:(ind.r+1)}
    indi=which(apply(cotrds.ary[ind.r,,i]>0,2,all))
  }
  mv.tmp.t2=mv.tmp[indi,]
  mv.tmp.t2=mv.tmp.t2%>%mutate(contribution=apply(matrix(cotrds.ary[ind.r,indi,i],nrow=length(ind.r)),2,sum))%>%
    filter(Ttype%in%c('A1','A2'))#%>%filter((turn_time+change_time)<=1)
  tarid2=mv.tmp.t2[,'ID']
  contri=mv.tmp.t2[,'contribution']
  
  if(length(tarid2)==0){}else{
    domino.mn=rbind(domino.mn, data.frame(Variable=rep(i),ID=rep(tarid),Ttype=mv.tmp.t2['Ttype'],Time=rep(cht.t),Domino=tarid2,Contribution=contri))
    mv.tmp.t2=NULL
  }
}


ch4.mn=filter(mv.tmp,ID%in%colnames(cotrdr0),Ttype%in%c('A1','A2'))
domino.mn=rbind(data.frame(Variable=0,ID='CH4',Ttype=ch4.mn[,'Ttype'],Time=tipp,Domino=ch4.mn[,'ID'],
                           Contribution=as.numeric(cotrdr0[1,ch4.mn[,'ID']])),domino.mn)

domisp.mn=intersect(colnames(crnotu3),unique(unlist(c(domino.mn[,c('Domino')])))) 

pathd2=pathd1[domino.mn[,'ID']];pathd2[is.na(pathd2)]=0
timm.mn=aggregate(cbind(Time=domino.mn[,'Time'],pathd2),list(domino.mn[,'Domino']),mean)

domi.net.mn=domi.netw.mn=matrix(0,length(domisp)+length(domisp.mn)+1,length(domisp)+length(domisp.mn)+1)
colnames(domi.net.mn)=rownames(domi.net.mn)=colnames(domi.netw.mn)=rownames(domi.netw.mn)=c('CH4',domisp,domisp.mn)
domi.link.mn=pastt(domino.mn[,c('Domino','ID')],'->')
domi.link2=c(domi.link,domi.link.mn)
for(i in 1:length(domi.link2)){
  link.i=unlist(strsplit(domi.link2[i],'->'))
  domino.i=filter(rbind(domino,domino.mn), Domino==link.i[1],ID==link.i[2])
  domi.net.mn[link.i[2],link.i[1]]=1
  domi.netw.mn[link.i[2],link.i[1]]=as.numeric(domino.i['Contribution'])
}


g2=graph.adjacency(t(domi.netw.mn), mode="directed",weighted =T)

mylay2 <- layout.auto(g2)
#mylay[,1]=c(tipp,timm,rep(0,length(domisp.mn)));mylay[,2]=c(0,pathd1, seq(min(pathd1),max(pathd1),length.out=length(domisp.mn))[order(timm.mn[,'pathd2'],decreasing = T)]);
mylay2[,1]=c(mylay[,1],rep(-2,length(domisp.mn)));
#mylay2[,2]=c(mylay[,2], seq(min(pathd1),max(pathd1),length.out=length(domisp.mn))[rank(-ind3.co[domisp.mn]+runif(length(domisp.mn),0,0.1))]);
mylay2[,2]=c(mylay[,2], seq(min(mylay[,2]),max(mylay[,2]),length.out=length(domisp.mn))[rank(-ind3.co[domisp.mn]+runif(length(domisp.mn),0,0.1))]);


names(ind3.co)=colnames(crnotu3)

vcol=ucol[c(1,ind3.co2[c(domisp,domisp.mn),'Abb']+1)]


domi.mvp.mn=filter(mv.tmp,ID%in%c(domisp.mn))
mvp.bi.mn=rep(0,length(domisp.mn));mvp.bi.mn[grep('1',domi.mvp.mn[,'Ttype'])]=1;mvp.bi.mn[grep('2',domi.mvp.mn[,'Ttype'])]=2
domi.taxa.mn= data.frame(filter(taxa.otu,OTU%in%domisp.mn),ind3.co2[c(domisp.mn[grep('OTU_',domisp.mn)]),])

domisp2=gsub('OTU_','',c(domisp,domisp.mn))
domisp3=c(domisp,domisp.mn)
domisp3[grep('OTU_',c(domisp,domisp.mn))]=rbind(domi.taxa,domi.taxa.mn)[,'fine.taxa']


V(g2)$color=vcol
V(g2)$label <- c('CH4',domisp3)
V(g2)$label.cex <- 0.7
V(g2)$size <- c(1,ccm.crnotu3[c(domisp,domisp.mn)])*8# 44
V(g2)$shape <- shape[c(2,mvp.bi, mvp.bi.mn)]


#E(g)$weight=abs(E(g)$weight)
E(g2)$arrow.size=0.3
E(g2)$width=c(rep(1.2, length(E(g)$weight)),rep(0.8, length(E(g2)$weight)-length(E(g)$weight)))

w.per=c(round(ecdf(E(g2)$weight)(E(g2)$weight)*100))
E(g2)$color=c(col.code[w.per[1:length(E(g)$weight)]+1],col.code2[w.per[-c(1:length(E(g)$weight))]+1])
#E(g2)$color=c(rep(colorRampPalette('white')(1),length(E(g)$weight)),col.code2[w.per[-c(1:length(E(g)$weight))]+1])


#windows(160,160);
postscript(paste('./Manuscript/Figure/Domino effect_mono_verc_',step.ch4,step.c,'.eps',sep=''), fonts=c("serif", "Palatino"))
plot(g2, layout=mylay2, vertex.label.family="serif",cex.lab=5,
     edge.label.family="Palatino", vertex.label.cex=0.7,edge.label.cex=2, xlab = "Warning time", ylab='')
legend('bottomright', indx, pch=rep(16,length(indx)), col=ucol, bty='n',cex=1)
axis(side=1,at=round(2*(seq(0,40,10)-min(mylay2[,1]))/(max(mylay2[,1])-min(mylay2[,1]))-1,2), labels =seq(0,40,10))
box()
dev.off()

domino.mls=data.frame(Type=as.matrix(ind3.co[domisp.mn]))
domino.mls[,1]=indx[-1][domino.mls[,1]]
arrange(domino.mls, Type)

