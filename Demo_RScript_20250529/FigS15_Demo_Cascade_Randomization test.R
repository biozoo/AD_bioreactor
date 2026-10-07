# R code for plotting Fig. 5
rm(list = ls());gc()
# Set path
OPath <- 'C:/data/Cooperation/Prodigio/Code/Demo_RScript_20250529'
setwd(OPath)

library(factoextra)
library(igraph)
library(mgcv)
library(gratia)
library(dplyr)
library(fields)
library(ggplot2)
library(patchwork)

seed <- 49563
set.seed(seed)
SaveFile <- F #T/F for saving files

tipp <- 38+1
tipp2 <- 44+1
nd <- 46+1
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

mag_scale <- function(x,mgf=0.5){
  if(prod(x)>=0){
    if(x[2]>0){y <- x;y[1] <- y[1]*(1.0-mgf);y[2]*(1.0+mgf)}else{
      y <- x;y[1] <- y[1]*(1.0+mgf);y[2]*(1.0-mgf)}
  }else{y <- x*(1.0+mgf)}
  return(y)
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
mco3.ins <- read.csv('./Output/Demo_CH4_otucrn_Perturbed_D47_NA_nin47_cvunit0.025_enet_jcof_Nmvx_Rallx_demo_NEW.csv',stringsAsFactors = F)[,-c(1:5)]

co1.ins <- read.csv('./Output/Control_otucrn_CH4_D47_Species_NA_cvunit0.025_enet_jcof_Nmvx_Rallx_demo.csv',stringsAsFactors = F)
co3.ins <- read.csv('./Output/Perturbed_otucrn_CH4_D47_Species_NA_cvunit0.025_enet_jcof_Nmvx_Rallx_demo.csv',stringsAsFactors = F)

co.ary3 <- array(NA,dim=c(nco3,nco3,nd-1))
for(i in 1:(nd-1)){
  ints.i31 <- as.matrix(filter(co3.ins,time==i)[,-c(1:4)])
  co.ary3[,,i] <- ints.i31
}

################################################################################
### CCM results
ccm.otu <- read.csv('./Output/Demo_ccm_OTU_CH4_Species_D47.csv', header=T)
ccm.crn <- read.csv('./Output/Demo_ccm_metabolite_CH4_D47.csv', header=T)

ccm.otu.c <- ccm.otu[,'CCM_Control']
ccm.crn.c <- ccm.crn[,'CCM_Control']

ccm.otu <- ccm.otu[,'CCM_Perturbed']
ccm.crn <- ccm.crn[,'CCM_Perturbed']

ccm.crnotu <- c(ccm.crn,ccm.otu)
names(ccm.crnotu) <- colnames(crnotu3)

################################################################################
# Contributions for each species to the CH4 critical transition at the tipping point
da.i <- data.frame(value=envi[,'CH4_O1'],time=1:nd)
g.t <- gam(value~s(time),data=da.i)
newd.i <- seq(1,nd,0.01)
yhat <- predict.gam(g.t,newdata=data.frame(time=newd.i))
yhatd1 <- derivatives(g.t, type = "central",data=data.frame(time=newd.i), order=1)$derivative

# Contribution of each node on the state changes in CH4 based on the dynamical interactome and time series of state variables 
impc <- mco3.ins*crnotu3s[-nd,]
cotr <- rbind(rep(NA),impc[-1,]-impc[-(nd-1),])
# Sign of the local trend (t, t+1) defined as the estimated first derivative evaluated at t
sdch4s <- sign(yhatd1[newd.i%in%(seq(1,nd-1,1))])
# Examine whether the sign of local trend equals to the sign of contribution  
drcotrs <- sign(cotr)*sdch4s
# Corrected contributions depending on whether the sign of contribution facilitate(+) or buffer (-) the observed trend in CH4
cotrds <- drcotrs*abs(cotr)

# Cumulative contribution
ind <- (tipp):(tipp2-1)#(tipp+1)#
rg.con <- 1:2#1#
nrg <- length(rg.con)
cotrdr0 <- cotrds[ind,];
ind.r <- which(apply(matrix(as.matrix(cotrdr0[rg.con,]),nrow=nrg)>0,2,all))#cottar(apply(cotrdr0[1:2,]>0,2,all))#ind.r <- which(cotrdr0[1,]>0)
cotrdr0 <- cotrdr0[,ind.r];rownames(cotrdr0)=paste('D',ind,sep='-')
cotr0 <- apply(matrix(as.matrix(cotrdr0[rg.con,]),nrow=nrg),2,sum)#cotr0 <- cotrdr0[1,]#
names(cotr0) <- colnames(cotrdr0)
# Examine the contributions for all the state changes among network nodes. 
# This is a similar process repeating the computation presented for CH4 

# smoothed trend estimate for all the variables 
smo.crnotu3 <- read.csv('./Output/smo_otucrn3_D47.csv', header = T)

crnotu3.smo <- matrix(smo.crnotu3[,'yhat'],nrow=tipp)
crnotu3.smod <- matrix(smo.crnotu3[,'est'],nrow=tipp)

cotrd.ary <- cotrds.ary <- array(NA, dim=c(NCOL(crnotu3),NCOL(crnotu3),tipp))
for(i in 1:NCOL(crnotu3)){
  impc.i <- t(co.ary3[i,,])*crnotu3s[-nd,]
  cotr.i <- rbind(rep(NA),impc.i[-1,]-impc.i[-(nd-1),])
  cotr.i <- cotr.i[1:tipp,]
  sdxis <- sign(crnotu3.smod[,i])
  drcotrs.i <- sign(cotr.i)*sdxis
  cotrds.i <- drcotrs.i*abs(cotr.i)
  cotrds.ary[i,,] <- t(cotrds.i)
}

mv.tmp <- read.csv('./Output/Localized_Transition_crnotu_D47.csv',header = T)
warn.crnotu <- read.csv('./Output/warn_time_CRNOTU_D47.csv',header = T)

rownames(mv.tmp) <- mv.tmp[,'ID']

domino <- NULL
for(i in 1:NCOL(crnotu3)){
  mv.tmp.t0 <- filter(mv.tmp, ID==colnames(crnotu3)[i])
  cht.t <- as.numeric(mv.tmp.t0[c('warn_time')])
  tarid <- as.character(mv.tmp.t0['ID'])
  
  if(mv.tmp.t0['Ttype']%in%c('D1','D2','B1','B2')){
    ind.r <- as.numeric(ceiling(cht.t)):as.numeric((nrg-1)+ceiling(cht.t)); 
    indi <- which(apply(matrix(as.matrix(cotrds.ary[i,,ind.r]),ncol=nrg)>0,1,all))
    mv.tmp.t2 <- mv.tmp[indi,]
    mv.tmp.t2 <- mv.tmp.t2%>%mutate(contribution=apply(matrix(as.matrix(cotrds.ary[i,indi,ind.r]),ncol=nrg),1,sum))%>%
      filter(Ttype%in%c('B1','B2','D1','D2'))%>%arrange(desc(warn_time))%>%filter((warn_time)<trunc(cht.t))
    tarid2 <- mv.tmp.t2[,'ID']
    contri <- mv.tmp.t2[,'contribution']
  }else{tarid2 <- NA; contri <- NA}
  if(length(tarid2)==0){tarid2 <- NA; contri <- NA}
  domino <- rbind(domino, data.frame(Variable=rep(i),ID=rep(tarid),Ttype=rep(mv.tmp.t0['Ttype']),Time=rep(cht.t),Domino=tarid2,Contribution=contri))
}

domino <- filter(domino, ID!=Domino)
domino <- filter(domino, ID%in%warn.crnotu[,1] & Domino%in%warn.crnotu[,1])
dim(domino)

ch4.contr <- mv.tmp%>%filter(ID%in%names(cotr0),ID%in%warn.crnotu[,1])%>%arrange(desc(warn_time))
domino <- rbind(data.frame(Variable=0,ID='CH4',Ttype='',Time=tipp,Domino=ch4.contr[,'ID'],
                           Contribution=as.numeric(cotr0[ch4.contr[,'ID']])),domino)

################################################################################
## Reconstruct various paths involved in the cascading mechanism
dp.ls <- list();dp.ls[[1]] <- cbind(ind=rep('CH4'), ID=ch4.contr[,'ID'])
for(i in 2:30){
  dpath.i <- dp.ls[[i-1]][,2]
  ud.i <- unique(dpath.i)
  ud.i <- ud.i[!is.na(ud.i)]
  dj2 <- NULL
  for(j in 1:length(ud.i)){
    dj <- filter(domino,ID==ud.i[j])[,'Domino']
    if(length(dj)!=0 & all(!is.na(dj))){
      dj2 <- rbind(dj2, cbind(ind=rep(ud.i[j],each=length(dj)),ID=rep(dj,length(ud.i[j]))))
    }else{
      dj2 <- rbind(dj2, cbind(ind=rep(ud.i[j],each=1),ID=rep(NA,1)))
    }
  }
  dp.ls[[i]] <- dj2;print(paste(i,NROW(dj2)))
  if(all(is.na(dj2[,2]))){break}else{rm(dj2);gc()}
}

# All individual cascading paths 
dpi1 <- dp.ls[[1]]
for(i in 2:(length(dp.ls))){
  dpi <- as.matrix(dp.ls[[i]])
  ud.t <- unique(dpi[,1])
  path.i <- NULL
  for(j in 1:length(ud.t)){
    dc.j <- dpi[dpi[,1]==ud.t[j],2]
    dc.j1 <- matrix(dpi1[which(dpi1[,NCOL(dpi1)]==ud.t[j]),],ncol=i)
    path.i <- rbind(path.i, cbind(repmat(matrix(dc.j1,ncol=i),length(dc.j),1),as.character(rep(dc.j,each=NROW(dc.j1)))))
  }  
  colnames(path.i)=NULL
  na.ind <- which(is.na(dpi1[,NCOL(dpi1)]))
  
  if(length(na.ind)>0) {
    dpi1[na.ind,NCOL(dpi1)] <- NA
    path.na <- cbind(dpi1[na.ind,],rep(NA,length(na.ind)))
    colnames(path.na) <- colnames(path.i)
    path.i <- rbind(path.i, path.na)
  }
  path.i <- path.i
  dpi1 <- path.i
}
dpi1 <- dpi1[,-NCOL(dpi1)]

# Network nodes involved in the cascading process
domisp <- intersect(colnames(crnotu3),unique(c(dpi1[,-1])))

# Timing for localized transition during the cascading process
dptime <- matrix(NA,NROW(dpi1),NCOL(dpi1))
dptime[,1] <- tipp
timm <- NULL
for(i in 1:length(domisp)){
  mv.i <- filter(mv.tmp,ID==domisp[i])
  timm <- c(timm, as.numeric(mv.i['warn_time']))
  dptime[which(dpi1==domisp[i])] <- timm[i]
}
names(timm) <- domisp
(timm)

# Extract each link during the cascading process
domi.link <- NULL
domi.netw <- matrix(0,length(domisp)+1,length(domisp)+1);
colnames(domi.netw) <- rownames(domi.netw) <- c('CH4', domisp)
for(i in 2:NCOL(dpi1)){
  dpi.i <- dpi1[,c(i,i-1)]
  dpi.i <- dpi.i[apply(!is.na(dpi.i),1,all),]
  if(!is.matrix(dpi.i)){dpi.i=matrix(dpi.i,1,2)}
  ub.i <- apply(dpi.i,1,function(x){y=paste(x,collapse = '->');return(y)})
  domi.link <- c(domi.link, ub.i)
  for(j in 1:NROW(dpi.i)){
    domino.j <- filter(domino, Domino==dpi.i[j,1],ID==dpi.i[j,2])
    domi.netw[dpi.i[j,2], dpi.i[j,1]] <- as.numeric(domino.j['Contribution'])
  }
}
(domi.link <- unique(domi.link))

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

mylay <- cbind(c(tipp,timm)-1, c(0, cen.i))

ucol <- c(CH4='#FF0000', VFA='#DBB100', Fingerprint='#66A61E',
          Methanogen='#E7298A', Acetogen="#7570B3", Other_taxa="grey80")
vcol <- c(ucol[1], ucol[mv.tmp[domisp,'Type']])

domisp2 <- c('CH4',domisp);
domisp2[grep('OTU_',domisp2)] <- clasdo[domisp2[grep('OTU_',domisp2)],'Finest_taxon']

V(g)$color <- vcol
V(g)$label <- domisp2
V(g)$label.cex <- 0.7
V(g)$size <- c(1, ccm.crnotu[domisp])*10 # Node size is scaled according to the strength of causal influences on CH4
V(g)$shape <- shape.bi

E(g)$arrow.size=0.3
# Edge color presented the quantitative contribution of one node to the localized transition of the other nodes
colfunc <- colorRampPalette(c("grey90", "black"))
col.code=colfunc(102)[-1]
w.per=round(ecdf(E(g)$weight)(E(g)$weight)*100)+1
E(g)$color=col.code[w.per]

#win.graph(100,60)

#postscript(paste('D:/data/Cooperation/Prodigio/Manuscript/Figure/D47/Cascading diagram.eps',sep=''), fonts=c("serif", "Palatino"))
windows(160,160);
n <- 100
x_seq <- seq(0, tipp2-1, length.out = tipp2)
y_seq <- seq(0, tipp2-1, length.out = tipp2)
z <- repmat(c(envi[(1:tipp),'CODt_O1'],rep(NA,tipp2-tipp)),1,tipp2)   # 用 x 作為漸層依據
image(x_seq, y_seq, z,
      col = colorRampPalette(c("white", "#fed5b1"))(n),
      useRaster = TRUE,xlim=c(-4.5,43),axes=F,
      xlab = "", ylab = "")
par(new=T)

plot(g, layout=mylay, vertex.label.family="serif",cex.lab=5,edge.label.family="Palatino",
     vertex.label.cex=0.7,edge.label.cex=2, xlab = "Transition time", ylab='')
legend('topright', names(ucol), pch=rep(16,length(ucol)), col=ucol, bty='n',cex=1)
axis(side=1,at=round(2*(seq(0,40,10)-min(mylay[,1]))/(max(mylay[,1])-min(mylay[,1]))-1,2), labels =seq(0,40,10))

rescale <- function(x){y <- round(2*(x-min(mylay[,1]))/(max(mylay[,1])-min(mylay[,1]))-1,2);return(y)}
polygon(x=rescale(c(1,tipp,tipp2)-1)[c(2:3,3:2)],y=rep(c(-1.1,-0.9),each=2), col="orange",density=15, border="white")

#Display the graph and the legend.
ordered <- order(col.code)
color <- vector(length = length(ordered),mode="double")
for (i in 1:length(ordered)) color[ ordered[i] ] <- col.code [ i ]

image.plot(legend.only=T, zlim=range(E(g)$weight), col=col.code, 
           legend.shrink= 0.25, legend.cex=0.75, smallplot=c(0.025,0.05,0.35,0.5))
#dev.off()

########################################################
### Cumulative effect of cascading process
narrow <- apply(domi.netw,1,sum)

warn_cascade <- mv.tmp[names(narrow),"warn_time"]
warn_cascade[1] <- tipp

cumu.data <- data.frame(warn_cascade, narrow)%>%
  arrange(warn_cascade)%>% mutate (arrow_cum= cumsum(narrow))
cumu.data <- rbind(c(1,0,0,0),cumu.data)  
cumu.data[,-1] <- (cumu.data[,-1])/repmat(cumu.data[NROW(cumu.data),-1],NROW(cumu.data))

narrow <- apply(domi.netw!=0,1,sum)

warn_cascade <- mv.tmp[names(narrow),"warn_time"]
warn_cascade[1] <- tipp

# cumulative effect of cascading diagram (Fig. 5C)
win.graph(100,50)
par(mar=c(4,4,1,4))
plot(envi[1:(tipp+3),'COD_rm_O1']~c((1:(tipp+3))-1),ylim=range(envi[1:(tipp+2),c('COD_rm_O1','COD_rm_O2')]),
     xlab='', ylab='tCOD removal efficiency (%)', cex.lab=1.5,pch=1,cex=2,
     lwd=2,xlim=c(0,tipp2-1),type='o',axes=F)
lines(envi[1:(tipp+3),'COD_rm_O2']~c((1:(tipp+3))-1),
      col='grey',pch=1,cex=2,lwd=2,type='o')
axis(1);axis(2)

par(new=T)
plot(I(arrow_cum)~I(warn_cascade-1), cumu.data,xlim=c(0,tipp2-1),
     type='o',cex.lab=1.5,pch=16,cex=2,lty=3,axes=F, 
     lwd=2,xlab='Day',ylab='')
#lines(I(arrow_mn_cum)~I(warn_cascade-1), cumu.data,xlim=c(0,tipp2-1),
#      type='o',cex.lab=1.5,pch=16,cex=2,col='grey',lwd=2)
axis(4,col.axis = 1)
mtext('Cumulative effect',side=4, col=1,line = 3, cex=1.5)


cor.obs <- cor(cumu.data[,'arrow_cum'],envi[1:NROW(cumu.data),'COD_rm_O1'])

################################################################################
## Randomly rewire the edges in cascade diagram while preserving node transition times
## Node transition times (warn_cascade) remain unchanged.
## Only the influence values (edges) are randomly reassigned to new node pairs,
## breaking the original "who influences whom" relationships.
################################################################################

b.num <- 1000
cor.netw <- numeric(b.num)

domi.z <- matrix(0, nrow = nrow(domi.netw), ncol = ncol(domi.netw))
rownames(domi.z) <- rownames(domi.netw)
colnames(domi.z) <- colnames(domi.netw)
for (i in 1:b.num) {
  
  # Get positions and values of all non-zero edges
  non_zero_idx <- which(domi.netw != 0, arr.ind = TRUE)
  non_zero_vals <- domi.netw[non_zero_idx]
  
  # Shuffle edge values
  shuffled_vals <- sample(non_zero_vals)
  
  # Create zero matrix
  domi_shuffled <- domi.z

  # Rewire edges to random node pairs
  shuffled_pos <- non_zero_idx[sample(1:nrow(non_zero_idx)), ]
  for (j in 1:length(shuffled_vals)) {
    domi_shuffled[shuffled_pos[j, 1], shuffled_pos[j, 2]] <- shuffled_vals[j]
  }
  
  # Recalculate cumulative cascade (node transition times unchanged)
  narrow_shuff <- apply(domi_shuffled, 1, sum)
  
  cumu_shuff <- data.frame(warn = warn_cascade, narrow = narrow_shuff) %>%
    arrange(warn) %>%
    mutate(arrow_cum = cumsum(narrow_shuff))
  
  cumu_shuff <- rbind(data.frame(warn = 1, narrow = 0, arrow_cum = 0), cumu_shuff)
  cumu_shuff$arrow_cum <- cumu_shuff$arrow_cum / cumu_shuff$arrow_cum[nrow(cumu_shuff)]
  
  cascade_interp_shuff <- approx(x = cumu_shuff$warn - 1,
                                 y = cumu_shuff$arrow_cum,
                                 xout = 0:(tipp2 - 1),
                                 rule = 2)$y
  
  cor.netw[i] <- cor(cascade_interp_shuff, envi[1:tipp2, 'COD_rm_O1'], use = "complete.obs")
}


################################################################################
## Randomly rewire edges while preserving temporal directionality:
## Only nodes that transition earlier can influence nodes that transition later
##
## NOTE: warn_cascade corresponds to rownames(domi.netw) in the SAME ORDER
################################################################################

b.num <- 1000
cor.netw.temporal <- numeric(b.num)

domi.z <- matrix(0, nrow = nrow(domi.netw), ncol = ncol(domi.netw))
rownames(domi.z) <- rownames(domi.netw)
colnames(domi.z) <- colnames(domi.netw)

# warn_cascade corresponds to rownames(domi.netw) in the same order
node_times <- warn_cascade
names(node_times) <- rownames(domi.netw)

# CH4 is fixed at tipp (methane collapse time)
node_times[1] <- tipp

# Find valid edges: influencer (col) must transition BEFORE recipient (row)
non_zero_idx <- which(domi.netw != 0, arr.ind = TRUE)
valid_idx_list <- list()

for (k in 1:nrow(non_zero_idx)) {
  r <- non_zero_idx[k, 1]
  c <- non_zero_idx[k, 2]
  
  if (node_times[c] < node_times[r]) {
    valid_idx_list[[length(valid_idx_list) + 1]] <- c(r, c)
  }
}

valid_idx <- do.call(rbind, valid_idx_list)
colnames(valid_idx) <- c("row", "col")

cat("Valid edges:", nrow(valid_idx), "/", nrow(non_zero_idx), "\n")

if (nrow(valid_idx) == 0) {
  stop("No valid edges found! Check warn_cascade values.")
}

valid_vals <- domi.netw[valid_idx]

for (i in 1:b.num) {
  
  shuffled_vals <- sample(valid_vals)
  
  domi_shuffled <- domi.z
  shuffled_pos <- valid_idx[sample(1:nrow(valid_idx)), ]
  
  for (j in 1:length(shuffled_vals)) {
    domi_shuffled[shuffled_pos[j, 1], shuffled_pos[j, 2]] <- shuffled_vals[j]
  }
  
  narrow_shuff <- apply(domi_shuffled, 1, sum)
  
  cumu_shuff <- data.frame(warn = warn_cascade, narrow = narrow_shuff) %>%
    arrange(warn) %>%
    mutate(arrow_cum = cumsum(narrow_shuff))
  
  cumu_shuff <- rbind(data.frame(warn = 1, narrow = 0, arrow_cum = 0), cumu_shuff)
  cumu_shuff$arrow_cum <- cumu_shuff$arrow_cum / cumu_shuff$arrow_cum[nrow(cumu_shuff)]
  
  cor.netw.temporal[i] <- cor(cumu_shuff[, 'arrow_cum'],
                              envi[1:NROW(cumu_shuff), 'COD_rm_O1'],
                              use = "complete.obs")
}

cor.netw.temporal <- cor.netw.temporal[!is.na(cor.netw.temporal)]


################################################################################
## Combined visualization: Two null distributions
################################################################################

# Create data frames
df_null_rewire <- data.frame(correlation = cor.netw, type = "Edge rewiring")
df_null_temporal <- data.frame(correlation = cor.netw.temporal, type = "Temporally-constrained")

# Combine
df_all <- rbind(df_null_rewire, df_null_temporal)

# Compute p-values
p_rewire <- which(abs(cor.netw) < abs(cor.obs))[1] / length(cor.netw)
p_temporal <- which(abs(cor.netw.temporal) < abs(cor.obs))[1] / length(cor.netw.temporal)

# Create labels
label_rewire <- paste0("Observed r = ", round(cor.obs, 3), "\np = ", round(p_rewire, 4))
label_temporal <- paste0("Observed r = ", round(cor.obs, 3), "\np = ", round(p_temporal, 4))

# Plot 1: Edge rewiring
p1 <- ggplot(df_null_rewire, aes(x = correlation)) +
  geom_histogram(bins = 50, 
                 fill = "lightblue", 
                 color = "white",
                 alpha = 0.8) +
  geom_vline(xintercept = cor.obs, 
             color = "red", 
             linewidth = 1.2, 
             linetype = "dashed") +
  annotate("text", 
           x = -Inf, 
           y = Inf,
           label = label_rewire,
           color = "red",
           size = 3.5,
           hjust = -0.1,
           vjust = 1.5) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(min(cor.obs, cor.netw) - 0.05,
                           max(cor.obs, cor.netw) + 0.05)) +
  labs(title = "A",
       subtitle = "Edge rewiring",
       x = "Correlation",
       y = "Frequency") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.6))

# Plot 2: Temporally-constrained rewiring
p2 <- ggplot(df_null_temporal, aes(x = correlation)) +
  geom_histogram(bins = 50, 
                 fill = "lightgreen", 
                 color = "white",
                 alpha = 0.8) +
  geom_vline(xintercept = cor.obs, 
             color = "red", 
             linewidth = 1.2, 
             linetype = "dashed") +
  annotate("text", 
           x = -Inf, 
           y = Inf,
           label = label_temporal,
           color = "red",
           size = 3.5,
           hjust = -0.1,
           vjust = 1.5) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(min(cor.obs, cor.netw.temporal) - 0.05,
                           max(cor.obs, cor.netw.temporal) + 0.05)) +
  labs(title = "B",
       subtitle = "Temporally-constrained rewiring",
       x = "Correlation",
       y = "Frequency") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.line = element_line(color = "black", linewidth = 0.6),
        axis.ticks = element_line(color = "black", linewidth = 0.6))

# Combine side by side

combined_plot <- p1 + p2 + plot_annotation(
  title = "Null distributions",
  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
)
win.graph(100,50)
combined_plot



ggsave(
  filename = "C:/data/Cooperation/Prodigio/Manuscript/Figure/D47/Cascade_null.pdf",
  plot     = combined_plot,
  width    = 18,
  height   = 7,
  units    = "in",      # 
  device   = cairo_pdf   # 
)

