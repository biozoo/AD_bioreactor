#########################################topological comparison
wtopo=function(X,attri=NULL,mai=NULL, cluster=F){
  n=NROW(X)
  if(any(is.na(X))){
    topov=matrix(NA,NCOL(X),5)
    topo=rep(NA,3)
  }else{
    deig=eigen(X)$value
    deig=abs(deig[which.max(abs(deig))])
    diag(X)=0
    Xa=abs(X)
    x.d.uw=graph.adjacency(t(X), mode="directed",weighted =T)
    x.d.aw=graph.adjacency(t(Xa), mode="directed",weighted =T)
    
    if(cluster){
      if(any(X<0)){
        Xp=Xn=Xa=X
        Xp[X<0]=0
        Xn[X>0]=0
        kcx.a <- cluster_walktrap(x.d.aw)
        modu.a=modularity(kcx.a)
        modu=modu.a
      }else{
        kcx <- cluster_walktrap(x.d.uw)
        modu=modularity(kcx)
      }
    }else{modu=NA;ncc=NA}
    topov=cbind(strength.in=graph.strength(x.d.uw, mode = "in"),
                strength.out=graph.strength(x.d.uw, mode = "out"),
                strength.all=graph.strength(x.d.uw, mode = "all"),
                hub=hub.score(x.d.uw)$vector,
                auth=authority_score(x.d.uw)$vector
                )
    topo=c(
      sdstrength.all=sd(X,na.rm=T),
      modularity=modu,
      deig=deig
    )
    
  }
  
  return(list(topo=topo,topov=topov))
}

istype.node=function(ins,ous,node,exclude=NULL,lnumber=T){
  ins[node]=NA
  ous[node]=NA
  if(!is.null(exclude)){ins[exclude]=NA;ous[exclude]=NA}
  if(all(is.na(c(ins,ous)))){
    typeis <- rep(NA,ifelse(lnumber,17,9))
  }else{
    upis=ins
    lwis=ous
    pris=apply(abs(cbind(upis,lwis)), 1, mean)
    npair=length(upis)
    
    ind1=which(upis>0 & lwis>0)
    ind2=which(upis<0 & lwis<0)
    ind3=which(upis<0 & lwis>0)
    ind3i=which(upis>0 & lwis<0)
    ind4=which((upis>0 & lwis==0))
    ind4i=which((upis==0 & lwis>0))
    ind5=which((upis<0 & lwis==0))
    ind5i=which((upis==0 & lwis<0))
    #ind6=which(upis==0 & lwis==0)
    
    num.is=c(length(ind1),length(ind2),length(ind3),length(ind3i),length(ind4),length(ind4i),length(ind5),length(ind5i))#,length(ind6))
    num.isp=num.is/npair
    strength.is=c(sum(pris[ind1],na.rm=T),sum(pris[ind2],na.rm=T),sum(pris[ind3],na.rm=T),sum(pris[ind3i],na.rm=T),
                  sum(pris[ind4],na.rm=T),sum(pris[ind4i],na.rm=T),sum(pris[ind5],na.rm=T),sum(pris[ind5i],na.rm=T))/num.is
    strength.is[num.is==0]=0

    Interaction_Diversity=-sum(num.isp[num.isp>0]*log(num.isp[num.isp>0]))
    typeis=c(strength.is,Interaction_Diversity)
    names(typeis) <- c("Cooperation","Competition","Exploitated","Exploitator",
                       "Commensalism1","Commensalism0","Ammensalism-1","Ammensalism0","Interaction_Diversity")
  }
  return(typeis)
}


blockj=function(ism, member){
  umem=sort(unique(member))
  block0=matrix(0, length(umem), length(umem))
  colnames(block0)=rownames(block0)=umem
  block.d=block.da=block.dp=block.dn=block0
  ism.p=ism.n=ism.a=ism
  ism.a=abs(ism.a)
  ism.p[ism.p<0]=0
  ism.n[ism.n>0]=0
  
  for(i in 1:length(umem)){
    for(j in 1:length(umem)){
      block.d[i,j]=mean(ism[member==umem[i],member==umem[j]],na.rm=T)
      block.da[i,j]=mean(ism.a[member==umem[i],member==umem[j]],na.rm=T)
      block.dp[i,j]=mean(ism.p[member==umem[i],member==umem[j]],na.rm=T)
      block.dn[i,j]=mean(ism.n[member==umem[i],member==umem[j]],na.rm=T)
    }
  }
  block.ja=list(oj=block.d, aj=block.da, pj=block.dp, nj=block.dn)
  return(block.ja)
}


blockj.a=function(ism.a, member){
  n3=dim(ism.a)[3]
  nb=length(unique(member))
  ism.b=array(NA, dim=c(nb, nb, n3))
  ism.l=list();
  ism.l[['oj']]=ism.l[['aj']]=ism.l[['pj']]=ism.l[['nj']]=ism.b
  for(i in 1:n3){
    bk.i=blockj(ism.a[,,i],member)
    ism.l[['oj']][,,i]=bk.i$oj
    ism.l[['aj']][,,i]=bk.i$aj
    ism.l[['pj']][,,i]=bk.i$pj
    ism.l[['nj']][,,i]=bk.i$nj
  }
  return(ism.l)
}


block.prop=function(ism,member,exclude=NULL,itype=T,centrality=T){
  n2=dim(ism)[2]
  n3=dim(ism)[3]
  nam=sort(unique(member))
  nb=length(nam)
  topo.x=NULL
  for(kj in 1:n3){
    #print(kj)
    ism.at=ism[,,kj]
    if(all(is.na(ism.at))){}else{
      top.at=wtopo(ism.at,cluster = T)
      if(centrality){
        topov.j <- aggregate(top.at$topov,list(as.factor(member)),mean,na.rm=T)
      }else{topov.j <- NULL}
      if(itype){
        ist.k=NULL
        for(k in 1:n2){
          if(!is.null(exclude)){
            if(any(exclude==k)){ist.k=rbind(ist.k, rep(NA,9))}else{
              ist.k=rbind(ist.k, istype.node(ism.at[k,], ism.at[,k], k, exclude=exclude, lnumber=F))
            }
          }else{ist.k=rbind(ist.k, istype.node(ism.at[k,],ism.at[,k],k,lnumber=F))}
        }
        ist.j <- aggregate(ist.k,list(as.factor(member)),mean,na.rm=T)
      }else{ist.j <- NULL}
    }
    
    
    for(i in 1:nb){
      #print(i)
      ism.t=ism[member==nam[i],member==nam[i],kj]
      if(all(is.na(ism.t))){
        topo.x=rbind(topo.x,c(member=nam[i], time=kj, 
                              rep(NA,3+(NCOL(ist.j)-1)+(NCOL(topov.j)-1))))
      }else{
        top.t=wtopo(ism.t,cluster = T)
        x.t=c(member=nam[i], time=kj, 
              top.t[['topo']], 
              topov.j[i,-1], 
              ist.j[i,-1])
        topo.x=rbind(topo.x,unlist(x.t))
      }
    }# end of i
  }#end of kj
  colnames(topo.x)=c('member', 'time', names(top.t[['topo']]), 
                     colnames(topov.j)[-1],colnames(ist.j)[-1])
  return(topo.x)
}

blockins=function(bj.t){
  sign.t=c(-1,1,11,0)
  bkins=NULL
  for(i in 1:4){ 
    bj.i=bj.t[[i]]
    n.i=dim(bj.i)[1]
    for(j in 1:n.i){
      for(k in 1:n.i){
        bkins=rbind(bkins,data.frame(type='Guild_Int',sign=sign.t[i],cause=k,effect=j, time=tind, value=bj.i[j,k,]))
      }
    }
  }
  return(bkins)
}


FindCycles.node=function(v1, g, L){
  Cycles.i  = NULL
  if(degree(g, v1, mode="in") != 0) {
    GoodNeighbors = neighbors(g, v1, mode="out")
    GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
    for(v2 in GoodNeighbors) {
      TempCyc = lapply(all_simple_paths(g, v2, v1, mode="out", cutoff=L-1 ), function(p) c(v1,p))
      TempCyc = TempCyc[which(sapply(TempCyc, length) == L+1)]
      TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles.i  = c(Cycles.i, TempCyc)
    }
  }
  return(Cycles.i)
}

# Exploring all the feedback loop
FindCycles.lapply = function(g, L) {
  t.begin <- proc.time()
  if(!is.igraph(g)){
    diag(g)=0
    g=graph.adjacency(t(g), mode="directed",weighted =T)}
    vg=V(g)
  GenerateCycle=function(v1){
    Cycles.i=tibble(FindCycles.node(v1, g, L))
    return(Cycles.i)
  }
  
  ptm <- proc.time()  
  output <- lapply(vg, GenerateCycle)
  output <- do.call(bind_rows, output)
  output <- do.call(rbind, output[[1]])
  lapplyTime <- proc.time() - ptm  
  cat(lapplyTime)
  
  return(output)
}

# Feedback analysis
fdwf=function(ism, fdm, L=3, memm=F,ind.t=NULL){
  fdw=NULL
  if(memm){fdmm=fdm[ind.t==1,]}else{fdmm=fdm}
  for(i in 1:dim(ism)[3]){
    ism.i=ism[,,i]
    if(L==3){fdw.i=c(ism.i[cbind(fdmm[,2],fdmm[,1])]*ism.i[cbind(fdmm[,3],fdmm[,2])]*ism.i[cbind(fdmm[,1],fdmm[,3])])}
    fdw.i=sign(fdw.i)*abs(fdw.i)^(1/L)
    if(all(is.na(fdw.i))){fdw=c(fdw,NA)}else{fdw=c(fdw,sum(fdw.i,na.rm=T))}
    rm(fdw.i);gc()
  }
  return(fdw)
}

datransf=function(topo.o, ind.t=c(1:2)){
  topp.n=topo.o[,ind.t]
  n.t=nrow(topo.o)
  topp.t=NULL
  for(i in (1:ncol(topo.o))[-ind.t]){
    topp.t=rbind(topp.t,data.frame(topp.n,variable=rep(colnames(topo.o)[i],n.t),value=topo.o[,i]))
  }
  rownames(topp.t)=NULL
  return(topp.t)
}

# bin data vector according to the membership
databin=function(dada,member,FUN='mean'){
  fun.t=match.fun(FUN)
  umem=sort(unique(member))
  nmem=length(umem)
  dabin=NULL
  for(i in 1:length(umem)){
    dabin=cbind(dabin,apply(dada[,member==umem[i]],1,fun.t,na.rm=T))
  }
  colnames(dabin)=umem
  return(dabin)
}

repmat = function(X,m=1,n=1){
  ##R equivalent of repmat (matlab)
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

####################################
## interpolation by cubic spline
cinter <- function(x,extra=F,fit.out=F){
  x=as.matrix(x)
  n=nrow(x)
  y=x
  y.fit=NULL
  for(i in 1:NCOL(x)){
    x.t=x[,i]
    id=which(is.na(x.t))
    id2=nacap(x.t)
    xtt=cbind(1:n,x.t)
    xtt=xtt[!apply(is.na(xtt),1,any),]
    inp=fields::splint(xtt[,1],xtt[,2],1:n)
    inp2=inp;
    if(any(id2==10^9)){}else{inp2[id2]=NA}
    y.fit=cbind(y.fit,inp2)
    if(length(id)!=0){
      ####exclude extrapolation 
      if(extra==F){id=setdiff(id,id2)}
      if(length(id)!=0){x.t[id]=inp[id];y[,i]=x.t}      
    }else{y[,i]=x.t}
  }
  if(fit.out){list(y=y,y.fit=y.fit)}else{return(y)} 
}

####################################
## Find NA at begiening or end
nacap=function(x,lists=F){
  n=length(x)
  ind=which(is.na(x))
  bna=ena=NULL
  if(any(ind==1)){
    bna=1
    for(i in 2:n){
      if(any(ind==i)&any(ind==(i-1))){bna=c(bna,i)}else{break}
    }
  }else{bna=NULL}
  
  if(any(ind==n)){
    ena=n
    for(i in (n-1):1){
      if(any(ind==i)&any(ind==(i+1))){ena=c(ena,i)}else{break}
    }
  }else{ena=NULL}
  naind=list(bna,ena)
  naind.f=naind
  if(length(c(bna,ena))==0){naind.f=10^9}
  if(lists){return(naind)}else{return(sort(unlist(naind.f)))}
}
