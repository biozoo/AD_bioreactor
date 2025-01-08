
## Additional function tools
##Making a matrix with repeated column or row 
repmat <- function(X,m=1,n=1){
  if(!is.null(dim(X))){X=as.matrix(X)}else{X=as.numeric(X)}
  if(is.matrix(X)){
    mx <- dim(X)[1]
    nx <- dim(X)[2]
    y <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }else if(m>1&n==1) {
    y <- matrix(rep(X,each=m),m,length(X))
    
  }else if(m==1&n>1){
    y <- matrix(rep(X,n),length(X),n)
  }
  return(y)
}

nomz=function(x, normalization=T, dseason=T, season_sd=T, sea=12, dtrend=T, dTtype="linear"){
  x=as.numeric(x)
  xt=x
  # Detrend
  if(dtrend==T & dTtype=="first"){xt=diff(xt)} else if (dtrend==T & dTtype=="linear"){
    lm.t=lm(xt~c(1:length(xt)))
    xt=xt-(lm.t$coefficients[1]+lm.t$coefficients[2]*c(1:length(xt)))}
  # Deseason
  if(dseason==T){
    xs=as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,mean,na.rm=T))
    xsd=as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,sd,na.rm=T))
    xt=xt-c(rep(xs,1+length(xt)%/%sea))[1:length(xt)]
    if(season_sd==T){xt=xt/(c(rep(xsd,1+length(xt)%/%sea))[1:length(xt)])}}
  # Normalization (zero mean & unity variance)
  if(normalization==T){xt=(xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)}
  return(xt)
}


####### Create a series of lag series 
lgmm <- function(x,E,nm="y"){
  nam <- NULL
  nmt <- paste(nm,"_t",sep="")
  if(E==0){lagm <- x; names(lagm) <- nmt}else{
    n <- length(x)
    lagm <- matrix(NA,nrow=n,ncol=E)
    for(i in 1:E){
      lagm[(i:n),i] <- x[1:(n-i+1)]
    }
    nam <- paste(nmt,"-",c(0:(E-1)),sep="")
    nam[1] <- nmt
    colnames(lagm) <- nam
  }  
  return(lagm)
}

logspace <- function(d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n)) 

#######################
# split a parameter vector into parts with roughly equal size
eqsplit <- function(x,nub){
  n <- length(x)
  nl <- trunc(n/nub)
  exs <- n-nl*nub
  nlseq <- rep(nl,nub)
  if(exs>0){nlseq[1:exs] <- nlseq[1:exs]+1}
  nlcu <- c(0,cumsum(nlseq))
  spn <- NULL
  for(i in 2:(nub+1)){
    spn <- rbind(spn, c(nlcu[i-1]+1,nlcu[i]))
  }
  return(spn)
}


ccm.fast.inferM <- function(ds,inferM=NULL,Epair=T,cri='rmse',Emax=10, deseason=F, 
                           period=12, silience=F, singleY=NULL,Lmin=10,Lint=10){
  if(is.null(singleY)&is.null(inferM)){inferM=matrix(1,ncol(ds),ncol(ds))}
  if(!is.null(singleY)){inferM=matrix(1,1,ncol(ds));inferM[1,singleY]=0}
  if(cri=='rho'){jfun <- match.fun('which.max')}else{jfun <- match.fun('which.min')}
  if(deseason){
    apply(ds,2,nomz,dTtype="linear",season_sd=T,sea=period)
  }else{
    ds=as.matrix(apply(ds,2,scale))
  }
  np <- nrow(ds) # time series length
  ns <- ncol(ds) # number of nodes
  lib.s <- c(seq(Lmin,nrow(ds),Lint),nrow(ds)) # sequence of library size 
  crirho <- qt(0.95,np-2)/sqrt(np-2+qt(0.95,np-2)^2) # critical values with alpha=0.05
  if(is.null(singleY)){ccm.rho <- ccm.sig <- matrix(0,ns,ns);ns1=ns
  }else{ccm.rho <- ccm.sig <- matrix(0,1,ns);ns1=1}
  for(i in 1:ns1){
    t.begin <- proc.time()
    for(j in 1:ns){
      if(inferM[i,j]==0|all(is.na(ds[,i]+ds[,j]))){}else{
        # select the optimal E for CCM based on hindcast at time= t-1 (tp=-1)
        ccm.E <- NULL  
        for(E.t in 2:Emax){
          ccm.E <- rbind(ccm.E,ccm(cbind(x=ds[,i],y=ds[,j]), E = E.t, tp=-1,
                                   lib_column = "x", target_column = "y", 
                                   lib_sizes = nrow(ds),  random_libs =F))
        }
        Eop <- ccm.E[jfun(ccm.E[,cri]),'E'] # The optimal E for the cross-mapping from node i to node j 
        
        
        # Perform CCM at time t (tp=0)      
        ccm.out <- ccm(cbind(x=ds[,i],y=ds[,j]), E = Eop, tp=0, 
                       lib_column = "x", target_column = "y", 
                       lib_sizes = lib.s,  random_libs =F)
        if(all(is.na(ccm.out[,'rho']))){}else{
          # aggregate the results with respect to each library size
          ccm.seq <- aggregate(ccm.out[,'rho'],list(ccm.out[,'lib_size']),mean,na.rm=T)
          ccm.seq <- ccm.seq[!(is.na(ccm.seq[,2])|is.infinite(ccm.seq[,2])),]
          ccm.seq[ccm.seq[,2]<0,2] <- 0
          termrho <- ccm.seq[nrow(ccm.seq),2]  # rho at the maximal library size (terminal rho)
          if(nrow(ccm.seq)>=3){
            kend <- MannKendall(ccm.seq[,2]);  # Kendall's tau test for mononotic increase
            # Causation is significant only if both (1) Kendall's tau and (2) terminal rho are significantly larger than zero
            ccm.sig[i,j] <- (kend$tau[1]>0)*(kend$sl[1]<0.05)*(termrho>crirho) # ccm.sig records the significance of each CCM
          }else{ccm.sig[i,j] <- 0}
          ccm.rho[i,j] <- termrho 
        }
      }# end of else
                                            # ccm.rho records the terminal rho
    }#end of j
    time.used <- proc.time() - t.begin 
    if(!silience){cat("variable", i, "ccm completed:", time.used[3],"sec\n")}
  }#end of i
  return(list(ccm.rho=ccm.rho,ccm.sig=ccm.sig))  
}

###############################################################################
## Embedding multiview SSR based on causal variable selected by CCM 
esim.lag.demo2 <- function(ds,ccm.rho,ccm.sig,inferM=NULL,Ed,kmax=10000,kn=100,max_lag=3,Emax=Emax){
  # ds= dataset
  # ccm.rho= matrix of CCM terminal rho
  # ccm.sig= matrix recoding the significance of CCM between each node pair
  # Ed= the optimal embedding dimension for each variable
  # max_lag = the maximal time lag included in multiview embedding
  # kmax= The maximal number for generating multiview SSR
  # kn= Select the kn best multiview SSR
  # Emax=maximal embedding dimension
  
  # create dataset with lags 
  dst <- NULL
  for(k in 1:ncol(ds)){dst <- cbind(dst,lgmm(ds[,k],E=max_lag,nm=k))}
  mx.na <- matrix(NA,kn,Emax); colnames(mx.na) <- paste('X',1:Emax,sep='_') 
  mx.na2 <- matrix(NA,kn,Emax); colnames(mx.na2) <- paste('vlag',1:Emax,sep='_')
  if(is.null(inferM)){inferM=matrix(1,ncol(ds),ncol(ds))}
  esele <- NULL
  for(i in 1:ncol(ds)){   # Multiview embedding for each node
    t.begin <- proc.time()
    if(all(inferM[i,]==0)){}else{
      ccm.t <- (ccm.rho*ccm.sig)[i,]
      cind <- which(ccm.t>0)   # Select significant causal variables
      num_vars <- length(cind) # Number of causal variable
      E.i <- Ed[i]             # The optimal embedding dimension
      # Preparing subset index for all lag terms of causal variables in the dataset, dst
      lagind <- NULL;for(j in 1:max_lag){lagind <- c(lagind,(cind-1)*max_lag+j)} 
      lagind <- sort(lagind)
      tarind <- (i-1)*max_lag+1   # index of the target node with lag 0
      lagind_i <- lagind[-tarind]
      
      allscore <- NULL
      k <- 1
      de <- length(lagind_i)-(E.i-1)
      # Performing multiview embedding
      if(de>=0){
        # The condition for number of causal variables (including lags) larger than optimal embedding dimension 
        # In a large network, this condition is much more common 
        while(k<=kmax){ 
          # Generating one of multiview embedding from the random subsets of causal variables
          e.tt <- c(tarind,sample(lagind_i,E.i-1))
          block.t <- dst[,e.tt] 
          # Multivariate simplex projection for the one-step foreward forecast of node i at time= t+1 (tp=1)
          allscore <- rbind(allscore,
                            data.frame(block_lnlp(block.t,target_column = 1,
                                                  tp=1,method = "simplex"),matrix(e.tt,nrow=1)))
          k <- k+1
        }
      }else{
        # The under-embedding condition: Number of causal variables (including lags) is less than the optimal embedding dimension
        # This condition is very rare in large networks
        print(paste('variable',i,'has less dimensions',de))
        for(k in 1:kn){
          # Solve under-embedding by introducing more lag terms of target node
          block.t <- cbind(dst[,lagind_i],lgmm(ds[,i],E=max_lag+abs(de))[,-c(1:max_lag)])
          allscore <- rbind(allscore,
                            data.frame(block_lnlp(block.t,target_column = 1,
                                                  tp=1,method = "simplex"),matrix(c(lagind_i,-1:de),nrow=1)))
        }
      }
      
      # Rank multiview SSR according to their performance and then select the kn best multiview SSR  
      allscore <- allscore[order(allscore[,'rho'],decreasing = T),][1:kn,] 
      
      e.t2 <- as.matrix(allscore[,-c(1:15)])
      mx.t <- mx.na       # The list of kn sets of embedding variables
      mx.t2 <- mx.na2     # The list of time lags in mx.t2
      sunn <- as.matrix(e.t2%/%max_lag) # Find causal variable
      fisn <- as.matrix(e.t2%%max_lag)  # Assign lags of the causal variables
      fisn[fisn==0] <- max_lag
      sunn[fisn!=max_lag] <- sunn[fisn!=max_lag]+1
      # processing under-embedding cases
      if(any(e.t2<0)){sunn[which(e.t2<0)] <- i}
      if(any(e.t2<0)){fisn[which(e.t2<0)] <- max_lag+abs(e.t2[e.t2<0])}
      mx.t[,1:ncol(e.t2)] <- sunn
      mx.t2[,1:ncol(e.t2)] <- fisn
      
      ese.t <- cbind(sp=rep(i),rho=allscore[,'rho'],mx.t,mx.t2) # Compiled multiview results for node i
      esele <- rbind(esele, ese.t) # Compiled multiview results for all node
      
      t.end <- proc.time()
      t.use <- (t.end-t.begin)[3]
      cat('sample',i,'complete,',t.use,"sec","\n")
    }#end of else
  }# end of i 
  return(esele)
}

##############################################################################
## The computation of multiview distance based on the results of multiview SSR 
mvdist.demo2 <- function(ds,ds.all,esele,inferM=NULL){
  if(is.null(inferM)){inferM=matrix(1,ncol(ds),ncol(ds))}
  dmatrix.train.mvx <- dmatrix.all.mvx <- dmatrix.test.mvx <- list()
  for(i in 1:ncol(ds)){
    if(all(inferM[i,]==0)){
      dmatrix.train.mvx[[i]] <- dmatrix.all.mvx[[i]] <- dmatrix.test.mvx[[i]] <- 0
    }else{
      ecomb.t <- esele[esele[,'sp']==i,]                  #The multiview results of the i-th node
      e.tt <- as.matrix(ecomb.t[,-c(1:2)])
      e.tt <- as.matrix(e.tt[,!apply(is.na(e.tt),2,all)]) # Embedded variables & lags
      E.i <- ncol(e.tt)                                   # Embedding dimension
      
      ett.vb <- e.tt[,-agrep('vlag_',colnames(e.tt))]     # Extract embedded variables 
      ett.lag <- e.tt[,agrep('vlag_',colnames(e.tt))]     # Extract time lags 
      max_lag <- max(ett.lag,na.rm=T)
      
      # Create lag dataset 
      dst.all <- dst <- NULL
      for(k in 1:ncol(ds)){
        dst <- cbind(dst,lgmm(ds[,k],E=max_lag,nm=k))
        dst.all <- cbind(dst.all,lgmm(ds.all[,k],E=max_lag,nm=k))
      }
      e.tt <- ((ett.vb-1)*max_lag+1)+(ett.lag-1)
      ecomb.t[ecomb.t[,'rho']<0,'rho']=0
      ww <- ecomb.t[,'rho']/sum(ecomb.t[,'rho'])           # Weighs the multiview results based on their performance
      Dx.t <- matrix(0,nrow(ds),nrow(ds))
      Dx.t2 <- matrix(0,nrow(ds.all),nrow(ds.all))
      # Computation of multiview distance as the weighted average of multivariate distance
      if(is.null(dim(e.tt))){e.tt=matrix(e.tt,ncol=1)}
      for(j in 1:NROW(e.tt)){
        dX <- as.matrix(dst[,as.numeric(e.tt[j,])])        # Multivariate distance calculated from each multiview embedding
        Dx.t <- Dx.t+as.matrix(dist(dX))*ww[j]             # The comutation of weighted average  
        dX2 <- as.matrix(dst.all[,as.numeric(e.tt[j,])])   # Distance among all data points
        Dx.t2 <- Dx.t2+as.matrix(dist(dX2))*ww[j]          
      }
      dmatrix.train.mvx[[i]] <- Dx.t                       # Multiview distances among in-sample data points
      dmatrix.all.mvx[[i]] <- Dx.t2
      dmatrix.test.mvx[[i]] <- Dx.t2[-c(1:nrow(ds)),][,c(1:nrow(ds))] # Multiview distances of out-of-sample from in-sample data points
    }
  }# end of i sp
  return(list(dmatrix.all.mvx=dmatrix.all.mvx,
              dmatrix.train.mvx=dmatrix.train.mvx,
              dmatrix.test.mvx=dmatrix.test.mvx))
}

##############################################################################################
## Cross-validation for MDR S-map for selecting parameters: (theta, lambda, alpha)
cv.MDR.infM3 <- function(ds,ds_tp1,dmatrix.list, inferM=NULL, parall=T,ncore=3,keep_intra=T,
                         lambda.seq = sort(logspace(-3,0,15),decreasing = T),rm_zero=F,absence_m=NULL,na_m=NULL,
                         tht = c(0, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8),
                         alpha.seq=seq(0.1,0.9,0.1)
){
  if(is.null(inferM)){inferM=matrix(1,ncol(ds),ncol(ds))}
  # Parallel computation
  if(parall){
    cl  <-  makeCluster(ncore)
    registerDoParallel(cl)
  }
  
  # All combinations of parameters (theta, alpha)
  para.grid <- expand.grid(tht,alpha.seq) 
  colnames(para.grid) <- c('theta','alpha')
  nlamd <- length(lambda.seq)
  cind <- 1:ncol(ds)
  cv.ind <- NULL
  for(i in 1:ncol(ds)){
    t.begin <- proc.time() 
    if(all(inferM[i,]==0)){}else{
      cind.i <- which(inferM[i,]==1)
      dmatrix.na <- dmatrix.list[[i]] # Multiview distance matrix for node i
      pf <- rep(1,ncol(ds))
      if(keep_intra){pf[i] <- 0} # zero-penalty for diagonal coefficients
      pf=pf[inferM[i,]==1]
      ## In-sample prediction
      for(pind in 1:nrow(para.grid)){
        para.t <- as.numeric(para.grid[pind,])
        thet.i <- para.t[1]
        alpha.i <- para.t[2]
        lambda.i <- para.t[3]
        predicted <- NULL
        nonzero <- 0
        
        # Leave-one-out CV throughout the time series
        for(j in 1:nrow(ds)){
          d.j <-  dmatrix.na[,j]           # Multiview distnaces between time j and time k, k!=j
          dpar <- mean(d.j[-j],na.rm=T)    # Mean multiview distance
          naind <- as.numeric(1*(is.na(ds[j,i])|(all(is.na(d.j[-j])))|is.na(ds_tp1[j,i])))# Does the data points at time j contain NA?
          pf2=pf
          if(rm_zero){
            pf2[as.numeric(absence_m[j,inferM[i,]==1])==1]=Inf
            pf2[as.numeric(na_m[j,inferM[i,]==1])==1]=Inf
            }
          if(naind==1){
            predicted <- rbind(predicted,rep(NA,nlamd)) # The result becomes NA, if data contains NA
          }else if(rm_zero & absence_m[j,i]==1){
            predicted <- rbind(predicted,rep(NA,nlamd)) # The result becomes NA (not included in CV), if the data point includes continuous zero
          }else{
            d.j[j] <- NA                   # Excluded time j from the training sets
            w.j <- exp(-thet.i*d.j/dpar)   # Weights as the expoenentially decay function of multiview distance
            y <- sqrt(w.j)*ds_tp1[,i]      # Weighted data vector y = node i at time j+1 
            x <- sqrt(w.j)*ds[,cind.i]       # Weighted data matrix x = all the network node at time j 
            
            da.j <- data.frame(y=y,x)      # Compiled x and y
            da.j <- da.j[!(apply(is.na(x),1,all)|is.na(y)),] # excludes NA from the data frame
            da.j <- apply(da.j,2,as.numeric)
            if(any(is.na(da.j))){ da.j[is.na(da.j)]=0}
            # Solving sparse regression with elastic-net regularization
            fit0 <- glmnet(da.j[,-1], da.j[,'y'], alpha=alpha.i,lambda=lambda.seq,family="gaussian",penalty.factor=pf2)
            zeroi <- as.numeric(apply(fit0$beta!=0,2,sum))
            nonzero <- nonzero+zeroi
            
            # Make the prediction
            prex <- as.numeric(ds[j,cind.i]) # Predicted time j+1 from the data of time j 
            pred.t <- NULL
            if(parall){# Parallel condition
              # Prediction by different lambda
              pred.t <- foreach(lambda.t=lambda.seq, .combine=c, .packages = 'glmnet') %dopar% {
                as.numeric(predict(fit0,s=lambda.t,newx=matrix(prex,nrow=1)))}
            }else{     # Non-parallel condition
              for(lambda.t in lambda.seq){
                pred.t <- c(pred.t, as.numeric(predict(fit0,s=lambda.t,newx=matrix(prex,nrow=1))))
              }
            }
            if(is.null(pred.t)){print(j)}
            predicted <- rbind(predicted,pred.t)
          }# end of else
        }# end of j
        
        # Evaluate the skills (mse & rho) of one-step foreward forecast under a given parameter combination 
        cv.ind <- rbind(cv.ind,data.frame(variable=i,
                                          theta=thet.i,
                                          alpha=alpha.i,
                                          lambda=lambda.seq,
                                          mse=apply((predicted-ds_tp1[,i])^2,2,mean, na.rm=T),
                                          msesd=apply((predicted-ds_tp1[,i])^2,2,sd,na.rm=T),
                                          nonzero=nonzero/sum(!is.na(predicted[,1])), # mean nonzero jacobian
                                          rho=apply(predicted,2,cor,y=ds_tp1[,i],use='pairwise.complete.obs')))
      }# end of pind
      time.used <- proc.time() - t.begin
      cat("effect", i, "in sample completed:", time.used[3],"sec\n")
    }# end of else
  } # end of i species
  if(parall){stopCluster(cl)}# end of parallel computation
  return(cv.ind)
}
#############################################3
## select the parameter set with the minimal MSE 
secv.demo <- function(cv.ind,nonze=T){
  cv.see <- NULL
  if(nonze){nonzer <- 0}else{nonzer <- -1}
  for(i in 1:ncol(ds)){
    # discard the parameter sets leading to trivial solution including only intercept
    cv.ind.i <- filter(cv.ind,variable==i,nonzero>nonzer)
    # Find the minimal MSE
    mind <- which.min(cv.ind.i[,'mse'])[1]
    cv.seet <- cv.ind.i[mind,]
    min.mse <- cv.seet[,'mse']
    cv.see <- rbind(cv.see,cv.seet)
  } # end of node i
  rownames(cv.see) <- NULL
  return(cv.see)
}


##############################################################################################
## MDR S-map based on the parameters selected by cross-validation 
MDRsmap.infM3 <- function(paracv,ds,ds_tp1,ds.test=NULL,dst_tp1=NULL, inferM=NULL,
                          dmatrix.list=dmatrix.list.mvx, out.sample=T, keep_intra=T,
                          dmatrix.test.list=dmatrix.test.list.mvx,rm_zero=F,absence_m=NULL,absence_mt=NULL,na_m=NULL,na_mt=NULL,
                          ptype='enet'
){
  # parameter sets with the minimal MSE
  thet <- paracv[,'theta']
  alph <- paracv[,'alpha']
  lamb <- paracv[,'lambda']
  if(is.null(inferM)){inferM=matrix(1,ncol(ds),ncol(ds))}
  # dimensions of dataset
  nrds <- nrow(ds)
  ncds <- ncol(ds)
  if(is.null(colnames(ds))){colnames(ds)=paste('Node',1:ncol(ds),sep='_')}
  # Estimated jacobian in in-sample & out-of-sample
  jcof <- jcof.test <- NULL
  enet.in <- enet.out <- NULL
  jcof.b <- matrix(NA,nrds,ncds+1)
  cind <- 1:ncol(ds)
  if(out.sample){jcof.bo <- matrix(NA,nrow(ds.test),ncds+1) }
  naresult <- rep(NA,4);names(naresult) <- c('obs','pred','mae','mse')
  for(i in 1:ncol(ds)){ # node i
    t.begin <- proc.time()
    if(all(inferM[i,]==0)){}else{
      # Multiview diance
      dmatrix.na <- dmatrix.list[[i]]
      if(out.sample){dmatrix.test.na <- dmatrix.test.list[[i]]}
      cind.i <- which(inferM[i,]==1)
      i.ind=which(cind.i==i)
      pf <- rep(1,ncol(ds)); intra.i <- 0
      if(keep_intra){pf[i] <- 0; intra.i <- i.ind} # zero-penalty for diagonal elements
      pf=pf[inferM[i,]==1]
      # In-sample prediction
      #### The optimal parameter sets for node i
      thet.i <- thet[i]
      alpha.i <- alph[i]
      lambda.i <- lamb[i]
      jcof.i <- jcof.b 
      enet.t <- NULL
      for(j in 1:nrow(ds)){ # time j
        d.j <- dmatrix.na[,j]
        dpar <- mean(d.j[-j],na.rm=T)
        naind <- as.numeric(1*(is.na(ds[j,i])|(all(is.na(d.j[-j])))|is.na(ds_tp1[j,i]))) # Any NA?
        pf2=pf
        if(rm_zero){
          pf2[as.numeric(absence_m[j,inferM[i,]==1])==1]=Inf
          pf2[as.numeric(na_m[j,inferM[i,]==1])==1]=Inf
          }
        if(naind==1|(rm_zero & (sum(absence_m[j,]==0)==1))){
          enet.t <- rbind(enet.t,naresult)
        }else if(rm_zero & absence_m[j,i]==1){
          enet.t <- rbind(enet.t, naresult)
          jcof.i[j, c(1,cind+1)] <- 0
        }else{
          d.j[j] <- NA # excluded the target time j
          w.j <- exp(-thet.i*d.j/dpar) # weights based on the multiview diance 
          y <- sqrt(w.j)*ds_tp1[,i]
          x <- sqrt(w.j)*ds[,cind.i]
          # Prediction X is also based on in-sample dataset
          prex <- as.numeric(ds[j,cind.i])
          da.j <- data.frame(y=y,x)
          da.j <- da.j[!(apply(is.na(x),1,all)|is.na(y)),] # excludes NA from the data frame
          da.j <- apply(da.j,2,as.numeric)
          if(any(is.na(da.j))){ da.j[is.na(da.j)]=0}
          
          if(ptype=='enet'){ # regularization using elastic-net
            enet.fit <- glmnet(da.j[,-1], da.j[,'y'], alpha=alpha.i, penalty.factor = pf2,
                               lambda=lambda.i,family="gaussian")
            predicted <-predict(enet.fit,s=lambda.i, newx=matrix(prex,nrow=1))
          }else if(ptype=='aenet'){ # regularization using adaptive elastic-net
            enet.fit <- msaenet.ck(x=da.j[,-1], y=da.j[,'y'], alpha.ini=alpha.i,lambda.ini=lambda.i, penalty.factor.init=pf2,
                                   keep_intra_i=intra.i, newx=matrix(prex,nrow=1),obs=as.numeric(ds_tp1[j,i]),nstep=1L)
            predicted <-as.numeric(enet.fit$predicted)
          }
          
          # calculating the prediction skills
          mse <- as.numeric((ds_tp1[j,i] - predicted)^2)
          mae <- as.numeric(abs(ds_tp1[j,i] - predicted))
          # High-dimensional Jacobian matrix was estimated from the the fitting of elastic-net
          jcof.i[j,1] <- as.numeric(enet.fit$a0)
          jcof.i[j,cind.i+1] <- as.numeric(enet.fit$beta)
          jcof.i[j,cind.i+1][na_m[j,]==1]=NA
          enet.t <- rbind(enet.t,c(obs=as.numeric(ds_tp1[j,i]),pred=predicted,mae=mae,mse=mse))
        }# end else
      }# end of j
      colnames(jcof.i) <- c('j0',colnames(ds))
      # output compiling parameter sets and prediction skills
      enet.in <- rbind(enet.in,c(theta=thet.i,
                                 alpha=alpha.i,
                                 lambda=lambda.i,
                                 rho=cor(enet.t[,'obs'],enet.t[,'pred'],use="pairwise.complete.obs"),
                                 mae=mean(enet.t[,'mae'],na.rm=T),
                                 rmse=sqrt(mean(enet.t[,'mse'],na.rm=T))
      ))
      # save Jacobians
      jcof <- rbind(jcof,data.frame(Insample=rep(1),time=1:nrow(ds),variable=rep(i),jcof.i))
      time.used <- proc.time() - t.begin
      cat("effect", i, "in sample completed:", time.used[3],"sec\n")
      
      ## Out-of-sample prediction follows the same approach presented for in-sample prediction,
      ## but the out-of-sample data points were not involved in the fitting
      if(out.sample){
        jcof.i <- jcof.bo 
        enet.t <- NULL
        for(j in 1:nrow(ds.test)){
          d.j <-  dmatrix.test.na[j,] 
          dpar <- mean(d.j[-j],na.rm=T)
          naind <- as.numeric(1*(is.na(ds[j,i])|(all(is.na(d.j[-j])))|is.na(ds_tp1[j,i])))
          pf2=pf
          if(rm_zero){
            pf2[as.numeric(absence_mt[j,inferM[i,]==1])==1]=Inf
            pf2[as.numeric(na_mt[j,inferM[i,]==1])==1]=Inf
            }
          if(naind==1){
            enet.t <- rbind(enet.t,naresult)
          }else if(rm_zero & absence_mt[j,i]==1){
            enet.t <- rbind(enet.t, naresult)
            jcof.i[j, c(1,cind+1)] <- 0
          }else{
            w.j <- exp(-thet.i*d.j/dpar)
            # Both x & y included only in-samples but no out-of-samples
            x <- sqrt(w.j)*ds[,cind.i]
            y <- sqrt(w.j)*ds_tp1[,i]
            # Prediction X is based on out-of-sample
            prex <- as.numeric(ds.test[j,cind.i])
            
            da.j <- data.frame(y=y,x)
            da.j <- da.j[!(apply(is.na(x),1,all)|is.na(y)),] # excludes NA from the data frame
            da.j <- apply(da.j,2,as.numeric)
            if(any(is.na(da.j))){ da.j[is.na(da.j)]=0}
            
            if(ptype=='enet'){
              enet.fit <- glmnet(da.j[,-1], da.j[,'y'], alpha=alpha.i, penalty.factor = pf2,
                                 lambda=lambda.i,family="gaussian")
              predicted <-predict(enet.fit,s=lambda.i, newx=matrix(prex,nrow=1))
            }else if(ptype=='aenet'){
              enet.fit <- msaenet.ck(x=da.j[,-1], y=da.j[,'y'], alpha.ini=alpha.i,lambda.ini=lambda.i,penalty.factor.init=pf2,
                                     keep_intra_i=intra.i, newx=matrix(prex,nrow=1),obs=as.numeric(dst_tp1[j,i]),nstep=1L)
              predicted <-as.numeric(enet.fit$predicted)
            }
            
            mse <- as.numeric((dst_tp1[j,i] - predicted)^2)
            mae <- as.numeric(abs(dst_tp1[j,i] - predicted))
            jcof.i[j,1] <- as.numeric(enet.fit$a0)
            jcof.i[j,cind.i+1] <- as.numeric(enet.fit$beta)
            jcof.i[j,cind.i+1][na_mt[j,]==1]=NA
            enet.t <- rbind(enet.t,c(obs=as.numeric(dst_tp1[j,i]),pred=predicted,mae=mae,mse=mse))
            
          }# end else
        }# end of j  
        colnames(jcof.i) <- c('j0',colnames(ds))
        
        enet.out <- rbind(enet.out,c(theta=thet.i,
                                     alpha=alpha.i,
                                     lambda=lambda.i,
                                     rho=cor(enet.t[,'obs'],enet.t[,'pred'],use="pairwise.complete.obs"),
                                     mae=mean(enet.t[,'mae'],na.rm=T),
                                     rmse=sqrt(mean(enet.t[,'mse'],na.rm=T))
        ))
        jcof.test <- rbind(jcof.test,data.frame(Insample=rep(0),time=1:nrow(ds.test),variable=rep(i),jcof.i))
        time.used2 <- proc.time() - t.begin - time.used
        cat("effect", i, "out of sample completed:", time.used2[3],"sec\n")
      } # end of if (out-of-sample)
    }# end of else
    
  } # end of i species
  
  #  Arranging the Jacobian coefficients by time and variables--> time-varying Jacobian matrix
  colnames(jcof) <- c('Insample','time','variable','j0',colnames(ds))
  jcof <- arrange(data.frame(jcof),time,variable)
  colnames(jcof) <- c('Insample','time','variable','j0',colnames(ds))
  
  if(out.sample){ 
    colnames(jcof.test) <- c('Insample','time','variable','j0',colnames(ds))
    jcof.test <- arrange(data.frame(jcof.test),time,variable)
    colnames(jcof.test) <- c('Insample','time','variable','j0',colnames(ds))
    jcof <- rbind(jcof,jcof.test)
    nr.out <- rbind(cbind(In_sample=rep(1),enet.in),cbind(In_sample=rep(0),enet.out))
  }else{
    nr.out <- cbind(In_sample=rep(1),enet.in)
  }
  return(list(nr.out=nr.out,jcof=jcof))
}

#   The algorithm for adaptive elastic-net is modified from the R code provided in the following paper 
#   Nan Xiao and Qing-Song Xu. (2015). Multi-step adaptive elastic-net: reducing false positives in
#   high-dimensional variable selection. Journal of Statistical Computation and Simulation 85(18),3755-3765.
msaenet.ck <- function (x, y, family = "gaussian",#c("gaussian", "binomial","poisson", "cox"), #init = c("enet", "ridge"), 
                        lambda.ini = NULL, alpha.ini = NULL, intra=F, penalty.factor.init = rep(1, ncol(x)),keep_intra_i=0,
                        alphas = seq(0.05, 0.95, 0.05), tune = 'cv',#c("cv", "ebic","bic", "aic"),
                        nfolds = 5L, rule = "lambda.min", #c("lambda.min","lambda.1se"), 
                        ebic.gamma = 1, nsteps = 2L,tune.nsteps = "max",#c("max","ebic", "bic", "aic"), 
                        ebic.gamma.nsteps = 1,scale = 1, lower.limits = -Inf, upper.limits = Inf, 
                        seed = 1001, parallel = FALSE, newx=NULL, obs=NULL) 
{ 
  best.alphas <- rep(NA, nsteps)
  best.lambdas <- rep(NA, nsteps)
  a0.list <- rep(NA, nsteps)
  beta.list <- vector("list", nsteps)
  model.list <- vector("list", nsteps)
  adapen.list <- vector("list", nsteps)
  
  best.alpha.enet <- alpha.ini
  best.lambda.enet <- lambda.ini
  #pf <- rep(1, ncol(x))
  if(keep_intra_i!=0){penalty.factor.init[keep_intra_i] <- 0}
  enet.full <- glmnet(x = x, y = y, family = family, alpha = best.alpha.enet, 
                      lambda = best.lambda.enet, lower.limits = lower.limits, 
                      upper.limits = upper.limits, penalty.factor=penalty.factor.init)
  a0.enent <- enet.full$a0
  bhat <- as.matrix(enet.full$beta)
  if (all(bhat == 0)){bhat <- rep(.Machine$double.eps * 2, length(bhat))}
  adpen <- (pmax(abs(bhat), .Machine$double.eps))^(-scale)
  if(any(is.infinite(penalty.factor.init))){adpen[is.infinite(penalty.factor.init)]=Inf}
  if(keep_intra_i!=0){adpen[keep_intra_i] <- 0}
  model.cv <- msaenet.tune.glmnet.ck(x = x, y = y, family = family,
                                     alphas = alphas, tune = tune, nfolds = nfolds, rule = rule, 
                                     ebic.gamma = ebic.gamma, lower.limits = lower.limits, 
                                     upper.limits = upper.limits, seed = seed + 1L, parallel = F, 
                                     penalty.factor = adpen)
  
  best.alphas[[1L]] <- model.cv$best.alpha
  best.lambdas[[1L]] <- model.cv$best.lambda
  model.list[[1L]] <- glmnet(x = x, y = y, family = family, 
                             alpha = best.alphas[[1L]], lambda = best.lambdas[[1L]], 
                             lower.limits = lower.limits, upper.limits = upper.limits, 
                             penalty.factor = adpen)
  
  bhat <- as.matrix(model.list[[1L]][["beta"]])
  if (all(bhat == 0)) 
    bhat <- rep(.Machine$double.eps * 2, length(bhat))
  a0.list[1L] <- as.numeric(model.list[[1L]]$a0)
  beta.list[[1L]] <- bhat
  adapen.list[[1L]] <- adpen
  
  if(nsteps>1L){
    for (i in 1L:(nsteps-1L)) {
      adpen.raw <- (pmax(abs(beta.list[[i]]), .Machine$double.eps))^(-scale)
      if(keep_intra_i!=0){adpen.raw[keep_intra_i]=0}
      if(any(is.infinite(penalty.factor.init))){adpen.raw[is.infinite(penalty.factor.init)]=Inf}
      adapen.list[[i]] <- as.vector(adpen.raw)
      adpen.name <- rownames(adpen.raw)
      names(adapen.list[[i]]) <- adpen.name
      model.cv <- msaenet.tune.glmnet.ck(x = x, y = y, family = family, 
                                         alphas = alphas, tune = tune, nfolds = nfolds, rule = rule, 
                                         ebic.gamma = ebic.gamma, lower.limits = lower.limits, 
                                         upper.limits = upper.limits, seed = seed + i, parallel = parallel, 
                                         penalty.factor = adapen.list[[i]])
      best.alphas[[i + 1L]] <- model.cv$best.alpha
      best.lambdas[[i + 1L]] <- model.cv$best.lambda
      model.list[[i + 1L]] <- glmnet(x = x, y = y, family = family, 
                                     alpha = best.alphas[[i + 1L]], lambda = best.lambdas[[i +1L]], lower.limits = lower.limits, upper.limits = upper.limits, 
                                     penalty.factor = adapen.list[[i]])
      bhat <- as.matrix(model.list[[i + 1L]][["beta"]])
      if (all(bhat == 0)) 
        bhat <- rep(.Machine$double.eps * 2, length(bhat))
      a0.list[i + 1L] <- as.numeric(model.list[[i + 1L]]$a0)
      beta.list[[i + 1L]] <- bhat
      adapen.list[[i + 1L]] <- adpen
    }
  }
  
  if(nsteps>1L){
    best.step <- which.min(mse)[1]
    predicted <- unlist(lapply(model.list,predict,newx=newx))
    mse=(obs-predicted)^2
  }else{
    best.step=1L
    predicted=predict(model.list[[1L]],newx=newx)
    mse=(obs-predicted)^2
  }
  msaenet.model <- list(a0=a0.list[[best.step]], beta = Matrix(beta.list[[best.step]],sparse = TRUE), model = model.list[[best.step]], 
                        best.step = best.step,best.alphas = best.alphas, best.lambdas = best.lambdas,
                        mse=mse,predicted=predicted,nvb=unlist(lapply(beta.list,function(x){return(sum(x!=0))})),
                        beta.list = beta.list, model.list = model.list, adapen.list = adapen.list, 
                        seed = seed, call = call)
  return(msaenet.model)
}


#   The algorithm for adaptive elastic-net is modified from the R code provided in the following paper
#   Nan Xiao and Qing-Song Xu. (2015). Multi-step adaptive elastic-net: reducing false positives in
#   high-dimensional variable selection. Journal of Statistical Computation and Simulation 85(18),3755-3765.
msaenet.tune.glmnet.ck <- function(
    x, y, family,
    alphas,
    lambdas=NULL,
    tune,
    nfolds, rule,
    ebic.gamma,
    lower.limits, upper.limits,
    seed, parallel, ...) {
  
  if (tune == "cv") {
    if (!parallel) {
      model.list <- vector("list", length(alphas))
      for (i in 1L:length(alphas)) {
        set.seed(seed)
        model.list[[i]] <- cv.glmnet(
          x = x, y = y, family = family,
          nfolds = nfolds, alpha = alphas[i],lambda=lambdas,
          lower.limits = lower.limits,
          upper.limits = upper.limits, ...
        )
      }
    } else {
      model.list <- foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        cv.glmnet(
          x = x, y = y, family = family,
          nfolds = nfolds, alpha = alphas,lambda=lambdas,
          lower.limits = lower.limits,
          upper.limits = upper.limits, ...
        )
      }
    }
    
    errors <- unlist(lapply(model.list, function(x) min(sqrt(x$"cvm"))))
    errors.min.idx <- which.min(errors)
    best.model <- model.list[[errors.min.idx]]
    
    best.alpha <- alphas[errors.min.idx]
    
    if (rule == "lambda.min") best.lambda <- best.model$"lambda.min"
    if (rule == "lambda.1se") best.lambda <- best.model$"lambda.1se"
    
    step.criterion <- errors[errors.min.idx]
  } else {
    if (!parallel) {
      model.list <- vector("list", length(alphas))
      for (i in 1L:length(alphas)) {
        set.seed(seed)
        model.list[[i]] <- glmnet(
          x = x, y = y, family = family,
          alpha = alphas[i],lambda=lambdas,
          lower.limits = lower.limits,
          upper.limits = upper.limits, ...
        )
      }
    } else {
      model.list <- foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        glmnet(
          x = x, y = y, family = family,
          alpha = alphas,lambda=lambdas,
          lower.limits = lower.limits,
          upper.limits = upper.limits, ...
        )
      }
    }
    
    if (tune == "aic") {
      ics.list <- mapply(
        .aic,
        deviance = lapply(model.list, .deviance),
        df = lapply(model.list, .df),
        SIMPLIFY = FALSE
      )
    }
    
    if (tune == "bic") {
      ics.list <- mapply(
        .bic,
        deviance = lapply(model.list, .deviance),
        df = lapply(model.list, .df),
        nobs = lapply(model.list, .nobs),
        SIMPLIFY = FALSE
      )
    }
    
    if (tune == "ebic") {
      ics.list <- mapply(
        .ebic,
        deviance = lapply(model.list, .deviance),
        df = lapply(model.list, .df),
        nobs = lapply(model.list, .nobs),
        nvar = lapply(model.list, .nvar),
        gamma = ebic.gamma,
        SIMPLIFY = FALSE
      )
    }
    
    ics <- sapply(ics.list, function(x) min(x))
    ics.min.idx <- which.min(ics)
    best.model <- model.list[[ics.min.idx]]
    
    best.alpha <- alphas[ics.min.idx]
    
    best.ic.min.idx <- which.min(ics.list[[ics.min.idx]])
    best.lambda <- best.model$"lambda"[[best.ic.min.idx]]
    
    step.criterion <- ics.list[[ics.min.idx]][[best.ic.min.idx]]
  }
  
  list(
    "best.model" = best.model,
    "best.alpha" = best.alpha,
    "best.lambda" = best.lambda,
    "step.criterion" = step.criterion
  )
}
