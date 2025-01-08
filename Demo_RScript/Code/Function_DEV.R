
smap_ls <- function(x, E, theta, coeff=F, fit=F){
  n <- length(x)
  X <- embed(x, E)
  Y <- c(X[-1,1],NA)
  da.t <- data.frame(Y,X)
  naind <- apply(is.na(da.t),1,any)
  XD <- as.matrix(dist(X, method = "euclidean"))
  if(length(theta)>1){coeff=F; fit=F}
  #diag(XD) <- NA
  XDD <- XD[!naind,!naind]
  da.rna <- da.t[!naind,]
  nef <- NROW(da.rna)
  rho <- mae <- rmse <- c()
  Yn <- matrix(da.rna[,1],ncol=1)
  Xn <- as.matrix(cbind(rep(1,nef),da.rna[,-1]))
  if(NROW(Xn)>=(NCOL(Xn)+1)){
    for(thet.i in theta){
      smapc <- NULL
      Yhat <- NULL
      for(i in 1:NROW(da.rna)){
        di <- XDD[,i]
        dpar <- mean(di[-i],na.rm=T)
        wi <- diag(exp(-thet.i*di[-i]/dpar))
        #dpar <- mean(di,na.rm=T)
        #wi <- diag(exp(-thet.i*di/dpar))
        Yni <- Yn[-i]
        Xni <- Xn[-i,]
        smapc.i <- matrix(ginv(t(Xni)%*%wi%*%Xni)%*%t(Xni)%*%wi%*%Yni,nrow=1)
        Yhat.i <- sum(smapc.i*Xn[i,])
        if(coeff) {smapc <- rbind(smapc, smapc.i)}
        Yhat <- c(Yhat, Yhat.i)
      }
      rho <- c(rho, cor(Yhat, Yn, use="pairwise.complete.obs" ))
      mae <- c(mae, mean(abs(Yhat-Yn),na.rm=T))
      rmse <- c(rmse, sqrt(mean((Yhat-Yn)^2,na.rm=T)))
    }# end of thet.i
    if(coeff) {
      smapc.f <- matrix(NA, length(E:(n-1)), E+1)
      smapc.f[which(!naind),] <- smapc
      colnames(smapc.f) <- c('Intercept',paste('X',1:E,sep=''))
      smapc.f <- data.frame(time=E:(n-1),smapc.f)
    }else{smapc.f <- NULL}
    if(fit) {
      Yhat.f <- rep(NA, length(E:n))
      Yhat.f[which(!naind)+1] <- Yhat
      Yhat.f <- data.frame(time=E:n,pred=Yhat.f, obs=x[E:n])
    }else{Yhat.f <- NULL}
    output <- list(stat=data.frame(E=E, theta=theta, npred= nef, rho=rho,mae=mae,rmse=rmse), 
                   smap_cof= smapc.f,Yhat=Yhat.f)
  }else{    
    output <- list(stat=data.frame(E=E, theta=theta, npred= 0, rho=NA, mae=NA, rmse=NA), 
                           smap_cof= NULL,Yhat=NULL)
  }# end of if
  return(output)
}

xi <- rnorm(21); 
(smp.t <- smap_ls(x=xi, theta=1, E=10, fit=T, coeff=T))[[1]]
#s_map(xi, theta=1, E=10)


dev_cvls=function(raw_time_series,                
                  E = seq(1,10,by=1),
                  tau = seq(1,1,by=1),
                  theta = 0, # Starting from linear model
                  window_size = seq(15, 30, by=1),
                  step_size = 1, stat=F
){
  param=expand.grid(window_size,E,tau,theta)
  colnames(param)=c('w','E','tau','theta')
  param <- filter(param, w >(2*E+1))
  rho_out <- c()
  for(i in 1:nrow(param)){
    rho_i <- c()
    m <- 0
    param_i <- lapply(param[i,],as.numeric)
    while(m <= length(raw_time_series) - param_i[['w']])
    {
      raw_ts_part <- raw_time_series[(m + 1):(m + param_i[['w']])]
      time_series <- (raw_ts_part - mean(raw_ts_part, na.rm=TRUE))/sd(raw_ts_part, na.rm=TRUE)
      edm.m <- smap_ls(time_series, E=param_i[['E']], theta=param_i[['theta']])
      
      rho_i <- c(rho_i, as.numeric(edm.m[[1]]['rho']))
      m <- m + step_size
    }
    rho_out <- c(rho_out, mean(rho_i, na.rm=T))
  }
  matrix_results <- data.frame(param, rho=rho_out)
  
  ### select the optimal results by rho~w curve
  matrix_rho <- c()
  for(i in 1:length(window_size))
  {
    res.t=matrix_results[matrix_results$w == window_size[i],]
    matrix_rho <- c(matrix_rho, res.t[which.max(res.t$rho)[1],'rho'])
  }
  
  d.t <- diff(c(matrix_rho));
  # Find window size with rather saturated improvement 
  ind.w1 <- which((d.t<median(d.t))&c(matrix_rho[-length(matrix_rho)]>0))
  if(length(ind.w1)>0){ind.w2=ind.w1[1]}else{ind.w2=which.max(c(matrix_rho))}
  # Select the optimal moving window
  opt.w <- window_size[ind.w2]
  results.w <- filter(matrix_results,w==opt.w)
  # Select the optimal parameter sets
  param.opt=results.w[which.max(results.w[,'rho']),]
  
  if(stat){return(list(Param_opt1=param.opt,Stat=results))}else{
    return(list(Param_opt1=param.opt))}
}

DEV_perform_ls=function(
    raw_time_series,
    param.opt,
    theta = c(0, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8),
    step_size = 1
){
  window_size <- param.opt[1]
  Eopt <- param.opt[2]
  
  ### start algorithm
  window_indices <- seq(window_size, NROW(raw_time_series), step_size)
  matrix_result <- matrix(NaN, nrow = length(window_indices), ncol = 8)
  index <- 0
  
  for(j in window_indices)
  {
    index <- index + 1
    rolling_window <- raw_time_series[(j-window_size+1):j]
    
    norm_rolling_window <- (rolling_window - mean(rolling_window, na.rm=TRUE))/sd(rolling_window, na.rm=TRUE)
    
    # calculate best theta
    smap <- smap_ls(norm_rolling_window, E=Eopt, theta=theta)[[1]]
    theta_best <- theta[order(-smap[,'rho'])[1]]
    param.opt['theta'] <- theta_best
    
    # calculate eigenvalues for best theta
    smap <- smap_ls(norm_rolling_window, E=Eopt, theta=theta_best, coeff=TRUE)
    if(NROW(smap[[2]])>0){
      smap_co <- as.matrix(smap[[2]][,-c(1:2)])
      rownames(smap_co)=NULL
      
      matrix_eigen <- matrix(NA, nrow = NROW(smap_co), ncol = 3)
      for(k in 1:NROW(smap_co))
      {
        if(!is.na(smap_co[k,1]))
        {
          M <- rbind(c(smap_co[k, 1:Eopt]), cbind(diag(Eopt - 1), rep(0, Eopt - 1)))
          M_eigen <- eigen(M)$values
          lambda1 <- M_eigen[order(abs(M_eigen))[Eopt]]
          
          matrix_eigen[k,1] <- abs(lambda1)
          matrix_eigen[k,2] <- Re(lambda1)
          matrix_eigen[k,3] <- Im(lambda1)
        }
      }
      mean_eigen <- apply(matrix_eigen,2,mean,na.rm=T)
    }else{
      mean_eigen <- rep(NA, 3)
    }# end of if/else

    # save results
    matrix_result[index,1] <- j
    matrix_result[index,2] <- mean_eigen[1]
    matrix_result[index,3] <- mean_eigen[2]
    matrix_result[index,4] <- mean_eigen[3]
    matrix_result[index,5:8] <- as.numeric(smap[[1]][c('rho','mae','rmse','theta')])
    smap_co <- NULL; rm(smap); gc()
  }
  colnames(matrix_result)=c('t_max','Abs','Re','Im','rho','mae','rmse','theta')
  return(matrix_result)
}

EDM_NA_fill=function(x, E_d=2:8,max_iter=1, data_only=T){
  spx=simplex(x, E=E_d)
  optE=spx[which.max(spx$rho),'E']
  smp=s_map(x, E=optE)
  opth=smp[which.max(smp$rho),'theta']
  smp=s_map(x, E=optE, theta=opth)
  na.ind=which(is.na(x))
  i=0
  x.new=x
  while(length(na.ind)>0 & i<max_iter){
    pred.set=cbind(na.ind-optE,na.ind)
    pred.set[pred.set<=0]=1
    smp.o=s_map(x, E=optE, theta=opth,pred=pred.set, stats_only = FALSE)
    smp.pred=smp.o$model_output[[1]]
    pred.o=smp.pred[smp.pred[,'time']%in%na.ind,'pred']
    x.new=x
    if(any(na.ind<min(smp.pred[,'time'],na.rm=T))){pred.o=c(rep(NA,sum(na.ind<min(smp.pred[,'time']))),pred.o)}
    x.new[na.ind]=pred.o
    x=x.new
    na.ind=which(is.na(x))
    i=i+1
  }
  if(data_only){return(x.new)}else{list(data=x.new,max_tp=i,statistic=smp[1:9])}
}


zeropad=function(x){
  zint <- which(x==0)
  z1 <- zL <- zint[1]
  zc <- 0
  zlist=list()
  if(length(zint)>1){
    for(i in 2:length(zint)){
      if((zint[i]-z1)<=1){
        zL=c(zL,zint[i])
      }else{
        if(length(zL)>1){zc=zc+1; zlist[[zc]]=zL}
        zL=zint[i]
      }
      z1=zint[i]
    }
    zlist[[zc+1]]=zL
  }else{zlist=NULL}
  zlist <- unlist(lapply(zlist,length))
  if(length(zlist)==0){zout <- -Inf}else{zout <- max(zlist)}
  return(zout)
}


################################################################################
## Bifurcation analysis for microbial OTUs
periodf=function(x){
  ust=unique(x)
  period.a=NULL
  for(ui in 1:length(ust)){
    y=which(x==ust[ui])
    if(length(y)!=1){
      emb2=embed(y,2)
      period.ind=emb2[,1]-emb2[,2]
      indi=which(period.ind!=1)
      period.i=matrix(c(emb2[1,2], max(emb2)),nrow=1,ncol=2)
      if(length(indi)>0){
        for(ki in 1:length(indi)){
          period.i[ki,2]=emb2[indi[ki],2]
          period.i=rbind(period.i,c(emb2[indi[ki],1], NA))
        }
        period.i[NROW(period.i),2]=max(emb2)
      }
      colnames(period.i)=c('initial','final')
      period.a=rbind(period.a, data.frame(state=rep(ust[ui]),period.i)) 
    }
  }
  return(period.a)
}

################################################################################
## Automatic detection of bifurcation event  
find.bif2 <- function(da.t, ind.tip=NULL, ploty=F){
  ne2 <- sum(!is.na(da.t[,'value']))
  g.t <- gam(value~s(time, k = trunc(ne2/2)),data=da.t[1:ind.tip,])
  xmin <- da.t[which(!is.na(da.t[,'value']))[1],'time']
  newd.i <- seq(xmin, ind.tip, 0.01)
  nad <- length(newd.i)
  fd <- fderiv(g.t, newdata = data.frame(time=newd.i))
  # Calculate 95% confidence intervals for derivatives
  cid <- as.matrix(confint(fd, type = "confidence")[,c('lower','est','upper')])
  ci <- confint(g.t,parm='s(time)',data=data.frame(time=newd.i))[,c('lower','est','upper')]+g.t$coefficients["(Intercept)"]
  yhat <- predict.gam(g.t,newdata=data.frame(time=newd.i))
  yhatd1 <- cid[,'est']
  # Find DEV values increased crossing the theoretical threshold =1
  crxp <- sort(as.numeric(which((sign(yhat-1)[-length(yhat)]-sign((yhat-1))[-1])<0)))
  if(length(crxp)!=0){
    sig.cx <- warn.ind <- warn.max <- NULL
    for(i in 1:length(crxp)){
      crxp.i=crxp[i]
      ind1 <- (cid[1:crxp.i,'est']>0) # any increasing trend before the crossing point?
      prd <- periodf((ci[max(1, c(crxp[i-1]),na.rm=T):crxp.i,'lower']>0)*1)# Does this increasing trend start from some values significantly < 1
      prd <- prd[prd[,'state']==1 & prd[,'final']==crxp.i ,]
      if(NROW(prd)>0){
        ind2 <- (ci[as.numeric(prd['initial']):crxp.i,'upper']<1) 
      }else{ind2 <- F}
      sig.cx <- c(sig.cx, any(ind1)&any(ind2))
      
      if(sig.cx[i]){
        above.i <- which(ci[(crxp.i:min(c(crxp[i+1],nad),na.rm=T)),'lower']>1)
        if(length(above.i)==0){above.ind <- NA}else{above.ind <- min(above.i,na.rm=T)}
        warn.ind.i <- newd.i[(crxp.i:min(c(crxp[i+1],nad),na.rm=T))][above.ind]
        warn.max.i <- newd.i[(crxp.i:min(c(crxp[i+1],nad),na.rm=T))][which.max(ci[(crxp.i:min(c(crxp[i+1],nad),na.rm=T)),'est'])]
        
        if(length(warn.ind.i)==0) {warn.ind.i <- NA}else{warn.ind.i <- min(warn.ind.i)}
        warn.ind <- c(warn.ind , warn.ind.i)
        warn.max <- c(warn.max , warn.max.i)
      }
    }
    crxp.f <- crxp[sig.cx]
    crip.i=data.frame(BFcross=newd.i[crxp.f], BFwarn=warn.ind, BFmax=warn.max)
    crip.i=crip.i[NROW(crip.i),]
  }else{
    crip.i=NULL
  }
  
  if(ploty){
    plot(value~time, da.t,type='o')
    lines(ci[,'lower']~newd.i,col='blue')
    lines(ci[,'upper']~newd.i,col='red')
    lines(yhat~newd.i,col='orange')
    abline(h=1,col='grey')
    abline(v=unlist(c(crip.i)))
  }
  return(crip.i)
}
