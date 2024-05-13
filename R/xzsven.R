
sven.xz <- function(x1,x2, y, w1 = sqrt(nrow(x1))/ncol(x1),w2 = sqrt(nrow(x2))/ncol(x2),
                 lam1 = nrow(x1)/ncol(x1)^2,lam2 = nrow(x2)/ncol(x2)^2, Ntemp = 3,
                 Tmax = (log(log(ncol(x1)+ncol(x2)))+log(ncol(x1)+ncol(x2))), Miter = 50, wam.threshold = 0.5, 
                 log.eps = -16, L = 20, verbose = TRUE) {
  result <- list()
  
  n <- length(y)
  p1=ncol(x1);p2=ncol(x2)
  ncovar <- p1+p2
  ys = scale(y)
  x1bar =colMeans(x1);x2bar = colMeans(x2)
  xbar = c(x1bar, x2bar)
  
  stopifnot(class(x1)[1] %in% c("dgCMatrix","matrix"))
  if(class(x1)[1] == "dgCMatrix") {
    D1 = 1/sqrt(colMSD_dgc(x1,x1bar))
  }  else   {
    D1 = apply(x1,2,sd)
    D1 = 1/D1
  }
  
  stopifnot(class(x2)[1] %in% c("dgCMatrix","matrix"))
  if(class(x2)[1] == "dgCMatrix") {
    D2 = 1/sqrt(colMSD_dgc(x2,xbar))
  }  else   {
    D2 = apply(x2,2,sd)
    D2 = 1/D2
  }
  
  
  x1ty = D1*as.numeric(crossprod(x1,ys))
  x2ty = D2*as.numeric(crossprod(x2,ys))
  xty=c(x1ty,x2ty)
  
  logp <- numeric(Miter * (Ntemp))
  RSS <- numeric(Miter * (Ntemp))
  size <- integer(Miter * (Ntemp))
  indices <- integer(Miter*100)
  
  if (verbose) cat("temperature = 1\n")
  o <- xzsven.notemp(x1,x2,y,x1ty,x2ty,lam1,lam2,w1,w2,L,D1,D2,x1bar,x2bar,n,ncovar,Miter)
  logp.best <- o$bestlogp
  r.idx.best <- o$bestidx
  logp[1:Miter] <- o$currlogp
  RSS[1:Miter] <- o$currRSS
  size[1:Miter] <- o$modelsizes
  ed <- sum(size)
  indices[1:ed] <- o$curridx
  
  stepsize <- Tmax/(Ntemp-1)
  for (t in 1:(Ntemp-1)) {
    if (verbose) cat("temperature =", t*stepsize, "\n")
    o <- xzsven.temp(x1,x2,y,x1ty,x2ty,lam1,lam2,w1,w2,topKeep=L,D1,D2,x1bar,x2bar,t,stepsize= stepsize,logp.best,r.idx.best,n,ncovar,Miter) 
    logp.best <- o$bestlogp
    r.idx.best <- o$bestidx
    logp[(Miter*t+1):(Miter*(t+1))] <- o$currlogp
    RSS[(Miter*t+1):(Miter*(t+1))] <- o$currRSS
    size[(Miter*t+1):(Miter*(t+1))] <- o$modelsizes
    indices[(ed+1):(ed+sum(o$modelsizes))] <- o$curridx
    ed <- ed + sum(o$modelsizes)
  }
  indices <- indices[indices>0]
  cumsize <- cumsum(size)
  modelSparse <- sparseMatrix(i=indices,p = c(0,cumsize),index1 = T,dims = c(ncovar,length(logp)), x = T)
  
  logp.uniq1 <- unique(logp)
  for (i in 1:(length(logp.uniq1)-1)) {
    for (j in 2:length(logp.uniq1)) {
      if (abs(logp.uniq1[i]-logp.uniq1[j]) < 1e-10) {
        logp.uniq1[j] <- logp.uniq1[i]
      }
    }
  }
  logp.uniq <- unique(logp.uniq1)
  logp.top <- sort(logp.uniq[(logp.best-logp.uniq)<(-log.eps)], decreasing = T)
  cols.top <- unlist(lapply(logp.top, FUN=function(x){which(x==logp)[1]}))
  size.top <- size[cols.top]
  RSS.top <- RSS[cols.top] 
  model.top <- modelSparse[, cols.top, drop=F]
  #model1.top = model.top[1:p1,];model2.top = model.top[(p1+1):(p1+p2),]
  
  
  logp.top1 <- logp.top-logp.best
  weight <- exp(logp.top1)/sum(exp(logp.top1))
  
  # beta1.est <- matrix(0, (p1+1), length(cols.top))
  # for(i in 1:length(cols.top)){
  #  if (size.top[i]==0){
  #    beta1 <- mean(y)
  #    beta1.est[1, i] <- beta1
  #  } else { 
  #    m_i = model.top[, i][1:p1]
  #    emni = sum(model.top[, i][1:p1])
  #    x1.est <- cbind(rep(1, n), scale(x1[, m_i], center = x1bar[m_i], scale = 1/D1[m_i]))
  #    beta1 <- solve(crossprod(x1.est) + lam1*diag(c(0, rep(1, emni))), crossprod(x1.est, y))
  #    beta1.est[c(T, m_i), i] <- c(beta1[1]-sum(beta1[-1]*x1bar[m_i]*D1[m_i]), beta1[-1] * D1[m_i])
  #  }
  # }
  # 
  # beta1.MAP <- beta1.est[, 1]
  # beta1.WAM <- rowSums(beta1.est%*%diag(weight, nrow=length(size.top)))
  
  
  # beta2.est <- matrix(0, (p2+1), length(cols.top))
  # for(i in 1:length(cols.top)){
  #   if (size.top[i]==0){
  #     beta2 <- mean(y)
  #     beta2.est[1, i] <- beta2
  #   } else {
  #     m_i = model.top[, i][(p1+1):(p1+p2)]
  #     emni=sum(model.top[, i][(p1+1):(p1+p2)])
  #     x2.est <- cbind(rep(1, n), scale(x2[, m_i], center = x2bar[m_i], scale = 1/D2[m_i]))
  #     beta2 <- solve(crossprod(x2.est) + lam2*diag(c(0, rep(1,emni))), crossprod(x2.est, y))
  #     beta2.est[c(T, m_i), i] <- c(beta2[1]-sum(beta2[-1]*x2bar[m_i]*D2[m_i]), beta2[-1] * D2[m_i])
  #   }
  # }
  # 
  # beta2.MAP <- beta2.est[, 1]
  # beta2.WAM <- rowSums(beta2.est%*%diag(weight, nrow=length(size.top)))
  
  
  
  MIP = rowSums(model.top%*%Diagonal(length(weight), weight))
  model.WAM <- sort(which(MIP >= wam.threshold))
  model.MAP <- sort(r.idx.best)
  MIP.MAP <- MIP[model.MAP]
  MIP.WAM <- MIP[model.WAM]
  
  mtop.idx <- apply(model.top, 2, which)
  model.union <- Reduce(union, mtop.idx)
  model1.union = model.union[model.union<=p1]
  model2.union = model.union[model.union>p1]
  Xm1 <- sparseMatrix(i = integer(0), j = integer(0), dims = c(n, p1))
  Xm1 <- as(Xm1, "dMatrix")
  Xm1[, model1.union] <- x1[, model1.union, drop=F]
  Xm2 <- sparseMatrix(i = integer(0), j = integer(0), dims = c(n, p2))
  Xm2 <- as(Xm2, "dMatrix")
  Xm2[, model2.union-p1] <- x2[, model2.union-p1, drop=F]
  
  
  result$model.map <- model.MAP
  result$model.wam <- model.WAM
  result$model.top <- model.top
  # result$beta2.map <- beta2.MAP
  # result$beta1.wam <- beta1.WAM
  # result$beta1.map <- beta1.MAP
  # result$beta1.wam <- beta1.WAM
  result$mip.map <- MIP.MAP
  result$mip.wam <- MIP.WAM
  result$pprob.map <- logp.best
  result$pprob.top <- logp.top
  result$stats <- list(RSS.top = RSS.top, weights = weight, Xm = Xm1, Zm=Xm2, y = y, lam1 = lam1, lam2 = lam2)
  
  class(result) <- c(class(result), "sven")
  return(result)
} 

