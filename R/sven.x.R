
sven.x <- function(X, y, w = sqrt(nrow(X))/ncol(X), lam = nrow(X)/ncol(X)^2, Ntemp = 3,
                 Tmax = (log(log(ncol(X)))+log(ncol(X))), Miter = 50, wam.threshold = 0.5, 
                 log.eps = -16, L = 20, verbose = TRUE) {
  result <- list()

  n <- length(y)
  ncovar <- ncol(X)
  ys = scale(y)
  xbar = colMeans(X)

  stopifnot(class(X)[1] %in% c("dgCMatrix","matrix"))
  if(class(X)[1] == "dgCMatrix") {
    D = 1/sqrt(colMSD_dgc(X,xbar))
  }  else   {
    D = apply(X,2,sd)
    D = 1/D
  }
  Xty = D*as.numeric(crossprod(X,ys))

  logp <- numeric(Miter * (Ntemp))
  RSS <- numeric(Miter * (Ntemp))
  size <- integer(Miter * (Ntemp))
  indices <- integer(Miter*100)

  if (verbose) cat("temperature = 1\n")
  o <- sven.notemp(X, ys, Xty, lam, w, L, D, xbar, n, ncovar, Miter)
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
    o <- sven.temp(X, ys, Xty, lam, w, L, D, xbar, t, stepsize=stepsize, logp.best, r.idx.best, n, ncovar, Miter)
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

  logp.top1 <- logp.top-logp.best
  weight <- exp(logp.top1)/sum(exp(logp.top1))

  beta.est <- matrix(0, (ncovar+1), length(cols.top))
  for(i in 1:length(cols.top)){
   if (size.top[i]==0){
     beta <- mean(y)
     beta.est[1, i] <- beta
   } else {
     m_i = model.top[, i]
     x.est <- cbind(rep(1, n), scale(X[, m_i], center = xbar[m_i], scale = 1/D[m_i]))
     beta <- solve(crossprod(x.est) + lam*diag(c(0, rep(1, size.top[i]))), crossprod(x.est, y))
     beta.est[c(T, m_i), i] <- c(beta[1]-sum(beta[-1]*xbar[m_i]*D[m_i]), beta[-1] * D[m_i])
   }
  }

  beta.MAP <- beta.est[, 1]
  beta.WAM <- rowSums(beta.est%*%diag(weight, nrow=length(size.top)))

  MIP = rowSums(model.top%*%Diagonal(length(weight), weight))
  model.WAM <- sort(which(MIP >= wam.threshold))
  model.MAP <- sort(r.idx.best)
  MIP.MAP <- MIP[model.MAP]
  MIP.WAM <- MIP[model.WAM]

  mtop.idx <- apply(model.top, 2, which)
  model.union <- Reduce(union, mtop.idx)
  Xm <- sparseMatrix(i = integer(0), j = integer(0), dims = c(n, ncovar))
  Xm <- as(Xm, "dMatrix")
  Xm[, model.union] <- X[, model.union, drop=F]

  result$model.map <- model.MAP
  result$model.wam <- model.WAM
  result$model.top <- model.top
  result$beta.map <- beta.MAP
  result$beta.wam <- beta.WAM
  result$mip.map <- MIP.MAP
  result$mip.wam <- MIP.WAM
  result$pprob.map <- logp.best
  result$pprob.top <- logp.top
  result$stats <- list(RSS.top = RSS.top, weights = weight, Xm = Xm, y = y, lam = lam)

  class(result) <- c(class(result), "sven")
  return(result)
}


