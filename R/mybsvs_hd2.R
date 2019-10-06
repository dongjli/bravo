#' Bayesian Variable Selection (ultra-high, high or low dimensional).
#' @rdname mybsvs
#' @description Perform Bayesian Variable Selection in Gaussian regression models
#' @param X  An \eqn{n x p} matrix. Sparse matrices are supported and every
#' care is taken not to make copies of this (typically) giant matrix.
#' No need to center or scale.
#' @param y  The response vector of length \code{n}.
#' @param w The prior inclusion probability of each variable.
#' @param lam The slab precision parameter.
#' as suggested by the theory of Wang et al. (2019).
#' @param temp.multip The temperature multiple. Default: 3.
#' @param M The number of iteration. Default: 200.

mybsvs <- function(X, y, w, lam, temp.multip=3, M=200, trueidx) {
  result <- list()
  k=20

  n <- length(y)
  ncovar <- ncol(X)
  ys = scale(y)
  xbar = colMeans(X)

  stopifnot(class(X) %in% c("dgCMatrix","matrix"))
  if(class(X) == "dgCMatrix") {
    D = 1/sqrt(colMSD_dgc(X,xbar))
  }  else   {
    D = apply(X,2,sd)
    D = 1/D
  }
  Xty = D*as.numeric(crossprod(X,ys))

  logp <- numeric()
  size <- integer()
  indices <- integer(M*100)

  o <- bsvs1(X, ys, Xty, lam, w, k, D, xbar, n, ncovar)
  saveRDS(o, file = "./results2/o0.rds")
  logp.best <- o$bestlogp
  r.idx.best <- o$bestidx
  nmodel <- o$nmodel
  logp[1:nmodel] <- o$currlogp
  size[1:nmodel] <- o$modelsizes
  ed <- sum(size)
  indices[1:ed] <- o$curridx


  t0 = 1
  for (t in t0:9) {
    cat("t =", t, "\n")
    o <- bsvs2_temp(X, ys, Xty, lam, w, k, D, xbar, t, temp.multip, logp.best, r.idx.best, n, ncovar)
    saveRDS(o, paste0("./results2/o", t, ".rds"))
    logp.best <- o$bestlogp
    r.idx.best <- o$bestidx
    nmodel_curr <- o$nmodel
    logp[(nmodel+1):(nmodel+nmodel_curr)] <- o$currlogp
    size[(nmodel+1):(nmodel+nmodel_curr)] <- o$modelsizes
    indices[(ed+1):(ed+sum(o$modelsizes))] <- o$curridx
    ed <- ed + sum(o$modelsizes)
    nmodel <- nmodel + nmodel_curr
  }
  indices <- indices[indices>0]
  cumsize <- cumsum(size)
  print(length(size))
  print(length(logp))
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
  logp.top <- sort(logp.uniq[(logp.best-logp.uniq)<16], decreasing = T)
  cols.top <- unlist(lapply(logp.top, FUN=function(x){which(x==logp)[1]}))
  size.top <- size[cols.top]
  model.top <- modelSparse[, cols.top]

  logp.top <- logp.top-logp.best
  weight <- exp(logp.top)/sum(exp(logp.top))

  beta.est <- matrix(0, (ncovar+1), length(cols.top))
  for(i in 1:length(cols.top)){
   if (size.top[i]==0){
     beta <- mean(y)
     beta.est[1, i] <- beta
   } else {
     m_i = model.top[, i]
     x.est <- cbind(rep(1, n), scale(X[, m_i], center = F, scale = 1/D[m_i]))
     beta <- solve(crossprod(x.est) + lam*diag(c(0, rep(1, size.top[i]))), crossprod(x.est, y))
     beta.est[c(T, m_i), i] <- c(beta[1], beta[-1] * D[m_i])
   }
  }

  beta.mode <- beta.est[, 1]
  beta.weighted <- rowSums(beta.est%*%diag(weight, nrow=length(size.top)))

  model.BMA <- sort(which(rowSums(model.top%*%Diagonal(length(weight), weight)) >= 0.5))
  model.MAP <- sort(r.idx.best)
  models <- list(model.BMA=model.BMA, model.MAP=model.MAP)

  nvarsBMA <- length(model.BMA)
  includeBMA <- ifelse (all(trueidx %in% model.BMA), 1, 0)
  if (length(model.BMA) != 0){
    fdrBMA <- sum(!(model.BMA %in% trueidx), na.rm = T)/length(model.BMA)
  } else {
    fdrBMA <- 0
  }
  fnrBMA <- (length(trueidx) - length(intersect(model.BMA, trueidx)))/(length(trueidx))
  jaccard.weighted <- length(intersect(model.BMA, trueidx)) / length(union(model.BMA, trueidx))

  nvarsMAP <- length(model.MAP)
  includeMAP <- ifelse (all(trueidx %in% model.MAP), 1, 0)
  fdrMAP <- sum(!(model.MAP %in% trueidx), na.rm = T)/length(model.MAP)
  fnrMAP <- (length(trueidx) - length(intersect(model.MAP, trueidx)))/(length(trueidx))
  jaccardMAP <- length(intersect(model.MAP, trueidx)) / length(union(model.MAP, trueidx))

  result$models <- models
  result$betaMAP <- beta.mode
  result$betaWeighted <- beta.weighted
  result$BMA$nvars <- nvarsBMA
  result$MAP$nvars <- nvarsMAP
  result$BMA$fdr <- fdrBMA
  result$MAP$fdr <- fdrMAP
  result$BMA$fnr <- fnrBMA
  result$MAP$fnr <- fnrMAP
  result$BMA$include <- includeBMA
  result$MAP$include <- includeMAP
  result$BMA$Jaccard <- jaccard.weighted
  result$MAP$Jaccard <- jaccardMAP
  return(result)
}

mybsvs <- cmpfun(mybsvs)
