### A function to compute the log (unormalized) posterior probabilities for all models
### in the "swapped" set for the current model.
swapvar <- function(model, x, ys, xty, lam, w, D, xbar, swapOnly=F) {
  n <- nrow(x)
  p <- ncol(x)
  p0 <- length(model)
  xtx <- n - 1
  yty <- n-1
  logw <- log(w/(1-w))
  if (swapOnly == T) {
    logp <- matrix(0, nrow=p, ncol=p0)
    RSS.mat <- matrix(0, nrow=p, ncol=p0)
    if (p0 == 1) {
      r1 <- addvar(model=NULL, x=x, ys=ys, xty=xty, lam=lam, w=w, D=D, xbar=xbar)
      logp[, 1] <- r1$logp
      logp[model, 1] <- -Inf
      RSS.mat[, 1] <- r1$RSS
      RSS.mat[model, 1] <- -Inf
    } else {
      x0 <- scale(x[, model, drop=F])
      xgx <- crossprod(x0) + lam*diag(p0)
      x0x1 <- crossprod(x0, x)
      x0x <- x0x1 %*% Diagonal(p,x=D)
      for (j in 1:p0) {
        # delete one variable in the current model
        model.temp <- model[-j]
        R0 <- chol(xgx[-j, -j])
        logdetR0 <- sum(log(diag(R0)))
        if(is.nan(logdetR0)) logdetR0 = Inf
        v0 <- backsolve(R0, xty[model.temp], transpose = T)

        # add back another variable from the remaning variables
        S <- backsolve(R0, x0x[-j, , drop=F], transpose = T)
        if(length(model.temp) == 1) S = matrix(S, nrow = 1)
        if(class(S)[1] == "dgeMatrix") {
          sts <- colSumSq_dge(S@x,S@Dim)
        } else {
          sts <- colSumSq_matrix(S)
        }
        sts[model] <- 0;
        s0 <- sqrt({xtx+lam} - sts)

        u <- (xty-crossprod(S, v0))/s0
        u[model] = 0;
        logdetR1 <- sum(log(diag(R0))) + log(s0)
        RSS <- {yty - sum(v0^2)} - u^2
        RSS[model] = 1 # whatever, they are going to be set to -Inf
        logp1 <- 0.5*(p0)*log(lam) - logdetR1 - 0.5*(n-1)*log(RSS) + p0*logw
        logp[, j] = as.numeric(logp1)
        logp[model, j] <- -Inf
        RSS.mat[, j] <- RSS
        RSS.mat[model, j] <- -Inf
      }
    }
  } else {
    logp <- matrix(0, nrow=p, ncol=p0+1)
    RSS.mat <- matrix(0, nrow=p, ncol=p0+1)
    if (p0 == 1) {
      logp.del <- -(n-1)/2*log(yty)
      RSS.del <- 0
      r1 <- addvar(model=NULL, x=x, ys=ys, xty=xty, lam=lam, w=w, D=D, xbar=xbar)
      logp[, 2] <- r1$logp
      logp[model, 2] <- -Inf
      RSS.mat[, 2] <- r1$RSS
      RSS.mat[model, 2] <- -Inf
    } else {
      logp.del <- numeric(p0)
      RSS.del <- numeric(p0)
      x0 <- scale(x[, model, drop=F])
      xgx <- crossprod(x0) + lam*diag(p0)
      x0x1 <- crossprod(x0, x)
      x0x <- x0x1 %*% Diagonal(p,x=D)
      for (j in 1:p0) {
        # delete one variable in the current model
        model.temp <- model[-j]
        R0 <- chol(xgx[-j, -j])
        logdetR0 <- sum(log(diag(R0)))
        if(is.nan(logdetR0)) logdetR0 = Inf
        v0 <- backsolve(R0, xty[model.temp], transpose = T)
        RSS0 <- yty - sum(backsolve(R0, xty[model.temp], transpose = T)^2)
        if(RSS0 <= 0) RSS0 = .Machine$double.eps
        logp.del[j] <- 0.5*(p0-1)*log(lam) - logdetR0 - 0.5*(n-1)*log(RSS0) + (p0-1)*logw
        RSS.del[j] <- RSS0

        # add back another variable from the remaning variables
        S <- backsolve(R0, x0x[-j, , drop=F], transpose = T)
        if(length(model.temp) == 1) S = matrix(S, nrow = 1)
        if(class(S)[1] == "dgeMatrix") {
          sts <- colSumSq_dge(S@x,S@Dim)
        } else {
          sts <- colSumSq_matrix(S)
        }
        sts[model] <- 0;
        s0 <- sqrt({xtx+lam} - sts)

        u <- (xty-crossprod(S, v0))/s0
        u[model] = 0;
        logdetR1 <- sum(log(diag(R0))) + log(s0)
        RSS <- {yty - sum(v0^2)} - u^2
        RSS[model] = 1 # whatever, they are going to be set to -Inf
        logp1 <- 0.5*(p0)*log(lam) - logdetR1 - 0.5*(n-1)*log(RSS) + p0*logw
        logp[, j+1] = as.numeric(logp1)
        logp[model, j+1] <- -Inf
        RSS.mat[, j+1] <- RSS
        RSS.mat[model, j+1] <- -Inf
      }
    }
    logp[model, 1] <- logp.del
    logp[-model, 1] <- -Inf
    RSS.mat[model, 1] <- RSS.del
    RSS.mat[-model, 1] <- -Inf
  }
  return(list(logp=logp, RSS=RSS.mat))
}

