### A function to compute the log (unormalized) posterior probabilities for all models
### in the "added" set for the current model.
addvar <- function (model, x, ys, xty, lam, w, R0 = NULL, v0 = NULL, D,xbar,BugShikari=FALSE) {
  n <- nrow(x)
  p <- ncol(x)
  p0 = length(model)


  yty <- n - 1
  xtx <- n - 1
  logw <- log(w/(1-w))

  if (p0 == 0) { # There is no variable in the model
    R1 <- sqrt(xtx+lam)
    RSS <- yty - (xty/R1)^2
    logp <- 0.5*log(lam)-0.5*log(xtx+lam)-0.5*(n-1)*log(RSS) + logw
    j = which.max(logp)
    return(list(logp=logp,R=sparseMatrix(i=1,j=1,x=R1,triangular=T),
                v=xty[j]/R1,which.max=j, RSS=RSS))
  }

  # else there's at least one variable

  D0 <- D[model]
  x0bar <- xbar[model]
  x0t <- D0*{t(x[,model,drop=FALSE]) - x0bar}

  if(is.null(R0)) {
    R0 = chol(tcrossprod(x0t) + diag(x = lam,nrow = p0));
  }
  if(is.null(v0)) {
    v0 = backsolve(R0,xty[model],transpose = T)
  }


  S1 <- backsolve(R0,x0t,transpose = T);
  if(p0 == 1) S1 = matrix(S1,nrow = 1)

  S1 <-  S1 - rowMeans(S1); # For numerical stability.

  S  <- S1 %*% x;

  S <- S %*% Diagonal(p,x=D)

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
  logp <- 0.5*(p0+1)*log(lam) - logdetR1 - 0.5*(n-1)*log(RSS) + (p0+1)*logw
  logp = as.numeric(logp)
  logp[model] <- -Inf
  j = which.max(logp)

  # now update R1tinv and v1 and return as a list
  sj = S[,j]
  s0j = s0[j]
  R1 = rbind(cbind(R0,sj),c(rep(0,p0),s0j))
  R1 = as(R1,"dtCMatrix")
  v1 = c(v0,u[j])

  if(BugShikari) {
    x1 = scale(x[,c(model,j),drop=FALSE])
    R1exact = chol(crossprod(x1) + diag(x=lam,nrow=p0+1))
    v1exact = backsolve(R1exact,crossprod(x1,ys),transpose = T)
    cat("Error in R = ",norm(as.matrix(R1-R1exact)),"\n")
    cat("Error in v = ",sum(abs(v1-v1exact)),"\n")
  }

  RSS[model] <- -Inf
  return(list(logp=logp, R = R1, v=v1, which.max=j, RSS=RSS))
}

