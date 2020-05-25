### A function to compute the log (unormalized) posterior probabilities for all models
### in the "deleted" set for the current model.
delvar <- function(model, x, xty, lam, w, D, xbar) {
  n <- nrow(x)
  p <- ncol(x)
  p0 <- length(model)
  yty <- n-1
  logw <- log(w/(1-w))
  logp <- numeric(p)
  if (p0 == 1) {
    logp.del <- -(n-1)/2*log(yty)
  } else {
    logp.del <- numeric(p0)
    RSS.del <- numeric(p0)
    x0 <- scale(x[, model, drop=F])
    xgx <- crossprod(x0) + lam*diag(p0)
    for (j in 1:p0) {
      # delete one variable in the current model
      model.temp <- model[-j]
      R0 <- chol(xgx[-j, -j])
      logdetR0 <- sum(log(diag(R0)))
      if(is.nan(logdetR0)) logdetR0 = Inf
      RSS0 <- yty - sum(backsolve(R0, xty[model.temp], transpose = T)^2)
      if(RSS0 <= 0) RSS0 = .Machine$double.eps
      logp.del[j] <- 0.5*(p0-1)*log(lam) - logdetR0 - 0.5*(n-1)*log(RSS0) + (p0-1)*logw
      RSS.del[j] <- RSS0
    }
  }
  logp[model] <- logp.del
  logp[-model] <- -Inf
  RSS.del[model] <- RSS.del
  RSS.del[-model] <- -Inf
  return(list(logp=logp, RSS=RSS.del))
}

