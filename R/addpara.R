addpara <- function(x, xty, model, lam, D, xbar) {
  D0 <- D[model]
  x0bar <- xbar[model]
  x0t <- D0*{t(x[,model,drop=FALSE]) - x0bar}
  R1 <- chol(tcrossprod(x0t) + lam * diag(length(model)))
  v1 <- backsolve(R1, xty[model], transpose = T)
  return(list(R1=R1, v1=v1))
}

