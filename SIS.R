
whichpartrev <- function(x, n=30) {
  which(x >= -sort(-x, partial=n)[n])
}

sis <- function(X,y){
  p = ncol(X)
  n = nrow(X)
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
  max.var = n;
  
  return(whichpartrev(Xty, n))
}

holp <- function(X,y){
  p = ncol(X)
  n = nrow(X)
  ys = scale(y)
  
  xbar = colMeans(X)
  
  stopifnot(class(X) %in% c("dgCMatrix","matrix"))
  
  beta <- crossprod(X,solve(tcrossprod(X),y))
  return(whichpartrev(beta, n))
}

Hat_mat <- function(X){
  X%*% tcrossprod(ginv(crossprod(X)),X)
}

residual_FR <- function(x, X0, y){
  X <- cbind(x,X0)
  SSR <- crossprod(y,Hat_mat(X))%*%y
  return(SSR)
}

addvar_FR <- function(model, X, ys){
  n <- nrow(X)
  p <- ncol(X)
  p0 = length(model)
  
  if(p0==0){
    res <- which.max(crossprod(X,ys))
    return(res)
  }
  
  X0 <- X[,model]
  SSR <- apply(X, 2, residual_FR, X0=X0, y=ys)
  SSR[model] <- -Inf
  j <- which.max(SSR)
  return(j)
}

FR <- function(X,y){
  p = ncol(X)
  n = nrow(X)
  ys = scale(y)
  
  xbar = colMeans(X)
  stopifnot(class(X) %in% c("dgCMatrix","matrix"))
  model <- integer(0L)
  
  for(i in 1:n){
    this <- addvar_FR(model = model, X=X, ys = ys)
    model = c(model,this)
  }
  return(model)
}