### A function to compute the log (un normalized) posterior probabilities for all models
### in the "swapped" set for the current model.


xzswapvar <- function(model,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,D1,D2,x1bar,x2bar,swapOnly=F) { 
  n <- nrow(x1)
  p1 <- ncol(x1)
  p2 <- ncol(x2)
  p=p1+p2
  model=sort(model)
  model1 <- model[model<=p1]
  model2 <- model[model>p1]
  p01 <- length(model1)
  p02 <-length(model2)
  p0=p01+p02
  #x=as.matrix(cbind(x1,x2))
  xbar=c(x1bar,x2bar)
  D=c(D1,D2)
  #x=scale(x)
  x1tx1 <- n-1       
  x2tx2 <- n-1
  xtx <- n-1
  yty <- n-1
  xty=c(x1ty,x2ty)
  logw1 <- log(w1/(1-w1));logw2 <- log(w2/(1-w2))
  if (swapOnly == T) {
    logp <- matrix(0, nrow=p, ncol=p0)
    RSS.mat <- matrix(0, nrow=p, ncol=p0)
    if (p0 == 1) {
      r1 <- xzaddvar(model=NULL,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,R0=NULL,v0=NULL,D1,D2,x1bar,x2bar)
      logp[, 1] <- r1$logp
      logp[model, 1] <- -Inf
      RSS.mat[, 1] <- r1$RSS
      RSS.mat[model, 1] <- -Inf
    } else {
      x0 <- cbind(scale(x1[, model1, drop=F]),scale(x2[, model2-p1, drop=F]))
      xgx <- crossprod(x0) + diag(c(rep(lam1,p01),rep(lam2,p02)))
      x0x11 <- crossprod(x0, x1); x0x12 <- crossprod(x0, x2)
      x0x <- cbind(x0x11 %*% Diagonal(p1,x=D1),x0x12 %*% Diagonal(p2,x=D2))
      for (j in 1:p0) {
        # delete one variable in the current model
        model.temp <- model[-j]
        R0 <- chol(xgx[-j, -j])
        logdetR0 <- sum(log(diag(R0))) 
        if(is.nan(logdetR0)) logdetR0 = Inf
        v0 <- backsolve(R0, xty[model.temp], transpose = T)
        
        # add back another variable from the remaining variables
        S <- backsolve(R0, x0x[-j, , drop=F], transpose = T)
        if(length(model.temp) == 1) S = matrix(S, nrow = 1)
        if(class(S)[1] == "dgeMatrix") {
          sts <- colSumSq_dge(S@x,S@Dim)
        } else {
          sts <- colSumSq_matrix(S)
        }
        sts[model] <- 0;
        sts1=sts[1:p1]
        sts2=sts[(p1+1):p]
        s01 <- sqrt({x1tx1+lam1} - sts1)
        sts2[model2-p1] <- 0;
        s02 <- sqrt({x2tx2+lam2} - sts2)
        s0=c(s01,s02)
        
        u <- (xty-crossprod(S, v0))/s0
        u[model] = 0;
        logdetR1 <- sum(log(diag(R0))) + log(s0)
        RSS <- {yty - sum(v0^2)} - u^2
        RSS[model] = 1 # whatever, they are going to be set to -Inf
        logp1 <- 0.5*(p01)*log(lam1) - logdetR1 - 0.5*(n-1)*log(RSS) + p01*logw1 +0.5*(p02)*log(lam2) +p02*logw2
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
      r1 <- xzaddvar(model=NULL,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,R0=NULL,v0=NULL,D1,D2,x1bar,x2bar)
      logp[, 2] <- r1$logp
      logp[model, 2] <- -Inf
      RSS.mat[, 2] <- r1$RSS
      RSS.mat[model, 2] <- -Inf
    } else {
      
      logp.del <- numeric(p0)
      RSS.del <- numeric(p0)
      x0 <- cbind(scale(x1[, model1, drop=F]),scale(x2[, model2-p1, drop=F]))
      xgx <- crossprod(x0) + diag(c(rep(lam1,p01),rep(lam2,p02)))
      x0x11 <- crossprod(x0, x1); x0x12 <- crossprod(x0, x2)
      x0x <- cbind(x0x11 %*% Diagonal(p1,x=D1),x0x12 %*% Diagonal(p2,x=D2))
      for (j in 1:p0) {
        # delete one variable in the current model
        
        
        model.temp <- model[-j]
        R0 <- chol(xgx[-j, -j])
        logdetR0 <- sum(log(diag(R0)))
        if(is.nan(logdetR0)) logdetR0 = Inf
        v0 <- backsolve(R0, xty[model.temp], transpose = T)
        RSS0 <- yty - sum(backsolve(R0, xty[model.temp], transpose = T)^2)
        if(RSS0 <= 0) RSS0 = .Machine$double.eps
        logp_init <- 0.5*p01*log(lam1) - logdetR0 - 0.5*(n-1)*log(RSS0) + p01*logw1 + p02*logw2 +0.5*p02*log(lam2)
        if(model[j] <=p01){logp.del[j] = logp_init -0.5*log(lam1)-logw1}
        else{logp.del[j] = logp_init -0.5*log(lam2)-logw2}
        RSS.del[j] <- RSS0
        
        # add back another variable from the remaining variables
        
        rfinal=xzaddvar(model=model.temp,x1,x2,ys,x1ty,x2ty,lam1,lam2,
                        w1,w2,R0=NULL,v0=NULL,D1,D2,x1bar,x2bar)
        
        
        logp1 = rfinal$logp
        RSS = rfinal$RSS

        logp[, j+1] = as.numeric(logp1)
        logp[model, j+1] <- -Inf
        #print(length(RSS))
        RSS.mat[, j+1] <- as.numeric(RSS)
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

#model=sample(1:60,1,replace=F)
#test1=xzswapvar(model,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,D1,D2,x1bar,x2bar,swapOnly=F)
#test2=bravo:::swapvar(model,cbind(x1,x2),y,c(x1ty,x2ty),lam1,w1,c(D1,D2),c(x1bar,x2bar),swapOnly=F)
#test1$logp - test2$logp >0.0000001
#test1$RSS == test2$RSS
