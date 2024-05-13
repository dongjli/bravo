### A function to compute the log (un normalized) posterior probabilities for all models
### in the "deleted" set for the current model.

xzdelvar <- function(model,x1,x2,x1ty,x2ty,lam1,lam2,w1,w2,D1,D2,x1bar,x2bar) {
  n <- nrow(x1)
  p1 <- ncol(x1);p2 <- ncol(x2)
  model1 = model[model<=p1];model2= model[model>p1]
  p01 <- length(model1) ; p02 <- length(model2)
  p0=p01+p02 ; p=p1+p2
  yty <- n-1
  xty=c(x1ty,x2ty)
  logw1 <- log(w1/(1-w1));logw2 <- log(w2/(1-w2))
  RSS=numeric(p1+p2); logp <- numeric(p1+p2)   
  if (p0 == 1) {
    logp.del <- -(n-1)/2*log(yty)
    RSS.del=yty
  } else {
    logp.del <- numeric(p0)
    RSS.del <- numeric(p0)
    x01 <- scale(x1[, model1, drop=F])
    x02 <- scale(x2[, model2-p1, drop=F])
    x0 <-cbind(x01,x02)
    xgx <- crossprod(x0) + Diagonal(p01+p02,x=c(rep(lam1,p01),rep(lam2,p02)))
    for (j in 1:p0) {
      # delete one variable in the current model
      model.temp <- model[-j]
      R0 <- chol(xgx[-j, -j])
      logdetR0 <- sum(log(diag(R0)))
      if(is.nan(logdetR0)) logdetR0 = Inf
      RSS0 <- yty - sum(backsolve(R0, xty[model.temp], transpose = T)^2)
      if(RSS0 <= 0) RSS0 = .Machine$double.eps
      logp_init <- 0.5*p01*log(lam1) - logdetR0 - 0.5*(n-1)*log(RSS0) + p01*logw1 +p02*logw2 +0.5*p02*log(lam2)
      if(j <=p01){logp.del[j] = logp_init -0.5*log(lam1)-logw1}
      else{logp.del[j] = logp_init -0.5*log(lam2)-logw2}
      RSS.del[j] <- RSS0
    }
    
  }
  
  
  
  #logp.del <- c(logp1.del,logp2.del)
  logp[model] <- logp.del
  logp[-model] <- -Inf
  #RSS.del <- c(RSS.del1,RSS.del2)
  RSS[model] <- RSS.del
  RSS[-model] <- -Inf
  return(list(logp=logp, RSS=RSS))
}

#model=sample(1:60,5,replace=F)
#library(bravo)
#test1=xzdelvar(model,x1,x2,x1ty,x2ty,lam1,lam2,w1,w2,D1,D2,xbar,x2bar)
#test2=bravo:::delvar(model,cbind(x1,x2),c(x1ty,x2ty),lam1,w1,c(D1,D2),c(x1bar,x2bar))
#test1$logp == test2$logp 
#test1$logp - test2$logp >1e-8
#test1$RSS == test2$RSS

#bravo te problem hochhe ekta variable dile
