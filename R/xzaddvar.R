
xzaddvar <- function(model,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,R0=NULL,v0=NULL,D1,D2,x1bar,x2bar) {
  n <- nrow(x1)
  p1 <- ncol(x1);p2 <-ncol(x2) 
  model = sort(model) 
  model1=model[model<=p1];model2=model[model>p1]
  p01=length(model1);p02=length(model2) 
  p =p1+p2;p0=p01+p02
  yty <- n-1               
  x1tx1 <- n-1       
  x2tx2 <- n-1
  xtx <- n-1
  xty=c(x1ty,x2ty)
  logw1 <- log(w1/(1-w1))
  logw2 <- log(w2/(1-w2))
  
  if (p0 == 0) { # There is no variable in the model
    R11 <- sqrt(x1tx1+lam1)
    RSS1 <- yty - (x1ty/R11)^2
    R12 <- sqrt(x2tx2+lam2)
    RSS2 <- yty - (x2ty/R12)^2
    RSS <- c(RSS1,RSS2) 
    logp1 <- 0.5*log(lam1)-0.5*log(x1tx1+lam1)-0.5*(n-1)*log(RSS1) + logw1
    logp2 <- 0.5*log(lam2)-0.5*log(x2tx2+lam2)-0.5*(n-1)*log(RSS2) + logw2
    logp <- c(logp1,logp2) 
    j = which.max(logp)
    if(j<=p1)
    {return(list(logp=logp,R=sparseMatrix(i=1,j=1,x=R11,triangular=T),
                 v=xty[j]/R11,which.max=j, RSS=RSS))}
    else{return(list(logp=logp,R=sparseMatrix(i=1,j=1,x=R12,triangular=T),
                     v=xty[j]/R12,which.max=j, RSS=RSS))}
  }
  
  # else there's at least one variable
  
  D01 <- D1[model1]
  D02 <- D2[model2-p1]
  D0 = c(D01,D02)
  x01bar <- x1bar[model1]
  x02bar <- x2bar[model2-p1] 
  x0bar = c(x01bar,x02bar)
  x0t1 <- D01*{t(x1[,model1,drop=FALSE]) - x01bar} 
  x0t2 <- D02*{t(x2[,model2-p1,drop=FALSE]) - x02bar} 
  x0t=rbind(x0t1,x0t2)
  x0=t(x0t)
  #x0tnew <- D0*{t(x[,model,drop=FALSE]) - x0bar} 
  if(is.null(R0)) {
    d1=rep(lam1,p01);d2=rep(lam2,p02)
    R0 = chol(tcrossprod(x0t) + Diagonal(p01+p02,x=c(d1,d2)));
    #R0 = chol(crossprod(x0) + diag(c(d1,d2)))
  }
  if(is.null(v0)) {
    v0 = backsolve(R0,c(x1ty[model1],x2ty[model2-p1]),transpose = T)
    #v0 = backsolve(R0,xty[model],transpose = T)
  } 
  
  
  S1 <- backsolve(R0,x0t,transpose = T);
  if(p0 == 1) S1 = matrix(S1,nrow = 1)
  
  S1 <-  S1 - rowMeans(S1); # For numerical stability.
  
  Sa  <- S1 %*% x1;
  
  Sa <- Sa %*% Diagonal(p1,x=D1)
  
  Sb  <- S1 %*% x2;
  
  Sb <- Sb %*% Diagonal(p2,x=D2)
  
  S = cbind(Sa,Sb)
  
  if(class(S)[1] == "dgeMatrix") {
    sts <- colSumSq_dge(S@x,S@Dim) 
  } else {
    sts <- colSumSq_matrix(S)
  }
  
  #if(class(Sb)[1] == "dgeMatrix") {
  # sts2 <- colSumSq_dge(Sb@x,Sb@Dim) 
  #} else {
  #sts2 <- colSumSq_matrix(Sb)
  #}
  
  sts[model] <- 0
  s01 <- sqrt({x1tx1+lam1} - sts[1:p1]) 
  s02 <- sqrt({x2tx2+lam2} - sts[(p1+1):p]) 
  
  s0 = c(s01,s02)
  
  u <- (xty-crossprod(S, v0))/s0
  
  #u <- as.matrix(rbind(u1,u2)) 
  #u <- (xty-crossprod(S, v0))/s0
  u[model] = 0;
  logdetR1 <- sum(log(diag(R0))) + log(s0)
  RSS <- {yty - sum(v0^2)} - u^2 
  RSS[model] = 1 # whatever, they are going to be set to -Inf
  logp1 <- 0.5*(p01+1)*log(lam1)+ 0.5*(p02)*log(lam2) - logdetR1[1:p1] - 0.5*(n-1)*log(RSS[1:p1]) + (p01+1)*logw1 +(p02)*logw2
  logp2 <- 0.5*(p01)*log(lam1)+ 0.5*(p02+1)*log(lam2) - logdetR1[(p1+1):p] - 0.5*(n-1)*log(RSS[(p1+1):p]) + (p01)*logw1 +(p02+1)*logw2
  logp = c(logp1,logp2)
  logp[model] <- -Inf
  j = which.max(logp)
  
  # now update R1tinv and v1 and return as a list
  sj = S[,j]
  s0j = s0[j]
  R1 = rbind(cbind(R0,sj),c(rep(0,p0),s0j))
  R1 = as.matrix(R1)   ##as(R1,"dtCMatrix") depricated. simply keeping dense matrix format
  v1 = c(v0,u[j])
  
  RSS[model] <- -Inf
  return(list(logp=logp, R = R1, v=v1, which.max=j, RSS=RSS))
}




#xzaddvar(model=NULL,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,R0=NULL,v0=NULL,D1,D2,x1bar,x2bar,BugShikari=FALSE)

#model=NULL
#model=sample(1:60,5,replace=F)
#library(bravo)
#test1=xzaddvar(model,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,R0=NULL,v0=NULL,D1,D2,x1bar,x2bar,BugShikari=FALSE)
#test2=bravo:::addvar(model,cbind(x1,x2),ys,c(x1ty,x2ty),lam1,w1,R0=NULL,v0=NULL,c(D1,D2),c(x1bar,x2bar),BugShikari=FALSE)
#actual functions for double checking 

#sum(test1$logp==test2$logp)
#test1$logp - test2$logp > 0.000000000001
#sum(test1$RSS == test2$RSS) 
#test1$which.max == test2$which.max
#test1$RSS - test2$RSS
