xzaddpara <- function(x1,x2,x1ty,x2ty,model,lam1,lam2,D1,D2,x1bar,x2bar) {
  p1=ncol(x1);p2=ncol(x2)
  model=sort(model)
  model1=model[model<=p1]
  model2=model[model>p1]
  p01=length(model1);p02=length(model2)
  D01 <- D1[model1];D02 <- D2[model2-p1]
  x01bar <- x1bar[model1];x02bar <- x2bar[model2-p1];
  x0t1 <- D01*{t(x1[,model1,drop=FALSE]) - x01bar}
  x0t2 <- D02*{t(x2[,model2-p1,drop=FALSE]) - x02bar}
  x0t <-rbind(x0t1,x0t2)
  R1 <- chol(tcrossprod(x0t) + Diagonal(p01+p02,x=c(rep(lam1,p01),rep(lam2,p02))))
  v1 <- backsolve(R1,c(x1ty[model1],x2ty[model2-p1]), transpose = T)
  return(list(R1=R1, v1=v1))
}

#xzaddpara(x1,x2,x1ty,x2ty,model,lam1,lam2,D1,D2,x1bar,x2bar)
#library(bravo)
#bravo:::addpara(cbind(x1,x2),c(x1ty,x2ty),model,lam1,c(D1,D2),c(x1bar,x2bar))
