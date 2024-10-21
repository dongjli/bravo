xzsven.temp <- function(x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,topKeep,D1,D2,x1bar,x2bar,temp,stepsize,logp.best,r.idx.best,n,ncovar,Miter) {
  currlogp <- numeric(Miter)
  currRSS <- numeric(Miter)
  curridx <- integer()
  model.sizes <- integer(Miter)
  start <- 2
  end <- 1
  
  logp.curr <- as.numeric(-(n-1)/2*log(n-1))
  
  r <- xzaddvar(model=NULL,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,R0=NULL,v0=NULL,D1,D2,x1bar,x2bar)
  res1 <- r$logp
  res1.RSS <- r$RSS
  res1.top <- res1[max(res1)-res1 < 6]
  k <- max(topKeep, length(res1.top))
  rc.k.idx <- order(res1, decreasing = T)[1:k]
  rc.k <- res1[rc.k.idx]
  RSS.k <- res1.RSS[rc.k.idx]
  probs = exp(rc.k-max(rc.k, na.rm = T)) / (temp*stepsize)
  probs[is.na(probs)] = 0
  s <- sample.int(k, 1, prob = probs)
  # s <- sample.int(k, 1, prob = exp((rc.k-max(rc.k, na.rm = T))/(temp*stepsize)))
  
  logp.curr <- rc.k[s]
  RSS.curr <- RSS.k[s]
  rc.idx <- rc.k.idx[s]
  if (logp.curr > logp.best) {
    logp.best <- logp.curr
    r.idx.best <- rc.idx
  }
  
  para.add <- xzaddpara(x1,x2,x1ty,x2ty,model=rc.idx,lam1,lam2,D1,D2,x1bar,x2bar)
  
  currlogp[1]<-logp.curr
  currRSS[1] <- RSS.curr
  model.sizes[1] <- 1
  curridx[1] <- rc.idx
  
  for (m in 1:(Miter-1)) {
    r.idx.old <- sort(rc.idx)
    if (length(rc.idx) > 0) {
      para.add <- xzaddpara(x1,x2,x1ty,x2ty,model=rc.idx,lam1,lam2,D1,D2,x1bar,x2bar)
      r.add <- xzaddvar(model=rc.idx,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,R0=para.add$R1,v0=para.add$v1,D1,D2,x1bar,x2bar)
      res.add <- r.add$logp
      RSS.add <- r.add$RSS
      r.swap <-  xzswapvar(model=rc.idx,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,D1,D2,x1bar,x2bar,swapOnly=F)
      res.swap <- r.swap$logp
      RSS.swap <- r.swap$RSS
    } else {
      r.add <- xzaddvar(model=NULL,x1,x2,ys,x1ty,x2ty,lam1,lam2,w1,w2,R0=NULL,v0=NULL,D1,D2,x1bar,x2bar)
      res.add <- r.add$logp
      RSS.add <- r.add$RSS
      res.swap <- c()
      RSS.swap <- c()
    }
    rc <- cbind(res.add, res.swap)
    rc.RSS <- cbind(RSS.add, RSS.swap)
    rc.top <- rc[max(rc)-rc < 6]
    k <- max(topKeep, length(rc.top))
    rc.k.idx <- order(rc, decreasing = T)[1:k]
    rc.k <- rc[rc.k.idx]
    RSS.k <- rc.RSS[rc.k.idx]
    probs = exp(rc.k-max(rc.k, na.rm = T)) / (temp*stepsize)
    probs[is.na(probs)] = 0
    s <- sample.int(k, 1, prob = probs)
    # s <- sample.int(k, 1, prob = exp((rc.k-max(rc.k, na.rm = T))/(temp*stepsize)))
    logp.curr <- rc.k[s]
    RSS.curr <- RSS.k[s]
    idx.pick = rc.k.idx[s]
    if(idx.pick <= ncovar) {
      rc.idx <- sort(c(r.idx.old, idx.pick))
    } else {
      idx.temp <- idx.pick - ncovar
      idx.add <- idx.temp %% ncovar
      # if(idx.add == 0) idx.add = ncovar
      idx.del <- idx.temp %/% ncovar
      if(idx.add > 0) {
        if (idx.del == 0) {
          rc.idx <- r.idx.old[r.idx.old!=idx.add]
        } else {
          rc.idx <- c(r.idx.old[-idx.del], idx.add)
        }
      } else {
        if (idx.del == 1) {
          rc.idx <- r.idx.old[r.idx.old!=ncovar]
        } else {
          rc.idx <- c(r.idx.old[-(idx.del-1)], ncovar)
        }
      }
    }
    
    if (logp.curr > logp.best) {
      logp.best <- logp.curr
      r.idx.best <- rc.idx
    }
    
    currlogp[m+1] <- logp.curr
    currRSS[m+1] <- RSS.curr
    currlength <- length(rc.idx)
    model.sizes[m+1] <- currlength
    end <- end + currlength
    if (currlength != 0){
      curridx[start:end] <- rc.idx
      start <- start + currlength
    }
  }
  #print(unique(currlogp))
  return(list(bestlogp=logp.best, bestidx=r.idx.best, currlogp=currlogp,
              modelsizes=model.sizes, curridx=curridx, currRSS=currRSS))
}
#m1=sort(test1$bestidx);m2=test1$bestlogp
#set.seed(481)
#temptest1=xzsven.temp(x1,x2,y,x1ty,x2ty,lam1,lam2,w1,w2,topKeep=20,D1,D2,x1bar,x2bar,temp=100,stepsize=50,logp.best=m2,r.idx.best=m1,n=nrow(x1),ncovar=ncol(x1)+ncol(x2),Miter=100) 
#m3=sort(test2$bestidx);m4=test2$bestlogp
#set.seed(481)
#temptest2=bravo:::sven.temp(cbind(x1,x2),y,c(x1ty,x2ty),lam1,w1,topKeep = 20,c(D1,D2),c(x1bar,x2bar),temp=100,stepsize=50,logp.best=m4,r.idx.best=m3,n=nrow(x1),ncovar=ncol(x1)+ncol(x2),Miter=100)
#sort(temptest1$bestidx)
#sort(temptest2$bestidx)
#temptest1$bestlogp
#temptest2$bestlogp
