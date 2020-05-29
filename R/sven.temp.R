sven.temp <- function(xmat, ys, xty, lam, w, topKeep, D, xbar, temp, stepsize, logp.best, r.idx.best, n, ncovar, Miter) {
  currlogp <- numeric(Miter)
  currRSS <- numeric(Miter)
  curridx <- integer()
  model.sizes <- integer(Miter)
  start <- 2
  end <- 1

  logp.curr <- as.numeric(-(n-1)/2*log(n-1))

  r <- addvar(model=NULL, x=xmat, ys=ys, xty=xty, lam=lam, w=w, D=D, xbar=xbar)
  res1 <- r$logp
  res1.RSS <- r$RSS
  res1.top <- res1[max(res1)-res1 < 6]
  k <- max(topKeep, length(res1.top))
  rc.k.idx <- order(res1, decreasing = T)[1:k]
  rc.k <- res1[rc.k.idx]
  RSS.k <- res1.RSS[rc.k.idx]
  s <- sample.int(k, 1, prob = exp((rc.k-max(rc.k, na.rm = T))/(temp*stepsize)))

  logp.curr <- rc.k[s]
  RSS.curr <- RSS.k[s]
  rc.idx <- rc.k.idx[s]
  if (logp.curr > logp.best) {
    logp.best <- logp.curr
    r.idx.best <- rc.idx
  }

  para.add <- addpara(x=xmat, xty=xty, model=rc.idx, lam=lam, D=D, xbar=xbar)

  currlogp[1]<-logp.curr
  currRSS[1] <- RSS.curr
  model.sizes[1] <- 1
  curridx[1] <- rc.idx

  for (m in 1:(Miter-1)) {
    r.idx.old <- rc.idx
    if (length(rc.idx) > 0) {
      para.add <- addpara(x=xmat, xty=xty, model=rc.idx, lam=lam, D=D, xbar=xbar)
      r.add <- addvar(model=rc.idx, x=xmat, ys=ys, xty=xty, lam=lam, w=w, R0=para.add$R1, v0=para.add$v1, D=D, xbar=xbar)
      res.add <- r.add$logp
      RSS.add <- r.add$RSS
      r.swap <- swapvar(model=rc.idx, x=xmat, ys=ys, xty=xty, lam=lam, w=w, D=D, xbar=xbar)
      res.swap <- r.swap$logp
      RSS.swap <- r.swap$RSS
    } else {
      r.add <- addvar(model=NULL, x=xmat, ys=ys, xty=xty, lam=lam, w=w, D=D, xbar=xbar)
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
    s <- sample.int(k, 1, prob = exp((rc.k-max(rc.k, na.rm = T))/(temp*stepsize)))
    logp.curr <- rc.k[s]
    RSS.curr <- RSS.k[s]
    idx.pick = rc.k.idx[s]
    if(idx.pick <= ncovar) {
      rc.idx <- c(r.idx.old, idx.pick)
    } else {
      idx.temp <- idx.pick - ncovar
      idx.add <- idx.temp %% ncovar
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
  return(list(bestlogp=logp.best, bestidx=r.idx.best, currlogp=currlogp,
              modelsizes=model.sizes, curridx=curridx, currRSS=currRSS))
}

