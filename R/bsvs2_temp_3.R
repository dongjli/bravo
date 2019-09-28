#' Bayesian Iterated Screening (ultra-high, high or low dimensional).
#' @rdname bsvs2_temp
#' @description Perform Bayesian iterated screening in Gaussian regression models

bsvs2_temp <- function(xmat, ys, xty, lam, w, k, D, xbar, temp, temp.multip, logp.best, r.idx.best, n, ncovar) {
  Miter = 200
  currlogp <- numeric(Miter)
  curridx <- integer()
  model.sizes <- integer(Miter)
  start <- 2
  end <- 1

  logp.curr <- as.numeric(-(n-1)/2*log(n-1))

  r <- addvar(model=NULL, x=xmat, ys=ys, xty=xty, lam=lam, w=w, D=D, xbar=xbar)
  res1 <- r$logp
  rc.k.idx <- order(res1, decreasing = T)[1:k]
  rc.k <- res1[rc.k.idx]
  mh.step <- mh(logp.curr, rc.k, k, t=temp, multip=temp.multip, temperature = T)
  logp.curr <- mh.step$logp
  while(!mh.step$success) {
    mh.step <- mh(logp.curr, rc.k, k, t=temp, multip=temp.multip, temperature = T)
    logp.curr <- mh.step$logp
  }
  rc.idx <- rc.k.idx[mh.step$idx]

  para.add <- addpara(x=xmat, xty=xty, model=rc.idx, lam=lam, D=D, xbar=xbar)

  currlogp[1]<-logp.curr
  model.sizes[1] <- 1
  curridx[1] <- rc.idx
  count <- 0

  for (m in 1:(Miter-1)) {
    #logp.curr.old <- logp.curr
    r.idx.old <- rc.idx
    if (length(rc.idx) == 0) {
      r <- addvar(model=NULL, x=xmat, ys=ys, xty=xty, lam=lam, w=w, D=D, xbar=xbar)
      res <- r$logp
      rc.k.idx <- order(res, decreasing = T)[1:k]
      rc.k <- res[rc.k.idx]
      mh.step <- mh(logp.curr, rc.k, k, t=temp, multip=temp.multip, temperature = T)
      logp.curr <- mh.step$logp
      while(!mh.step$success) {
        mh.step <- mh(logp.curr, rc.k, k, t=temp, multip=temp.multip, temperature = T)
        logp.curr <- mh.step$logp
      }
      rc.idx <- rc.k.idx[mh.step$idx]
      para.add <- addpara(x=xmat, xty=xty, model=rc.idx, lam=lam, D=D, xbar=xbar)

    } else {
      die <- sample.int(3, 1)
      if (count == 1000 && die == die.prev) next
      die.prev <- die
      if (die == 1) {
        r <- addvar(model=rc.idx, x=xmat, ys=ys, xty=xty, lam=lam, w=w, R0=para.add$R1, v0=para.add$v1, D=D, xbar=xbar)
        res <- r$logp
        rc.k.idx <- order(res, decreasing = T)[1:k]
        rc.k <- res[rc.k.idx]
        mh.step <- mh(logp.curr, rc.k, k, t=temp, multip=temp.multip, temperature = T)
        logp.curr <- mh.step$logp
        count <- 0
        while(!mh.step$success && count < 1000) {
          mh.step <- mh(logp.curr, rc.k, k, t=temp, multip=temp.multip, temperature = T)
          logp.curr <- mh.step$logp
          count <- count + 1
        }
        if (count == 1000) next
        idx.add = rc.k.idx[mh.step$idx]
        rc.idx <- c(r.idx.old, idx.add)
        para.add <- addpara(x=xmat, xty=xty, model=rc.idx, lam=lam, D=D, xbar=xbar)

      } else if (die == 2) {
        res <- delvar(model=rc.idx, x=xmat, xty=xty, lam=lam, w=w, D=D, xbar=xbar)
        l <- length(rc.idx)
        rc.l.idx <- order(res, decreasing = T)[1:l]
        rc.l <- res[rc.l.idx]
        mh.step <- mh(logp.curr, rc.l, l, t=temp, multip=temp.multip, temperature = T)
        logp.curr <- mh.step$logp
        count <- 0
        while(!mh.step$success && count < 1000) {
          mh.step <- mh(logp.curr, rc.l, l, t=temp, multip=temp.multip, temperature = T)
          logp.curr <- mh.step$logp
          count <- count + 1
        }
        if (count == 1000) next
        idx.del <- rc.l.idx[mh.step$idx]
        rc.idx <- r.idx.old[r.idx.old != idx.del]
        if (length(rc.idx)>0) {
          para.add <- addpara(x=xmat, xty=xty, model=rc.idx, lam=lam, D=D, xbar=xbar)
        }

      } else {
        res <- swapvar(model=rc.idx, x=xmat, ys=ys, xty=xty, lam=lam, w=w, D=D, xbar=xbar, swapOnly = T)
        rc.k.idx <- order(res, decreasing = T)[1:k]
        rc.k <- res[rc.k.idx]
        mh.step <- mh(logp.curr, rc.k, k, t=temp, multip=temp.multip, temperature = T)
        logp.curr <- mh.step$logp
        count <- 0
        while(!mh.step$success && count < 1000) {
          mh.step <- mh(logp.curr, rc.k, k, t=temp, multip=temp.multip, temperature = T)
          logp.curr <- mh.step$logp
          count <- count + 1
        }
        if (count == 1000) next
        idx.temp <- rc.k.idx[mh.step$idx]
        idx.add <- idx.temp %% ncovar
        idx.del <- idx.temp %/% ncovar
        if(idx.add > 0) {
          rc.idx <- c(r.idx.old[-(idx.del+1)], idx.add)
        } else {
          rc.idx <- c(r.idx.old[-idx.del], ncovar)
        }
        para.add <- addpara(x=xmat, xty=xty, model=rc.idx, lam=lam, D=D, xbar=xbar)
      }
    }

    if (logp.curr > logp.best) {
      logp.best <- logp.curr
      r.idx.best <- rc.idx
    }

    print(sort(rc.idx))
    print(sort(r.idx.best))
    print(c(logp.best, logp.curr))

    currlogp[m+1] <- logp.curr
    currlength <- length(rc.idx)
    model.sizes[m+1] <- currlength
    end <- end + currlength
    if (currlength != 0){
      curridx[start:end] <- rc.idx
      start <- start + currlength
    }
  }
  currlogp <- currlogp[currlogp!=0]
  model.sizes <- model.sizes[model.sizes!=0]
  return(list(bestlogp=logp.best, bestidx=r.idx.best, currlogp=currlogp,
              modelsizes=model.sizes, curridx=curridx, temp=temp))
}

bsvs2_temp <- cmpfun(bsvs2_temp)
