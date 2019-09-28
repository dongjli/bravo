mh <- function(logp.curr, rc.k, k, t=1, multip=3, temperature) {
  s <- sample.int(k, 1, prob = exp(rc.k-max(rc.k, na.rm = T)))
  logp.prop <- rc.k[s]
  if (!temperature) {
    alpha <- logp.prop-logp.curr/t
  } else {
    alpha <- (logp.prop-logp.curr)/(t*multip)
  }
  if (log(runif(1)) <= alpha) {
    logp.next <- logp.prop
    success = T
  } else {
    logp.next <- logp.curr
    s = NULL
    success = F
  }
  logp.curr <- logp.next
  return(list(logp=logp.curr, idx=s, success=success))
}
mh <- cmpfun(mh)


