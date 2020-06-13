### Unnormalized log-posterior probability of a model.
logp.u <- function(model,X,y,lam,w)
{

  if(is.logical(model)) model = which(model);

  model = as.integer(model)

  n = nrow(X)
  logw = log(w/(1-w))

  p0 = length(model)
  if(p0 == 0)
    return( -0.5*(n-1)*log(n-1) );

  x.g = scale(X[,model,drop=FALSE])

  ys = scale(y)
  xtx = crossprod(x.g) + diag(x = lam,nrow = p0)

  R = chol(xtx)

  z = backsolve(R,crossprod(x.g,ys),transpose = T)

  logp = 0.5*p0*log(lam) - sum(log(diag(R))) - 0.5*(n-1)*log(sum(ys^2) - sum(z^2)) + p0*logw

  return(logp)

}
