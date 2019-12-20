#' Unnormalized log-posterior probability of a model.
#' @rdname logp.u
#' @description Unnormalized log-posterior of a model.
#' @param model An integer vector containing the model
#' @param X  An \eqn{n x p} matrix. Sparse matrices are supported and every
#' care is taken not to make copies of this (typically) giant matrix.
#' No need to center or scale.
#' @param y  The response vector of length \code{n}.
#' @param lam The slab precision parameter.
#' @param w The prior inclusion probability of each variable.
#' @return The unnormalized probability of the model (see detail).
#'
#' @section Details:
#' The unnormalized posterior probability of the model \code{model} is given by
#' \deqn{0.5p_{0}\log(lam) - 0.5*\log\det(A) - 0.5*(n-1)*\log(RSS) + p_{0}*\log(w/(1-w))}
#' where
#' \enumerate{
#' \item \eqn{p_{0}} = \code{length(model)}, the model size,
#' \item \eqn{Xs} = \code{scale(X[,model])},
#' \item \code{A = crossprod(Xs) + lam I},
#' \item \code{ys = scale(y)}, so that \code{sum(ys^2) = n-1} and
#' \deqn{RSS = n-1 - ys^{T}Xs A^{-1}Xs^{T}ys},
#' }
#' Obviously the computations are not as bruteforce as in pseudocode above.
#'@examples
#'\dontrun{
#'set.seed(1234)
#'X <- matrix(runif(100*200),100,200)
#'y <- 1 + X[,1:5]%*%rep(5,5) + 0.4*rnorm(100)
#' logp.u(model = 1:5,X = X,y = y,lam = 0.001,w = 1/200) # -20.2828
#'}
#'
#' @export
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
