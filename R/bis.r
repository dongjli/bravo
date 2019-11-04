#' Bayesian Iterated Screening (ultra-high, high or low dimensional).
#' @rdname bis
#' @description Perform Bayesian iterated screening in Gaussian regression models
#' @param X  An \eqn{n x p} matrix. Sparse matrices are supported and every
#' care is taken not to make copies of this (typically) giant matrix.
#' No need to center or scale.
#' @param y  The response vector of length \code{n}.
#' @param lam The slab precision parameter. Default: \code{n/p^2}.
#' @param w The prior inclusion probability of each variable. Default: \code{sqrt(n)/p}
#' as suggested by the theory of Wang et al. (2019).
#' @param criteria The stopping criteria. Could be "PP" for
#' posterior probability stopping rule, or "eBIC" for extended BIC stopping rule,
#' or "both" (default). Note that for "eBIC" the value of \code{w} is not used.
#'
#' @return A list with components
#' \item{model.pp}{An integer vector of screened model under posterior probability stopping rule.
#'  This will be null if only "eBIC" stopping criterion was used.}
#' \item{mdoel.ebic}{An integer vector of screened model under eBIC criterion. This will be NULL if
#'  only "PP" stopping criterion was used.}
#' \item{postprobs}{The sequence of posterior probabilities until the last included variable.
#'  This will be null if only "eBIC" stopping criterion was used. Here the last included variable
#'  is the last one included by either "PP" or "eBIC" if criteria="both" was selected}
#' \item{ebics}{The sequence of eBIC values until the last included variable.
#'  This will be null if only "PP" stopping criterion was used.  Here the last included variable
#'  is the last one included by either "PP" or "eBIC" if criteria="both" was selected}
#' @export
bis <- function(X,y,lam=nrow(X)/ncol(X)^2,w = sqrt(nrow(X))/ncol(X),criteria="PP")
{
  p = ncol(X)
  n = nrow(X)
  ys = scale(y)

  xbar = colMeans(X)

  stopifnot(class(X) %in% c("dgCMatrix","matrix"))

  if(class(X) == "dgCMatrix") {
    D = 1/sqrt(colMSD_dgc(X,xbar))
  }  else   {
    D = apply(X,2,sd)
    D = 1/D
  }
  Xty = D*as.numeric(crossprod(X,ys))


  max.var = n; # Intially allocate for maximum of n variables.

  model = integer(0L)
  postprob = numeric(max.var+1)
  R0 = NULL
  v0 = NULL

  postprob[1] = -0.5*(n-1)*log(n) # The posterior probability of the null model
  cat("\n Including: ")
  for(ii in 1:n)
  {
    this <- addvar(model = model,x = X, ys = ys, xty = Xty, lam = lam, w = w,
                   R0 = R0, v0 = v0,D = D,xbar = xbar)
    j = this$which.max
    if(this$logp[j] < postprob[ii])      break;
    cat(j,", ",sep = "")
    postprob[ii+1] <- this$logp[j]
    model = c(model,j)
    R0 = this$R
    v0 = this$v
  }
  cat(" Done.\n")

  return(list(model.pp = model, postprobs=postprob[1:ii]))
}
