#' Bayesian Iterated Screening (ultra-high, high or low dimensional).
#' @rdname bits
#' @description Perform Bayesian iterated screening in Gaussian regression models
#' @param X  An \eqn{n\times p} matrix. Sparse matrices are supported and every
#' care is taken not to make copies of this (typically) giant matrix.
#' No need to center or scale.
#' @param y  The response vector of length \code{n}.
#' @param lam The slab precision parameter. Default: \code{1}.
#' @param w The prior inclusion probability of each variable. Default: \code{1/2}.
#' @param pp Boolean: If \code{FALSE} (default) the algorithm stops after including \code{max.var} many variables.
#' If true, the posterior probability stopping rule is used.
#' @param max.var The maximum number of variables to be included.
#' @param verbose If \code{TRUE} (default) will show the variable index included in each iteration.
#' @return A list with components
#' \item{model.pp}{An integer vector of the screened model.}
#' \item{postprobs}{The sequence of posterior probabilities until the last included variable.}
#' \item{lam}{The value of lam, the slab precision parameter.}
#' \item{w}{The value of w, the prior inclusion probability.}
#' @references Wang, R., Dutta, S., Roy, V. (2021) Bayesian iterative screening in ultra-high dimensional 
#' settings. https://arxiv.org/abs/2107.10175
#' @examples
#' n=50; p=100;
#' TrueBeta <- c(rep(5,3),rep(0,p-3))
#' 
#' rho <- 0.6
#' x1 <- matrix(rnorm(n*p), n, p)
#' X <- sqrt(1-rho)*x1 + sqrt(rho)*rnorm(n)
#' y <- 0.5 + X %*% TrueBeta + rnorm(n)
#' res<-bits(X,y, pp=TRUE)
#' res$model.pp # the vector of screened model
#' res$postprobs # the log (unnormalized) posterior probabilities corresponding to the model.pp.
#' @export

bits <- function(X,y,lam=1, w=0.5, pp = FALSE,max.var = nrow(X), verbose = TRUE)
{
  p = ncol(X)
  n = nrow(X)
  ys = scale(y)

  pp = (pp != 0);



  stopifnot(class(X)[1] %in% c("dgCMatrix","matrix"))
  stopifnot(max.var <= n && max.var < p)
  if(!is.null(w)) stopifnot(w > 0 && w < 1)
  logw = ifelse(is.null(w), 0, log(w/{1-w}))

  xbar = colMeans(X)

  if(class(X)[1] == "dgCMatrix") {
    D = 1/sqrt(colMSD_dgc(X,xbar))
  }  else   {
    D = apply(X,2,sd)
    D = 1/D
  }

  xty = D*as.numeric(crossprod(X,ys))


  yty <- n - 1
  xtx <- n - 1


  model = integer(n)
  postprob = numeric(max.var+1)


  R = matrix(NA,max.var,max.var)
  sumv2 = 0
  logdetR = 0;
  z = numeric(p)
  u   = numeric(p)
  v = numeric(max.var)


  postprob[1] = -0.5*(n-1)*log(yty) # The posterior probability of the null model
  if(verbose) cat("\n Including: ")

  # First variable
  b0 = sqrt(xtx + lam)
  logdetR = log(b0)

  logp <- as.numeric(0.5*log(lam)-logdetR - 0.5*(n-1)*log(yty - (xty/b0)^2) + logw)


  j = which.max(logp)
  if(verbose) cat(j)
  model[1] = j;
  postprob[2] = logp[j]
  if(postprob[2]<postprob[1] && pp){
    if(verbose) cat(" Done.\n")
    return(list(model.pp = NULL, postprobs=postprob[1],lam=lam))
  }

  # Need to do the second variable by hand
  if(max.var >= 2)
  {
    R[1,1] = b0;
    xjc = (X[,j] - xbar[j])*D[j]
    v[1] = xty[j]/b0
    sumv2 = v[1]^2

    S = D*crossprod(X,xjc)/b0
    z = S^2

    w1 = sqrt(xtx+lam - z)
    u = {xty - v[1]*S}/w1;

    RSS = yty - sumv2 - u^2
    RSS[j] = Inf

    logp = as.numeric(0.5*2*log(lam) - logdetR - log(w1) - 0.5*{n-1}*log(RSS) + 2*logw)

    j = which.max(logp)
    if(verbose) cat(", ",j)
    model[2] = j
    postprob[3] = logp[j]
    if(postprob[3]<postprob[2] && pp){
      if(verbose) cat(" Done.\n")
      return(list(model.pp = model[1:2], postprobs=postprob[1],lam=lam))
    }

  }



  if(max.var > 3)
  {
    for(ii in 3:max.var)
    {

      model.prev = model[1:{ii-2}]


      xjc = (X[,j] - xbar[j])*D[j]


      X1 = X[,model.prev,drop=FALSE]
      D1 = D[model.prev]
      Xbar1 = xbar[model.prev]
      a1 = backsolve(R,D1*crossprod(X1,xjc),k = ii-2,transpose = T)
      # print(a1)
      b1 = w1[j]
      logdetR = logdetR + log(b1)

      v[ii-1] = u[j];
      sumv2 = sumv2 + u[j]^2


      temp1 = D1*backsolve(R,a1,transpose = FALSE,k = ii-2)
      temp2 = as.numeric(xjc - X1 %*% temp1);
      temp2 = temp2 - mean(temp2)

      eta = D*crossprod(X,temp2)
      eta = eta/b1


      z = z + eta^2

      temp3  = xtx + lam - z;
      temp3[model[1:{ii-1}]] = 1;
      w2 = sqrt(temp3)
      u = {u*w1 - u[j]*eta}/w2;


      RSS = yty - sumv2 - u^2
      RSS[model[1:{ii-1}]] = Inf

      logp = as.numeric(0.5*ii*log(lam) - logdetR - log(w2) - 0.5*{n-1}*log(RSS) + ii*logw)

      j = which.max(logp)

      if(verbose) cat(", ",j)
      postprob[ii+1] <- logp[j]
      model[ii] = j
      if(postprob[ii+1]<postprob[ii] && pp){
        if(verbose) cat(" Done.\n")
        return(list(model.pp = model[1:(ii-1)], postprobs=postprob[1:ii],lam=lam,w=w))
      }

      w1 = w2;
      if(ii <= max.var)
      {
        R[1:{ii-2},ii-1] = a1;
        R[ii-1,ii-1] = b1;
      }

    }
  }
  if(verbose) cat(" Done.\n")

  return(list(model.pp = model[1:ii], postprobs=postprob,lam=lam,w=w))
}

