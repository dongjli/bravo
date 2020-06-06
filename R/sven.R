#' Selection of variables with embedded screening using Bayesian methods (SVEN) 
#' in Gaussian linear models (ultra-high, high or low dimensional).
#' @rdname sven
#' @description SVEN is an approach to selecting variables with embedded screening 
#' using Bayesian method. It is also a variable selection method in the spirit of 
#' the stochastic shotgun search algorithm. However, by embedding a unique model 
#' based screening and using fast Cholesky updates, SVEN produces a highly scalable 
#' algorithm to explore gigantic model spaces and rapidly identify the regions of 
#' high posterior probabilities. It outputs the log (unnormalized) posterior 
#' probability of a set of best (highest probability) models. 
#' For more details, see Li et al. (2020).
#' 
#' @param X The \eqn{\code{n} \times \code{p}} covariate matrix. Sparse matrices
#' are supported and every care is taken not to make copies of this (typically)
#' giant matrix. No need to center or scale.
#' @param y The response vector of length \code{n}.
#' @param w The prior inclusion probability of each variable. Default: \code{sqrt(n)/p}.
#' @param lam The slab precision parameter. Default: \code{n/p^2}
#' as suggested by the theory of Li et al. (2020).
#' @param Ntemp The number of temperatures. Default: 3.
#' @param Tmax The maximum temperature. Default: log(log(\code{p}))+log(\code{p}).
#' @param Miter The number of iteration. Default: 50.
#' @param wam.threshold The threshold probability to select the covariates for WAM.
#' A covariate will be included in WAM if its corresponding marginal inclusion
#' probability is greater than the threshold. Default: 0.5.
#' @param log.eps The tolerance to choose the number of top models. Default: -16.
#'
#' @details
#' SVEN is developed based on a hierarchical Gaussian linear model with priors placed 
#' on the regression coefficients as well as on the model space. Degenerate spike priors 
#' on inactive variables and Gaussian slab priors on active covariates makes the posterior 
#' probability (up to a normalizing constant) of a model \eqn{P(\gamma|Y)} available in 
#' explicit form. The variable selection starts from an empty model and update the model 
#' according to the posterior probability of its neighboring models for some pre-specified 
#' number of iterations. In each iteration, the models with large probabilities are screened 
#' out in order to quickly identify the regions of high posterior probabilities. A temperature 
#' schedule is used in order to recover models with large posterior probabilities and mitigate 
#' posterior multimodality associated with variable selection models. The default maximum 
#' temperature is guided by the posterior model selection consistency asymptotics in
#' Li et al. (2020).
#' 
#' SVEN provides the maximum a posteriori (MAP) model as well as the weighted average model 
#' (WAM). WAM is obtained in the following way: (1) keep the best (highest probability) \eqn{K} 
#' models where \eqn{K} is chosen so that 
#' \eqn{\log \left(P(\gamma^{(K)}|y)/P(\gamma^{(1)}|y)\right) > \code{log.eps}};
#' (2) assign the weights \deqn{w_i = P(\gamma^{(i)}|y)/\sum_{k=1}^K P(\gamma^{(k)}|y)}
#' to the model \eqn{\gamma^{(i)}}; (3) define the approximate marginal inclusion probabilities 
#' for the \eqn{j}th variable as \deqn{\hat\pi_j = \sum_{k=1}^K w_k I(\gamma^{(k)}_j = 1).} 
#' Then, the WAM is defined as the model containing variables \eqn{j} with 
#' \eqn{\hat\pi_j > \code{wam.threshold}}. SVEN also provides all the top \eqn{K} models which
#' are stored in an \eqn{p \times K} sparse matrix, along with their corresponding log (unnormalized) 
#' posterior probabilities. 
#' 
#' 
#' @return A list with components
#' \item{model.map}{A vector of indices corresponding to the selected variables
#' in the MAP model.}
#' \item{model.wam}{A vector of indices corresponding to the selected variables
#' in the WAM.}
#' \item{model.top}{A sparse matrix storing the top models.}
#' \item{beta.map}{The ridge estimator of regression coefficients in the MAP model.}
#' \item{beta.wam}{The ridge estimator of regression coefficients in the WAM.}
#' \item{mip.map}{The marginal inclusion probabilities of the variables in the MAP model.}
#' \item{mip.wam}{The marginal inclusion probabilities of the variables in the WAM.}
#' \item{pprob.map}{The log (unnormalized) posterior probability corresponding
#' to the MAP model.}
#' \item{pprob.top}{A vector of the log (unnormalized) posterior probabilities
#' corresponding to the top models.}
#' \item{stats}{Additional statistics.}
#'
#' @author Dongjin Li and Somak Dutta\cr Maintainer:
#' Dongjin Li <dongjl@@iastate.edu>
#' @references
#' @examples
#' n = 400; p = 20000; nonzero = 15
#' trueidx <- 1:15
#' nonzero.value <- 3
#' TrueBeta = numeric(ncovar)
#' TrueBeta[trueidx] <- nonzero.value
#' 
#' # Generate a matrix with group structure
#' z0 <- matrix(rnorm(n*(p-nonzero)), n, p-nonzero)
#' zz = matrix(rnorm(3*n), n, 3)
#' z = cbind(matrix(rep(zz, nonzero/3), nrow = n), z0)
#' e = cbind(matrix(rnorm(n*nonzero, mean = 0, sd = 0.1), nrow = n), matrix(0, n, p-nonzero))
#' X = z + e
#' y <- 0.5 + X %*% TrueBeta + rnorm(n)
#' res <- sven(X = X, y = y, lam = 200, w = 0.02, Ntemp = 9, Miter = 200)
#' res$model.map # the MAP model
#' res$model.wam # the WAM
#' res$mip.map # the marginal inclusion probabilities of the variables in the MAP model
#' res$mip.wam # the marginal inclusion probabilities of the variables in the WAM
#' res$pprob.top # the log (unnormalized) posterior probabilities corresponding to the top models.
#' 
#' res$beta.map # the ridge estimator of regression coefficients in the MAP model 
#' res$beta.wam # the ridge estimator of regression coefficients in the WAM
#' @export
sven <- function(X, y, w = sqrt(nrow(X))/ncol(X), lam = nrow(X)/ncol(X)^2, Ntemp = 3,
                 Tmax = (log(log(ncol(X)))+log(ncol(X))), Miter = 50, wam.threshold = 0.5, log.eps = -16) {
  result <- list()
  k=20

  n <- length(y)
  ncovar <- ncol(X)
  ys = scale(y)
  xbar = colMeans(X)

  stopifnot(class(X)[1] %in% c("dgCMatrix","matrix"))
  if(class(X)[1] == "dgCMatrix") {
    D = 1/sqrt(colMSD_dgc(X,xbar))
  }  else   {
    D = apply(X,2,sd)
    D = 1/D
  }
  Xty = D*as.numeric(crossprod(X,ys))

  logp <- numeric(Miter * (Ntemp))
  RSS <- numeric(Miter * (Ntemp))
  size <- integer(Miter * (Ntemp))
  indices <- integer(Miter*100)

  cat("temperature = 1\n")
  o <- sven.notemp(X, ys, Xty, lam, w, k, D, xbar, n, ncovar, Miter)
  logp.best <- o$bestlogp
  r.idx.best <- o$bestidx
  logp[1:Miter] <- o$currlogp
  RSS[1:Miter] <- o$currRSS
  size[1:Miter] <- o$modelsizes
  ed <- sum(size)
  indices[1:ed] <- o$curridx

  stepsize <- Tmax/(Ntemp-1)
  for (t in 1:(Ntemp-1)) {
    cat("temperature =", t*stepsize, "\n")
    o <- sven.temp(X, ys, Xty, lam, w, k, D, xbar, t, stepsize=stepsize, logp.best, r.idx.best, n, ncovar, Miter)
    logp.best <- o$bestlogp
    r.idx.best <- o$bestidx
    logp[(Miter*t+1):(Miter*(t+1))] <- o$currlogp
    RSS[(Miter*t+1):(Miter*(t+1))] <- o$currRSS
    size[(Miter*t+1):(Miter*(t+1))] <- o$modelsizes
    indices[(ed+1):(ed+sum(o$modelsizes))] <- o$curridx
    ed <- ed + sum(o$modelsizes)
  }
  indices <- indices[indices>0]
  cumsize <- cumsum(size)
  modelSparse <- sparseMatrix(i=indices,p = c(0,cumsize),index1 = T,dims = c(ncovar,length(logp)), x = T)

  logp.uniq1 <- unique(logp)
  for (i in 1:(length(logp.uniq1)-1)) {
    for (j in 2:length(logp.uniq1)) {
      if (abs(logp.uniq1[i]-logp.uniq1[j]) < 1e-10) {
        logp.uniq1[j] <- logp.uniq1[i]
      }
    }
  }
  logp.uniq <- unique(logp.uniq1)
  logp.top <- sort(logp.uniq[(logp.best-logp.uniq)<(-log.eps)], decreasing = T)
  cols.top <- unlist(lapply(logp.top, FUN=function(x){which(x==logp)[1]}))
  size.top <- size[cols.top]
  RSS.top <- RSS[cols.top]
  model.top <- modelSparse[, cols.top, drop=F]

  logp.top1 <- logp.top-logp.best
  weight <- exp(logp.top1)/sum(exp(logp.top1))

  beta.est <- matrix(0, (ncovar+1), length(cols.top))
  for(i in 1:length(cols.top)){
   if (size.top[i]==0){
     beta <- mean(y)
     beta.est[1, i] <- beta
   } else {
     m_i = model.top[, i]
     x.est <- cbind(rep(1, n), scale(X[, m_i], center = xbar[m_i], scale = 1/D[m_i]))
     beta <- solve(crossprod(x.est) + lam*diag(c(0, rep(1, size.top[i]))), crossprod(x.est, y))
     beta.est[c(T, m_i), i] <- c(beta[1]-sum(beta[-1]*xbar[m_i]*D[m_i]), beta[-1] * D[m_i])
   }
  }

  beta.MAP <- beta.est[, 1]
  beta.WAM <- rowSums(beta.est%*%diag(weight, nrow=length(size.top)))

  MIP = rowSums(model.top%*%Diagonal(length(weight), weight))
  model.WAM <- sort(which(MIP >= wam.threshold))
  model.MAP <- sort(r.idx.best)
  MIP.MAP <- MIP[model.MAP]
  MIP.WAM <- MIP[model.WAM]

  mtop.idx <- apply(model.top, 2, which)
  model.union <- Reduce(union, mtop.idx)
  Xm <- sparseMatrix(i = integer(0), j = integer(0), dims = c(n, ncovar))
  Xm <- as(Xm, "dgCMatrix")
  Xm[, model.union] <- X[, model.union, drop=F]

  result$model.map <- model.MAP
  result$model.wam <- model.WAM
  result$model.top <- model.top
  result$beta.map <- beta.MAP
  result$beta.wam <- beta.WAM
  result$mip.map <- MIP.MAP
  result$mip.wam <- MIP.WAM
  result$pprob.map <- logp.best
  result$pprob.top <- logp.top
  result$stats <- list(RSS.top = RSS.top, weights = weight, Xm = Xm, y = y, lam = lam)

  return(result)
}

