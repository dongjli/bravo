#' Selection of variables with embedded screening using Bayesian methods (SVEN) 
#' in Gaussian linear models (ultra-high, high or low dimensional).
#' @rdname sven
#' @description SVEN is an approach to selecting variables with embedded screening 
#' using a Bayesian hierarchical model. It is also a variable selection method in the spirit of 
#' the stochastic shotgun search algorithm. However, by embedding a unique model 
#' based screening and using fast Cholesky updates, SVEN produces a highly scalable 
#' algorithm to explore gigantic model spaces and rapidly identify the regions of 
#' high posterior probabilities. It outputs the log (unnormalized) posterior 
#' probability of a set of best (highest probability) models. 
#' For more details, see Li et al. (2023, https://doi.org/10.1080/10618600.2022.2074428)
#' 
#' @param X The \eqn{n\times p} covariate matrix or list of two matrices without intercept. 
#' The following classes are supported: \code{matrix} and \code{dgCMatrix}. Every care is taken not to make copies of these (typically)
#' giant matrices. No need to center or scale these matrices manually. Scaling is performed implicitly and 
#' regression coefficient are returned on the original scale. Typically, in a combined GWAS-TWAS type 
#' analysis, \code{X[[1]]} should be a sparse matrix and \code{X[[2]]} should be a dense matrix.
#' @param y The response vector of length \eqn{n}. No need to center or scale.
#' @param w The prior inclusion probability of each variable. Default: NULL, whence it is set as
#' \eqn{\sqrt{n}/p} if \eqn{X} is a matrix. Or \eqn{(\sqrt{n}/p_1,\sqrt{n}/p_2)} if $X$ is a list of 
#' two matrices with \eqn{p_1} and \eqn{p_2} columns.
#' @param lam The slab precision parameter. Default: NULL, whence it is set as \eqn{n/p^2} for 
#' as suggested by the theory of Li et al. (2023). Similarly, it's a vector of length two with values
#' \eqn{\sqrt{n}/P_1^2} and \eqn{\sqrt{n}/p_2^2} when \code{X} is a list.
#' @param Ntemp The number of temperatures. Default: 10.
#' @param Tmax The maximum temperature. Default: \eqn{\log\log p+\log p}.
#' @param Miter The number of iterations per temperature. Default: \code{50}.
#' @param wam.threshold The threshold probability to select the covariates for WAM.
#' A covariate will be included in WAM if its corresponding marginal inclusion
#' probability is greater than the threshold. Default: 0.5.
#' @param log.eps The tolerance to choose the number of top models. See detail. Default: -16.
#' @param L The minimum number of neighboring models screened. Default: 20.
#' @param verbose If \code{FALSE}, the function prints the current temperature SVEN is at; the default is TRUE. 
#'
#' @details
#' SVEN is developed based on a hierarchical Gaussian linear model with priors placed 
#' on the regression coefficients as well as on the model space as follows:
#' \deqn{y | X, \beta_0,\beta,\gamma,\sigma^2,w,\lambda \sim N(\beta_01 + X_\gamma\beta_\gamma,\sigma^2I_n)}
#' \deqn{\beta_i|\beta_0,\gamma,\sigma^2,w,\lambda \stackrel{indep.}{\sim} N(0, \gamma_i\sigma^2/\lambda),~i=1,\ldots,p,}
#' \deqn{(\beta_0,\sigma^2)|\gamma,w,p \sim p(\beta_0,\sigma^2) \propto 1/\sigma^2}
#' \deqn{\gamma_i|w,\lambda \stackrel{iid}{\sim} Bernoulli(w)}
#' where \eqn{X_\gamma} is the \eqn{n \times |\gamma|} submatrix of \eqn{X} consisting of 
#' those columns of \eqn{X} for which \eqn{\gamma_i=1} and similarly, \eqn{\beta_\gamma} is the 
#' \eqn{|\gamma|} subvector of \eqn{\beta} corresponding to \eqn{\gamma}.
#' Degenerate spike priors on inactive variables and Gaussian slab priors on active 
#' covariates makes the posterior 
#' probability (up to a normalizing constant) of a model \eqn{P(\gamma|y)} available in 
#' explicit form (Li et al., 2020).
#' 
#' The variable selection starts from an empty model and updates the model 
#' according to the posterior probability of its neighboring models for some pre-specified 
#' number of iterations. In each iteration, the models with small probabilities are screened 
#' out in order to quickly identify the regions of high posterior probabilities. A temperature 
#' schedule is used to facilitate exploration of models separated by valleys in the posterior 
#' probability function, thus mitigate posterior multimodality associated with variable selection models.
#' The default maximum temperature is guided by the asymptotic posterior model selection consistency results
#' in Li et al. (2020).
#' 
#' SVEN provides the maximum a posteriori (MAP) model as well as the weighted average model 
#' (WAM). WAM is obtained in the following way: (1) keep the best (highest probability) \eqn{K} 
#' distinct models \eqn{\gamma^{(1)},\ldots,\gamma^{(K)}} with 
#' \deqn{\log P\left(\gamma^{(1)}|y\right) \ge \cdots \ge \log P\left(\gamma^{(K)}|y\right)}
#' where \eqn{K} is chosen so that 
#' \eqn{\log \left\{P\left(\gamma^{(K)}|y\right)/P\left(\gamma^{(1)}|y\right)\right\} > \code{log.eps}};
#' (2) assign the weights \deqn{w_i = P(\gamma^{(i)}|y)/\sum_{k=1}^K P(\gamma^{(k)}|y)}
#' to the model \eqn{\gamma^{(i)}}; (3) define the approximate marginal inclusion probabilities 
#' for the \eqn{j}th variable as \deqn{\hat\pi_j = \sum_{k=1}^K w_k I(\gamma^{(k)}_j = 1).} 
#' Then, the WAM is defined as the model containing variables \eqn{j} with 
#' \eqn{\hat\pi_j > \code{wam.threshold}}. SVEN also provides all the top \eqn{K} models which
#' are stored in an \eqn{p \times K} sparse matrix, along with their corresponding log (unnormalized) 
#' posterior probabilities. 
#' 
#' When \code{X} is a list with two matrices, say, \code{W} and \code{Z}, the above method is extended 
#' to \code{ncol(W)+ncol(Z)} dimensional regression. However, the hyperparameters \code{lam} and \code{w}
#' are chosen separately for the two matrices, the default values being  \code{nrow(W)/ncol(W)^2}
#' and \code{nrow(Z)/ncol(Z)^2} for \code{lam} and \code{sqrt(nrow(W))/ncol(W)} and
#'  \code{sqrt(nrow(Z))/ncol(Z)} for \code{w}.
#'  
#'  The marginal inclusion probabities can be extracted by using the function \code{mip}.
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
#' @author Dongjin Li, Debarshi Chakraborty, and Somak Dutta\cr Maintainer:
#' Dongjin Li <liyangxiaobei@@gmail.com>
#' 
#' @seealso [mip.sven()] for marginal inclusion probabilities, [predict.sven()](via [predict()]) for prediction for .
#' 
#' 
#' @references Li, D., Dutta, S., and Roy, V. (2023). Model based screening embedded Bayesian variable 
#' selection for ultra-high dimensional settings. Journal of Computational and Graphical Statistics, 
#' 32(1), 61-73.
#' @examples
#' \donttest{
#' n <- 50; p <- 100; nonzero <- 3
#' trueidx <- 1:3
#' truebeta <- c(4,5,6)
#' X <- matrix(rnorm(n*p), n, p) # n x p covariate matrix
#' y <- 0.5 + X[,trueidx] %*% truebeta + rnorm(n)
#' res <- sven(X=X, y=y)
#' res$model.map # the MAP model
#' 
#' 
#' Z <- matrix(rnorm(n*p), n, p) # another covariate matrix
#' y2 = 0.5 + X[,trueidx] %*% truebeta  + Z[,1:2] %*% c(-2,-2) + rnorm(n)
#' res2 <- sven(X=list(X,Z), y=y2)
#' }
#' 
#' @export
sven <- function(X, y, w = NULL, lam = NULL, Ntemp = 10,
                 Tmax = NULL, Miter = 50, wam.threshold = 0.5, 
                 log.eps = -16, L = 20, verbose = FALSE) {
  if(is.list(X))
  {
    if(length(X) != 2)
      stop("Input X must be a matrix or a list of two matrices")
      X1 = X[[1]]
      X2 = X[[2]]
      stopifnot(class(X1)[1] %in% c("dgCMatrix","matrix"))
      stopifnot(class(X2)[1] %in% c("dgCMatrix","matrix"))
      
      if(is.null(w))
      {
        w1 = sqrt(nrow(X1))/ncol(X1)
        w2 = sqrt(nrow(X2))/ncol(X2)
      } else {
        w1 = w[1]
        w2 = w[2]
        stopifnot(is.numeric(w) && length(w) == 2 && all(w > 0) && all(w < 1))
      }
      if(is.null(lam)){
        lam1 = nrow(X1)/ncol(X1)^2
        lam2 = nrow(X2)/ncol(X2)^2
      } else {
        lam1 = lam[1]
        lam2 = lam[2]
        stopifnot(is.numeric(lam) && length(lam) == 2 && all(lam > 0))
      }
      Tmax = log(log(ncol(X1)+ncol(X2))) + log(ncol(X1) + ncol(X2))
      result = sven.xz(x1 = X1,x2 = X2,y = y,w1 = w1,w2 = w2,lam1 = lam1,lam2 = lam2,
                        Ntemp = Ntemp, Tmax = Tmax, Miter = Miter, wam.threshold = wam.threshold,
                        log.eps = log.eps, L = L, verbose = verbose)
  } else   {
    stopifnot(class(X)[1] %in% c("dgCMatrix","matrix"))
    if(is.null(w)) 
      w = sqrt(nrow(X))/ncol(X)
    
    stopifnot(is.numeric(w) && {w > 0} && {w < 1})
    if(is.null(lam))
      lam = nrow(X)/ncol(X)^2
    stopifnot(is.numeric(lam) && {lam > 0})
    Tmax = log(log(ncol(X))) + log(ncol(X))
    
    result = sven.x(X,y,w = w, lam = lam, Ntemp = Ntemp, Tmax = Tmax, Miter = Miter,
                    wam.threshold = wam.threshold,
                    log.eps = log.eps, L = L, verbose = verbose)
  }
  
  result$stats$islist = is.list(X)
  if(is.list(X))
  {
    result$stats$colnamesX = colnames(X[[1]])
    result$stats$colnamesZ = colnames(X[[2]])
  } else   {
    result$stats$colnamesX = colnames(X)
  }
  return(result)
}


