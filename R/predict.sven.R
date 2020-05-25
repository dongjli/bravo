#' Bayesian Variable Selection (ultra-high, high or low dimensional).
#' @rdname predict.sven
#' @description This function makes point predictions and computes prediction intervals
#' from a fitted "sven" object.
#' @param object A fitted "sven" object
#' @param Xnew Matrix of new values for \code{X} at which predictions are to be made. Must be a matrix;
#' can be sparse as in Matrix package
#' @param model The model to be used to make point predictions. Model "MAP" gives the predictions calculated
#' using coefficients estimated from the MAP model; model "WAM" gives the predictions calculated using
#' coefficients estimated from the WAM. Default: "MAP".
#' @param interval Logical; if true, the function gives the prediction interval in addition to the point
#' prediction. Default: TRUE.
#' @param MC only requeired if interval = TRUE. if TRUE, the function gives the Monte Carlo prediction
#' interval. Default: TRUE.
#' @param Nsim only required if MC = TRUE. The Monte Carlo sample size. Default: 10000.
#' @param conf.level Confidence level of the interval. Default: 0.95.
#' @param alpha Type one error rate. Default: 1-conf.level.
#' @return The object returned depends on "interval" argumemt. If interval = TRUE, the object is an
#' \eqn{\code{ncol(X)}\times 3} matrix with the point predicitons in the first column and the prediciton
#' in the second and third columns; otherwise, the object is an \eqn{\code{ncol(X)}\times 1} vector of the
#' point predicitons.
#' @author Dongjin Li and Somak Dutta\cr Maintainer:
#' Dongjin Li <dongjl@@iastate.edu>
#' @references
#' @examples
#' n = 400; p = 20000; nonzero = 5
#' trueidx <- 1:5
#' nonzero.value <- c(0.50, 0.75, 1.00, 1.25, 1.50)
#' TrueBeta = numeric(p)
#' TrueBeta[trueidx] <- nonzero.value
#'
#' X <- matrix(rnorm(n*p), n, p)
#' y <- 0.5 + X %*% TrueBeta + rnorm(n)
#' res <- sven(X=X, y=y)
#' newx <- matrix(rnorm(n*p), 200, p)
#' # predicted values at a new data matrix using MAP model
#' yhat <- predict(object = res, Xnew = newx, model = "MAP")
#' # 95% Monte Carlo prediction interval
#' MC.interval <- predict(object = res, Xnew = newx, interval = TRUE, MC = TRUE, conf.level = 0.95)
#' # 95% Z-prediction interval
#' Z.interval <- predict(object = res, Xnew = newx, interval = TRUE, conf.level = 0.95)
#' @method predict sven
#' @export

predict.sven <- function(object, Xnew, model="MAP", interval = FALSE, MC = FALSE,
                         Nsim = 10000, conf.level = 0.95, alpha = 1-conf.level) {
  yhat <- predPoint(object = object, newx = Xnew, model = model)
  if (!interval){
    return(yhat)
  } else {
    CI <- predIntv(object = object, Xnew = Xnew, MC = MC, Nsim = Nsim, alpha = alpha)
    return(cbind())
  }
}
