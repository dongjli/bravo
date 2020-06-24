#' Make predictions from a fitted "sven" object.
#' @rdname predict.sven
#' @description This function makes point predictions and computes prediction intervals
#' from a fitted "sven" object.
#' @param object A fitted "sven" object
#' @param Xnew Matrix of new values for \code{X} at which predictions are to be made. Must be a matrix;
#' can be sparse as in Matrix package.
#' @param model The model to be used to make point predictions. Model "MAP" gives the predictions calculated
#' using coefficients estimated from the MAP model; model "WAM" gives the predictions calculated using
#' coefficients estimated from the WAM. Default: "MAP".
#' @param interval Logical; if \code{TRUE}, the function gives the prediction interval in addition to the point
#' prediction. Default: \code{FALSE}.
#' @param MC only required if \code{interval} = \code{TRUE}. if \code{TRUE}, the function gives the Monte Carlo prediction
#' interval. Default: \code{FALSE}.
#' @param Nsim only required if \code{MC} = \code{TRUE}. The Monte Carlo sample size. Default: 10000.
#' @param return.draws only required if \code{MC} = \code{TRUE}. if \code{TRUE}, the Monte Carlo draws is returned.
#' Default: \code{FALSE}.
#' @param conf.level Confidence level of the interval. Default: 0.95.
#' @param alpha Type one error rate. Default: 1-\code{conf.level}.
#' 
#' @details
#' 
#' @return The object returned depends on "interval" argument. If \code{interval} = \code{TRUE}, the object is an
#' \eqn{\code{ncol(Xnew)}\times 3} matrix with the point predictions in the first column and the prediction
#' in the second and third columns; otherwise, the object is an \eqn{\code{ncol(Xnew)}\times 1} vector of the
#' point predictions.
#' 
#' if return.draws is \code{TRUE}, a list with the following components is returned:
#' \item{prediction}{vector or matrix as above}
#' \item{mc.draws}{an \eqn{\code{ncol(Xnew)} \times \code{Nsim}} matrix of the Monte Carlo samples}
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

predict.sven <- function(object, Xnew, model="WAM", interval = FALSE, MC = FALSE, return.draws = FALSE,
                         Nsim = 10000, conf.level = 0.95, alpha = 1-conf.level) {
  y.pred <- predPoint(object = object, newx = Xnew, model = model)
  if (!interval){
    return(y.pred)
  } else {
    intv <- predIntv(object = object, model = model, Xnew = Xnew, MC = MC, Nsim = Nsim, 
                     conf.level = conf.level, alpha = alpha)
    if (return.draws){
      prediction <- cbind(y.pred, intv$ci)
      colnames(prediction) <- c("y.pred", colnames(intv$ci))
      return(list(prediction = prediction, mc.draws = intv$mc.draws))
    } else {
      prediction <- cbind(y.pred, intv)
      colnames(prediction) <- c("y.pred", colnames(intv))
      return(prediction)
    }
  }
}
