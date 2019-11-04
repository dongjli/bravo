#' make predictions from a "bsvs" object.
#' @rdname predict.bsvs
#' @description Similar to other predict methods, this functions predicts fitted
#' values from a fitted "bsvs" object.
#' @param bsvs.obj  fitted bsvs model object
#' @param newx  Matrix of new values for x at which predictions are to be made. Must be a
#' matrix; can be sparse as in Matrix package.
#' @param model The model to be used to make predictions. Model "MAP" gives the fitted
#' values calculated using coefficients estimated from the MAP model; model "WAM" gives
#' the predicted values calculated using coefficients estimated from the WAM model. Default is "MAP".

#' @return a vector of fitted values
#'
#' @export
predict.bsvs <- function(bsvs.obj, newx, model="MAP") {
  if(model=="MAP") {
    beta <- bsvs.obj$betaMAP
  } else if (model=="WAM") {
    beta <- bsvs.obj$betaWAM
  } else {
    cat("model must be \"MAP\" or \"WAM\".")
  }
  yhat <- beta[1] + newx %*% beta[-1]
  return(yhat)
}
