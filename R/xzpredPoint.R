predPoint.xz <- function(object, x1new, x2new, model) {
  if(model=="MAP") {
    beta <- object$beta.map
  } else if (model=="WAM") {
    beta <- object$beta.wam
  } else {
    cat("model must be \"MAP\" or \"WAM\".")
  }
  intercept <- beta[1]
  coeffs <- beta[-1]
  beta1 <- coeffs[1:ncol(x1new)] ; beta2 <- coeffs[(ncol(x1new)+1):(ncol(x1new)+ncol(x2new))]
  y.pred <- beta[1] + x1new %*% beta1 + x2new %*% beta2
  return(y.pred)
}

#verified, code is perfect
