predPoint <- function(object, Xnew, model) {
  if(model=="MAP") {
    beta <- object$beta.map
  } else if (model=="WAM") {
    beta <- object$beta.wam
  } else {
    cat("model must be \"MAP\" or \"WAM\".")
  }
  y.pred <- beta[1] + Xnew %*% beta[-1]
  return(y.pred)
}
