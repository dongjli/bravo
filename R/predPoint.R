predPoint <- function(object, newx, model="MAP") {
  if(model=="MAP") {
    beta <- object$beta.map
  } else if (model=="WAM") {
    beta <- object$beta.wam
  } else {
    cat("model must be \"MAP\" or \"WAM\".")
  }
  y.pred <- beta[1] + newx %*% beta[-1]
  return(y.pred)
}
