#roxygen removed

xzpredict.sven <- function(object, newdata1, newdata2, model = c("WAM", "MAP"), interval = c("none", "MC", "Z"),
                         return.draws = FALSE, Nsim = 10000, level = 0.95, alpha = 1-level, ...) {
  if (alpha < 0 | alpha > 1) {
    stop("alpha has to be between 0 and 1")
  }
  if (level < 0 | level > 1) {
    stop("level has to be between 0 and 1")
  }
  interval <- match.arg(interval)
  model <- match.arg(model)
  y.pred <- predPoint.xz(object = object, x1new = newdata1, x2new=newdata2, model = model)
  if (interval == "none"){
    return(y.pred)
  } else {
    intv <- predIntv.xz(object = object, model = model, x1new = newdata1, x2new=newdata2,
                     interval = interval, Nsim = Nsim,
                     return.draws = return.draws, conf.level = level, alpha = alpha)
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
