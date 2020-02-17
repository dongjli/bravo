calsigma <- function(Rsquared, trueidx, truebeta, type, rho=NULL, 
                     factor.mat=NULL, ngroup=NULL, delta=NULL) {
  p <- length(trueidx)
  truebeta <- truebeta[trueidx] 
  var.x <- matrix(0, nrow = p, ncol=p)
  
  if (type == "compSym"){
    diag(var.x) <- 1
    var.x[upper.tri(var.x) | lower.tri(var.x)] <- rho
    var.xb <- truebeta %*% var.x %*% truebeta
  } 
  
  if (type == "ar") {
    for (i in 1:p){
      for (j in i:p){
        var.x[i, j] <- rho^abs(trueidx[i]-trueidx[j])
      }
    }
    var.x[lower.tri(var.x)] <- (t(var.x))[lower.tri(var.x)]
    var.xb <- truebeta %*% var.x %*% truebeta
  }
  
  if (type == "factor") {
    var.xb = sum((factor.mat[,1:p] %*% truebeta)^2) + sum(truebeta^2) 
  }
  
  if (type == "group") {
    temp = 0
    for (i in 1:(p-ngroup)){
      if (i %% ngroup == 1) temp <- temp+1
      for (m in 0:(p/ngroup-temp)){
        var.x[i, i+3*m] <- 1
      }
    }
    diag(var.x) <- 1 + delta^2
    var.x[lower.tri(var.x)] <- (t(var.x))[lower.tri(var.x)]
    var.xb <- truebeta %*% var.x %*% truebeta
  }
  
  if (type == "extrmCor" | type == "iid") {
    var.xb <- sum(truebeta^2)
  }
  
  if (type == "extrmCorgroup") {
    temp = 0
    for (i in 1:(p-ngroup)){
      if (i %% ngroup == 1) temp <- temp+1
      for (m in 0:(p/ngroup-temp)){
        var.x[i, i+3*m] <- 1
      }
    }
    diag(var.x) <- 1 + delta^2
    var.x[lower.tri(var.x)] <- (t(var.x))[lower.tri(var.x)]
    var.xb <- truebeta %*% var.x %*% truebeta
    #var.xb <- 1/2 * var.xb + 1/2 * sum(truebeta^2)
    var.xb <- 1/(2+delta^2) * var.xb + 1/(2+delta^2) * sum(truebeta^2)
  }

    #var.xb = (9*2/3 + 27 * 0.1) * truebeta[1]^2
  
  
  if (type == "extrmCorgroup2") {
    # temp = 0
    # for (i in 1:(p-ngroup)){
    #   if (i %% ngroup == 1) temp <- temp+1
    #   for (m in 0:(p/ngroup-temp)){
    #     var.x[i, i+3*m] <- 1
    #   }
    # }
    # diag(var.x) <- 1 + delta^2
    # var.x[lower.tri(var.x)] <- (t(var.x))[lower.tri(var.x)]
    # var.xb <- truebeta %*% var.x %*% truebeta
    # #var.xb <- 1/2 * var.xb + 1/2 * sum(truebeta^2)
    # var.xb <- 1/(1+delta^2) * var.xb 
    var.xb=(27+9*delta^2)/(1+delta^2) * truebeta[1]^2
  }
  
  #var.xb = (9*2/3 + 27 * 0.1) * truebeta[1]^2
  sig2 <- var.xb/Rsquared - var.xb
  return(as.numeric(sqrt(sig2)))
}

calsigma <- cmpfun(calsigma)

gen_y <- function(x, sigma, truebeta) {
  y = 0.5 + x %*% truebeta + sigma * rnorm(dim(x)[1])
}

gen_y <- cmpfun(gen_y)

