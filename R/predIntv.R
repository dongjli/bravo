predIntv <- function(object, model, Xnew, interval, Nsim, return.draws, conf.level, alpha) {
  X <- object$stats$Xm
  y <- object$stats$y
  lam <- object$stats$lam

  Models <- object$model.top
  weights <- object$stats$weights
  if (model == "MAP") weights[-1] <- 0
  RSS <- object$stats$RSS.top

  n <- nrow(X)
  p <- ncol(X)
  n.new <- nrow(Xnew)
  Xbar <- colMeans(X)
  if(class(X)[1] == "dgCMatrix") {
    D = 1/sqrt(colMSD_dgc(X, Xbar))
  }  else   {
    D = apply(X,2,sd)
    D = 1/D
  }
  var_y <- as.vector(var(y))
  mean_y <- mean(y)
  
  if(interval == "MC") {
    model.sim <- sample(1:length(weights), Nsim, prob = weights, replace = T)
    uni.model <- unique(model.sim)
    getR_a <- function(x, y, n, model.idx, lam) {
      model <- which(Models[, model.idx])
      model.size <- length(model)
      xmat <- scale(x[, model, drop=F])
      A <- as.matrix(crossprod(xmat) + lam*diag(model.size))
      R = chol(A)
      a = backsolve(R, crossprod(xmat, y), transpose = T)
      return(list(R=R, a=a))
    }

    gen_newy <- function(nnew, xnew, RSS, mean_y, var_y, xbar, Ra, model.idx) {
      ### simulate sigma^2
      sig2 <- 1/rgamma(1, shape=(n-1)/2, rate=var_y*RSS[model.idx]/2)
      ### simulate regression coefficients
      model <- which(Models[, model.idx])
      model.size <- length(model)
      R <- Ra$R
      a <- Ra$a
      sig <- sqrt(sig2)
      beta.sc <- backsolve(R, a+rnorm(model.size, 0, sig))
      beta <- D[model]*beta.sc
      beta1 <- numeric(p)
      beta1[model] <- beta
      beta0 <- rnorm(1, mean_y - sum(xbar*beta1), sig/sqrt(n))
      ynew <- beta0 + xnew %*% beta1 + rnorm(nnew, 0, sig)
      return(ynew)
    }

    Ra_list <- lapply(sort(uni.model), FUN=getR_a, x=X, y=y, n=n, lam=lam)
    model.sim1 <- table(model.sim)
    ystar <- list()
    for (i in 1:length(uni.model)) {
      n_models <- model.sim1[i]
      n_vec <- rep(n.new, n_models)
      ynew <- lapply(n_vec, FUN=gen_newy, xnew=Xnew, RSS=RSS, mean_y=mean_y, var_y=var_y,
                     xbar=Xbar, Ra=Ra_list[[i]], model.idx=as.numeric(names(n_models)))
      ystar <- append(ystar, ynew)
    }

    ystar1 <- lapply(ystar, FUN = as.matrix)
    ynew_mat <- matrix(unlist(ystar1), ncol = Nsim)
    ci <- t(apply(ynew_mat, 1, FUN = quantile, probs=c(alpha/2, 1-alpha/2)))
    colnames(ci) <- c(paste0("lower ", conf.level*100, "%"), paste0("upper ", conf.level*100, "%"))
    if(return.draws) {
      return(list(ci=ci, mc.draws=ynew_mat))
    } else {
      return(ci)
    }
  } else if(interval == "Z") {
    E = var_y * RSS / (n - 3)
    V = matrix(0, nrow = n.new, ncol = length(weights))
    U = matrix(0, nrow = n.new, ncol = length(weights))
    for (i in 1:length(weights)) {
      model <- Models[, i]
      Xnew_s <- scale(Xnew[, model, drop=F], center = Xbar[model], scale = 1/D[model])
      X_s <- scale(X[, model, drop=F], center = Xbar[model], scale = 1/D[model])
      E_m <- RSS[i] * var_y/(n-3)
      A <- crossprod(X_s) + lam * diag(sum(model))
      R <- chol(A)
      V[, i] <- as.numeric(Xnew_s %*% solve(A, crossprod(X_s, y)))
      U[, i] <- 1/n + rowSums((Xnew_s %*% backsolve(R, diag(sum(model))))^2)
    }
    var_new = sum(E * weights) + rowSums((V - rowSums(V %*% diag(weights)))^2 %*% diag(weights)) + rowSums(U %*% diag(E * weights))
    mean_new = rowSums((mean_y + V) %*%  diag(weights))
    prec <- qnorm(1-alpha/2) * sqrt(var_new)
    ci <- cbind(mean_new-prec, mean_new+prec)
    colnames(ci) <- c(paste0("lower ", conf.level*100, "%"), paste0("upper ", conf.level*100, "%"))
    return(ci)
  }
}
