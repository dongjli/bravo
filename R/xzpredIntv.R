library(Matrix)

#MC and Z both intervals verified, code is perfect

predIntv.xz <- function(object, model, x1new, x2new, interval, Nsim, return.draws, conf.level, alpha) {
 
  
  x1 <- object$stats$Xm; x2 <- object$stats$Zm
  y <- object$stats$y
  lam1 <- object$stats$lam1;lam2 <- object$stats$lam2

  Models <- object$model.top
  weights <- object$stats$weights
  if (model == "MAP") weights[-1] <- 0
  RSS <- object$stats$RSS.top

  n <- nrow(x1)
  p1 <- ncol(x1);p2 <- ncol(x2); p <- p1+p2
  n.new <- nrow(x1new)
  x1bar <- colMeans(x1); x2bar <- colMeans(x2); xbar <- c(x1bar,x2bar)
  
  
  
  if(class(x1)[1] == "dgCMatrix") {
    D1 = 1/sqrt(colMSD_dgc(x1, x1bar))
  }  else   {
    D1 = apply(x1,2,sd)
    D1 = 1/D1
  }
  
  if(class(x2)[1] == "dgCMatrix") {
    D2 = 1/sqrt(colMSD_dgc(x2, x2bar))
  }  else   {
    D2 = apply(x2,2,sd)
    D2 = 1/D2
  }
  
  
  var_y <- as.vector(var(y))
  mean_y <- mean(y)
  
  if(interval == "MC") {
    model.sim <- sample(1:length(weights), Nsim, prob = weights, replace = T)
    
    
    uni.model <- unique(model.sim)
    getR_a <- function(x1,x2, y, n, model.idx, lam1,lam2)  {
      model <- which(Models[, model.idx])
      model1 <- model[model<=p1] ; model2 <- model[model>p1] 
      model1.size <- length(model1) ; model2.size <- length(model2)
      x1mat <- scale(x1[, model1, drop=F]) ; x2mat <- scale(x2[, model2-p1, drop=F])
      umat <- cbind(x1mat,x2mat)
      B1 <- lam1*diag(rep(1,model1.size)) ; B2 <- lam2*diag(rep(1,model2.size))
      B <- bdiag(B1,B2)
      A <- as.matrix(crossprod(umat) + B)
      R = chol(A)
      a = backsolve(R, crossprod(umat, y), transpose = T)
      return(list(R=R, a=a))
    }

    gen_newy <- function(nnew, x1new,x2new, RSS, mean_y, var_y, x1bar,x2bar, Ra, model.idx)  {
      ### simulate sigma^2
      sig2 <- 1/rgamma(1, shape=(n-1)/2, rate=var_y*RSS[model.idx]/2)
      ### simulate regression coefficients
      model <- which(Models[, model.idx])
      model1 <- model[model<=p1] ;model2 <- model[model>p1]
      model.size <- length(model)
      R <- Ra$R
      a <- Ra$a
      sig <- sqrt(sig2)
      beta.sc <- backsolve(R, a+rnorm(model.size, 0, sig))
      beta <- c(D1[model1],D2[model2-p1])*beta.sc
      beta1 <- numeric(p)
      beta1[model] <-beta
      beta0 <- rnorm(1, mean_y - sum(xbar*beta1), sig/sqrt(n))
      unew <- cbind(x1new,x2new)
      ynew <- beta0 + unew %*% beta1 + rnorm(nnew, 0, sig)
      return(ynew)
    }

    Ra_list <- lapply(sort(uni.model), FUN=getR_a, x1=x1,x2=x2, y=y, n=n, lam1=lam1,lam2=lam2)
    model.sim1 <- table(model.sim)
    ystar <- list()
    for (i in 1:length(uni.model)) {
      n_models <- model.sim1[i]
      n_vec <- rep(n.new, n_models)
      ynew <- lapply(n_vec, FUN=gen_newy, x1new=x1new, x2=x2new, RSS=RSS, mean_y=mean_y, var_y=var_y,
                     x1bar=x1bar,x2bar=x2bar, Ra=Ra_list[[i]], model.idx=as.numeric(names(n_models))) 
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
      model <- which(Models[, i])
      model1 <- model[model<=p1] ;model2 <- model[model>p1] 
      x1new_s <- scale(x1new[, model1, drop=F], center = x1bar[model1], scale = 1/D1[model1])
      x2new_s <- scale(x2new[, model2-p1, drop=F], center = x2bar[model2-p1], scale = 1/D2[model2-p1])
      x1_s <- scale(x1[, model1, drop=F], center = x1bar[model1], scale = 1/D1[model1])
      x2_s <- scale(x2[, model2-p1, drop=F], center = x2bar[model2-p1], scale = 1/D2[model2-p1])
      u_s = cbind(x1_s,x2_s) ; unew_s = cbind(x1new_s,x2new_s)
      E_m <- RSS[i] * var_y/(n-3)
      mult = c(rep(lam1,length(model1)),rep(lam2,length(model2)))
      A <- crossprod(u_s) + mult * diag(length(model))
      R <- chol(A)
      V[, i] <- as.numeric(unew_s %*% solve(A, crossprod(u_s, y)))
      U[, i] <- 1/n + rowSums((unew_s %*% backsolve(R, diag(length(model))))^2) 
    }
    var_new = sum(E * weights) + rowSums((V - rowSums(V %*% diag(weights)))^2 %*% diag(weights)) + rowSums(U %*% diag(E * weights))
    mean_new = rowSums((mean_y + V) %*%  diag(weights))
    prec <- qnorm(1-alpha/2) * sqrt(var_new)
    ci <- cbind(mean_new-prec, mean_new+prec)
    colnames(ci) <- c(paste0("lower ", conf.level*100, "%"), paste0("upper ", conf.level*100, "%"))
    return(ci)
  }
}
