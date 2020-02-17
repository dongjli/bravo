
#------------------------------- Compound symmetry ----------------------------------#
dtgen.compSym <- function(n, p, rho) {
  x <- matrix(rnorm(n*p), n, p)
  xmat <- sqrt(1-rho)*x + sqrt(rho)*rnorm(n)
  return(xmat)
}

#--------------------------------------- AR -----------------------------------------#
dtgen.ar <- function(n, p, rho) {
  x0 <- rnorm(n)
  x2 = matrix(0, n, p)
  sig2 = 1-rho^2
  x2[, 1] = rho*x0 + rnorm(n, 0, sqrt(sig2))
  for (i in 2:p) {
    z = rnorm(n, 0, sqrt(sig2))
    x2[, i] = rho*x2[, i-1]+z
  }
  return(x2)
}


#--------------------------------- Factor modeles -----------------------------------#
dtgen.factor <- function(n, p, factor, k) {
  phi <- matrix(rnorm(n*k), n, k)
  x <- phi %*% factor + matrix(rnorm(n*p),n,p)
  return(x)
}
dtgen.factor <- cmpfun(dtgen.factor)

#---------------------------------- Group Structure ---------------------------------#
# m = number of groups
dtgen.group <- function(n, p, p.true, m, delta) {
  z0 <- matrix(rnorm(n*(p-p.true)), n, p-p.true)
  zz = matrix(rnorm(m*n), n, m)
  z = cbind(matrix(rep(zz, p.true/m), nrow = n), z0)
  e = cbind(matrix(rnorm(n*p.true, mean = 0, sd = delta), nrow = n), matrix(0, n, p-p.true))
  x4 = z + e
  return(x4)
}

#-------------------------------- Extreme Correlation -------------------------------#
dtgen.extrmCor <- function(n, p, p.true) {
  w = matrix(rnorm(n*p), n, p.true)
  z = matrix(rnorm(n*p), n, p)
  x5 = matrix(0, n, p)
  x5[, 1:p.true] = (z[, 1:p.true] + w) / sqrt(2)
  x5[, (p.true+1):p] = (z[, (p.true+1):p] + rowSums(w)) / 2
  return(x5)
} 

#---------------------------------- Extreme Group Structure ---------------------------------#
# m = number of groups
dtgen.extrmgroup <- function(n, p, p.true, m, delta) {
  z0 <- matrix(rnorm(n*(p-p.true)), n, p-p.true)
  zz = matrix(rnorm(m*n), n, m)
  z = cbind(matrix(rep(zz, p.true/m), nrow = n), z0)
  e = cbind(matrix(rnorm(n*p.true, mean = 0, sd = delta), nrow = n), matrix(0, n, p-p.true))
  x6 = z + e
  #w = matrix(rnorm(n*p.true), n, p.true/m)
  w = matrix(rnorm(n*p.true), n, p.true)
  #w = matrix(rep(w, p.true/m), nrow = n)
  x6[, 1:p.true] = (z[, 1:p.true] + w) / sqrt(2+delta^2)
  x6[, (p.true+1):p] = (z[, (p.true+1):p] + rowSums(w)) / sqrt(1+p.true)
  return(x6)
}

#---------------------------------- Extreme Group Structure 2 ---------------------------------#
# m = number of groups
dtgen.extrmgroup2 <- function(n, p, p.true, m, delta) {
  z0 <- matrix(rnorm(n*(p-p.true)), n, p-p.true)
  zz = matrix(rnorm(m*n), n, m)
  z = cbind(matrix(rep(zz, p.true/m), nrow = n), z0)
  e = cbind(matrix(rnorm(n*p.true, mean = 0, sd = delta), nrow = n), matrix(0, n, p-p.true))
  x6 = z + e
  #w = matrix(rnorm(n*p.true), n, p.true/m)
  w = rowSums(zz)
  #w = matrix(rep(w, p.true/m), nrow = n)
  #x6[, 1:p.true] = (z[, 1:p.true] + w) / sqrt(2+delta^2)
  x6[, 1:p.true]  = x6[, 1:p.true]/ (1+delta^2)
  x6[, (p.true+1):p] = (z[, (p.true+1):p] + (w)) / sqrt(1+m)
  return(x6)
}


