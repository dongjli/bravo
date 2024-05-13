#' Compute marginal inclusion probabilities from a fitted "sven" object.
#' @rdname mip.sven
#' @description This function computes the marginal inclusion probabilities of all variables from
#' a fitted "sven" object.
#' @param object A fitted "sven" object
#' @param threshold marginal inclusion probabilities above this threshold are stored. Default 0.
#' 
#' @return The object returned is a data frame if the \code{sven} was run with a single matrix, 
#' or a list of two data frames if \code{sven} was run with a list of two matrices.
#' The first column are the variable names (or numbers if column names of were absent). 
#' Only the nonzero marginal inclusion probabilities are stored.
#' 
#' 
#' @author Somak Dutta\cr Maintainer:
#' Somak Dutta <somakd@@iastate.edu>
#' @examples
#' \donttest{
#' n <- 50; p <- 100; nonzero <- 3
#' trueidx <- 1:3
#' truebeta <- c(4,5,6)
#' X <- matrix(rnorm(n*p), n, p) # n x p covariate matrix
#' y <- 0.5 + X[,trueidx] %*% truebeta + rnorm(n)
#' res <- sven(X=X, y=y)
#' res$model.map # the MAP model
#' mip.sven(res)
#' 
#' Z <- matrix(rnorm(n*p), n, p) # another covariate matrix
#' y2 = 0.5 + X[,trueidx] %*% truebeta  + Z[,1:2] %*% c(-2,-2) + rnorm(n)
#' res2 <- sven(X=list(X,Z), y=y2)
#' mip.sven(res2) # two data frames, one for X and another for Z
#' }
##' @export
##' @export mip.sven
mip.sven <- function(object, threshold = 0) {
  pp <- exp(object$pprob.top - max(object$pprob.top))
  pp <- pp/sum(pp)
  
  mips <- as.numeric(object$model.top %*% pp)
  
  if(!object$stats$islist) {
    id <- mips > threshold
    if(!is.null(object$stats$colnamesX)) {
      Hits <- object$stats$colnamesX[id]
    } else {
      Hits <- which(id)
    }
    return(data.frame(Hits=Hits,MIP = mips[id]))
  }
  
  # else
  p1 <- ncol(object$stats$Xm)
  p2 <- ncol(object$stats$Zm)
  id <- which(mips > threshold)
  idX <- id[id <= p1]
  idZ <- id[id > p1] - p1
  
  if(!is.null(object$stats$colnamesX)) {
    Hits1 <- object$stats$colnamesX[idX]
  } else {
    Hits1 <- idX
  }
  
  if(!is.null(object$stats$colnamesZ)) {
    Hits2 <- object$stats$colnamesZ[idZ]
  } else {
    Hits2 <- idZ
  }
  
  return(list(HitsX = data.frame(Hits = Hits1, MIP = mips[idX]),
              HitsZ = data.frame(Hits = Hits2, MIP = mips[idZ+p1])
             ))
}

