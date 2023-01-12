
#' Auxiliary function for ARI
#'
#' @noRd
#' @param x vector
#' @param y vector
#'
#' @return ARI
ARI2vectors <- function(x, y) {
  if (is.matrix(x))
    x <- apply(x, 2, which.max)
  if (is.matrix(y))
    y <- apply(y, 2, which.max)

  # first, get crosstabs
  ctab <- table(x,y)
  # now calculate 4 intermediary sums
  cellsum <- sum(ctab*(ctab-1)/2)
  totsum <- sum(ctab)*(sum(ctab)-1)/2
  # use matrix multiplication to get row and column marginal sums
  rows <- ctab %*% rep(1, ncol(ctab))
  rowsum <- sum(rows*(rows-1)/2)
  cols <- rep(1, nrow(ctab)) %*% ctab
  colsum <- sum(cols*(cols-1)/2)
  # now put them together
  denom <- (rowsum +colsum)/2 - rowsum*colsum/totsum
  if (denom!=0){
    adj.rand <- (cellsum - (rowsum*colsum/totsum))/denom
  }else{
    d <- x[1] - y[1]
    adj.rand <- if (sum(x==(y+d))==length(x)) 1 else NA
  }
  return(adj.rand)
}


#' Adjusted Rand index
#'
#' ARI to compare two clusterings or to compare two entire lists of clusterings
#'
#' @param x vector with clustering, matrix with hot-one-encoding of the
#' clustering, or a list of clusterings (in vector or matrix form)
#' @param y as x
#'
#' @return ARI (scalar of vector)
#' @export
#'
#' @examples
#' x <- c(1,1,2,2,3,3)
#' y <- c(1,1,1,2,2,2)
#' ARI(x,y)
#'
#' x <- matrix(0, 3, 6)
#' x[1,1] <- x[1,2] <- x[2,3] <- x[2,4] <- x[3,5] <- x[3,6] <- 1
#' y <- matrix(0, 2, 6)
#' y[1,1] <- y[1,2] <- y[1,3] <- y[2,4] <- y[2,5] <- y[2,6] <- 1
#' ARI(x,y)
#'
#' X <- list(c(1,1,2,2,3,3), rep(1,10))
#' Y <- list(c(1,1,1,2,2,2), rep(1:2,each=5))
#' ARI(X,Y)
ARI <- function(x, y) {
  if ( (is.matrix(x)|is.vector(x, mode='numeric'))&(is.matrix(y)|is.vector(y, mode='numeric')))
    adj.rand <- ARI2vectors(x,y)

  if (is.list(x)){
    L <- length(x)
    adj.rand <- sapply(1:L, function(l) ARI2vectors(x[[l]], y[[l]]))
  }

  return(adj.rand)
}
