#' Sort stochastic block model parameter in a unique way using its graphon
#'
#' @param thetaInit stochastic block model parameter to be sorted
#' @param outTheta if TRUE returns the sorted stochastic block model parameter
#' @param outPerm if TRUE returns the permutation of the blocks of the
#' stochastic block model to provide the sorted stochastic block model
#' parameter
#'
#' @return according to the values of outTheta and outPerm the function returns
#' the sorted stochastic block model parameter or the associated permutation of
#' the blocks of the stochastic block model or a list with both of them
#' @export
#'
#' @examples
#' theta1 <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' degreeSort(theta1)
#' theta2 <- list(pi=c(.5,.5), gamma=matrix(4:1/8,2,2))
#' degreeSort(theta2)
degreeSort<- function(thetaInit, outTheta=TRUE, outPerm=FALSE){
  K <- length(thetaInit$pi)
  graphonMat <- thetaInit$gamma * matrix(thetaInit$pi, K, K)
  margin1 <- colSums(graphonMat)
  permut <- order(margin1, decreasing = TRUE)
  theta <- permutParam(thetaInit, permut)
  if (length(unique(margin1)) < K){
    margin1 <- margin1[permut]
    ind <- (1:(K-1))[margin1[1:(K-1)] == margin1[2:K]]
    graphonMat <- theta$gamma * matrix(theta$pi, K, K)
    margin2 <- rowSums(graphonMat)
    for (i in ind){
      if (margin2[i]<margin2[i+1]){
        permut[c(i, i+1)] <- permut[c(i+1, i)]
        theta <- permutParam(thetaInit, permut)
      }
    }
  }

  output <- theta
  if (outPerm){
    if(!outTheta)
      output <- permut
    else
      output <- list(theta=theta, permut=permut)
  }

  return(output)
}


#' (squared) L2-norm of the graphons associated with two stochastic block model parameters
#'
#' @param theta1 a stochastic block model parameter
#' @param theta2 a stochastic block model parameter
#'
#' @return (squared) L2-norm of the graphons associated with two stochastic block model parameters
#' @export
#'
#' @examples
#' theta1 <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' theta2 <- list(pi=c(.5,.5), gamma=matrix(4:1/8,2,2))
#' graphonL2norm(theta1, theta2)
graphonL2norm <- function(theta1, theta2){
  if (sum(theta1$pi==0)>0){
    ind <- theta1$pi==0
    theta1$pi <- theta1$pi[!ind]
    theta1$gamma <- theta1$gamma[(!ind)%*%(t(!ind))==1]
  }
  if (sum(theta2$pi==0)>0){
    ind <- theta2$pi==0
    theta2$pi <- theta2$pi[!ind]
    theta2$gamma <- theta2$gamma[(!ind)%*%(t(!ind))==1]
  }
  K1 <- length(theta1$pi)
  K2 <- length(theta2$pi)

  s <- c(cumsum(theta1$pi), cumsum(theta2$pi))
  s.order <- order(s)
  Ks <- length(s) - 1

  gInd <- matrix(1, Ks, 2)

  if (Ks>1){
    for (j in 2:Ks){
      if (s.order[j-1]<=K1){
        gInd[j, 1] <- gInd[j-1, 1] + 1
        gInd[j, 2] <- gInd[j-1, 2]
      } else {
        gInd[j, 1] <- gInd[j-1, 1]
        gInd[j, 2] <- gInd[j-1, 2] + 1
      }
    }
  }
  s <- sort(s[1:Ks])
  length.rect.s <- if (Ks>1) s - c(0,s[1:(Ks-1)]) else 1
  area.s <- length.rect.s %*% t(length.rect.s)

  indGamma1 <- matrix(c(rep(gInd[,1], Ks), rep(gInd[,1], each=Ks)), ncol=2)
  indGamma2 <- matrix(c(rep(gInd[,2], Ks), rep(gInd[,2], each=Ks)), ncol=2)

  L2norm <- sum((matrix(theta1$gamma[indGamma1], Ks, Ks) - matrix(theta2$gamma[indGamma2], Ks, Ks))^2 * area.s)

  return(L2norm)
}



#' Permute block labels of a stochastic block model parameter
#'
#' @param theta a SBM parameter with say K blocks
#' @param permut a permutation of the block labels 1,2,...,K
#'
#' @return stochastic block model parameter with permuted block labels
#' @export
#'
#' @examples
#' theta1 <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' theta2 <- list(pi=c(.5,.5), gamma=matrix(4:1/8,2,2))
#' permutParam(theta1, 2:1)
#' permutParam(theta2, 2:1)
permutParam <- function(theta, permut){
  theta$pi <- theta$pi[permut]
  K <- length(permut)
  ind <- matrix(c(rep(permut, K), rep(permut, each=K)), ncol=2)
  theta$gamma <- matrix(theta$gamma[ind], K, K)
  return(theta)
}


#' (squared) norm between two stochastic block models
#'
#' the norm is the minimal graphon distance between two stochastic block model
#' parameters obtained with the best permutations of the parameters
#'
#' @param theta1 a stochastic block model parameter
#' @param theta2 a stochastic block model parameter
#'
#' @return (squared) norm between two stochastic block models
#' @export
#'
#' @examples
#' theta1 <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' theta2 <- list(pi=c(.5,.5), gamma=matrix(4:1/8,2,2))
#' theta3 <- list(pi=c(.5,.5), gamma=matrix(1:4/4,2,2))
#' sbmNorm(theta1, theta2)
#' sbmNorm(theta1, theta3)
#' sbmNorm(theta2, theta3)
sbmNorm <- function(theta1, theta2){
  return(graphonL2norm(degreeSort(theta1), degreeSort(theta2)))
}
