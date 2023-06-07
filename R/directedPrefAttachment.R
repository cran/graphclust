 #

#
#' generation of a network of the directed preferential attachment (DPA) model
#'
#' @param param  list with the following elements: $R (= number of iterations),
#' $alpha, $beta , $deltaIn, $deltaOut (parameters of the DPA model)
#'
#' @return adjacency matrix of generated network
#' @export
#'
#' @examples
#' param <- list(R = 500, alpha = .04, beta = .02, deltaIn = 100, deltaOut = 100)
#' A <- sampleDPA(param)
#' A
sampleDPA <- function(param){
  R <- param$R
  alpha <- param$alpha
  beta <- param$beta
  deltaIn <- param$deltaIn
  deltaOut <- param$deltaOut
  n <- 1
  A <- matrix(0, n, n) # current adjacency matrix
  if ((alpha<=0) | (beta<=0)){
    return(warning("alpha and beta must be positive"))
  }
  if ((alpha+beta) > 1){
    s <- alpha + beta
    alpha <- alpha/s
    beta <- beta/s
  }
  for (k in 1:R){
    u <- sample(1:3, 1, prob = c(alpha, beta, 1-alpha-beta))
    if (u==3){ # a new directed edge between existing nodes, but no new edge
      outDegrees <- rowSums(A)
      inDegrees <- colSums(A)
      i <- sample(1:n, 1, prob = outDegrees + deltaOut)
      j <- sample(1:n, 1, prob = inDegrees + deltaIn)
      if (i!=j)
        A[i, j] <- 1
    }else{
      A <- rbind(A, rep(0, n))
      n <- n + 1
      A <- cbind(A, rep(0, n))
      if (u==1){ # add a new node with an outgoing edge
        inDegrees <- colSums(A)[1:(n-1)]
        j <- sample(1:(n-1), 1, prob = inDegrees + deltaIn)
        A[n, j] <- 1
      }else{ # add a new node with an outgoing edge
        outDegrees <- rowSums(A)[1:(n-1)]
        i <- sample(1:(n-1), 1, prob = outDegrees + deltaOut)
        A[i, n] <- 1
      }
    }
  }
  return(A)
}


#' Generation of a mixture of directed preferential attachment (DPA) models
#'
#' @param M number of desired networks
#' @param param list of list of parameters of the DPA models. Each element of
#' param is a list with the following elements: $prop (weight of the mixture
#' component), $R (= number of iterations),
#' $alpha, $beta , $deltaIn, $deltaOut (parameters of the DPA model)
#'
#' @return list of 2 lists : the first ($listAdj) is a list of M adjacency matrices, the
#' second a list ($graphGroups) contains the true cluster labels
#' @export
#'
#' @examples
#' param <- vector('list', 3)
#' param[[1]] <- list(prop = 1/3, # component 1 :  alpha > beta
#'                    alpha = .04,
#'                    beta = .02,
#'                    deltaIn = 100,
#'                    deltaOut = 100,
#'                    R = 500
#' )
#' param[[2]] <- list(prop = 1/3, # component 2 : just permute alpha and beta ;
#'                    alpha = .01,
#'                    beta = .02,
#'                    deltaIn = 100,
#'                    deltaOut = .1,
#'                    R = 1000
#' )
#' param[[3]] <- list(prop = 1/3, # component 3 : alpha=beta
#'                    alpha = .015,
#'                    beta = .015,
#'                    deltaIn = .1,
#'                    deltaOut = .1,
#'                    R = 1000
#' )
#' obs <- sampleDPAMixture(M=20, param)
sampleDPAMixture <- function(M, param){
  nbComp <- length(param)
  prop <- sapply(param, function(comp) comp$prop)
  graphGroups <- sample(1:nbComp, M, replace=TRUE, prob=prop)
  listAdj <- lapply(1:M, function(m) sampleDPA(param[[graphGroups[m]]]))
  res <- list(graphGroups=graphGroups, listAdj=listAdj)
  return(res)
}

