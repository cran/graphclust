#' Simulate a network of a stochastic block model
#'
#' @param n number of vertices
#' @param theta stochastic block model parameter with latent group probabilities $pi and
#' connectivy parameters $gamma
#' @param directed directed network (TRUE by default) or undirected (FALSE)
#'
#' @return list with simulated adjacency matrix ($adj) and node labels ($Z)
#' @export
#'
#' @examples
#' theta1 <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' rsbm(10, theta1)
rsbm <- function(n, theta, directed=TRUE){
  Q <- length(theta$pi)
  Z <- sample(1:Q, n, replace=TRUE, prob=theta$pi)
  # adjacency matrix
  A <- matrix(0, n, n)
  for (i in 1:n){
    A[i, -i] <- rbinom(n-1, 1, theta$gamma[Z[i],Z[(1:n)[-i]]])
  }
  return(list(adj=A, Z=Z))
}


#' Simulate a sample of networks of a stochastic block model
#'
#' @param vec_n vector with number of vertices
#' @param theta stochastic block model parameter with latent group probabilities $pi and
#' connectivy parameters $gamma
#' @param directed directed networks (TRUE by default) or undirected (FALSE)
#'
#' @return list with a list of adjacency matrices ($listGraphs) and a list of
#' node labels ($listLatentZ)
#' @export
#'
#' @examples
#' theta1 <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' rCollectSBM(2:4, theta1)
rCollectSBM <- function(vec_n, theta, directed=TRUE){
  listRes <-lapply(vec_n, function(n) rsbm(n, theta, directed))
  listGraphs <- lapply(listRes, function(el) el$adj)
  listLatentZ <- lapply(listRes, function(el) el$Z)
  return(list(listGraphs=listGraphs, listLatentZ=listLatentZ))
}

# list_theta : includes prop = cluster proportions
#' Simulate a collection of networks of a mixture of stochastic block models
#'
#' @param vec_n vector with number of vertices
#' @param thetaMixSBM K-list for a mixture with K components. Each field
#' is a list with the stochastic block model parameter ($pi and $gamma) and a cluster
#' proportion ($prop)
#' @param directed directed networks (TRUE by default) or undirected (FALSE)
#'
#' @return list with a list of adjacency matrices ($listGraphs), a list of
#' node labels ($listLatentZ) and a vector with the graph clustering ($label)
#' @export
#' @importFrom stats rbinom
#' @examples
#' theta1 <- list(prop=.2, pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' theta2 <- list(prop=.8, pi=c(.5,.5), gamma=matrix(4:1/8,2,2))
#' thetaMixSBM <- list(NULL)
#' thetaMixSBM[[1]] <- theta1
#' thetaMixSBM[[2]] <- theta2
#' obs <- rMixSBM(vec_n=rep(10,3), thetaMixSBM)
rMixSBM <- function(vec_n, thetaMixSBM, directed=TRUE){
  M <- length(vec_n)
  G <- length(thetaMixSBM)
  prop <- sapply(1:G, function(g) thetaMixSBM[[g]]$prop)
  graphGroups <- sample(1:G, M, replace=TRUE, prob=prop)
  listNetworks <- list(label = graphGroups,
                       listGraphs = vector('list', M),
                       listLatentZ = vector('list', M))
  for (g in 1:G){
    groupSize <- sum(graphGroups==g)
    if (groupSize>0){
      networkGroup <- rCollectSBM(vec_n[graphGroups==g], thetaMixSBM[[g]])
      networkIndices <- (1:M)[graphGroups==g]
      listNetworks$listGraphs[networkIndices] <- networkGroup$listGraphs
      listNetworks$listLatentZ[networkIndices] <- networkGroup$listLatentZ
    }
  }
  return(listNetworks)
}
