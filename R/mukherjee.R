
#' Computation of graph moments of a network
#'
#' @param A adjacency matrix
#' @param k order of the largest graph moments to be considered
#'
#' @return vector with the first k (normalized) graph moments of the network A
#' @export
#'
#' @examples
#' param <- list(R = 500, alpha = .04, beta = .02, deltaIn = 100, deltaOut = 100)
#' A <- sampleDPA(param)
#' moments(A)
moments <- function(A, k=3){
  # A = adjacency matrix of the desired graph
  # k = no. of moments required
  # method = 'exact' or 'approx': 'exact' returns exact counts and
  # 'approx' returns normalized traces of powers of adjacency matrix
  y <- matrix(NA, k)
  n <- nrow(A)
  for (t in 1:k){
    if (t==1)
      temp <- A
    else
      temp <- temp %*% temp
    y[t] <- sum(diag(temp))/ exp((1+t/2)*log(n)); # this normalization is for dense graphs
  }
  return(y)
}



#' Graph clustering method using graph moments
#'
#' Graph clustering method based on graph moments by Mukherjee et al. (2017)
#'
#' @param Networks list of adjacency matrices
#' @param nbMoments order of the largest graph moments to be considered
#' @param nbClusters desired number of clusters
#'
#' @return vector with the clustering of the networks
#' @export
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
#' res <- graphMomentsClustering(obs$listAdj, 3, 3)
#' table(res, obs$graphGroups)
graphMomentsClustering <- function(Networks, nbMoments=3, nbClusters){
  # Networks: list of adjacency matrices
  M <- length(Networks)
  countstat <- sapply(Networks, function(net) moments(net, nbMoments))

  # Euclidean distance matrix between estimated graphons (Frobenius distance)
  DM_counts <- matrix(0, M, M)
  for (i in 1:(M-1)){
    for (j in (i+1):M){
      DM_counts[i, j] <- DM_counts[j, i] <- sqrt(sum((countstat[ , i] - countstat[ , j])^2))
    }
  }

  # Now we do spectral clustering on distance matrix (ce n'est pas du vrai spectral clustering car pas de Laplacien !!)
  spec <- eigen(DM_counts, symmetric=TRUE)
  eigvecMat_U <- spec$vectors[ , M:(M-nbClusters+1)]

  clustering <- stats::kmeans(eigvecMat_U, centers = nbClusters, nstart = 100)$cluster

  return(clustering)
}
