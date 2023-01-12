# 30/11/2022
# graphon based clustering


#' Graph clustering using the pairwise graphon distances and spectral clustering
#'
#' @param allAdj list of adjacency matrices
#' @param nbClusters number of clusters to be found
#' @param nbCores number of cores for parallelization. Default: detectCores().
#'
#' @return list with the obtained graph clusteirng ($clust) and the matrix with the
#' pairwise graphon distances between all pairs of networks
#' @export
#' @importFrom sClust spectralPAM
#' @examples
#' theta <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' obs <- rCollectSBM(rep(10,4), theta)$listGraphs
#' res <- graphonSpectralClustering(obs, 2, nbCores=1)
graphonSpectralClustering <- function(allAdj, nbClusters, nbCores=detectCores()){
  M <- length(allAdj)
  message("initialization of individual SBMs...\n")
  listTheta <- fitSimpleSBM(allAdj, nbCores = nbCores, outCountStat = FALSE)

  message("pairwise graphon distances...\n")
  graphonDist <- matrix(0, M, M)
  for (k in 2:M){
    for (l in 1:(k-1)){
      graphonDist[k, l] <- graphonDist[l, k] <- sbmNorm(listTheta[[k]], listTheta[[l]])
    }
  }
  clustering <- spectralPAM(graphonDist, K = nbClusters, flagDiagZero = TRUE)$cluster

  return(list(clust=clustering, graphonDist=graphonDist))
}





