#' Plot dendrogram to visualize the clustering obtained by the hierarchical
#' clustering algorithm
#'
#' @param res output of graphClustering()
#' @param labels network labels, default (NULL) network number.
#' @param labcex size of labels in the figure
#'
#' @return dendrogram
#' @export
#' @importFrom graphics abline
#' @examples
#' theta <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' obs <- rCollectSBM(rep(10,4), theta)$listGraphs
#' res <- graphClustering(obs, nbCores=2)
#' plotDendrogram(res)
plotDendrogram <- function(res, labels=NULL, labcex=0.5){ # include cluster labels

  resclust <- gethclust(res, labels)
  P <- plot(resclust, xlab=' ', ylab='Delta ICL', sub=' ', labels=labels,
            cex=labcex)
  cst <- 1.2
  abline(h=cst * sum(res$histDeltaICL) , col= "black")
  return(P)
}


#' Auxiliary function for plotDendrogram()
#'
#' @param res output of graphClustering()
#' @param labels network labels, default (NULL) network number.
#'
#' @return hclust object
#' @importFrom stats as.dist hclust
#' @noRd
gethclust <- function(res, labels=NULL){
  cst <- 1.2
  M <- ncol(res$histGraphGroups)
  I <- nrow(res$histGraphGroups) - 1

  cumDeltaICL <- cumsum(res$histDeltaICL)

  D <- matrix(0, M, M)
  for (it in 1:I){
    label_m <- res$histFusedClusters[it, 1]
    label_l <- res$histFusedClusters[it, 2]
    ind_m <- (1:M)[res$histGraphGroups[it, ] == label_m]
    ind_l <- (1:M)[res$histGraphGroups[it, ] == label_l]
    nb_m <- length(ind_m)
    nb_l <- length(ind_l)

    ind <- matrix(c(rep(ind_m, times=nb_l), rep(ind_l, each=nb_m)), ncol=2)
    D[ind] <- cumDeltaICL[it]
  }

  D <- D + t(D)
  D[D==0] <- cst * cumDeltaICL[I]
  diag(D) <- 0
  if(!is.null(labels))
    colnames(D) <-rownames(D) <- labels

  D <- as.dist(D)
  resclust <- hclust(D)
  return(resclust)
}


#' Plot the metagraph of the parameter of the stochastic block model associated
#' with one of the estimated graph clusters
#'
#' @param nb number of the cluster we are interested in
#' @param res output of graphClustering()
#' @param title title of the figure
#' @param edge.width.cst width of edges in the metagraph
#'
#' @return none
#' @export
#' @importFrom igraph graph_from_adjacency_matrix
#' @examples
#' theta <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' obs <- rCollectSBM(rep(10,4), theta)$listGraphs
#' res <- graphClustering(obs, nbCores=2)
#' metagraph(1, res)
metagraph <- function(nb, res, title=NULL, edge.width.cst=10){
  adj <- res$thetaMixSBM[[nb]]$theta$gamma
  G <- graph_from_adjacency_matrix(adj, "directed", weighted=TRUE)
  plot(G, edge.arrow.size=.5,
       edge.width=adj * edge.width.cst,
       vertex.size=res$thetaMixSBM[[nb]]$theta$pi*100,
       vertex.label=round(res$thetaMixSBM[[nb]]$theta$pi, digits=2),
       edge.color='black',
       main=title
       )
  return()
}

