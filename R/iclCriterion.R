

#' Evaluation of ICL criterion for a sample of networks of a common SBM
#' @noRd
#' @param countStat list of count statistics
#' @param hyperParam hyperparameters of prior distributions
#'
#' @return value of the ICL criterion
iclCriterionSBMiid <- function(countStat, hyperParam){

  if(is.list(countStat$S)){
    M <- length(countStat$S)
    K <- length(countStat$S[[1]])
    sum_skm <- colSums(matrix(unlist(countStat$S), M, K, byrow=TRUE))
    sum_rklm <- rowSums(matrix(unlist(countStat$R), nrow=K^2))
    sum_aklm <- rowSums(matrix(unlist(countStat$A), nrow=K^2))

  }else{
    K <- ncol(countStat$S)
    M <- nrow(countStat$S)
    sum_skm <- rowSums(countStat$S)
    sum_rklm <- c(apply(countStat$R, 1:2, sum))
    sum_aklm <- c(apply(countStat$A, 1:2, sum))
  }

  sum_nm <- sum(sum_skm)

  icl <- -K^2*lbeta(hyperParam$eta, hyperParam$zeta)
  icl <- icl + lgamma(K*hyperParam$alpha) - K*lgamma(hyperParam$alpha)
  icl <- icl - lgamma(K*hyperParam$alpha + sum_nm)
  icl <- icl + sum(lgamma(hyperParam$alpha + sum_skm))
  icl <- icl + sum(lbeta(hyperParam$eta + sum_aklm, hyperParam$zeta + sum_rklm - sum_aklm))

  return(icl)
}


#' Evaluation of ICL criterion for a collection of networks and a given graph clustering
#'
#' @param countStat list of count statistics
#' @param graphGroups vector with graph clustering
#' @param hyperParam hyperparameters of prior distributions
#'
#' @return value of the ICL criterion
#' @noRd
iclCriterionSBMmix <- function(countStat, graphGroups, hyperParam){
  iclSBMgroups <- sapply(unique(graphGroups), function(l){
    iclCriterionSBMiid(extractCountStat(countStat, ind = (graphGroups==l)), hyperParam)})

  L <- length(unique(graphGroups))
  M <- length(graphGroups)
  m_c <- sapply(1:L, function(c) sum(graphGroups==c))
  penalty <- lgamma(L*hyperParam$lambda) - L*lgamma(hyperParam$lambda) -
    lgamma(L*hyperParam$lambda + M) + sum(lgamma(m_c + hyperParam$lambda))

  icl <- sum(iclSBMgroups) + penalty
  return(icl)
}




#' Compute the variations of the ICL criterion by all possible pairwise cluster
#' aggregations
#'
#' @param allAdj list of adjacency matrices
#' @param countStat list of count statistics
#' @param hyperParam hyperparameters of prior distributions
#'
#' @return matrix with variations of the ICL criterion for all possible pairwise cluster
#' aggregations
#' @noRd
getAllDeltaICL <- function(allAdj, countStat, hyperParam){
    M <- length(allAdj)
    # constPenalty <- globalPenalty(hyperParam$lambda), M, M)
    # m1 <- m2 <- 1
    diffU_partial <-  lgamma(hyperParam$lambda + 2) -
      lgamma(hyperParam$lambda + 1) - lgamma(hyperParam$lambda - 1)
    iclPerNetwork <- sapply(1:M, function(m) iclCriterionSBMiid(extractCountStat(countStat, m), hyperParam))
    matDeltaICL <- matrix(NA, M, M)
    colnames(matDeltaICL) <- rownames(matDeltaICL) <- 1:M
    for (ind1 in 1:(M-1)){
      for (ind2 in (ind1+1):M){
        resMerge <- merge2graphGroups(
          allAdj[ind1], extractCountStat(countStat, ind1),
          allAdj[ind2], extractCountStat(countStat, ind2),
          hyperParam)
        iclNew <- iclCriterionSBMiid(mergeCountStat(resMerge$countStat1, resMerge$countStat2), hyperParam)
        iclClust1 <- iclPerNetwork[ind1] # iclCriterionSBMiid(extractCountStat(countStat, ind1), hyperParam)
        iclClust2 <- iclPerNetwork[ind2] # iclCriterionSBMiid(extractCountStat(countStat, ind2), hyperParam)
        matDeltaICL[ind1, ind2] <- iclNew - iclClust1 - iclClust2 + diffU_partial
      }
    }
    return(matDeltaICL)
  }




#' Evaluation of the variation of the ICL criterion by aggregating two given
#' clusters
#'
#' @param resMerge result of merge2graphGroups() for the merge of clusters ind1 and ind2
#' @param countStat list of count statistics
#' @param ind1 cluster label of first cluster to be aggregated
#' @param ind2 cluster label of second cluster to be aggregated
#' @param hyperParam hyperparameters of prior distributions
#'
#' @return value of the variation of the ICL criterion by aggregating the two
#' clusters
#' @noRd
iclIncreaseCrit <- function(resMerge, countStat, ind1, ind2, hyperParam){
  iclNew <- iclCriterionSBMiid(mergeCountStat(resMerge$countStat1, resMerge$countStat2), hyperParam)
  iclClust1 <- iclCriterionSBMiid(extractCountStat(countStat, ind1), hyperParam)
  iclClust2 <- iclCriterionSBMiid(extractCountStat(countStat, ind2), hyperParam)
  m1 <- sum(ind1)
  m2 <- sum(ind2)
  diffU_partial <-  lgamma(hyperParam$lambda + m1 + m2) -
    lgamma(hyperParam$lambda + m1) - lgamma(hyperParam$lambda - m2)
  delta <- iclNew - iclClust1 - iclClust2 + diffU_partial
  return(delta)
}

#' Evaluation of the global penalty term in the ICL criterion
#'
#' @param lambda hyperparameter of Dirichlet prior
#' @param M number of networks in the collection
#' @param L number of current clusters
#'
#' @return global penalty term in the ICL criterion
#' @noRd
globalPenalty <- function(lambda, M, L){
  pen <- lgamma((L-1)*lambda) + lgamma(lambda) +
    lgamma(L*lambda + M) - lgamma((L-1)*lambda + M)-
    lgamma(L*lambda)
  return(pen)
}


# brauchen wir
#' Update of the matrix with variations of the ICL criterion of all possible
#' cluster aggregations after the aggregation of two clusters
#'
#' @param deltaICLmatrix old matrix with variations of the ICL criterion of all possible
#' cluster aggregations
#' @param allAdj list of adjacency matrices
#' @param countStat list of count statistics after aggregation
#' @param ind_ml clusters that have just been merged
#' @param graphGroups vector with graph clustering after aggregation
#' @param hyperParam hyperparameters of prior distributions
#'
#' @return updated matrix with variations of the ICL criterion of all possible
#' cluster aggregations
#' @noRd
updateDeltaICLmatrix <- function(deltaICLmatrix, allAdj, countStat, ind_ml,
                                 graphGroups, hyperParam){
  ind_m <- min(ind_ml)
  ind_l <- max(ind_ml)
  deltaICLmatrix <- deltaICLmatrix[-ind_l, ]
  deltaICLmatrix <- deltaICLmatrix[, -ind_l]

  # update deltaICL involving ind_m
  ind_row <- if (ind_m>1) (1:(ind_m-1)) else NULL
  L <- nrow(deltaICLmatrix)
  ind_col <- if (ind_m<L) ((ind_m+1):L) else NULL

  ind1 <- (graphGroups == rownames(deltaICLmatrix)[ind_m])
  countStat1 <- extractCountStat(countStat, ind1)
  for (k in ind_row){
    ind2 <- (graphGroups == rownames(deltaICLmatrix)[k])
    resMerge <- merge2graphGroups(
      allAdj[ind1], countStat1,
      allAdj[ind2], extractCountStat(countStat, ind2),
      hyperParam)

    deltaICLmatrix[k, ind_m] <- iclIncreaseCrit(resMerge, countStat, ind1, ind2, hyperParam)
  }
  for (k in ind_col){
    ind2 <- (graphGroups == rownames(deltaICLmatrix)[k])
    resMerge <- merge2graphGroups(
      allAdj[ind1], countStat1,
      allAdj[ind2], extractCountStat(countStat, ind2),
      hyperParam)

    deltaICLmatrix[ind_m, k] <- iclIncreaseCrit(resMerge, countStat, ind1, ind2, hyperParam)
  }
  return(deltaICLmatrix)
}




