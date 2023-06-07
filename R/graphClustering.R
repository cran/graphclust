#' Hierarchical graph clustering algorithm
#'
#' Applies the hierarchical graph clustering algorithm to a collection of networks
#' and fits a finite mixture model of stochastic block models to the data
#'
#' @param allAdj list of adjacency matrices
#' @param hyperParam hyperparameters of prior distributions
#' @param nbCores number of cores for parallelization
#' @param returnInitial Boolean. Return SBM parameters from initialization or not. Default is FALSE.
#' @param nbClust desired number of clusters. Default NULL, which means that the
#' number of clusters is chosen automatically via the ICL criterion
#' @param nbSBMBlocks upper bound for the number of blocks in the SBMs of the mixture components. Default is Inf
#' @param initCountStat initial count statistics may be provided to the method. Default is NULL.
#' @param initDeltaICL initial deltaICL-matrix may be provided to the method. Default is NULL.
#'
#' @return list with the following fields: $graphGroups is the graph clustering,
#' $nodeClusterings is a list with the node labels for each networks,
#' $thetaMixSBM contains the estimated parameter of the mixture of SBMs,
#' $ICL is the value of the ICL criterion of the final clustering,
#' $histGraphGroups traces the history of the cluster aggregations,
#' $histDeltaICL traces the evolution of the deltaICL value,
#' $histFusedClusters traces the history of the aggregated cluster numbers
#' @export
#'
#' @examples
#' theta <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' obs <- rCollectSBM(rep(10,4), theta)$listGraphs
#' res <- graphClustering(obs, nbCores=1)
graphClustering <- function(allAdj,
                            hyperParam = list(alpha = .5, eta = .5, zeta = .5, lambda = 0.5),
                            returnInitial = FALSE, nbClust = NULL, nbSBMBlocks = Inf,
                            initCountStat = NULL, initDeltaICL = NULL,
                            nbCores=1){
  directed <- TRUE
  M <- length(allAdj)
  if (M<2){
    return(warning("There are not enough networks. At least two networks are required for clustering."))
  }

  # exploreSplit <- rep(c(FALSE, (nbSBMBlocks==Inf)), M)
  exploreSplit <- rep(c(FALSE, TRUE), M)
  nbExploreSplit <- 1

  # initialiaztion
  message("initialization of individual SBMs...\n")
  if (is.null(initCountStat)){
    countStat <- fitSimpleSBM(allAdj, directed, nbSBMBlocks, nbCores=nbCores)
  }else{
    countStat <- initCountStat
  }
  if (returnInitial){
    initCountStat <- countStat
    initTheta <- estimMAPmixSBM(allAdj, graphGroups = 1:M, initCountStat,
                                hyperParam, directed=TRUE)

  }
  message("initialization of DelatICL...\n")
  if(is.null(initDeltaICL)){
    deltaICLmatrix <- getAllDeltaICL(allAdj, countStat, hyperParam)
  }else{
    deltaICLmatrix <- initDeltaICL
  }
  if (returnInitial){
    initDeltaICLmatrix <- deltaICLmatrix
  }
  message("clustering algorithm...\n")
  graphGroups <- 1:M
  it <- 0
  L <- nrow(deltaICLmatrix)

  saveGraphGroups <- saveFusedClusters <- list(NULL)
  saveDeltaICL <- NULL

  stopCrit <- FALSE
  if (!is.null(nbClust))
    if (nbClust >= M)
      stopCrit <- TRUE

  while (!stopCrit){
    it <- it + 1
    # best merge from deltaICLmatrix (if it still increases ICL)
    maxDeltICL <- max(deltaICLmatrix, na.rm=TRUE)
    penaltyTerm <- globalPenalty(hyperParam$lambda, M, L)
    if( (!is.null(nbClust)) | ((maxDeltICL > -penaltyTerm) & (is.null(nbClust))) ){
      ind_ml <- which(deltaICLmatrix == maxDeltICL, arr.ind=TRUE)
      twoGroups <- range(as.numeric(rownames(deltaICLmatrix)[ind_ml]))
      saveFusedClusters[[it]] <- twoGroups
      # merge
      ind1 <- (graphGroups == twoGroups[1])
      ind2 <- (graphGroups == twoGroups[2])

      resMerge <- merge2graphGroups(
        allAdj[ind1], extractCountStat(countStat, ind1),
        allAdj[ind2], extractCountStat(countStat, ind2),
        hyperParam)

      saveDeltaICL <- c(saveDeltaICL, maxDeltICL)

      # update countStat and graphGroups
      countStat <- updateCountStat(countStat, ind1, resMerge$countStat1)
      countStat <- updateCountStat(countStat, ind2, resMerge$countStat2)
      labelGraphGroup <- min(twoGroups)
      graphGroups[ind1] <- labelGraphGroup
      graphGroups[ind2] <- labelGraphGroup
      saveGraphGroups[[it]] <- graphGroups

      currentNbBlocks <- max(sapply(resMerge$countStat1$S, length))
      if(exploreSplit[it]&(currentNbBlocks<nbSBMBlocks)){
        ind <- rep(1:M, 2)[c(ind1, ind2)]
        for (rep in 1:nbExploreSplit){
          newCountStat <- splitBlock(allAdj[ind],
                                     extractCountStat(countStat, ind),
                                     S = 2, applyICL = TRUE,
                                     hyperParam = hyperParam)

          countStat <- updateCountStat(countStat, ind, newCountStat)
        }
      }

      # update  deltaICLmatrix
      if (L>2){
        deltaICLmatrix <- updateDeltaICLmatrix(deltaICLmatrix, allAdj, countStat, ind_ml, graphGroups, hyperParam)
        L <- nrow(deltaICLmatrix)
        if (is.null(nbClust)){
          # stopCrit <- (sum(deltaICLmatrix>0, na.rm=TRUE)==0)
          penaltyTerm <- globalPenalty(hyperParam$lambda, M, L)
          stopCrit <- (sum(deltaICLmatrix>-penaltyTerm, na.rm=TRUE)==0)
        }else{
          stopCrit <- (length(unique(graphGroups))==nbClust)
        }
      }else{
        if (is.null(nbClust)){
          stopCrit <- TRUE
        }else{
          stopCrit <- (length(unique(graphGroups))==nbClust)
        }
      }

    }
  }

  bestICL <- iclCriterionSBMmix(countStat, graphGroups, hyperParam)
  graphGroups <- relabelGraphGroups(graphGroups)
  thetaMixSBM <- estimMAPmixSBM(allAdj, graphGroups, countStat, hyperParam)

  saveGraphGroups <- matrix(c(1:M, unlist(saveGraphGroups)), ncol=M, byrow=TRUE)
  saveFusedClusters <- matrix(unlist(saveFusedClusters), ncol=2, byrow=TRUE)
  res <- list(graphGroups = graphGroups,
              nodeClusterings = countStat$Z,
              thetaMixSBM = thetaMixSBM,
              ICL = bestICL,
              histGraphGroups = saveGraphGroups,
              histDeltaICL = saveDeltaICL,
              histFusedClusters = saveFusedClusters
  )

  if (returnInitial){
    res$initCountStat <- initCountStat
    res$initTheta <- initTheta
    res$initDeltaICL <- initDeltaICLmatrix
  }

  return(res)
}


#' Fit a unique stochastic block model to a collection of networks
#'
#' fitSBMcollection() is a subversion of graphClustering() where
#' no stopping criterion is applied. So all networks are ultimately merged to
#' a single cluster and considered as i.i.d realisations of a single
#' stochastic block model.


#' @param allAdj list of adjacency matrices
#' @param hyperParam hyperparameters of prior distributions
#' @param nbCores number of cores for parallelization
#'
#' @return list with the following fields:
#' $nodeClusterings is a list with the node labels for each networks,
#' $theta contains the estimated SBM parameter,
#' $ICL is the value of the ICL criterion of the final clustering
#' @export
#'
#' @examples
#' theta <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' obs <- rCollectSBM(rep(10,4), theta)$listGraphs
#' res <- fitSBMcollection(obs, nbCores=1)
fitSBMcollection <- function(allAdj,
                                 hyperParam = list(alpha = .5, eta = .5, zeta = .5, lambda = 0.5),
                                 nbCores=1){
  directed <- TRUE
  M <- length(allAdj)
  exploreSplit <- rep(c(FALSE,TRUE), M)#TRUE# FALSE  # for mixture model: splits increase ICL and diminish the number of graph clusters
  nbExploreSplit <- 1 # can improve results, but not always

  # initialiaztion
  message("initialization...",'\n')
  countStat <- fitSimpleSBM(allAdj, directed, nbCores=nbCores)
  message("start merge algorithm...",'\n')

  graphGroups <- 1:M

  stopCrit <- (M<2)
  it <- 0
  while (!stopCrit){   # construction of an entire tree
    it <- it + 1
    twoGroups <- sample(unique(graphGroups), 2)
    # merge
    ind1 <- (graphGroups == twoGroups[1])
    ind2 <- (graphGroups == twoGroups[2])

    resMerge <- merge2graphGroups(
      allAdj[ind1], extractCountStat(countStat, ind1),  # les deux graphGroups choisis n'ont pas forcement le meme nombre de blocs
      allAdj[ind2], extractCountStat(countStat, ind2),
      hyperParam)

    # update countStat and graphGroups
    countStat <- updateCountStat(countStat, ind1, resMerge$countStat1)
    countStat <- updateCountStat(countStat, ind2, resMerge$countStat2)
    labelGraphGroup <- min(twoGroups)
    graphGroups[ind1] <- labelGraphGroup
    graphGroups[ind2] <- labelGraphGroup

    if(exploreSplit[it]){
      ind <- rep(1:M, 2)[c(ind1, ind2)]
      for (rep in 1:nbExploreSplit){
        newCountStat <- splitBlock(allAdj[ind],
                                   extractCountStat(countStat, ind),
                                   S = 2, applyICL = TRUE,
                                   hyperParam = hyperParam)

        countStat <- updateCountStat(countStat, ind, newCountStat)
      }
    }

    stopCrit <- length(unique(graphGroups))==1
  }

  bestICL <- iclCriterionSBMmix(countStat, graphGroups, hyperParam)
  graphGroups <- relabelGraphGroups(graphGroups)
  thetaMixSBM <- estimMAPmixSBM(allAdj, graphGroups, countStat, hyperParam)

  res <- list(
    nodeClusterings = countStat$Z,
    theta = thetaMixSBM[[1]]$theta,
    ICL = bestICL
  )

  return(res)
}




#' Merge two graph groups
#'
#' Merge two graph groups to a single cluster and update the count statistics
#' accordingly.
#'
#' @param Adj1 list of adjacency matrices of first group
#' @param countStat1 count statistics associated with Adj1
#' @param Adj2 list of adjacency matrices of second group
#' @param countStat2 count statistics associated with Adj2
#' @param hyperParam hyperparameters of prior distributions
#'
#' @return list with updated count statistics for Adj1 and Adj2 after
#' cluster aggregation
#' @noRd
merge2graphGroups <- function(Adj1, countStat1, Adj2, countStat2, hyperParam){
  M1 <- length(countStat1$S)
  M2 <- length(countStat2$S)
  K1 <- length(countStat1$S[[1]])
  K2 <- length(countStat2$S[[1]])

  theta1 <- estimMAPCollectSBM(Adj1, countStat1, hyperParam)
  theta2 <- estimMAPCollectSBM(Adj2, countStat2, hyperParam)

  # normally, order should be ok, but check to make sure
  thetaSort1 <- degreeSort(theta1, outPerm = TRUE)
  if( sum(thetaSort1$permut == (1:K1)) < K1){
    countStat1 <- permutCountStat(countStat1, thetaSort1$permut)
    theta1 <- thetaSort1$theta
  }
  thetaSort2 <- degreeSort(theta2, outPerm = TRUE)
  if( sum(thetaSort2$permut==(1:K2)) < K2){
    countStat2 <- permutCountStat(countStat2, thetaSort2$permut)
    theta2 <- thetaSort2$theta
  }

  if (K1!=K2){
    if (K1<K2){
      countStat1 <- mergeUnequalBlocks(countStat1, countStat2$Z, Adj1, Adj2,
                                       theta1, theta2, hyperParam)
    }else{
      countStat2 <- mergeUnequalBlocks(countStat2, countStat1$Z, Adj2, Adj1,
                                       theta2, theta1, hyperParam)
    }
  }

  # apply ICL to improve clustering
  resICL_countStat <- maxICL(c(Adj1, Adj2),
                             mergeCountStat(countStat1, countStat2),
                             hyperParam)

  countStat1 <- extractCountStat(resICL_countStat, 1:M1)
  countStat2 <- extractCountStat(resICL_countStat, (M1+1):(M1+M2))

  return(list(countStat1 = countStat1, countStat2 = countStat2))
}


#' Merge two blocks of a SBM to a single one
#'
#' @param Z vector of node labels
#' @param blocksToFuse vector of length 2 containing indices of blocks to be merged
#'
#' @return vector with new node labels
#' @noRd
fuse2blocks <- function(Z, blocksToFuse){# Z is a list
  Znew <- lapply(Z, function(el){
    el[el %in% blocksToFuse] <- blocksToFuse[1] # ???
    return(el)
  })

  # rename block labels such that no blocks are empty
  blockLabels <- unique(unlist(Znew))
  finalK <- length(blockLabels)
  blockOrder <- order(blockLabels)
  Zfinal <- Znew
  M <- length(Zfinal)
  for (j in 1:finalK){
    Zfinal <- lapply(1:M, function(m){
      Zfinal[[m]][Znew[[m]]==blockLabels[j]] <- (1:finalK)[blockOrder[j]]
      return(Zfinal[[m]])
    })
  }

  return(Zfinal)
}


#' Update of count statistics when mergin two graph groups whose associated
#' stochastic block models do not have
#' the same number of blocks
#'
#' Merge two graph groups to a single cluster and update the count statistics
#' accordingly.
#'
#' @param countStat.a count statistics associated with Adj.a
#' @param Z.b vector of node labels of the other graph group
#' @param Adj.a list of adjacency matrices of the graph group that has the smaller number of blocks in the stochastic block model
#' @param Adj.b list of adjacency matrices of the graph group that has the larger number of blocks in the stochastic block model
#' @param theta.a parameter of the stochastic block model associated with Adj.a
#' @param theta.b parameter of the stochastic block model associated with Adj.b
#' @param hyperParam hyperparameters of prior distributions
#'
#' @return updated count statistics associated with Adj.a
#' @noRd
mergeUnequalBlocks <- function(countStat.a, Z.b, Adj.a, Adj.b, theta.a, theta.b,
                                hyperParam){
  K.a <- length(theta.a$pi)
  K.b <- length(theta.b$pi)
  bestFusionL2 <- Inf
  diffK <- K.b - K.a
  for (s in 1:K.a){
    blocksToFuse <- s:(s+diffK)
    Z2fused <- fuse2blocks(Z.b, blocksToFuse)
    countsFused <- getListCounts(Adj.b, Z2fused, K.a) # output is automatically of size K.a
    theta2fused <- estimMAPCollectSBM(Adj.b, countsFused, hyperParam)
    L2fused <- graphonL2norm(theta.a, theta2fused)
    if (L2fused<bestFusionL2){
      bestFusionL2 <- L2fused
      bestBlocksToFuse <- blocksToFuse
    }
  }
  # split bestBlocksToFuse[1] into two blocks randomly
  countStat.a <- splitBlock(Adj.a, countStat.a, bestBlocksToFuse[1],
                            length(bestBlocksToFuse), FALSE, hyperParam)

  return(countStat.a)
}


#' Increase the number of blocks of a stochastic block model associated with a
#' single network
#'
#' for a single network and associated node labels, split a designated block randomly
#' into several smaller blocks
#'
#' @param Z vector of node labels
#' @param k block k is to be split
#' @param S split block k into S blocks
#' @param theta parameter of the stochastic block model
#' @param adj adjacency matrix
#'
#' @return vector with new node labels
#' @noRd
splitBlock_1graph <- function(Z, k, S, K, adj, rowClustering=TRUE){  ## we work
  # with a single network : Z is a vector, adj is a single matrix
  newK <- K + S - 1
  if (k<K){# relabel last blocks
    for (l in K:(k+1))
      Z[Z==l] <- l + S - 1
  }

  toBeSplit <- Z==k
  if (sum(toBeSplit)>1){
    subAdj <- if(rowClustering) adj[toBeSplit, ] else t(adj[, toBeSplit])
    smallerBlocks <- try(stats::kmeans(subAdj, S)$cluster, silent = TRUE)
    if (!is.numeric(smallerBlocks)){
      smallerBlocks <- sample(1:S, sum(toBeSplit), replace = TRUE)
    }
    Z[toBeSplit] <- smallerBlocks + k - 1
  }

  return(Z)
}


#' Split block
#'
#' @param listAdj list of adjacency matrices
#' @param countStat list of associated count statistics
#' @param k_star block k_star is to be split. If NULL (default), then a block is chosen at random.
#' @param S block k_star is split into S smaller blocks. Default: S=2.
#' @param applyICL If TRUE, apply a round of ICL maximization after the split,
#' which may undone the proposed split if irrelevant. Default: FALSE.
#' @param hyperParam hyperparameters of prior distributions
#'
#' @return list of new count statistics
#' @noRd
splitBlock <- function(listAdj, countStat, k_star=NULL, S=2, applyICL=FALSE,
                       hyperParam){
  if (applyICL){
    oldICL <- iclCriterionSBMiid(countStat, hyperParam)
  }

  theta <- estimMAPCollectSBM(listAdj, countStat, hyperParam)
  K <- length(theta$pi)

  if(is.null(k_star)){
    k_star <- sample(1:K, 1, prob=theta$pi)
  }

  M <- length(listAdj)
  newZ <- lapply(1:M, function(m)
    splitBlock_1graph(countStat$Z[[m]], k_star, S, K,
                      listAdj[[m]], rowClustering=sample(c(TRUE, FALSE), 1)))

  countStatNew <- getListCounts(listAdj, newZ, K=K+S-1)

  if (applyICL){
    countStatNewICL <- maxICL(listAdj, countStat, hyperParam)
    newICL <- iclCriterionSBMiid(countStatNewICL, hyperParam)
    if ((newICL-oldICL) > 0){
      countStat <- countStatNewICL
    }
  }else{
    countStat <- countStatNew
  }

  # sort sbm blocks for identifiability
  theta <- estimMAPCollectSBM(listAdj, countStat, hyperParam)
  thetaSort <- degreeSort(theta, outPerm = TRUE)
  K <- length(theta$pi)
  if( sum(thetaSort$permut == (1:K)) < K){
    countStat <- permutCountStat(countStat, thetaSort$permut)
    theta <- thetaSort$theta
  }

  return(countStat)
}

#' relabel clusters to avoid cluster labels referring to empty clusters
#'
#' @param graphGroups vector with a graph clustering
#'
#' @return clustering vector with relabeled clusters
#' @noRd
relabelGraphGroups <- function(graphGroups){
  groupNames <- sort(unique(graphGroups))
  G <- length(groupNames)
  for (g in 1:G){
    graphGroups[graphGroups==groupNames[g]] <- g
  }
  return(graphGroups)
}


