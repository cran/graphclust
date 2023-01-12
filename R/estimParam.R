#' MAP estimator of parameters of a stochastic block model for a collection of
#' i.i.d. networks
#' @noRd
#' @param listAdj list of adjacency matrices
#' @param countStat list of count statistics
#' @param hyperParam hyperparameters of prior distributions
#'
#' @return list with MAP estimators of parameters of the stochastic block model
estimMAPCollectSBM <- function(listAdj, countStat, hyperParam){
  # only appropriate for directed networks
  eps <- 1e-6

  K <- length(countStat$S[[1]])
  M <- length(listAdj)
  n_m <- sapply(countStat$S, sum) # nb of nodes per network (M-vector)

  sum_s_k <- rowSums(matrix(unlist(countStat$S), nrow=K))  # K-vec
  nenner <- (sum(n_m)+K*(hyperParam$alpha-1))
  if (nenner>0){
    pi <- (sum_s_k+hyperParam$alpha-1)/nenner
    if (sum(pi<0)>0){
      pi[pi<0] <- eps
      pi <- pi/sum(pi)
    }
  }else{
    nenner <- sum(n_m)
    pi <- sum_s_k/nenner
  }

  sum_a_kl <- apply(array(unlist(countStat$A), dim=c(K, K, M)), c(1, 2), sum)  # KxK

  sum_r_kl <- apply(array(unlist(countStat$R), dim=c(K, K, M)), c(1, 2), sum)  # KxK

  zahler <- sum_a_kl+hyperParam$eta-1
  nenner <- sum_r_kl+hyperParam$eta+hyperParam$zeta-2

  gamma <- zahler/nenner

  ind <- (gamma<=1) & (gamma>=eps)
  if(sum(!ind)>0){
    gamma[!ind] <- eps
  }
  theta <- degreeSort(list(pi=pi, gamma=gamma))
  return(theta)
}

#' MAP estimator of parameters of a mixture of stochastic block models for a
#' collection of networks
#' @noRd
#' @param listAdj list of adjacency matrices
#' @param graphGroups vector with graph clustering
#' @param countStat list of count statistics
#' @param hyperParam hyperparameters of prior distributions
#' @param directed Default TRUE. Networks must be directed.
#'
#' @return list with MAP estimators of parameters of a mixture of stochastic block models
estimMAPmixSBM <- function(listAdj, graphGroups, countStat, hyperParam,
                           directed=TRUE){
  M <- length(listAdj)
  groupNames <- unique(graphGroups)
  propGroups <- as.numeric(table(graphGroups))/M
  thetaMixSBM <- lapply(groupNames,
                        function (gr) list(label = gr,
                                           prop = propGroups[groupNames==gr],
                                           theta = degreeSort(estimMAPCollectSBM(listAdj[graphGroups==gr],
                                                                                 extractCountStat(countStat, (1:M)[graphGroups==gr]),
                                                                                 hyperParam) )))

  return(thetaMixSBM)
}

#' Fit a stochastic block model to every network in a collection of networks.
#'
#' Applies the variational EM-algorithm implemented in the package blockmodels to every network.
#'
#' @param allAdj list of adjacency matrices
#' @param directed Networks are directed (TRUE by default) or undirected (FALSE).
#' @param nbCores number of cores for parallelization. Default: detectCores().
#' @param outCountStat If TRUE (default), the output is a list of count
#' statistics for every network. If FALSE, the output is a list of parameters of
#' the stochastic block models fitted to every network.
#'
#' @return list of count statistics for every network or list of parameters of
#' the stochastic block models fitted to every network.
#' @export
#'
#' @examples
#' theta <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
#' obs <- rCollectSBM(rep(10,4), theta)$listGraphs
#' res <- fitSimpleSBM(obs, outCountStat=FALSE, nbCores=2)
fitSimpleSBM <- function (allAdj, directed=TRUE, nbCores=detectCores(), outCountStat=TRUE){
  M <- length(allAdj)
  if ((nbCores==1) | (M < 10)){
    res <- fitSimpleSBM_sequential(allAdj, directed, outCountStat, nbCores=nbCores)
  }else{
    res <- fitSimpleSBM_parallel(allAdj, nbCores, directed, outCountStat)
  }
  return(res)
}


#' Fit a stochastic block model to every network in a collection of networks without parallelization.
#'
#' @noRd
#' @param allAdj list of adjacency matrices
#' @param directed Networks are directed (TRUE by default) or undirected (FALSE).
#' @param outCountStat If TRUE (default), the output is a list of count
#' statistics for every network. If FALSE, the output is a list of parameters of
#' the stochastic block models fitted to every network.
#' @param silent If FALSE (default), prints a message when the computation starts.
#'
#' @return list of count statistics for every network or list of parameters of
#' the stochastic block models fitted to every network.
#' @importFrom blockmodels BM_bernoulli
fitSimpleSBM_sequential <- function (allAdj, directed=TRUE, outCountStat=TRUE, silent=FALSE, nbCores){
  M <- length(allAdj)
  clustering <- vector("list", M)
  if (!outCountStat)
    listTheta <- vector("list", M)
  for (m in 1:M){
    if (!silent){
      message("initialization of network", m)
    }
    if (directed){
      myModel <- BM_bernoulli("SBM", allAdj[[m]], verbosity=0, plotting='', ncores=nbCores)
    }
    else{
      myModel <- BM_bernoulli("SBM_sym", allAdj[[m]], verbosity=0, plotting='', ncores=nbCores)
    }
    myModel$estimate()
    selectedModel <- which.max(myModel$ICL) + 1

    theta <- list(
      pi = colMeans(myModel$memberships[[selectedModel]]$Z),
      gamma = myModel$model_parameters[[selectedModel]]$pi)
    permut <- degreeSort(theta, outTheta=FALSE, outPerm=TRUE)

    clustering[[m]] <- permutLabels(myModel$memberships[[selectedModel]]$map()$C,
                                    permut)

    # avoid empty blocks
    if(length(unique(clustering[[m]])) < max(unique(clustering[[m]]))){
      theta <- list(
        pi = colMeans(myModel$memberships[[selectedModel-1]]$Z),
        gamma = myModel$model_parameters[[selectedModel-1]]$pi)
      permut <- degreeSort(theta, outTheta=FALSE, outPerm=TRUE)

      clustering[[m]] <- permutLabels(myModel$memberships[[selectedModel-1]]$map()$C,
                                      permut)
    }

    if (!outCountStat){
      listTheta[[m]] <- theta
    }
  }

  if (outCountStat){
    countStat <- getListCounts(allAdj, clustering)
    res <- countStat
  }else{
    res <- listTheta
  }

  return(res)
}




#' Fit a stochastic block model to every network in a collection of networks with parallelization.
#'
#' @noRd
#' @param allAdj list of adjacency matrices
#' @param nbCores number of cores for parallelization. Default: 2.
#' @param directed Networks are directed (TRUE by default) or undirected (FALSE).
#' @param outCountStat If TRUE (default), the output is a list of count
#' statistics for every network. If FALSE, the output is a list of parameters of
#' the stochastic block models fitted to every network.
#'
#' @return list of count statistics for every network or list of parameters of
#' the stochastic block models fitted to every network.
#' @importFrom parallel mclapply detectCores
fitSimpleSBM_parallel <- function (allAdj, nbCores=2, directed=TRUE, outCountStat=TRUE){
  M <- length(allAdj)
  nbCores <- min(c(ceiling(M/10), nbCores))
  Mpart <- ceiling(M/nbCores)
  init.values <- vector('list', nbCores)
  ind <- 1
  for (k in 1:nbCores){
    indEnd <- if (k<nbCores) ind + Mpart-1 else M
    init.values[[k]] <- allAdj[ind:indEnd]
    ind <- ind + Mpart
  }

  res_parallel <- mclapply(init.values,
                           function(el){
                             fitSimpleSBM_sequential(allAdj = el, directed, outCountStat, silent=TRUE, nbCores=nbCores)
                           },
                           mc.cores = nbCores
  )

  if (outCountStat){
    countStat <- res_parallel[[1]]
    if(nbCores>1){
      for (k in 2:nbCores){
        countStat$S <- c(countStat$S, res_parallel[[k]]$S)
        countStat$A <- c(countStat$A, res_parallel[[k]]$A)
        countStat$R <- c(countStat$R, res_parallel[[k]]$R)
        countStat$Z <- c(countStat$Z, res_parallel[[k]]$Z)
      }
    }
    res <- countStat
  }else{
    listTheta <- res_parallel[[1]]
    if(nbCores>1){
      for (k in 2:nbCores){
        listTheta <- c(listTheta, res_parallel[[k]])
      }
    }
    res <- listTheta
  }

  return(res)
}



