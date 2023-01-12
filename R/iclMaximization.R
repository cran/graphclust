
#' Auxiliary function
#'
#' compute log(gamma(a+z)/gamma(a))
#'
#' @param a real number
#' @param z real number
#'
#' @return log(gamma(a+z)/gamma(a))
#' @noRd
logGamGam <- function(a, z){
  lgamma(a+z)-lgamma(a)
}

#' ICL maximization algorithm for a single stochastic block model for a sample of networks
#'
#' @param Adj list of adjacency matrices
#' @param countStat list of count statistics
#' @param hyperParam hyperparameters of prior distributions
#' @param maxNbEpochs maximal number of epochs
#'
#' @return final list of count statistics
#' @noRd
maxICL <- function(Adj, countStat, hyperParam, maxNbEpochs=3){
  M <- length(Adj)
  K <- max(unlist(countStat$Z))
  S <- matrix(unlist(countStat$S), K, M)
  A <- array(unlist(countStat$A), dim = c(K, K, M))
  R <- array(unlist(countStat$R), dim = c(K, K, M))

  nNodes <- colSums(S) # M-vec

  nbAllNodes <- sum(nNodes)
  allNodes <- matrix(NA, nrow=nbAllNodes, ncol=2)
  allNodes[ ,1] <- unlist(sapply(1:M, function(m) 1:nNodes[m]))
  allNodes[ ,2] <- unlist(sapply(1:M, function(m) rep(m, nNodes[m])))

  notConverged <- TRUE
  countIter <- 0
  while(notConverged&(countIter<maxNbEpochs)){
    # select a random order to visit all nodes
    all_iStars_ind <- sample(1:nbAllNodes)
    notConverged <- FALSE
    countIter <- countIter + 1

    for (t in 1:nbAllNodes){
      iStar <- allNodes[all_iStars_ind[t], 1]
      mStar <- allNodes[all_iStars_ind[t], 2]
      g <- countStat$Z[[mStar]][iStar]
      n_mStar <- nNodes[mStar]

      Zmask <- matrix(0, K, n_mStar)
      Zmask[matrix(c(countStat$Z[[mStar]], 1:n_mStar), ncol=2)] <- 1
      delta_mStar <- matrix(NA, 2, K)
      delta_mStar[1, ] <- c(Zmask%*%Adj[[mStar]][iStar, ]) # iStar dot
      delta_mStar[2, ] <- c(Zmask%*%Adj[[mStar]][ , iStar]) # dot iStar

      sum_a_m <- apply(A, 1:2, sum)  # KxK matrix
      sum_r_m <- apply(R, 1:2, sum)  # KxK matrix
      sum_s_m <- rowSums(S)

      # compute delta icl
      resDelta <- DeltaICL_iStar(g, S_mStar=S[ , mStar], delta_mStar, sum_s_m,
                                 sum_a_m, sum_r_m, hyperParam)

      bestMove_h <- which.max(resDelta$DeltaICL)

      if (bestMove_h!=g){ # perform swap
        notConverged <- TRUE

        # update Z, S, R, A
        countStat$Z[[mStar]][iStar] <- bestMove_h

        S[c(g, bestMove_h), mStar] <- S[c(g, bestMove_h) , mStar] + c(-1, 1)
        A[ , , mStar] <- A[ , , mStar] + resDelta$diff_a_mStar[ , , bestMove_h]
        R[ , , mStar] <- R[ , , mStar] + resDelta$diff_r_mStar[ , , bestMove_h]

        if (sum_s_m[g]==1){# diminish K
          K <- K-1
          S <- S[-g ,]
          if(!is.matrix(S)){
            S <- matrix(S, K, M)
          }
          A <- A[-g, -g, ]
          if(!is.array(A)){
            A <- array(A, c(K,K,M))
          }
          R <- R[-g, -g, ]
          if(!is.array(R)){
            R <- array(R, c(K,K,M))
          }
          if (g<=K){ # adapt labels of Z
            countStat$Z <- lapply(countStat$Z, function(el){
              el[el>g] <- el[el>g] - 1
              return(el)
            }
            )
          }
        }
      }
    }
  }

  # sort sbm blocks for identifiability
  countStat <- buildCountStat(countStat$Z, S, A, R)

  theta <- estimMAPCollectSBM(Adj, countStat, hyperParam)
  thetaSort <- degreeSort(theta, outPerm = TRUE)
  K <- length(theta$pi)
  if( sum(thetaSort$permut == (1:K)) < K){
    countStat <- permutCountStat(countStat, thetaSort$permut)
    theta <- thetaSort$theta
  }

  return(countStat)
}

#' Evaluate variation of ICL criterion of changing the block label of a given node
#'
#' @param g indice of node to be changed
#' @param S_mStar some count statistics
#' @param delta_mStar some count statistics
#' @param sum_s_m some count statistics
#' @param sum_a_m some count statistics
#' @param sum_r_m some count statistics
#' @param hyperParam hyperparameters of prior distributions
#'
#' @return a list with multiple information on the variation of ICL criterion
#' resulting from a change of the block label of a given node
#' @noRd
DeltaICL_iStar <- function(g, S_mStar, delta_mStar, sum_s_m, sum_a_m, sum_r_m,
                           hyperParam){
  K <- length(S_mStar)

  DeltaICL <- rep(0, K)
  DeltaICL[-g] <- log((sum_s_m[-g] + hyperParam$alpha)/(sum_s_m[g] + hyperParam$alpha - 1))

  sum_b_m <- sum_r_m - sum_a_m
  diff_r_mStar <- diff_a_mStar <- array(0, dim=c(K, K, K))
  for (h in setdiff(1:K, g)){
    diff_r_mStar[g, , h] <- -S_mStar
    diff_r_mStar[h, , h] <- S_mStar
    diff_r_mStar[ , g, h] <- diff_r_mStar[ , g, h] - S_mStar
    diff_r_mStar[ , h, h] <- diff_r_mStar[ , h, h] + S_mStar
    diff_r_mStar[g ,g, h] <- diff_r_mStar[g ,g, h] + 2
    diff_r_mStar[g ,h, h] <- diff_r_mStar[g ,h, h] - 1
    diff_r_mStar[h ,g, h] <- diff_r_mStar[h ,g, h] - 1

    DeltaICL[h] <- DeltaICL[h] -
      sum(logGamGam(hyperParam$eta+hyperParam$zeta+sum_r_m[g, ], diff_r_mStar[g, , h])) -
      sum(logGamGam(hyperParam$eta+hyperParam$zeta+sum_r_m[h, ], diff_r_mStar[h, , h])) -
      sum(logGamGam(hyperParam$eta+hyperParam$zeta+sum_r_m[-c(g,h), g], diff_r_mStar[-c(g,h), g, h])) -
      sum(logGamGam(hyperParam$eta+hyperParam$zeta+sum_r_m[-c(g,h), h], diff_r_mStar[-c(g,h), h, h]))

    diff_a_mStar[g, , h] <- -delta_mStar[1, ]
    diff_a_mStar[h, , h] <- delta_mStar[1, ]
    diff_a_mStar[ , g, h] <- diff_a_mStar[ , g, h] - delta_mStar[2, ]
    diff_a_mStar[ , h, h] <- diff_a_mStar[ , h, h] + delta_mStar[2, ]

    DeltaICL[h] <- DeltaICL[h] +
      sum(logGamGam(hyperParam$eta+sum_a_m[g, ], diff_a_mStar[g, , h])) +
      sum(logGamGam(hyperParam$eta+sum_a_m[h, ], diff_a_mStar[h, , h])) +
      sum(logGamGam(hyperParam$eta+sum_a_m[-c(g,h), g], diff_a_mStar[-c(g,h), g, h])) +
      sum(logGamGam(hyperParam$eta+sum_a_m[-c(g,h), h], diff_a_mStar[-c(g,h), h, h]))

    diff_b_mStar <- diff_r_mStar[ , , h] - diff_a_mStar[ , , h]
    DeltaICL[h] <- DeltaICL[h] +
      sum(logGamGam(hyperParam$zeta+sum_b_m[g, ], diff_b_mStar[g, ])) +
      sum(logGamGam(hyperParam$zeta+sum_b_m[h, ], diff_b_mStar[h, ])) +
      sum(logGamGam(hyperParam$zeta+sum_b_m[-c(g,h), g], diff_b_mStar[-c(g,h), g])) +
      sum(logGamGam(hyperParam$zeta+sum_b_m[-c(g,h), h], diff_b_mStar[-c(g,h), h]))
  }

  if (sum_s_m[g]==1){# add further term
    extraTerm <- lgamma((K-1)*hyperParam$alpha) - lgamma(K*hyperParam$alpha) +
      lgamma(K*hyperParam$alpha + sum(sum_s_m)) - lgamma((K-1)*hyperParam$alpha + sum(sum_s_m))
    DeltaICL[-g] <- DeltaICL[-g] + extraTerm
  }

  return(list(DeltaICL=DeltaICL, diff_r_mStar=diff_r_mStar, diff_a_mStar=diff_a_mStar))
}








