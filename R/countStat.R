#' Compute the count statistics of a network and given node labels
#'
#' @param adj adjacency matrix
#' @param Z vector with node labels
#' @param K Number of blocks of the associated stochastic block model. If NULL
#'  (by default), the number of empty blocks present in the node label vector Z is used.
#'  Otherwise, a stochastic block model with K blocks is consider implying the presence
#'  of empty blocks.
#' @param directed directed network (TRUE by default) or undirected (FALSE).
#'
#' @return list with the number of nodes per block ($S), the number of edges
#'  for all pairs of blocks ($A) and the maximal number of edges
#'  for all pairs of blocks ($R)
#' @noRd
getCounts <- function(adj, Z, K=NULL, directed=TRUE){
  n_m <- nrow(adj)
  labNames <- as.numeric(names(table(Z)))

  if (is.null(K)) {
    s_k <- as.numeric(table(Z))
    K <- J <- length(s_k)
    s <- NULL
    # relabel Z:
    if (max(Z)>K){
      for (k in 1:K){
        Z[Z==labNames[k]] <- k
      }
    }
  }else{
    s_k <- rep(0, K)
    s <- as.numeric(table(Z))
    J <- length(labNames)
    for (j in 1:J){
      s_k[labNames[j]] <- s[j]
    }
  }

  Zmask <- matrix(0, K, n_m)
  Zmask[matrix(c(Z, 1:n_m), ncol=2)] <- 1
  a_kl  <-  Zmask %*% adj %*% t(Zmask)

  if (K>1){
    r_kl <- s_k %*% t(s_k) - diag(s_k)
  }else{
    r_kl <- s_k * (s_k-1)
  }
  return(list(S=s_k, A=a_kl, R=r_kl))
}


#' Compute the count statistics for all networks in a collection and given node labels
#'
#' @param listAdj list of adjacency matrices
#' @param listZ list of node labels
#' @param K Number of blocks of the associated stochastic block model. If NULL
#'  (by default), the number of empty blocks present in the node label vector Z is used.
#'  Otherwise, a stochastic block model with K blocks is consider implying the presence
#'  of empty blocks.
#' @param directed directed network (TRUE by default) or undirected (FALSE).
#'
#' @return list of lists: with the number of nodes per block ($S), the number of edges
#'  for all pairs of blocks ($A), the maximal number of edges
#'  for all pairs of blocks ($R) and the node labels ($Z)
#' @noRd
getListCounts <- function(listAdj, listZ, K=NULL, directed=TRUE){
  M <- length(listAdj)
  listCountsPerGraph <- lapply(1:M,
                               function(m) getCounts(listAdj[[m]], listZ[[m]], K)
  )

  listCountsPerStat <- list(S=NULL, A=NULL, R=NULL)
  listCountsPerStat$S <- lapply(listCountsPerGraph, function(el) el$S)
  listCountsPerStat$A <- lapply(listCountsPerGraph, function(el) el$A)
  listCountsPerStat$R <- lapply(listCountsPerGraph, function(el) el$R)

  listCountsPerStat$Z <- listZ
  return(listCountsPerStat)
}

#' Extract count statistics for a given set of networks
#'
#' @param countStat list of count statistics of a collection of networks
#' @param ind vector of network indices to be extracted
#'
#' @return list of count statistics of the desired networks
#' @noRd
extractCountStat <- function(countStat, ind){
  return(list(Z=countStat$Z[ind], S=countStat$S[ind],
              R=countStat$R[ind], A=countStat$A[ind]))
}


#' Update entries of the list of count statistics of the collection of networks
#'
#' @param countStat list of count statistics of a collection of networks
#' @param ind vector of network indices to be updated
#' @param newCountStat list of new count statistics
#'
#' @return updated list of count statistics of the collection
#' @noRd
updateCountStat <- function(countStat, ind, newCountStat){
  countStat$Z[ind] <- newCountStat$Z
  countStat$S[ind] <- newCountStat$S
  countStat$A[ind] <- newCountStat$A
  countStat$R[ind] <- newCountStat$R
  return(countStat)
}


#' Merge two lists of count statistics
#'
#' @param countStat1 list of count statistics
#' @param countStat2 list of count statistics
#'
#' @return merged list of count statistics
#' @noRd
mergeCountStat <- function(countStat1, countStat2){
  return(list(Z=c(countStat1$Z, countStat2$Z),
              S=c(countStat1$S, countStat2$S),
              R=c(countStat1$R, countStat2$R),
              A=c(countStat1$A, countStat2$A)))
}


#' Create list of count statistics from arrays of count statistics
#'
#' @param Z list of node labels
#' @param S matrix of count statistics
#' @param A array of count statistics
#' @param R array of count statistics
#'
#' @return list of count statistics
#' @noRd
buildCountStat <- function(Z, S, A, R){
  countStat <- list(Z=Z)
  M <- ncol(S)
  countStat$S <- lapply(1:M, function(m) S[ , m])
  countStat$A <- lapply(1:M, function(m) A[ , , m])
  countStat$R <- lapply(1:M, function(m) R[ , , m])
  return(countStat)
}

#' Update list of count statistics when permuting block labels
#'
#' @param countStat list of count statistics
#' @param permut permutation of block labels
#'
#' @return updated list of count statistics
#' @noRd
permutCountStat <- function(countStat, permut){
    countStat$Z <- permutListZ(countStat$Z, permut)
    countStat$S <- permutListOfVectors(countStat$S, permut)
    countStat$R <- permutListOfMatrices(countStat$R, permut)
    countStat$A <- permutListOfMatrices(countStat$A, permut)
    return(countStat)
  }

#' Permute labels of a vector of node labels
#'
#' @param Z vector of node labels
#' @param permut permutation of block labels
#'
#' @return vector of permuted node labels
#' @noRd
permutLabels <- function(Z, permut){
  Znew <- Z
  K <- length(permut)
  for(k in 1:K)
    Znew[Z==permut[k]] <- k

  return(Znew)
}

#' Permute labels of a list of vectors of node labels
#'
#' @param Z vector of node labels
#' @param permut permutation of block labels
#'
#' @return updated list of vectors of node labels
#' @noRd
permutListZ <- function(Z, permut){
  Z <- lapply(Z, function(z) permutLabels(z, permut))
  return(Z)
}


#' Permute rows and columns of a list of matrices
#'
#' @param A list of square matrices of same size
#' @param permut permutation of rows and columns of matrices
#'
#' @return permuted list of matrices
#' @noRd
permutListOfMatrices <- function(A, permut){
  A <- lapply(A, function(mat) permutMatrice(mat, permut))
  return(A)
}

#' Permute rows and columns of a matrix
#'
#' @param A square matrix
#' @param permut permutation of rows and columns of A
#'
#' @return permuted matrix
#' @noRd
permutMatrice <- function(A, permut){
  K <- length(permut)
  ind <- matrix(c(rep(permut, K), rep(permut, each=K)), ncol=2)
  A <- matrix(A[ind], K, K)
  return(A)
}


#' Permute entries of a list of vectors of same size
#'
#' @param A list of vectors of same size
#' @param permut permutation of vector entries
#'
#' @return list of permuted vectors
#' @noRd
permutListOfVectors <- function(A, permut){
  A <- lapply(A, function(vec) vec[permut])
  return(A)
}

