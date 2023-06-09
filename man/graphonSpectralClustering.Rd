% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphonSpectralClustering.R
\name{graphonSpectralClustering}
\alias{graphonSpectralClustering}
\title{Graph clustering using the pairwise graphon distances and spectral clustering}
\usage{
graphonSpectralClustering(allAdj, nbClusters, sig = 0.1, nbCores = 1)
}
\arguments{
\item{allAdj}{list of adjacency matrices}

\item{nbClusters}{number of clusters to be found}

\item{sig}{parameter for Gaussian kernel used for the similarity matrix}

\item{nbCores}{number of cores for parallelization.}
}
\value{
list with the obtained graph clusteirng ($clust) and the matrix with the
pairwise graphon distances between all pairs of networks
}
\description{
Graph clustering using the pairwise graphon distances and spectral clustering
}
\examples{
theta <- list(pi=c(.5,.5), gamma=matrix((1:4)/8,2,2))
obs <- rCollectSBM(rep(10,4), theta)$listGraphs
res <- graphonSpectralClustering(obs, 2, nbCores=1)
}
