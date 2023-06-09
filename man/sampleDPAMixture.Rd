% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/directedPrefAttachment.R
\name{sampleDPAMixture}
\alias{sampleDPAMixture}
\title{Generation of a mixture of directed preferential attachment (DPA) models}
\usage{
sampleDPAMixture(M, param)
}
\arguments{
\item{M}{number of desired networks}

\item{param}{list of list of parameters of the DPA models. Each element of
param is a list with the following elements: $prop (weight of the mixture
component), $R (= number of iterations),
$alpha, $beta , $deltaIn, $deltaOut (parameters of the DPA model)}
}
\value{
list of 2 lists : the first ($listAdj) is a list of M adjacency matrices, the
second a list ($graphGroups) contains the true cluster labels
}
\description{
Generation of a mixture of directed preferential attachment (DPA) models
}
\examples{
param <- vector('list', 3)
param[[1]] <- list(prop = 1/3, # component 1 :  alpha > beta
                   alpha = .04,
                   beta = .02,
                   deltaIn = 100,
                   deltaOut = 100,
                   R = 500
)
param[[2]] <- list(prop = 1/3, # component 2 : just permute alpha and beta ;
                   alpha = .01,
                   beta = .02,
                   deltaIn = 100,
                   deltaOut = .1,
                   R = 1000
)
param[[3]] <- list(prop = 1/3, # component 3 : alpha=beta
                   alpha = .015,
                   beta = .015,
                   deltaIn = .1,
                   deltaOut = .1,
                   R = 1000
)
obs <- sampleDPAMixture(M=20, param)
}
