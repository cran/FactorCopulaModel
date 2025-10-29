\name{residDep}
\alias{residDep}
\title{correlation matrix for 1-factor plus 1-truncated vine (for residual dependence)}
\description{
correlation matrix for 1-factor plus 1-truncated vine (for residual dependence)
}

\usage{
residDep(cormat,loading)
}
\arguments{
\item{cormat}{dxd correlation matrix}
\item{loading}{d-dimensional loading vector (for latent factor), -1<loading[j]<1}
}
\value{
 list with R = correlation matrix for structure of 1-factor+Markov tree residual dependence;incl = d*(d-1)/2 binary vector: indicator of edges [1,2], [1,3], [2,3], [1,4], ...[d-1,d] edges in tree with d-1 edges;
partcor = conditional correlation matrix given the latent variable.
}
\details{
MST algorithm with weights log(1-rho^2), rho's are partial correlations
fiven the latent variable.
not exported
}
