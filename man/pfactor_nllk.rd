\name{pfactor_nllk}
\alias{pfactor_nllk}
\title{log-likelihood Gaussian p-factor structure correlation matrix}
\description{
log-likelihood Gaussian p-factor structure correlation matrix with gradient
}

\usage{
pfactor_nllk(rhvec,Robs,nsize)
}
\arguments{
\item{rhvec}{vector of length d*p with partial corr representation of loadings}
\item{Robs}{dxd empirical correlation matrix}
\item{nsize}{sample size}
}
\value{
 negative log-likelihood and gradient for Gaussian p-factor model}
