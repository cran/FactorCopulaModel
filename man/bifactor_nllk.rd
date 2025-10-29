\name{bifactor_nllk}
\alias{bifactor_nllk}
\title{log-likelihood Gaussian bi-factor structure correlation matrix}
\description{
log-likelihood Gaussian bi-factor structure correlation matrix
}

\usage{
bifactor_nllk(rhvec,grsize,Robs,nsize)
}
\arguments{
\item{rhvec}{vector of length d*2 for partial correlation representation of loadings,}
\item{grsize}{vector of group sizes for bi-factor model}
\item{Robs}{dxd empirical correlation matrix}
\item{nsize}{sample size}
}
\value{
 negative log-likelihood and gradient for Gaussian bi-factor model}
