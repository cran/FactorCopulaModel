\name{oblique_grad_nllk}
\alias{oblique_grad_nllk}
\title{log-likelihood Gaussian oblique factor structure correlation matrix}
\description{
log-likelihood Gaussian oblique factor structure correlation matrix with gradient
for d variables and m groups, 
}

\usage{
oblique_grad_nllk(theta, grsize, Robs, nsize=100)
}
\arguments{
\item{theta}{vector of length d + m*(m-1)/2; d loading parameters followed by m*(m-1)/2 entries in correlation matrix of latent variables (lower triangle by row)}
\item{grsize}{vector of group sizes (variables ordered by group)}
\item{Robs}{dxd empirical correlation matrix}
\item{nsize}{sample size if available}
}
\value{
 negative log-likelihood and gradient for Gaussian p-factor model}
