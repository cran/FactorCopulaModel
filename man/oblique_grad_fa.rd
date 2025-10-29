\name{oblique_grad_fa}
\alias{oblique_grad_fa}
\title{Gaussian oblique factor structure correlation matrix}
\description{
MLE of parameters in the Gaussian oblique factor model 
for d variables and m groups, 
}

\usage{
oblique_grad_fa(grsize, start, data=1, cormat=NULL, 
                            n=100, prlevel=0, mxiter=100)   

}
\arguments{
\item{grsize}{vector of group sizes (variables ordered by group)}
\item{start}{starting vector of length d + m*(m-1)/2; d loading parameters followed by m*(m-1)/2 entries in correlation matrix of latent variables (lower triangle by row)}
\item{data}{n x d data set to compute the correlation matrix if correlation matrix (cormat) not given}
\item{cormat}{dxd (empirical) correlation matrix of normal scores}
\item{n}{sample size, if available}
\item{prlevel}{print.level for nlm()}
\item{mxiter}{maximum number of iterations for nlm()}
}
\value{
  list with nllk: negative log-likeilihood;rhovec: the estimated mle;
loadings: loading matrix;
cor_lat: correlation matrix of the latent variables; 
Rmod: the correlation matrix with optimized parameters.
}
