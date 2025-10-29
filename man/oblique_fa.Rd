\name{oblique_fa}
\alias{oblique_fa}
\title{Gaussian oblique factor structure correlation matrix}
\description{
MLE of parameters in the Gaussian oblique factor model
for d variables and m groups, 
}

\usage{
oblique_fa(grsize, start, data=1, cormat=NULL, 
                       n=100, prlevel=0, mxiter=100)   

}
\arguments{
\item{grsize}{vector of group sizes (variables ordered by group)}
\item{start}{starting point should have dimension d+m*(m-1)/2}
\item{data}{nxd data set to compute the correlation matrix if cormat not given}
\item{cormat}{dxd empirical correlation matrix (of normal scores)}
\item{n}{sample size, if available}
\item{prlevel}{print.level for nlm()}
\item{mxiter}{maximum number of iterations for nlm()}
}
\value{
  list with nllk: negative log-likelihood; rhovec: the estimated mle;
loadings: loading matrix;
cor_lat: correlation matrix of the latent variables;
Rmod: the correlation matrix with optimized parameters.
}
\examples{
 # See example in bifact_fa() for a comparison with a data example
# Simpler example below
rhpar = c(0.81,0.84,0.84, 0.54,0.57,0.49, 0.51,0.54,0.55,0.70,  0.53,0.56,0.53,0.67,0.70)
cormat = corvec2mat(rhpar)
grsize = c(3,3)
mgrp = length(grsize)
d = sum(grsize)
start = rep(0.7,d+mgrp*(mgrp-1)/2)
res = oblique_fa(grsize, start, cormat=cormat, n=100, prlevel=1)
# iteration = 18
# Parameter:
# [1] 0.9005171 0.9064707 0.9270258 0.8202611 0.8480912 0.8243349 0.7039589
# Function Value
# [1] 622.5533
# Gradient:
# [1]  1.796252e-05  1.464286e-04  9.424639e-05 -8.344614e-05 -9.640644e-05
# [6] -9.799805e-05 -1.705303e-05
#
# Relative gradient close to zero.
# Current iterate is probably solution.
}
