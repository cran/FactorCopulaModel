\name{oblique_nllk}
\alias{oblique_nllk}
\title{log-likelihood Gaussian oblique factor structure correlation matrix}
\description{
negative log-likelihood of the Gaussian oblique factor model 
for d variables and m groups, 
}

\usage{
oblique_nllk(theta, grsize, Robs, nsize=100)
}
\arguments{
\item{theta}{vector of length d + m*(m-1)/2; d loading parameters followed by m*(m-1)/2 entries in correlation matrix of latent variables (lower triangle by row)}
\item{grsize}{vector of group sizes (variables ordered by group)}
\item{Robs}{dxd (empirical) correlation matrix of normal scores}
\item{nsize}{sample size used to get Robs if available}
}
\value{
 negative log-likelihood value of the oblique Gaussian factor modelwith fixed group size at MLE 
}
\examples{
 rhpar = c(0.81,0.84,0.84, 0.54,0.57,0.49, 0.51,0.54,0.55,0.70, 0.53,0.56,0.53,0.67,0.70)
cormat = corvec2mat(rhpar)
print(cormat)
#    [,1] [,2] [,3] [,4] [,5] [,6]
#[1,] 1.00 0.81 0.84 0.54 0.51 0.53
#[2,] 0.81 1.00 0.84 0.57 0.54 0.56
#[3,] 0.84 0.84 1.00 0.49 0.55 0.53
#[4,] 0.54 0.57 0.49 1.00 0.70 0.67
#[5,] 0.51 0.54 0.55 0.70 1.00 0.70
#[6,] 0.53 0.56 0.53 0.67 0.70 1.00
grsize = c(3,3)
mgrp = length(grsize)
d = sum(grsize)
theta = c(rep(0.3,d+mgrp*(mgrp-1)/2))
ml_obl = oblique_nllk(theta=theta, grsize, Robs=cormat)
print(ml_obl)
# 806.7432
}
