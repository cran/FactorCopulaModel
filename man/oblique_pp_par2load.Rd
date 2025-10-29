\name{oblique_pp_par2load}
\alias{oblique_pp_par2load}
\title{oblique factor correlation structure for d variables and m groups  include determinant and inverse 
}
\description{
For oblique factor correlation structure for d variables and m groups, 
convert the vector of parameters theta into
the loading matrix and correlation matrix of latent variables
The variables are assumed ordered by group.
}

\usage{
oblique_pp_par2load(theta,grsize, icheck=FALSE)
}
\arguments{
\item{theta}{vector of length d + m*(m-1)/2; d loading parameters followed by m*(m-1)/2 entries in correlation matrix of latent variables (lower triangle by row)}
\item{grsize}{vector of group sizes (variables ordered by group)}
\item{icheck}{flag, if TRUE checks are made }
}
\value{
 loadings: loading matrix; cor_lat: correlation matrix of latent variables; Rmod: correlation matrix based on theta for oblique factor model
}
\examples{
 theta = c(0.6,0.7,0.8,0.7,0.6,0.5,0.5)
oblique_pp_par2load(theta,grsize=c(3,3))
#$loadings
#[,1] [,2]
#[1,]  0.6  0.0
#[2,]  0.7  0.0
#[3,]  0.8  0.0
#[4,]  0.0  0.7
#[5,]  0.0  0.6
#[6,]  0.0  0.5
#
#$cor_lat
#[,1] [,2]
#[1,]  1.0  0.5
#[2,]  0.5  1.0
# 
#$Rmod
#[,1]  [,2] [,3]  [,4] [,5]  [,6]
#[1,] 1.00 0.420 0.48 0.210 0.18 0.150
#[2,] 0.42 1.000 0.56 0.245 0.21 0.175
#[3,] 0.48 0.560 1.00 0.280 0.24 0.200
#[4,] 0.21 0.245 0.28 1.000 0.42 0.350
#[5,] 0.18 0.210 0.24 0.420 1.00 0.300
#[6,] 0.15 0.175 0.20 0.350 0.30 1.000
}
