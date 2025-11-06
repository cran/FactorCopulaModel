\name{bifactor2cor_v2}
\alias{bifactor2cor_v2}
\title{Bi-factor partial correlations to correlation matrix  version 2, using the inverse and determinant of a smaller matrix 
}
\description{
Bi-factor partial correlations to correlation matrix, determinant, inverse
}

\usage{
bifactor2cor_v2(grsize,rh1,rh2)
}
\arguments{
\item{grsize}{vector with group sizes: d_1,d_2,...,d_G for G groups}
\item{rh1}{vector of length sum(grsize) of correlation with global latent variable, ordered by group index}
\item{rh2}{vector of length sum(grsize) of partial correlation with group latent variable given global}
}
\value{
 list with Rmat = correlation matrix; det = det(Rmat); Rinv = solve(Rmat)
}
\examples{
# see examples for bifactor2cor()
}
