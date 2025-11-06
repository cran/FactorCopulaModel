\name{bifactor2cor}
\alias{bifactor2cor}
\title{Bi-factor partial correlations to correlation matrix }
\description{
Bi-factor partial correlations to correlation matrix, determinant, inverse
}

\usage{
bifactor2cor(grsize,rh1,rh2)
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
grsize = c(5,5,3) 
d = sum(grsize)
bifpar = c(0.84,0.63,0.58,0.78,0.79, 0.87,0.80,0.74,0.71,0.57, 0.83,0.77,0.80,
0.67,0.58,0.15,0.70,0.47,   0.32,0.27,0.73,0.19,0.12,   0.35,0.23,0.53)
bifobj = bifactor2cor(grsize,bifpar[1:d],bifpar[(d+1):(2*d)])
rmat = bifobj$Rmat
print(det(rmat)-bifobj$det)
print(max(abs(solve(rmat)-bifobj$Rinv)))
bifobj2 = bifactor2cor_v2(grsize,bifpar[1:d],bifpar[(d+1):(2*d)])
rmat2 = bifobj2$Rmat
print(det(rmat2)-bifobj2$det)
print(max(abs(solve(rmat2)-bifobj2$Rinv)))
}
