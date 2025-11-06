\name{semiCorTable}
\alias{semiCorTable}
\title{Semi-correlation table for a multivariate data set}
\description{
Semi-correlation table for several variables
}

\usage{
semiCorTable(mdat, varnames, inscore=FALSE)
}
\arguments{
\item{mdat}{nxd multivariate data set with d>=2 columns}
\item{varnames}{d-vector of (abbreviated) variable names }
\item{inscore}{TRUE if mdat has already been converted to normal scores, default FALSE}
}
\value{
 d*(d-1)/2 by 8-column dataframe with columns j1,j2,ncor,lcor,ucor,bvnsemic, varnames[j1] varnames[j2]for 2 variable indices, correlation of normal scores, lower semi-correlation,
upper semi-correlation and BVN semi-correlation assuming Gaussian copula
Stronger than Gaussian dependence in upper tail if ucor is larger than bvn semicor
}
\examples{
rmat = toeplitz(c(1,.7,.4))
print(rmat)
set.seed(1234)
zdat = rmvn(n=500,rmat)
set.seed(12345)
tdat = rmvt(n=500,rmat,nu=5)
vnames=c("V1","V2","V3")
zsemi = semiCorTable(zdat,vnames)
tsemi = semiCorTable(tdat,vnames)
cat("trivariate normal\n")
print(zsemi)
cat("trivariate t(5)\n")
print(tsemi)
}
