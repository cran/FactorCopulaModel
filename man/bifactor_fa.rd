\name{bifactor_fa}
\alias{bifactor_fa}
\title{Gaussian bi-factor structure correlation matrix}
\description{
Gaussian bi-factor structure correlation matrix with quasi-Newton
}

\usage{
bifactor_fa(grsize,start,data=1,cormat=NULL,n=100,prlevel=0,mxiter=100)
}
\arguments{
\item{grsize}{vector of group sizes for bi-factor model}
\item{start}{starting point should have dimension 2*d}
\item{data}{nsize x d data set to compute the correlation matrix if correlation matrix (cormat) not given}
\item{cormat}{dxd empirical correlation matrix}
\item{n}{sample size}
\item{prlevel}{print.level for nlm()}
\item{mxiter}{maximum number of iterations for nlm()}
}
\value{
 a list with$nllk, $parmat= dx2 matrix of correlations and partial correlations
}
\examples{
data(rainstorm)
rmat = rainstorm$cormat
n = nrow(rainstorm$zprecip)
d = ncol(rmat)
grsize = rainstorm$grsize
fa1 = pfactor_fa(1,start=rep(0.8,d),cormat=rmat,n=n,prlevel=1)
fa2 = pfactor_fa(2,start=rep(0.8,2*d),cormat=rmat,n=n,prlevel=1)
st3 = c(rep(0.7,grsize[1]),rep(0.1,d),rep(0.7,grsize[2]),rep(0.1,d),rep(0.7,grsize[3]))
fa3 = pfactor_fa(3,start=st3,cormat=rmat,n=n,prlevel=1)
names(fa3)
loadmat_rotated = fa3$loading%*%fa3$rotmat # more values closer to 0
# compare factanal
fa3b = factanal(factors=3,covmat=rmat)
compare = cbind(loadmat_rotated,fa3b$loadings)
print(round(compare,3)) # order of factors is different but interpretation similar
#
bifa = bifactor_fa(grsize,start=c(rep(0.8,d),rep(0.2,d)),cormat=rmat,n=n,prlevel=1)
mgrp = length(grsize)
# oblique factor model is much more parsimonious than bi-factor
obfa = oblique_fa(grsize,start=rep(0.7,d+mgrp), cormat=rmat, n=n, prlevel=1)
#
cat(fa1$nllk, fa2$nllk, fa3$nllk, bifa$nllk, obfa$nllk,"\n")
}
