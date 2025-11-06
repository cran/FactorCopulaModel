\name{pfactor_fa}
\alias{pfactor_fa}
\title{Gaussian p-factor structure correlation matrix}
\description{
Gaussian p-factor structure correlation matrix with quasi-Newton
}

\usage{
pfactor_fa(factors,start,data=1,cormat=NULL,n=100,prlevel=0,mxiter=100)
}
\arguments{
\item{factors}{p = #factors}
\item{start}{starting point should have dimension 2*d}
\item{data}{nsize x d data set to compute the correlation matrix if correlation matrix (cormat) not given}
\item{cormat}{dxd empirical correlation matrix}
\item{n}{sample size}
\item{prlevel}{print.level for nlm()}
\item{mxiter}{maximum number of iterations for nlm()}
}
\value{
a list with $nllk, $rhmat = dxp matrix of partial correlations,
$loading = dxp loading matrix after varimax,
$rotmat = pxp rotation matrix used by varimax
}
\examples{
# See example in bifactor_fa()
}
