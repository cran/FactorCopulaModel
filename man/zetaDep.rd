\name{zetaDep}
\alias{zetaDep}
\title{Empirical version of zeta(alpha) tail-weighted dependence measure}
\description{
Empirical version of zeta(alpha) tail-weighted dependence measure
}

\usage{
zetaDep(dat,alpha,rank=TRUE,lowertail=FALSE)
}
\arguments{
\item{dat}{nx2 data matrix with values (in (0,1) if rank=FALSE)}
\item{alpha}{vector of alpha>0 for zeta measure}
\item{rank}{TRUE (default) if to convert data matrix to uniform scores in (0,1)}
\item{lowertail}{TRUE if lower tail-weighted dependence measure, default is FALSE}
}
\value{
Dependence measure zeta(alpha)}
\examples{
data(euro07gf)
udat = euro07gf$uscore
euro07names = colnames(udat)
d = ncol(udat)
for(j2 in 2:d)
{ for(j1 in 1:(j2-1))
  { zetaU = zetaDep(udat[,c(j1,j2)],alpha=15,rank=FALSE,lowertail=FALSE)
    zetaL = zetaDep(udat[,c(j1,j2)],alpha=15,rank=FALSE,lowertail=TRUE)
    zeta1 = zetaDep(udat[,c(j1,j2)],alpha=1,rank=FALSE,lowertail=FALSE)
    cat(j1,j2,round(zeta1,3),round(zetaL,3),round(zetaU,3),euro07names[j1],euro07names[j2],"\n")
  }
}
}
\details{
This is a central dependence measure if alpha =1 and
upper tail-weighted is alpha>>1
}
\references{
Lee D, Joe H, Krupskii P (2018). J Nonparametric Statistics, 30(2), 262-290
}
