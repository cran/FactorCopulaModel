\name{tailDep}
\alias{tailDep}
\title{Tail dependence parameter estimation}
\description{
Tail dependence parameter estimate based on extrapolating zeta(alpha)
}

\usage{
tailDep(u1,u2, lowertail=FALSE, eps=0.1, semictol=0.1, rank=TRUE,
  iprint=FALSE)

}
\arguments{
\item{u1}{nx1 vector with values (in (0,1) if rank=FALSE)}
\item{u2}{nx1 vector with values (in (0,1) if rank=FALSE)}
\item{lowertail}{TRUE if lower tail-weighted dependence measure, default is FALSE}
\item{eps}{tolerance (default 0.1) for rate parameter in method 2; non-linear (use method 2) if rate < 1-eps}
\item{semictol}{tolerance (default 0.1) for exceedance of normal semicorrelation to treat as tail dependent;  use something like semicol=-0.5 if normal scores plot suggest tail dependence }
\item{rank}{TRUE (default) if data matrix needs to be converted to uniform scores in (0,1)}
\item{iprint}{TRUE for intermediate prints}
}
\value{
 upper tail dependence parameter based on extrapolation of zeta(alpha) for large alpha}
\examples{
mytest = function(qcond,cpar,n=500,seed=123)
{ set.seed(seed)
  u1 = runif(n)
  u2 = qcond(runif(n),u1,cpar)
  #convert to uniform scores (marginals are usually not known)
  u1 = (rank(u1)-0.5)/n
  u2 = (rank(u2)-0.5)/n
  alp = c(1,5,10:20)
  zetaL = zetaDep(cbind(u1,u2),alp,rank=FALSE,lowertail=FALSE)
  zetaU = zetaDep(cbind(u1,u2),alp,rank=FALSE,lowertail=TRUE)
  print(cbind(alp,zetaL,zetaU))
  utd = tailDep(u1,u2, lowertail=FALSE, eps=0.1, semictol=0.1, rank=FALSE, iprint=TRUE)
  ltd = tailDep(u1,u2, lowertail=TRUE, eps=0.1, semictol=0.1, rank=FALSE, iprint=TRUE)
  cat(ltd,utd,"\n")
  utd = tailDep(u1,u2, lowertail=FALSE, eps=0.1, semictol=0.1, rank=FALSE, iprint=FALSE)
  ltd = tailDep(u1,u2, lowertail=TRUE, eps=0.1, semictol=0.1, rank=FALSE, iprint=FALSE)
  cat(ltd,utd,"\n")
  #par(mfrow=c(2,1))
  #zetaPlot(cbind(u1,u2),alp,ylim=c(0,1),inverse=FALSE)
  #zetaPlot(cbind(u1,u2),alp,ylim=c(0,1),inverse=TRUE)
  0
}
mytest(qcondFrank,3)
mytest(qcondbvtcop,c(0.6,5))
}
\references{
Lee D, Joe H, Krupskii P (2018). J Nonparametric Statistics, 30(2), 262-290
}
