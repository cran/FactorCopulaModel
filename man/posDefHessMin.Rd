\name{posDefHessMin}
\alias{posDefHessMin}
\title{Minimization with modified Newton-Raphson iterations, Hessian is modified to be positive definite at each step. 
 Algorithm and code produced by Pavel Krupskii (2013) 
 see PhD thesis Krupskii (2014), UBC and 
 Section 6.2 of # Joe (2014) Dependence Models with Copulas. Chapman&Hall/CRC 
}
\description{
modified Newton-Raphson minimization with positive Hessian 
}

\usage{
posDefHessMin(param,objfn,dstruct,LB,UB,mxiter=30,eps=1.e-6,bdd=5,iprint=FALSE)
}
\arguments{
\item{param}{starting point for minimization }
\item{objfn}{function to be minimized with gradient and Hessian}
\item{dstruct}{list with  data set and other variables used by objfn}
\item{LB}{lower bound vector}
\item{UB}{upper bound vector}
\item{mxiter}{max number of iterations}
\item{eps}{tolerance for Newton-Raphson iterations}
\item{bdd}{bound on difference of 2 consecutive iterations (useful is starting point is far from solution and func is far from convex)}
\item{iprint}{control on amount of printing, FALSE for no printing of iterations and TRUE for printing x^(k) on each iteration.}
}
\value{
 list withfnval = function value at minimum;
parmin = param for minimum;
invh = inverse Hessian;
iconv = 1 if converged, -1 for a boundary point, 0 otherwise;
iter = number of iterations.
}
\examples{
# See examples in onefactorcop_nllk(), bifactorcop_nllk(), nestfactorcop_nllk()
}
