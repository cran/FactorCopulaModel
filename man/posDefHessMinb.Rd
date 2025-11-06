\name{posDefHessMinb}
\alias{posDefHessMinb}
\title{Version with ifixed as argument}
\description{
modified Newton-Raphson minimization with positive Hessian 
}

\usage{
posDefHessMinb(param,objfn,ifixed,dstruct,LB,UB,mxiter=30,eps=1.e-6,
  bdd=5,iprint=FALSE) 

}
\arguments{
\item{param}{starting point for minimization }
\item{objfn}{function to be minimized with gradient and Hessian}
\item{ifixed}{vector of length(param) of TRUE/FALSE, such that ifixed[i]=TRUE iff param[i] is fixed at the given value}
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
