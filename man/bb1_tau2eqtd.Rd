\name{bb1_tau2eqtd}
\alias{bb1_tau2eqtd}
\title{BB1, given 0<tau<1, find theta and delta with lower tail dependence equal upper tail dependence 
}
\description{
BB1, given 0<tau<1, find theta and delta with equal lower/upper tail dependence
}

\usage{
bb1_tau2eqtd(tau,destart=1.5,mxiter=30,eps=1.e-6,iprint=FALSE)
}
\arguments{
\item{tau}{Kendall tau value}
\item{destart}{starting point for delta}
\item{mxiter}{maximum number of iterations}
\item{eps}{tolerance for convergence}
\item{iprint}{print flag for iterations}
}
\value{
 copula parameter (theta,delta) with ltd=utd given tau}
\examples{
bb1_tau2eqtd(c(0.1,0.2,0.5))
}
