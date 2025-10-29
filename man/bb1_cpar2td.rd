\name{bb1_cpar2td}
\alias{bb1_cpar2td}
\title{BB1 copula parameter (theta,delta) to tail dependence parameters }
\description{
BB1 copula parameter (theta,delta) to tail dependence parameters 
}

\usage{
bb1_cpar2td(cpar)
}
\arguments{
\item{cpar}{copula parameter with theta>0, delta>1 (vector of length 2) or mx2 matrix with columns for theta and delta}
}
\value{
 vector or matrix with lower and upper tail dependence}
\examples{
cpar = matrix(c(0.5,1.5,0.8,1.2),byrow=TRUE,ncol=2)
bb1_cpar2td(cpar)
}
