\name{ml1factor}
\alias{ml1factor}
\title{max likelihood (min negative log-likelihood) for 1-factor copula model }
\description{
min negative log-likelihood for 1-factor copula model 
}

\usage{
ml1factor(nq,start,udata,dcop,LB=0,UB=1.e2,prlevel=0,mxiter=100)
}
\arguments{
\item{nq}{number of quadrature points}
\item{start}{starting point (d-vector or m*d vector, e.g. 2*d vector for BB1)}
\item{udata}{nxd matrix of uniform scores}
\item{dcop}{name of function for a bivariate copula density (common for all variables)}
\item{LB}{lower bound on parameters (scalar or same dimension as start) }
\item{UB}{upper bound on parameters (scalar or same dimension as start)}
\item{prlevel}{printlevel for nlm()}
\item{mxiter}{maximum number of iteration for nlm()}
}
\value{
 nlm object with minimum, estimate, hessian at MLE}
\examples{
# See example in r1factor() 
}
