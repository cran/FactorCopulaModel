\name{ml1factor_f90}
\alias{ml1factor_f90}
\title{min negative log-likelihood for 1-factor copula with nlm()}
\description{
min negative log-likelihood for 1-factor copula with nlm()
}

\usage{
ml1factor_f90(nq,start,udata,copname,LB=0,UB=40,ihess=FALSE,prlevel=0,
  mxiter=100,nu=3)

}
\arguments{
\item{nq}{number of quadrature points}
\item{start}{starting point (d-vector or m*d vector, e.g. 2*d vector for BB1)}
\item{udata}{nxd matrix of uniform scores}
\item{copname}{name of copula family such as "gumbel", "frank", "bb1", "t" (copname common for all variables)}
\item{LB}{lower bound on parameters (scalar or same dimension as start) }
\item{UB}{upper bound on parameters (scalar or same dimension as start)}
\item{ihess}{flag for hessian option in nlm()}
\item{prlevel}{printlevel for nlm()}
\item{mxiter}{max number of iterations for nlm()}
\item{nu}{degree of freedom parameter if copname ="t"}
}
\value{
 MLE as nlm object (estimate, Hessian, SEs, nllk)}
\examples{
# See example in r1factor() 
}
