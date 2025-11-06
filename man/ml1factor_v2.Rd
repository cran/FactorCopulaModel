\name{ml1factor_v2}
\alias{ml1factor_v2}
\title{min negative log-likelihood for 1-factor copula model (some parameters can be fixed)}
\description{
min negative log-likelihood (nllk) for 1-factor copula model (some parameters can be fixed)
}

\usage{
ml1factor_v2(nq,start,ifixed,udata,dcop,LB=0,UB=1.e2,prlevel=0,mxiter=100)
}
\arguments{
\item{nq}{number of quadrature points}
\item{start}{starting point (d-vector or m*d vector, e.g. 2*d vector for BB1)}
\item{ifixed}{vector of length(param) of True/False, such that ifixed[i]=TRUE iff param[i] is fixed at the given value start[j]}
\item{udata}{nxd matrix of uniform scores}
\item{dcop}{name of function for a bivariate copula density (common for all variables)}
\item{LB}{lower bound on parameters (scalar or same dimension as start) }
\item{UB}{upper bound on parameters (scalar or same dimension as start)}
\item{prlevel}{printlevel for nlm()}
\item{mxiter}{maximum number of iteration for nlm()}
}
\value{
 nlm object with nllk value, estimate, hessian at MLE}
