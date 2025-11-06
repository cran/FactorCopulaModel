\name{latentUpdate1factor}
\alias{latentUpdate1factor}
\title{Compute new proxies for 1-factor copula based on the mean of observations}
\description{
Compute new proxies for 1-factor copula for 1-parameter or 2-parameter linking copulas
}

\usage{
latentUpdate1factor(cpar_est,udata,nq,family) 
}
\arguments{
\item{cpar_est}{estimated parameters (based on complete likelihood with latent variables known or estimated)}
\item{udata}{nxd matrix of data in (0,1)}
\item{nq}{number of nodes for Gaussian-Legendre quadrature }
\item{family}{vector of code for d linking copula families (choices 1,2,4,5,7,10,14,17,20)}
}
\value{
 latent_est: proxies as estimates of latent variables}
\examples{
# See examples in onefactorEstWithProxy()
}
