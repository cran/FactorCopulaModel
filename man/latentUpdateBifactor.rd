\name{latentUpdateBifactor}
\alias{latentUpdateBifactor}
\title{Conditional expectation proxies for bi-factor copula models   with linking copulas in different copula families 
}
\description{
Conditional expectation proxies for bi-factor copula models
with linking copulas in different copula families
}

\usage{
latentUpdateBifactor(udata,cparvec,grsize,family, nq)
}
\arguments{
\item{udata}{nxd matrix with valies in (0,1)}
\item{cparvec}{parameters for linking copulas;  order is global_par1, global_par2, local_par1, local_par2}
\item{grsize}{group size vector of length mgrp}
\item{family}{codes for linking copula (VineCopula)}
\item{nq}{number of Gaussian-Legendre points}
}
\value{
 v0: proxies of the global latent variable andvg: proxies of the local latent variables
}
\examples{
# See example in bifactorEstWithProxy()
}
