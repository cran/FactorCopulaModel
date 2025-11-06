\name{d1factcop}
\alias{d1factcop}
\title{Integrand for 1-factor copula with 1-parameter bivariate linking copula families; or for m-parameter bivariate linking copulas  
}
\description{
Integrand for 1-factor copula 
}

\usage{
d1factcop(u0,uvec,dcop,param)
}
\arguments{
\item{u0}{latent variable for integrand}
\item{uvec}{vector of length d, components in (0,1)}
\item{dcop}{name of function of bivariate copula density (common for all variables), dcop accepts input of form  of d-vector or dxm matrix}
\item{param}{d-vector or mxd matrix, parameter of dcop }
}
\value{
 integrand for 1-factor copula density}
