\name{frank_rhoS2cpar}
\alias{frank_rhoS2cpar}
\title{Frank: Spearman rho to copula parameter}
\description{
Frank: Spearman rho to copula parameter
}

\usage{
frank_rhoS2cpar(rho)
}
\arguments{
\item{rho}{vector of Spearman values, -1<rho<1}
}
\value{
 vector of Frank copula parameters with the given rho}
\examples{
rho = seq(-0.2,0.6,0.1)
frank_rhoS2cpar(rho)
}
