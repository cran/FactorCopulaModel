\name{gumbel_rhoS2cpar}
\alias{gumbel_rhoS2cpar}
\title{Gumbel: Spearman rho to copula parameter}
\description{
Gumbel: Spearman rho to copula parameter
}

\usage{
gumbel_rhoS2cpar(rho)
}
\arguments{
\item{rho}{vector of Spearman values, 0<rho<1}
}
\value{
 vector of Gumbel copula parameters with the given rho}
\examples{
rho = seq(0.1,0.5,0.1)
gumbel_rhoS2cpar(rho)
}
