\name{bvnSemiCor}
\alias{bvnSemiCor}
\title{Semi-correlation for bivariate normal/Gaussian distribution}
\description{
semicorrelation assuming bivariate normal/Gaussian copula 
}

\usage{
bvnSemiCor(rho)
}
\arguments{
\item{rho}{correlation in (-1,1)}
}
\value{
 Cor(Z1,Z2| Z1>0,Z2>0) when (Z1,Z2)~bivariate standard normal(rho)}
\references{
 Joe (2014), Dependence Modeling with Copulas, Chapman&Hall/CRC; p 71
}
