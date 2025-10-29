\name{gaussLegendre}
\alias{gaussLegendre}
\title{R interface for Gauss-Legendre quadrature}
\description{Gauss-Legendre quadrature nodes and weights
}

\usage{
gaussLegendre(nq)
}
\arguments{
\item{nq}{number of quadrature points}
}
\value{
 structures with}
\examples{
out = gaussLegendre(15)
# same as statmod::gauss.quad.prob(15,dist="uniform") in library(statmod)
print(sum(out$weights)) # should be 1
print(sum(out$weights*out$nodes)) # should be 0.5  = E(U), U~Uniform(0,1)
print(sum(out$weights*out$nodes^2)) # should be 1/3 = E(U^2)
}
\details{
links to C code translation of jacobi.f in Stroud and Secrest (1966)
}
