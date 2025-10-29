\name{gumbel_beta2cpar}
\alias{gumbel_beta2cpar}
\title{Gumbel: Blomqvist's beta to copula parameter}
\description{
Gumbel: Blomqvist's beta to copula parameter, vectorized
}

\usage{
gumbel_beta2cpar(beta)
}
\arguments{
\item{beta}{vector of Blomqvist's beta values, 0<beta<1}
}
\value{
 vector of Gumbel copula parameters with the given betas}
\examples{
b = seq(0.1,0.5,0.1)
gumbel_beta2cpar(b)
}
