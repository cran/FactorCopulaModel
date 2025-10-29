\name{qcondFrank}
\alias{qcondFrank}
\title{C_[2|1]^[-1](p|u) for bivariate Frank copula}
\description{
Frank bivariate copula conditional quantile
}

\usage{
qcondFrank(p,u,cpar)
}
\arguments{
\item{p}{0<p<1, could be a vector}
\item{u}{0<u<1, could be a vector }
\item{cpar}{copula parameter: cpar>0 or cpar<0; cpar=0 input will not work}
}
\value{
 conditional quantiles of bivariate Frank copula}
\details{
1-exp(-cpar) becomes 1 in double precision for cpar>37.4;
any argument can be a vector, but all vectors must have same length.
Form of inputs not checked (for readability of code).
}
