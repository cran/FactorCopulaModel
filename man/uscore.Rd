\name{uscore}
\alias{uscore}
\title{Rank-based uniform scores transform}
\description{
Rank-based uniform scores transform
}

\usage{
uscore(data,aunif=-0.5)
}
\arguments{
\item{data}{dataframe or matrix, or vector, of reals}
\item{aunif}{adjustment 'a' for (rank+a)/(n+1+2*a) as scores. n=sample size , default is -0,5}
}
\value{
 matrix or vector of uniform scores}
