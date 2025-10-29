\name{DJ20142016gf}
\alias{DJ20142016gf}
\alias{dateindex}
\alias{lab}
\alias{dj1416gf}
\title{GARCH-filtered log returns for Dow Jones stocks 2014-2016}
\usage{
data(DJ20142016gf)
}
\description{
R workspace file with
GARCH-filtered log returns for Dow Jones stocks 2014-2016

Three objects:

1. dj1416gf is a list with 
$filter, $uscore   $zscore   $uscmodel $zscmodel $sigmat $coefficient;
dimensions of matrices $uscore, $zscore are 764 x 30 (n x d).

filter: GARCH filtered data nxd, before transform

uscore: empirical uniform scores (nxd)

zscore: empirical normal scores (nxd)

uscmodel: model-based uniform scores (nxd) from the GARCH fit

zscmodel: model-based normal scores (nxd) from the GARCH fit

sigmat: matrix of estimated volatilities (nxd)

coef: matrix of GARCH parameters (6xd or 5xd depending on where AR(1) was
used for GARCH model;
the parameters are mu, (ar1), omega, alpha1, beta1, shape.


2. dateindex is an object with the dates for the rows of $uscore, $zscore

3. lab is an object with the ticker names for the 30 stocks in dj1416gf
}
\keyword{datasets}
