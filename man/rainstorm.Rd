\name{rainstorm}
\alias{rainstorm}
\title{Precipitation by rainstorm at 28 stations}
\usage{
data(rainstorm)
}
\description{
R workspace file with
transformed precipitation by rainstorm at 28 stations

There are 4 components:

1. $uprecip: 256x28 matrix with total precipitation in mm by storm 
at 28 stations, after rank transform to U(0,1);

2. $zprecip: 256x28 matrix with total precipitation in mm by storm 
at 28 stations, after rank transform to N(0,1);

3. $cormat: the correlation matrix of $zprecip;

4. $grsize: vector (6,12,10) for sizes of 3 goups of stations
found by a variable clustering method.
}

\keyword{datasets}
