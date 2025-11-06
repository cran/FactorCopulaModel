\name{semiCor}
\alias{semiCor}
\title{Semi-correlations for two variables}
\description{
semi-correlations (lower and upper) applied to data after normal scores transform
}

\usage{
semiCor(bivdat,inscore=FALSE)
}
\arguments{
\item{bivdat}{nx2 data set}
\item{inscore}{TRUE if bivdat has already been converted to normal scores, default FALSE}
}
\value{
 3-vector with rhoN = correlation of normal scores (vander Waerden correlation) and lower/ upper semi-correlations}
\examples{
# See example in semiCorTable()
}
