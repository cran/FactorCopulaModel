\name{zetaPlot}
\alias{zetaPlot}
\title{Plot zeta(alpha) against alpha}
\description{
Plot zeta(alpha) against alpha 
}

\usage{
zetaPlot(dat,alpha,ylim=c(0,1),inverse=FALSE)
}
\arguments{
\item{dat}{nx2 data matrix with values (u-data in (0,1))}
\item{alpha}{vector of alpha>0 for zeta measure}
\item{ylim}{limits for yaxis to pass to plot}
\item{inverse}{if TRUE, plot zeta against 1/alpha}
}
\value{
 nothing is returned, but a plot is produced}
