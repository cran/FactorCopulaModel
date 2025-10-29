\name{frank_beta2cpar}
\alias{frank_beta2cpar}
\title{Frank: Blomqvist's beta to copula parameter}
\description{
Frank: Blomqvist's beta to copula parameter, vectorized
}

\usage{
frank_beta2cpar(beta, cpar0=0,mxiter=20,eps=1.e-8,iprint=FALSE)
}
\arguments{
\item{beta}{vector of Blomqvist's beta values, -1<beta<1}
\item{cpar0}{starting point for Newton-Raphson iterations}
\item{mxiter}{maximum number of iterations, default 20}
\item{eps}{tolerance for convergence, default 1.e-8}
\item{iprint}{print flag for iterations, default FALSE}
}
\value{
 vector of Frank copula parameters with the given betas}
\examples{
b = seq(-0.2,0.5,0.1) 
frank_beta2cpar(b)
frank_beta2cpar(b,iprint=TRUE)
}
\details{
Solve equation to get cpar given Blomqvist's beta, Newton-Raphson iterations;
vectorized input beta is OK, beta=0 fails
}
