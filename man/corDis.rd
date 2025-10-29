\name{corDis}
\alias{corDis}
\title{Discrepancy of model-based and observed correlation matrices based on Gaussian log-likelihood}
\description{
Discrepancy of model-based and observed correlation matrices 
}

\usage{
corDis(Rmodel,Rdata,n=0,npar=0)
}
\arguments{
\item{Rmodel}{model-based correlation matrix}
\item{Rdata}{empirical correlation matrix (could be observed or polychoric)}
\item{n}{sample size (if positive integer)}
\item{npar}{#parameters in the correlation structure}
}
\value{
 vector with discrepancy Dfit, and also nllk2 (wice negative log-likelihood), BIC, AIC if n and npar are inputted
}
\examples{
Rmodel = matrix(c(1,.3,.4,.4,.3,1,.5,.6,.4,.5,1,.7,.4,.6,.7,1),4,4)
print(Rmodel); print(chol(Rmodel))
Rdata = matrix(c(1,.32,.38,.41,.32,1,.53,.61,.38,.53,1,.67,.41,.61,.67,1),4,4)
print(corDis(Rmodel,Rdata))
print(corDis(Rmodel,Rdata,n=400,npar=3))
}
