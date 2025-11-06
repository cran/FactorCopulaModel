\name{cparBounds}
\alias{cparBounds}
\title{lower and upper bounds for copula parameters (1-parameter, 2-parameter families) }
\description{
lower and upper bounds for copula parameters for use in min negative log-likelihood 
}

\usage{
cparBounds(familyvec)
}
\arguments{
\item{familyvec}{vector of family codes linking copula families}
}
\value{
 lower bound LB1/LB2 and upper bound LB2/UB2 for par1 and par2}
\examples{
famvec = c(1,4,5,14,2,7,17,10,20, 4,5,1)
out = cparBounds(famvec)
print(out)
print(out$LB1)
}
