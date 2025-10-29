\name{isPosDef}
\alias{isPosDef}
\title{Check if a square symmetric matrix is positive definite}
\description{
Check if a square symmetric matrix is positive definite
}

\usage{
isPosDef(amat)
}
\arguments{
\item{amat}{symmetric matrix}
}
\value{
 TRUE if amat is positive definite, FALSE otherwise}
\examples{
a1 = matrix(c(1,.5,.5,1),2,2)
a2 = matrix(c(1,1.5,1.5,1),2,2)
t1 = try(chol(a1))
t2 = try(chol(a2))
print(isPosDef(a1))
print(isPosDef(a2))
}
