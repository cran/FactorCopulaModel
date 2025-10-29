\name{corvec2mat}
\alias{corvec2mat}
\title{Convert from correlations in vector form to a correlation matrix}
\description{
Convert from correlations in vector form to a correlation matrix
}

\usage{
corvec2mat(rvec)
}
\arguments{
\item{rvec}{correlations in vector form of length d*(d-1)/2 in the order r12,r13,r23,r14,... r[d-1,d]}
}
\value{
 dxd correlation matrix}
\examples{
rvec = c(0.3,0.4,0.5,0.4,0.6,0.7)
Rmat = corvec2mat(rvec)
print(Rmat) # column 1 has 1, 0.3, 0.4, 0.4
}
