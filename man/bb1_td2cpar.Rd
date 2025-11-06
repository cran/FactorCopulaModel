\name{bb1_td2cpar}
\alias{bb1_td2cpar}
\title{BB1 tail dependence parameters to copula parameter (theta,delta) }
\description{
BB1 map (lower,upper) tail dependence to copula parameter vector
}

\usage{
bb1_td2cpar(taildep)  
}
\arguments{
\item{taildep}{tail dependence parameter in mx2 matrix, by row (ltd,utd) in (0,1)^2 }
}
\value{
 matrix of copula parameters, by row (theta,delta), theta>0, delta>1}
\examples{
cpar = bb1_td2cpar(c(0.4,0.6))
print(cpar)
#         theta    delta
#[1,] 0.3672112 2.060043
print(bb1_cpar2td(cpar))
}
