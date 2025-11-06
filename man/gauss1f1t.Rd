\name{gauss1f1t}
\alias{gauss1f1t}
\title{Compute correlation matrix according to 1-factor + 1-truncated vine (residual dependence) model }
\description{
Compute correlation matrix according to a 1-factor + 1-truncated vine (for residual dependence) model 
}

\usage{
gauss1f1t(cormat, start_loading, iter=10, est="mle", plots=TRUE, trace=TRUE)
}
\arguments{
\item{cormat}{dxd correlation matrix}
\item{start_loading}{dx1 loading vector (for latent factor)}
\item{iter}{number of iterations for modified EM}
\item{est}{"mle" or "mom"}
\item{plots}{flag that is TRUE to show plots of EM steps}
\item{trace}{flag that is TRUE to print every 100th integer for iter}
}
\value{
 components:loading = final estimate for loading vector; 
R = correlation matrix (from MLE for 1F1T structure); 
Psi = vector of residual variances; 
loadings = matrix where ith row has the ith iteration; 
Rmats = list of correlation matrices; Rmats[[i]] has the ith iteration; 
Rstart = starting value of R based on starting values; 
dists = vector of distance measures as GOF criterion, ith entry for ith iteration; 
incls = iter x d*(d-1)/2 matrix, ith row for ith iteration
columns are whether edge 12, 13, 23, 14, .... (d-1,d) are in residual tree; 
partcor = dxd matrix of partial correlations given latent variable;  
loglik = vector of loglik values, ith entry for ith iteration.
}
\examples{
\donttest{
library(igraph)
data(DJ20142016gf)
zdat = dj1416gf$zscore # GARCH-filtered returns that have been transformed to N(0,1)
rzmat = cor(zdat)
d = ncol(zdat)
cat("\n1-factor start for 1f1t\n")
fa = factanal(factors=1,covmat=rzmat)
start = c(fa$loading)
# fitting 1Factor 1Truncated vine residual dependence structure
out1f1t = gauss1f1t(rzmat,start_loading=start,iter=20,plots=FALSE)
print(out1f1t$loading)
cat("\nedges for tree of residual dependence\n")
i = 1:d
nn = i*(i-1)/2
niter = 21  # above iteration bound +1
incls = out1f1t$incls[niter,]
for(j in 2:d)
{ for(k in 1:(j-1)) 
  { if(incls[nn[j-1]+k]) 
    { cat("edge ", k,j," "); cat(out1f1t$partcor[k,j],"\n") } 
  }
}
# extract three columns of this output to use with factor1trvine_nllk
# See example in cop1f1t()
}
}
\details{
A modified EM algorithm is used
-- first step of the M-step performed either by MLE or by method of moments
-- the second step assumes the moment estimator has been used
}
\references{
Brechmann EC and Joe H (2014). 
Parsimonious parameterization of correlation matrices using truncated
vines and factor analysis.
Computational Statistics and Data Analysis, 77, 233-251.
}
