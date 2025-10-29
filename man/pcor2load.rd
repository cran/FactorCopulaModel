\name{pcor2load}
\alias{pcor2load}
\title{Partial correlation representation to loadings for p-factor}
\description{
Partial correlation representation to loadings for p-factor
}

\usage{
pcor2load(rhomat)
}
\arguments{
\item{rhomat}{dxp matrix with correlations for factor 1 in column 1, and partial correlations with factor k given previous factors in column k}
}
\value{
 loading matrix}
\examples{
grsize = c(5,4,3) 
# bi-factor parameters: 13 correlations with global latent and then 
# 5 partial correlations for group1 latent given global,
# 4 partial correlations for group2 latent given global,
# 3 partial correlations for group3 latent given global
par_bifact = c(0.84,0.63,0.58,0.78,0.79,  0.87,0.80,0.74,0.71,  0.83,0.77,0.80,
0.67,0.58,0.15,0.70,0.47,  0.32,0.27,0.73,0.19,  0.35,0.23,0.53)
mgrp = length(grsize)
d = sum(grsize)
pcmat = matrix(0,d,mgrp+1) # bi-factor structure has mgrp+1 factors
pcmat[,1] = par_bifact[1:d]
iend = cumsum(grsize)
ibeg = iend+1; ibeg = c(1,1+iend[-mgrp])
for(g in 1:mgrp)
{ pcmat[ibeg[g]:iend[g],g+1] = par_bifact[d+ibeg[g]:iend[g]] }
print(pcmat)
aload = pcor2load(pcmat)
print(aload)
}
\details{
Partial correlation representation to loadings for p-factor
}
