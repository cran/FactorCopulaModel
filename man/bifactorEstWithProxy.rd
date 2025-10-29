\name{bifactorEstWithProxy}
\alias{bifactorEstWithProxy}
\title{Sequential parameter estimation for bi-factor copula with estimated latent variables using VineCopula::BiCopSelect 
}
\description{
Sequential parameter estimation for bi-factor copula with estimated latent variables
}

\usage{
bifactorEstWithProxy(udata,vglobal,vgroup,grsize, 
  famset_global, famset_group, iprint=FALSE)

}
\arguments{
\item{udata}{nxd matrix with values in (0,1)}
\item{vglobal}{n-vector is estimated global latent variables (or test with known values)}
\item{vgroup}{n*mgrp matrix with estimated group-based latent variables}
\item{grsize}{G-vector with group sizes with mgrp=G=length(grsize)groups}
\item{famset_global}{codes for allowable copula families for d global linking copulas}
\item{famset_group}{codes for allowable copula families for d global linking copulas VineCopula: current choices to cover a range of tail behavior are: 1 = Gaussian/normal; 2 = t; 4 = Gumbel; 5 = Frank; 7 = BB1; 10 = BB8; 14 = survival Gumbel; 17 = survival BB1; 20 = survival BB8.}
\item{iprint}{if TRUE print intermediate results}
}
\value{
 list with fam = 2*d vector of family codes chosen via BiCopSelect;dx2 matrix global_par; dx2 matrix group_par;
these contain par,par2 for the selected copula families 
in the 2-truncated vine rooted at the latent variables.
}
\examples{
\donttest{
# BB1/Frank bi-factor copula
set.seed(2024)
th1_range = c(0.3,1)
th2_range = c(1.1,2.5)
th3_range = c(8.5,18.5)
grsize = rep(10,3)
mgrp = length(grsize)
d = sum(grsize)
parbi = c(runif(d,th1_range[1],th1_range[2]),  # BB1 theta
runif(d,th2_range[1],th2_range[2]),  # BB1 delta
runif(d,th3_range[1],th3_range[2]))  # Frank
n = 500
data = rbifactor(n,grsize=grsize,cop=7,parbi)
udata = data$data
vlat = cbind(data$v0,data$vg)
fam_true = c(rep(7,d),rep(5,d))
#
guess = c(rep(0.7,d),rep(0.5,d)) 
bif_obj = bifactorScore(udata, start=guess, grsize, prlev=1)
proxy_init = bif_obj$proxies
# selection of linking copula families and estimation of their parameters 
select1 = bifactorEstWithProxy(udata,proxy_init[,1],proxy_init[,-1], grsize, 
famset_global=c(1,4,5,7,10,17,20), famset_group=c(1,2,4,5))
fam1 = select1$fam
parglo1 = select1$global_par
pargrp1 = select1$group_par
print(fam1)  # 7,17,1 for global; 5 for group
print(parglo1)
print(pargrp1)
plot(parbi[(2*d+1):(3*d)],pargrp1[,1])
cor(parbi[(2*d+1):(3*d)],pargrp1[,1])
#
condExpProxy = latentUpdateBifactor(udata=udata, cparvec=c(parglo1,pargrp1),
grsize=grsize, family=fam1, nq=25)
proxy_improved = cbind(condExpProxy$v0, condExpProxy$vg)
par(mfrow=c(2,2))
for(j in 1:(mgrp+1))
{ plot(proxy_improved[,j],vlat[,j])
  print(cor(proxy_improved[,j],vlat[,j]))
}
rmse_values = sapply(1:ncol(proxy_improved),
function(i) sqrt(mean((proxy_improved[,i] - vlat[,i])^2)))
round(rmse_values,3)
#
# With improved proxies,
# selection of linking copula families and estimation of their parameters 
select2 = bifactorEstWithProxy(udata,proxy_improved[,1],proxy_improved[,-1],
grsize, famset_global=c(1,4,5,7,10,17,20), famset_group=c(1,2,4,5))
parglo2=select2$global_par
pargrp2=select2$group_par
fam2 = select2$fam  # 7,17 for global; 5 for local
cbind(fam_true,fam1,fam2)  
plot(parbi[(2*d+1):(3*d)],pargrp2[,1])
cor(parbi[(2*d+1):(3*d)],pargrp2[,1]) # higher correlation than pargrp1
}
}
\references{
1. Krupskii P and Joe H (2013).
Factor copula models for multivariate data.
Journal of Multivariate Analysis, 120, 85-101.
2. Fan X and Joe H (2024).
High-dimensional factor copula models with estimation of latent variables
Journal of Multivariate Analysis, 201, 105263.
}
\details{
It is best if variables have been oriented to be positively related
to the latent variable
}
