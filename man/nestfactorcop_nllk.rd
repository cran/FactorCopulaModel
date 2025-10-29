\name{nestfactorcop_nllk}
\alias{nestfactorcop_nllk}
\title{negative log-likelihoods of nested factor structured factor copula and derivatives computed in f90  for input to posDefHessMinb 
}
\description{
negative log-likelihoods of nested factor structured factor copula and derivatives
}

\usage{
nestfactorcop_nllk(param,dstruct,iprfn=FALSE)
}
\arguments{
\item{param}{parameter vector; parameters for copulas linking U_ij and V_j go *at the end* (i's with j=1 then j=2 etc) parameters for copulas linking V_j and V_0 go *first* (j=1,2 etc). For BB1 linking copulas for the global latent, the order is theta1,..,theta[d],delta1, ...,delta[d]; V_0 is the global/common latent variable that loads on all variables; V_j is a latent variable that loads only for variables in group j (by group g=1,2,..,mgrp etc).}
\item{dstruct}{list with  data set $data, copula name $copname, $quad is a Gauss-Legendre quadrature object, $repar is a flag for reparametrization (for Gumbel, BB1), $nu is a scalar or 2-vector for degree of freedom parameter(s), $grsize is a vector with group sizes; if dstruct$pdf == 1 the function evaluates nllk only  (and returns zero gradient and hessian). Options for copname are: frank, gumbel, frankgumbel, frankbb1, gumbelbb1, tbb1, t. For tbbb1, nu is a scalar. For t, nu is a 2-vector.}
\item{iprfn}{indicator for printing of function and gradient (within Newton-Raphson iterations)}
}
\value{
  nllk, grad, hess (gradient and hessian included)}
\examples{
\donttest{
grsize = c(4,4,3)
d = sum(grsize)
n = 500
mgrp = length(grsize)
set.seed(222)
par_nest = c(rep(1.7,3),seq(1.7,3.7,0.2))
udat_obj = rnestfactor(n,grsize,cop=4,par_nest)
udat = udat_obj$data
summary(udat_obj$v0)
summary(udat_obj$vg)
zdat = qnorm(udat)
rmat  = cor(zdat)
round(cor(zdat),3)
# run oblique_fa to get oblqiue factor correlation matrix
obfa = oblique_fa(grsize,start=c(rep(0.8,d),rep(0.5,mgrp)),cormat=rmat,n=n,prlevel=0)
loading1 = rowSums(obfa$loadings)
corlat = obfa$cor_lat # correlations of group latent variables
fa1 = factanal(covmat=corlat,factors=1)
loadlat = c(fa1$loadings)
print(loadlat)  
# starting values for different cases
# convert loading/latcor to Frank, Gumbel and BB1 parameters etc
start_frk1 = frank_rhoS2cpar(loading1)
start_frk0 = frank_rhoS2cpar(loadlat)
start_frk = c(start_frk0,start_frk1)
start_gum1 = gumbel_rhoS2cpar(loading1)
start_gum0 = gumbel_rhoS2cpar(loadlat)
start_gum = c(start_gum0,start_gum1)
start_frkgum = c(start_frk0,start_gum1)
start_tnu = c(loadlat,loading1)
tau = bvn_cpar2tau(loading1)
start_bb1 = bb1_tau2eqtd(tau)
start_bb1 = c(start_bb1[,1:2]) # all thetas and then all deltas
start_frkbb1 = c(start_frk0,start_bb1)
start_gumbb1 = c(start_gum0,start_bb1)
start_tnubb1 = c(loadlat,start_bb1)
#
gl = gaussLegendre(25)
npar = mgrp+d
dstrfrk = list(data=udat,copname="frank",quad=gl,repar=0,grsize=grsize)
out = nestfactorcop_nllk(start_frk, dstrfrk)
print(out$fnval)
print(out$grad)
ml_frk = posDefHessMinb(rep(3,npar),nestfactorcop_nllk, ifixed=rep(FALSE,npar), 
dstrfrk, LB=rep(-20,npar), UB=rep(30,npar), mxiter=30, eps=5.e-5, iprint=TRUE)
dstrgum = list(data=udat,copname="gumbel",quad=gl,repar=0,grsize=grsize)
ml_gum = posDefHessMinb(start_gum,nestfactorcop_nllk, ifixed=rep(FALSE,npar), 
dstrgum, LB=rep(1,npar), UB=rep(20,npar), mxiter=30, eps=5.e-5, iprint=TRUE)
dstrfrkgum = list(data=udat,copname="frankgumbel",quad=gl,repar=0,grsize=grsize)
ml_frkgum = posDefHessMinb(start_frkgum,nestfactorcop_nllk, ifixed=rep(FALSE,npar), 
dstrfrkgum, LB=c(rep(-20,mgrp),rep(1,d)), UB=rep(25,npar), mxiter=30, 
eps=5.e-5, iprint=TRUE)
dstrtnu = list(data=udat,copname="t",quad=gl,repar=0,grsize=grsize, nu=c(10,20))
ml_tnu = posDefHessMinb(start_tnu,nestfactorcop_nllk, ifixed=rep(FALSE,npar), 
dstrtnu, LB=c(rep(-1,npar)), UB=rep(1,npar), mxiter=30, eps=5.e-5, iprint=TRUE)
# diverges with parameters approaching 1
#
npar2 = mgrp+2*d
dstrfrkbb1 = list(data=udat,copname="frankbb1",quad=gl,repar=0,grsize=grsize)
out = nestfactorcop_nllk(start_frkbb1,dstrfrkbb1)
print(out$fnval)
print(out$grad)
ml_frkbb1 = posDefHessMinb(start_frkbb1,nestfactorcop_nllk, ifixed=rep(FALSE,npar2), 
dstrfrkbb1, LB=c(rep(-20,mgrp),rep(0,d),rep(1,d)), UB=rep(20,npar2), 
mxiter=30, eps=5.e-5, iprint=TRUE)
dstrgumbb1 = list(data=udat,copname="gumbelbb1",quad=gl,repar=0,grsize=grsize)
ml_gumbb1 = posDefHessMinb(start_gumbb1,nestfactorcop_nllk, ifixed=rep(FALSE,npar2), 
dstrgumbb1, LB=c(rep(1,mgrp),rep(0,d),rep(1,d)), UB=rep(20,npar2),
mxiter=30, eps=5.e-5, iprint=TRUE)
dstrtnubb1 = list(data=udat,copname="tbb1",quad=gl,repar=0,grsize=grsize, nu=20)
ml_tnubb1 = posDefHessMinb(start_tnubb1,nestfactorcop_nllk, ifixed=rep(FALSE,npar2), 
dstrtnubb1, LB=c(rep(-1,mgrp),rep(0,d),rep(1,d)), UB=c(rep(1,mgrp),rep(20,2*d)), 
mxiter=30, eps=5.e-5, iprint=TRUE)
#
# compare nllk and number of iterations
cat(ml_frk$fnval, ml_gum$fnval, ml_frkgum$fnval, ml_tnu$fnval,
ml_frkbb1$fnval, ml_gumbb1$fnval, ml_tnubb1$fnval, "\n")
# -1438.187 -1760.851 -1725.286 -5555.964 -1729.274 -1764.83 -1746.629
cat(ml_frk$iter, ml_gum$iter, ml_frkgum$iter, ml_tnu$iter,
ml_frkbb1$iter, ml_gumbb1$iter, ml_tnubb1$iter, "\n")
# 7 8 6 23 15 16 16 
# nested-factor t(10)/t(20) failed because some parameters approached the
# upper bound of 1 in which case the numerical integration is inaccurate.
}
}
