\name{onefactorcop_nllk}
\alias{onefactorcop_nllk}
\title{negative log-likelihood of 1-factor copula for input to posDefHessMin and posDefHessMinb }
\description{
negative log-likelihood (nllk) of 1-factor copula for input to posDefHessMin and posDefHessMinb
}

\usage{
onefactorcop_nllk(param,dstruct,iprfn=FALSE)
}
\arguments{
\item{param}{parameter vector}
\item{dstruct}{list with  data set $data, copula name $copname, $quad is list with quadrature weights and nodes,  $repar is code for reparametrization (for Gumbel, BB1), $nu is positive degree of freedom parameter (for t)  (linking copula is common for all variables). Options for copname are: frank, gumbel, bb1, t. For reflected gumbel or bb1, use something like dstruct$dat = 1-udata}
\item{iprfn}{print flag for function and gradient (within Newton-Raphson) iterations) for BB1, param is 2*d-vector with th1,de1,th2,de2,...}
}
\value{
 nllk, grad (gradient), hess (hessian) at MLE}
\examples{
\donttest{
cpar_gum = seq(1.9,3.7,0.2)
d = length(cpar_gum)
n = 300
param = c(rbind(cpar_gum,rep(0,d)))  # second par2 is 0 for VineCopula
set.seed(111)
gum_obj = r1factor(n,d,param,famvec=rep(4,d)) # uses VineCopula
udat = gum_obj$udata
zdat = qnorm(udat)
rmat = cor(zdat)
print(round(rmat,3))
# run factanal to get loading (rho in normal scale close to Spearman rho)
fa1 = factanal(covmat=rmat,factors=1)
loadings = c(fa1$loading)
# convert loadings to Frank, Gumbel and BB1 parameters 
start_frk = frank_rhoS2cpar(loadings)
start_gum = gumbel_rhoS2cpar(loadings)
tau = bvn_cpar2tau(loadings)
start_bb1 = bb1_tau2eqtd(tau)
start_bb1 = c(t(start_bb1[,1:2]))
gl = gaussLegendre(25)
dstrfrk1 = list(copname="frank",data=udat,quad=gl,repar=0, pdf=1)
dstrfrk = list(copname="frank",data=udat,quad=gl,repar=0, pdf=0)
obj1 = onefactorcop_nllk(start_frk,dstrfrk1) #nllk only
obj = onefactorcop_nllk(start_frk,dstrfrk) # nllk, grad, hess
print(obj1$fnval)
print(obj$grad)
#
ml_frk = posDefHessMinb(start_frk,onefactorcop_nllk,ifixed=rep(FALSE,d),
  dstruct=dstrfrk, LB=rep(-30,d),UB=rep(30,d),iprint=TRUE,eps=1.e-5)
dstrgum = list(copname="gumbel",data=udat,quad=gl,repar=0)
ml_gum = posDefHessMinb(start_gum,onefactorcop_nllk,ifixed=rep(FALSE,d),
  dstruct=dstrgum, LB=rep(-30,d),UB=rep(30,d),iprint=TRUE,eps=1.e-5)
dstrgumr = list(copname="gumbel",data=1-udat,quad=gl,repar=0)
ml_gumr = posDefHessMinb(start_gum,onefactorcop_nllk,ifixed=rep(FALSE,d),
  dstruct=dstrgumr, LB=rep(-30,d),UB=rep(30,d),iprint=TRUE,eps=1.e-5)
dstrtnu = list(copname="t",data=udat,quad=gl,repar=0,nu=10)
ml_tnu = posDefHessMinb(loadings,onefactorcop_nllk,ifixed=rep(FALSE,d),
  dstruct=dstrtnu, LB=rep(-1,d),UB=rep(1,d),iprint=TRUE,eps=1.e-5)
dstrbb1 = list(copname="bb1",data=udat,quad=gl,repar=0)
ml_bb1 = posDefHessMinb(start_bb1,onefactorcop_nllk,ifixed=rep(FALSE,2*d),
  dstruct=dstrbb1, LB=rep(c(0,1),d),UB=rep(20,2*d),iprint=TRUE,eps=1.e-5)
dstrbb1r = list(copname="bb1",data=1-udat,quad=gl,repar=0)
ml_bb1r = posDefHessMinb(start_bb1,onefactorcop_nllk,ifixed=rep(FALSE,2*d),
  dstruct=dstrbb1r, LB=rep(c(0,1),d),UB=rep(20,2*d),iprint=TRUE,eps=1.e-5)
cat(ml_frk$fnval, ml_gum$fnval, ml_gumr$fnval, ml_tnu$fnval, ml_bb1$fnval, ml_bb1r$fnval, "\n")
# -1342.936 -1560.391 -1198.837 -1454.907 -1560.45 -1549.935
cat(ml_frk$iter, ml_gum$iter, ml_gumr$iter, ml_tnu$iter, ml_bb1$iter, ml_bb1r$iter, "\n")
# 4 4 5 4 16 6
}
}
\details{
linked to Fortran 90 code for speed
}
