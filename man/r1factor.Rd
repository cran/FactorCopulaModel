\name{r1factor}
\alias{r1factor}
\title{simulate from 1-factor copula model with different linking copula families}
\description{
simulate from 1-factor copula model and include corresponding latent variables
}

\usage{
r1factor(n,d,param,famvec,vlatent=NULL)
}
\arguments{
\item{n}{sample size}
\item{d}{dimension}
\item{param}{copula parameter vector of dimension 2*d (par[j],par2[j], j=1,...,d) for d linking copulas, for one-parameter families set par2=0}
\item{famvec}{family vector for d linking copulas same index with VineCopula package, for a range of tail properties, select  1:Gaussian; 2:t; 4:Gumbel; 14:survival Gumbel;  5:Frank;   7:BB1; 17:survival BB1;   10:BB8; 20:survival BB8;  }
\item{vlatent}{given n-vector of latent variavles in U(0,1);  use for simulating from a group for oblique factor copula (in this case, apply this function G times for G groups, extracting dependent latent variables from a nxG matrix)}
}
\value{
 list with udata: nxd matrix in (0,1) andvlatent n-vector of corresponding latent variables in (0,1).
}
\examples{
 # Example 1
cpar_frk = c(12.2,3.45,4.47,4.47,5.82)
d = length(cpar_frk)
cpar2_frk = rep(0,d)
param = c(rbind(cpar_frk,cpar2_frk))
n = 300
set.seed(123)
frk_obj = r1factor(n,d,param,famvec=rep(5,d)) # uses VineCopula
frkdat = frk_obj$udata
print(cor(frkdat))
print(summary(frk_obj$vlatent))
dfrank = function(u,v,cpar)
{ t1 = 1.-exp(-cpar); tem1 = exp(-cpar*u); tem2 = exp(-cpar*v);
  pdf = cpar*tem1*tem2*t1; tem = t1-(1.-tem1)*(1.-tem2);
  pdf = pdf/(tem*tem);
  pdf
}
cat("\nFrank 1-factor MLE: standalone R\n")
out1_frk = ml1factor(nq=21,cpar_frk,frkdat,dfrank,LB=-30,UB=30,prlevel=1,mxiter=100)
cat("\nFrank 1-factor MLE: nlm with f90 code\n")
out2_frk = ml1factor_f90(nq=21,cpar_frk,frkdat,copname="frank",LB=-30,UB=30,prlevel=1,mxiter=100)
cat("\nFrank 1-factor MLE: posDefhessMinb with f90 code\n")
gl21 = gaussLegendre(21)
dstrfrk = list(copname="frank",data=frkdat,quad=gl21,repar=0)
out3_frk = posDefHessMinb(cpar_frk,onefactorcop_nllk,ifixed=rep(FALSE,d),
dstruct=dstrfrk, LB=rep(-30,d),UB=rep(30,d),iprint=TRUE,eps=1.e-5)
cat(out1_frk$minimum, out2_frk$minimum, out3_frk$fnval,"\n")
print(cbind(out1_frk$estimate, out2_frk$estimate, out3_frk$parmin))
print(sqrt(diag(out3_frk$invh))) # SEs
#
# Example 2 (oblique factor with 3 groups)
n = 500
d = 10
ltd1 = c(0.3,0.4,0.6,0.7,0.6); utd1 = c(0.5,0.6,0.7,0.5,0.4) 
ltd2 = c(0.3,0.4,0.6,0.5,0.4); utd2 = c(0.6,0.3,0.5,0.7,0.6)
ltd3 = c(0.5,0.4,0.5,0.5); utd3 = c(0.5,0.4,0.5,0.5)
grsize = c(5,5,4)
rmat = toeplitz(c(1,0.5,0.5))  # for Gaussian copula parameter of  latent
cpar1 = bb1_td2cpar(cbind(ltd1,utd1))
cpar2 = bb1_td2cpar(cbind(ltd2,utd2))
cpar3 = bb1_td2cpar(cbind(ltd3,utd3))
set.seed(205)
zmat = rmvn(n,rmat); vmat = pnorm(zmat)
# Any vine copula could be used for latent variables, besides multivariate normal
data1g = r1factor(n=500,grsize[1],param=c(t(cpar1)),
famvec=rep(7,grsize[1]), vlatent=vmat[,1])
data2g = r1factor(n=500,grsize[2],param=c(t(cpar2)),
famvec=rep(7,grsize[2]), vlatent=vmat[,2])
data3g = r1factor(n=500,grsize[3],param=c(t(cpar3)),
famvec=rep(7,grsize[3]), vlatent=vmat[,3])
udata_oblf = cbind(data1g$udata,data2g$udata,data3g$udata)
rr_oblf = cor(udata_oblf,method="spearman")
print(round(rr_oblf,2))
}
