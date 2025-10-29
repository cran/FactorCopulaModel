\name{factor1trvine_nllk}
\alias{factor1trvine_nllk}
\title{negative log-likelihood with gradient and Hessian computed in f90 for copula from 1-factor/1-truncated vine (tree for residual dependence 
   conditional on a latent variable); models included are BB1 for latent with Frank or Gaussian(bvncop) for truncated vine residual dependence 
}
\description{
negative log-likelihood for 1-factor with tree for residual dependence
}

\usage{
factor1trvine_nllk(param,dstruct,iprint=FALSE)
}
\arguments{
\item{param}{parameter vector (length 2*d+(d-1)); BB1 parameters for links of observed variable U_j to latent V, j=1,...,d;     Frank or Gaussian parameters for edge_k=(node_k[1],node_k[2]) observed variables nodek[1],nodek[2] given V, k=1,...,d-1, where the d-1 edges form a tree. 2*d BB1 parameters (theta_1,delta_1,...,theta_d,delta_d); d-1 Frank or BVN parameters for residual dependence.}
\item{dstruct}{list with nxd data set $data on (0,1)^d ; $quad is list with quadrature weights and nodes;  $edg1, $edg2 (d-1)-vectors for node labels of edges 1,...,d-1; for the tree for residual dependence; $fam for copula family label, where $fam=1 for Frank copula for residual dependence conditional on latent; $fam=2 for bivariate Gaussian/normal copula.}
\item{iprint}{indicator for printing of function and gradient}
}
\value{
 nllk, grad, hess (negative log-likelihood, gradient, Hessian at MLE)}
\references{
Joe (2018). Parsimonious graphical dependence models constructed from vines.
Canadian Journal of Statistics 46(4), 532-555.
}
\examples{
\donttest{
data(DJ20142016gf)
udat = dj1416gf$uscore
d = ncol(udat)
loadings = c(0.53,0.68,0.68,0.62,0.69,0.58,0.60,0.67,0.73,0.72,
0.64,0.70,0.65,0.69,0.75,0.56,0.56,0.79,0.59,0.66,
0.57,0.60,0.60,0.71,0.57,0.73,0.70,0.56, 0.52,0.61) 
# starting values for BB1 that have roughly dependence of bivariate Gaussian
tau = bvn_cpar2tau(loadings)
par0 = bb1_tau2eqtd(tau)
par0 = par0[,c(1,2)]
par0=c(t(par0))
# from gauss1f1t()
edg1 = c(4,4,4,2,5,10,10,16,9,14,1,3,12,13,8,11,
17,19,4,10,16,16,22,3,4,21,16,23,6)
edg2 = c(6,7,9,10,13,14,15,17,18,19,20,20,20,20,21,21,
21,22,23,23,23,24,25,26,26,27,28,29,30)
pc = c(0.2986594,0.1660604,0.213082,0.2975893,0.2534793,-0.1485211,
0.6179934,0.2207242,0.1877172,0.2625087,0.1414915,-0.1171917,
0.1396517,0.2632831,0.1964068,0.2850865,0.1679997,0.3683564,
-0.1666632,-0.2291848,0.3600637,0.1737616,0.2022859,0.2136862,
0.1948507,0.1731044,0.186075,0.2177455,0.6606208)
# convert above 29 rho parameters from above to Frank parameters with roughly
# the same Spearman rhos, e.g. frank_rhoS2cpar()
print(frank_rhoS2cpar(pc))
parfrk = c(1.87, 1.00, 1.30, 1.86, 1.57, -0.90, 4.67, 1.35, 1.14, 1.63, 
0.85, -0.70, 0.84, 1.63, 1.20, 1.78, 1.02, 2.37, -1.01, -1.41, 
2.30, 1.05, 1.23, 1.31, 1.19, 1.05, 1.13, 1.33, 5.23)
# parameters for tree 1 (latent) and tree 2 (residual dependence)
gl = gaussLegendre(21)
bffxvar = rep(FALSE,d-1) 
ifixed = c(rep(FALSE,2*d),bffxvar)
LB = c(rep(c(0.01,1.01),d),rep(-15,d-1)) 
UB = c(rep(c(10,10),d),rep(20,d-1))
st_bb1frk = c(par0,parfrk)
dstrbb1frk = list(data=udat,quad=gl, edg1=edg1, edg2=edg2, fam=1)
tem_frk = factor1trvine_nllk(st_bb1frk,dstrbb1frk)
print(tem_frk$fnval)
# -6600.042
ml_bb1frk = posDefHessMinb(st_bb1frk, factor1trvine_nllk, ifixed= ifixed, 
dstrbb1frk, LB, UB, mxiter=30, eps=5.e-5, bdd=5, iprint=TRUE)
#1 -6600.042 0.537193 
#2 -6645.547 0.06360556 
#3 -6645.968 0.0008600067 
#4 -6645.968 6.065411e-07
cat("BB1frk: parmin (mle for udata), SEs\n")
print(ml_bb1frk$fnval)
#  -6645.968
print(ml_bb1frk$parmin)
print(sqrt(diag(ml_bb1frk$invh)))
#
dstrbb1gau = list(data=udat,quad=gl, edg1=edg1, edg2=edg2, fam=2)
LB = c(rep(c(0.01,1.01),d),rep(-1,d-1)) 
UB = c(rep(c(10,10),d),rep(1,d-1))
st_bb1gau = c(par0,pc)
tem_gau = factor1trvine_nllk(st_bb1gau,dstrbb1gau)
print(tem_gau$fnval)
# -6587.514
ml_bb1gau = posDefHessMinb(st_bb1gau, factor1trvine_nllk, ifixed= ifixed, 
dstrbb1gau, LB, UB, mxiter=30, eps=5.e-5, bdd=5, iprint=TRUE)
#1 -6587.514 0.1732503 
#2 -6621.988 0.03526498 
#3 -6622.375 0.000806789 
#4 -6622.375 7.05569e-07 
cat("BB1gau: parmin (mle for udata), SEs\n")
print(ml_bb1gau$fnval)
# -6622.375
print(ml_bb1gau$parmin)
print(sqrt(diag(ml_bb1gau$invh)))
}
}