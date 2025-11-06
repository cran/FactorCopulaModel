\name{rhoS}
\alias{rhoS}
\title{Spearman's rho for bivariate copula with parameter cpar}
\description{
Spearman's rho for bivariate copula with parameter cpar
}

\usage{
rhoS(cpar,cop, zero=0, icond=FALSE, tol=0.0001)
}
\arguments{
\item{cpar}{copula parameter}
\item{cop}{function name of joint cdf or conditional cdf C_2|1}
\item{zero}{0 or something like 1.e-6 (to avoid endpoint problems)}
\item{icond}{icond flag for using condition cdf of copula;
if icond = TRUE, cop = conditional cdf pcondcop 
if icond = FALSE, cop = joint cdf pcop 
  default is to integrate on [zero,1-zero]^2 for icond=TRUE.}
\item{tol}{accuracy for 2-dimensional integration}
}
\value{
 Spearman rho value for copula}
\examples{
# Bivariate margin of 1-factor copula
# using conditional cdf via VineCopula::BiCopHfunc2
# cpar1 = (par,par2) for fam1 for first variable
# cpar2 = (par,par2) for fam2 for second variable
# family codes 1:Gaussian, 2:t, 4:Gumbel, 5:Frank, 7:BB1, 14:survGumbel, 17:survBB1
p1factbiv = function(u1,u2,cpar1,cpar2,fam1,fam2,nq)
{ if(length(u1)==1) u1 = rep(u1,nq)
  if(length(u2)==1) u2 = rep(u2,nq)
  gl = gaussLegendre(nq)
  wl = gl$weights; vl = gl$nodes
  a1 = VineCopula::BiCopHfunc2(u1,vl,family=fam1,par=cpar1[1], par2=cpar1[2])
  a2 = VineCopula::BiCopHfunc2(u2,vl,family=fam2,par=cpar2[1], par2=cpar2[2])
  sum(wl*a1*a2)
}
#
# Spearman rho
rhoS_1factor = function(cpar1,cpar2,fam1,fam2,nq)
{ gl = gaussLegendre(nq)
  wl = gl$weights; xl = gl$nodes
  pcint1 = rep(0,nq); pcint2 = rep(0,nq)
  for(iq in 1:nq)
  { a1 = VineCopula::BiCopHfunc2(xl,rep(xl[iq],nq),family=fam1,par=cpar1[1],par2=cpar1[2])
    a2 = VineCopula::BiCopHfunc2(xl,rep(xl[iq],nq),family=fam2,par=cpar2[1],par2=cpar2[2])
    pcint1[iq] = sum(wl*a1)
    pcint2[iq] = sum(wl*a2)
  }
  tem = sum(wl*pcint1*pcint2)
  12*tem-3
}
#
# Tests for Spearman rho with BB1 linking copulas to latent variable
param = matrix(c(0.5,1.5,0.6,1.2),2,2,byrow=TRUE)
# Gauss-Legendre quadrature
rho_1factbb1_gl = rhoS_1factor(param[1,],param[2,],fam1=7,fam2=7,nq=21)
# reflected/survival BB1
rho_1factbb1r_gl = rhoS_1factor(param[1,],param[2,],fam1=17,fam2=17,nq=21)
cat(rho_1factbb1_gl, rho_1factbb1r_gl, "\n")
# 0.3401764 0.339885 
#
pcop1fact_bb1 = function(u1,u2,param)
{ p1factbiv(u1,u2,param[c(1,3)],param[c(2,4)],fam1=7,fam2=7,nq=21) }
#
pcop1fact_bb1r = function(u1,u2,param)
{ p1factbiv(u1,u2,param[c(1,3)],param[c(2,4)],fam1=17,fam2=17,nq=21) }
#
# Using function rhoS() based on adaptive integration
library(cubature)
rho_1factbb1_ai = rhoS(c(param),pcop1fact_bb1,zero=0.00001,tol=0.00001)
rho_1factbb1r_ai = rhoS(c(param),pcop1fact_bb1r,zero=0.00001,tol=0.00001)
cat(rho_1factbb1_ai, rho_1factbb1r_ai, "\n")
# 0.3400871 0.3397855
# same to 3 decimal places
}
