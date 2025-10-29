\name{rnestfactor}
\alias{rnestfactor}
\title{Simulate data from nested copula or Gaussian model}
\description{
Simulate data from nested copula or Gaussian model
}

\usage{
rnestfactor(n,grsize,cop,param)
}
\arguments{
\item{n}{sample size }
\item{grsize}{G-vector of group sizes for G groups}
\item{cop}{number code: 1: Gaussian; 2: t; 4: Gumbel; 14: survival Gumbel; 5: Frank; 7: Gumbel+BB1}
\item{param}{vector of parameters, length is d+mgrp(+1 for cop==2): }
}
\value{
 d-dimensional random sample with U(0,1) margins or N(0,1) or t(df) margin; v0: n-vector of corresponding global latent variables; and
vg: nxG matrix of corresponding local (group) latent variables
}
\examples{
grsize = c(3,3,3)
cop = 4; param4 = c(1.1,1.1,1.1,  seq(1.5,2.3,0.1))
cop = 14; param14 = c(1.1,1.1,1.1,  seq(1.5,2.3,0.1))
cop = 5; param5 = c(1.1,1.1,1.1,  seq(1.5,2.3,0.1))
cop = 1; param1 = c(0.4,0.4,0.4, 0.5,0.6,0.7,0.8,0.9,0.4,0.5,0.6,0.7 ) 
cop = 2; param2 = c(0.4,0.4,0.4, 0.5,0.6,0.7,0.8,0.9,0.4,0.5,0.6,0.7, 7 ) 
cop = 7; param7 = c(1.5,1.5,1.5, seq(0.5,1.3,0.1), 1.5,1.6,1.7,1.8,1.9,1.4,1.5,1.6,1.7) 
set.seed(123)
gumdat = rnestfactor(n=10, grsize=grsize, cop=4, param=param4)
gumrdat = rnestfactor(n=10, grsize=grsize, cop=14, param=param14)
frkdat = rnestfactor(n=10, grsize=grsize, cop=5, param=param5)
gaudat = rnestfactor(n=10, grsize=grsize, cop=1, param=param1)
bvtdat = rnestfactor(n=10, grsize=grsize, cop=2, param=param2)
gumbb1dat = rnestfactor(n=10, grsize=grsize, cop=7, param=param7)
summary(gumbb1dat$data)
round(cor(gumbb1dat$data),2)
summary(gumbb1dat$v0)
summary(gumbb1dat$vg)
}
\details{
The user can modify this code to get other linking copulas.
}
