# Functions for simulation from structured copula models
# Code written by Pavel Krupskii
# Modified in February 2025 to include latent variable in output list.
# VineCopula:BiCopHinv2 is used for C_{1|2}^{-1}

#' Simulate data from nested copula or Gaussian model
#'
#' @description
#' Simulate data from nested copula or Gaussian model
#'
#' @param n sample size 
#' @param grsize G-vector of group sizes for G groups
#' @param cop number code: 1: Gaussian; 2: t; 4: Gumbel; 14: survival Gumbel; 5: Frank; 7: Gumbel+BB1
#  if cop = 1, data have standard normal marginals
#  if cop = 2, data have standard t(df) marginals
#  if cop > 2, data have uniform marginals
#' @param param vector of parameters, length is d+mgrp(+1 for cop==2): 
#     starting with mgrp
#     length(param) = d+d+mgrp if cop=7 (in this case thetas, then deltas.
# The order in param is the same as in start for mvtBifactor(full=FALSE).
#'
#' @return d-dimensional random sample with U(0,1) margins or N(0,1) or t(df) margin; 
#'   v0: n-vector of corresponding global latent variables; and
#'   vg: nxG matrix of corresponding local (group) latent variables
#'
#' @examples
#' grsize = c(3,3,3)
#' cop = 4; param4 = c(1.1,1.1,1.1,  1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3)
#' cop = 14; param14 = c(1.1,1.1,1.1,  1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3)
#' cop = 5; param5 = c(1.1,1.1,1.1,  1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3)
#' cop = 1; param1 = c(0.4,0.4,0.4, 0.5,0.6,0.7,0.8,0.9,0.4,0.5,0.6,0.7 ) 
#' cop = 2; param2 = c(0.4,0.4,0.4, 0.5,0.6,0.7,0.8,0.9,0.4,0.5,0.6,0.7, 7 ) 
#' cop = 7; param7 = c(1.5,1.5,1.5, 0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3, 1.5,1.6,1.7,1.8,1.9,1.4,1.5,1.6,1.7) 
#' set.seed(123)
#' gumdat = rnestfactor(n=10, grsize=grsize, cop=4, param=param4)
#' gumrdat = rnestfactor(n=10, grsize=grsize, cop=14, param=param14)
#' frkdat = rnestfactor(n=10, grsize=grsize, cop=5, param=param5)
#' gaudat = rnestfactor(n=10, grsize=grsize, cop=1, param=param1)
#' bvtdat = rnestfactor(n=10, grsize=grsize, cop=2, param=param2)
#' gumbb1dat = rnestfactor(n=10, grsize=grsize, cop=7, param=param7)
#' summary(gumbb1dat$data)
#' round(cor(gumbb1dat$data),2)
#' summary(gumbb1dat$v0)
#' summary(gumbb1dat$vg)
#'
#' @details
#' The user can modify this code to get other linking copulas.
#'
#' @export
#'
rnestfactor = function(n,grsize,cop,param)
{ d = sum(grsize)
  mgrp = length(grsize)
  if(cop==2) { df0 = param[mgrp+d+1] }
  if(cop==1 || cop==2) { z0 = rnorm(n) } else { z0 = runif(n) }
  z = matrix(0,nrow=n,ncol=mgrp)
  zdata = matrix(0,nrow=n,ncol=d)
  if(cop==1 || cop==2) # Gaussian or t
  { for (jg in 1:mgrp) { z[,jg] = z0*param[jg]+ sqrt(1-param[jg]^2)*rnorm(n) }}
  else if(cop==5) # Frank 
  #{ for (jg in 1:mgrp) { z[,jg]= qcondfrk(runif(n),z0,param[jg]); } } 
  { for (jg in 1:mgrp) { z[,jg] = VineCopula::BiCopHinv2(runif(n),z0,family=5,par=param[jg]) } } 
  else if(cop==4 || cop==14 || cop==7)  # Gumbel/Gum survGumbel/survGum or Gumbel/BB1 
  { coplat = cop
    if(coplat==7) coplat = 4
    for (jg in 1:mgrp)  
    { z[,jg] = VineCopula::BiCopHinv2(runif(n),z0,family=coplat,par=param[jg]) } 
  }
  else { message("this copula family is not implemented"); return(NA) }
  
  ind = 0
  for (jg in 1:mgrp)  
  { ind1 = ind+1  
    ind2 = ind+grsize[jg]
    ind = ind+grsize[jg]  
    for(ij in ind1:ind2)  
    { ijm = ij+mgrp
      if(cop==1 || cop==2)  # Gaussian or t
      { zdata[,ij] = z[,jg]*param[ijm]+sqrt(1-param[ijm]^2)*rnorm(n); }
      else if(cop==5) { zdata[,ij] = VineCopula::BiCopHinv2(runif(n),z[,jg],family=5,par=param[ijm]) }
      else if(cop==4 || cop==14) # Gumbel or survGumbel
      { zdata[,ij] = VineCopula::BiCopHinv2(runif(n),z[,jg],family=cop,par=param[ijm]) } 
      else if(cop==7) # Gumbel/BB1
      #{ for(i in 1:n) { zdata[i,ij] = qcondbb1(runif(1),z[i,jg],param[c(ijm,ijm+d)]); } } 
      { zdata[,ij] = VineCopula::BiCopHinv2(runif(n),z[,jg],family=7,par=param[ijm],par2=param[ijm+d]) }  
    }
  }
  if(cop==2) # t
  { for(i in 1:n)
    { zdata[i,] = zdata[i,]/sqrt(rchisq(1,df=df0)/df0); 
    }
  }
  list(data=zdata, v0=z0, vg=z) 
}

