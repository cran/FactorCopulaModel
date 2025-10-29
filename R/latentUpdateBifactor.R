

#' Conditional expectation proxies for bi-factor copula models
#'   with linking copulas in different copula families
#'
#' @description 
#' Conditional expectation proxies for bi-factor copula models
#'   with linking copulas in different copula families
#'
#' @param udata nxd matrix with valies in (0,1)
#' @param cparvec parameters for linking copulas; 
#'   order is global_par1, global_par2, local_par1, local_par2
#' @param grsize group size vector of length mgrp
#' @param family codes for linking copula (VineCopula)
#' @param nq number of Gaussian-Legendre points
#'
#' @return v0: proxies of the global latent variable and
#'   vg: proxies of the local latent variables
#'
#' @examples
#' # See example in bifactorEstWithProxy()
#'
#' @export
#'
latentUpdateBifactor = function(udata,cparvec,grsize,family, nq)
{ mgrp = length(grsize)
  npar = length(cparvec)
  n = nrow(udata)
  npar = length(cparvec)
  mgrp = length(grsize)
  dvar = sum(grsize)
  gl = gaussLegendre(nq)
  out = .Fortran("latupdatebifact3",
     as.integer(npar),as.double(cparvec),
     as.integer(mgrp), as.integer(family),
     as.integer(dvar),as.integer(n),
     as.integer(grsize),as.double(udata),
     as.integer(nq),as.double(gl$nodes),as.double(gl$weights),
     v0mat=rep(0,n), vgmat=matrix(0.0,nrow=n,ncol=mgrp))
  return(list(v0=out$v0mat,vg=out$vgmat))
}
