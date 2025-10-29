
# Currently allowable copula familes using code of VineCopula
# 1 = Gaussian/normal
# 2 = t
# 4 = Gumbel
# 5 = Frank
# 7 = BB1
# 10 = BB8
# 14 = survival Gumbel
# 17 = survival BB1
# 20 = survival BB8
# this covers a range of different tail behaviors, other families could be added later

#======================================================================

#' Compute new proxies for 1-factor copula based on the mean of observations
#'
#' @description 
#' Compute new proxies for 1-factor copula for 1-parameter linking copulas
#'
#' @param cpar_est estimated parameters (based on complete likelihood with latent variables known or estimated)
#' @param udata nxd matrix of data in (0,1)
#' @param nq number of nodes for Gaussian-Legendre quadrature 
#' @param family vector of code for d linking copula families (choices 1,4,5,14)
#'
#' @return latent_est: proxies as estimates of latent variables
#'
#' @examples
#' # See examples in onefactorEstWithProxy()
#'
#' @export
#'
latentUpdate1factor1 = function(cpar_est,udata,nq,family) 
{ n = dim(udata)[1]
  d = dim(udata)[2]
  gl = gaussLegendre(nq)
  out = .Fortran("latupdate",as.double(cpar_est),as.integer(n),as.integer(d),
     as.double(udata),as.integer(nq),as.double(gl$nodes),as.double(gl$weights),
     as.integer(family), latent=rep(0,n))
  latent_est = out$latent
  return(latent_est)
}

#' Compute new proxies for 1-factor copula based on the mean of observations
#'
#' @description 
#' Compute new proxies for 1-factor copula for 1-parameter or 2-parameter linking copulas
#'
#' @param cpar_est estimated parameters (based on complete likelihood with latent variables known or estimated)
#' @param udata nxd matrix of data in (0,1)
#' @param nq number of nodes for Gaussian-Legendre quadrature 
#' @param family vector of code for d linking copula families (choices 1,2,4,5,7,10,14,17,20)
#'
#' @return latent_est: proxies as estimates of latent variables
#'
#' @examples
#' # See examples in onefactorEstWithProxy()
#'
#' @export
#'
latentUpdate1factor = function(cpar_est,udata,nq,family) 
{ n = dim(udata)[1]
  d = dim(udata)[2]
  gl = gaussLegendre(nq)
  out = .Fortran("latupdate3",as.double(cpar_est),as.integer(n),as.integer(d),
     as.double(udata),as.integer(nq),as.double(gl$nodes),as.double(gl$weights),
     as.integer(family), latent=rep(0,n))
  latent_est = out$latent
  return(latent_est)
}

