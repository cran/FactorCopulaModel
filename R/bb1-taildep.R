
#' Kendall's tau for bivariate normal 
#'
#' @description
#' Kendall's tau for bivariate normal 
#'
#' @param rho in (-1,1)
#'
#' @return Kendall's tau = 2*arcsin(rho)/pi
#'
#' @export
#'
bvn_cpar2tau = function(rho) { (2/pi)*asin(rho) }


#' BB1, given 0<tau<1, find theta and delta with lower tail dependence
#' equal upper tail dependence
#'
#' @description
#' BB1, given 0<tau<1, find theta and delta with equal lower/upper tail dependence
#'
#' @param tau Kendall tau value
#' @param destart starting point for delta
#' @param mxiter maximum number of iterations
#' @param eps tolerance for convergence
#' @param iprint print flag for iterations
#'
#' @return copula parameter (theta,delta) with ltd=utd given tau
#'
#' @examples
#' bb1_tau2eqtd(c(0.1,0.2,0.5))
#'
#' @export
#'
bb1_tau2eqtd = function(tau,destart=1.5,mxiter=30,eps=1.e-6,iprint=FALSE)
{ iter = 0 
  diff = 1.
  ln2 = log(2.)
  de = destart
  rhs = 2/(1-tau)
  while(iter<mxiter & max(abs(diff))>eps)
  { tem = 2-2^(1/de) 
    ltem = -log(tem)
    g = ln2/ltem +2*de -rhs
    gp = (ln2/ltem/de)^2 * (2-tem)/tem +2
    iter = iter+1
    diff = g/gp
    de = de-diff
    while(min(de)<=1. | max(abs(diff))>5.) { diff = diff/2.; de = de+diff; }
    if(iprint) print(c(iter,de,diff))
  }
  if(iter>=mxiter) message("did not converge\n")
  th = ln2/(de*ltem)
  ltd = 2^(-1./(th*de))
  utd = 2.-2^(1./de)
  if(length(tau)==1)
  { out = c(th,de,ltd,utd)
    names(out)=c("theta","delta","ltd","utd")
  }
  else
  { out = cbind(th,de,ltd,utd) 
    names(out)=c("theta","delta","ltd","utd")
  }
  out
}


#' BB1 copula parameter (theta,delta) to tail dependence parameters 
#'
#' @description
#' BB1 copula parameter (theta,delta) to tail dependence parameters 
#'
#' @param cpar copula parameter with theta>0, delta>1 (vector of length 2)
#'           or mx2 matrix with columns for theta and delta
#'
#' @return vector or matrix with lower and upper tail dependence
#'
#' @examples
#' cpar = matrix(c(0.5,1.5,0.8,1.2),byrow=TRUE,ncol=2)
#' bb1_cpar2td(cpar)
#'
#' @export
#'
bb1_cpar2td = function(cpar)
{ if(is.matrix(cpar)) { th = cpar[,1]; de = cpar[,2] }
  else { th = cpar[1]; de = cpar[2] }
  ltd = 2^(-1/(th*de)); utd = 2-2^(1/de)
  out = cbind(ltd,utd)
  colnames(out) = c("lowertaildep", "uppertaildep")
  out
}

