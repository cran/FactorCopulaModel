# Density and log-likelihood functions for 1-factor copula
# fortran 90 with gradient and Hessian (for large dimensions) and R standalone versions

#' Integrand for 1-factor copula with 1-parameter bivariate linking copula families;
#' or for m-parameter bivariate linking copulas 
#'
#' @description
#' Integrand for 1-factor copula 
#'
#' @param u0 latent variable for integrand
#' @param uvec vector of length d, components in (0,1)
#' @param dcop name of function of bivariate copula density (common for all variables),
#'     dcop accepts input of form  of d-vector or dxm matrix
#' @param param d-vector or mxd matrix, parameter of dcop 
#'
#' @return integrand for 1-factor copula density
#'
#' @export
#'
d1factcop = function(u0,uvec,dcop,param)
{ # dcop(u0,uvec,param) is a vector
  #if(is.matrix(param)) param=t(param)  # for BB1?
  tem = dcop(u0,uvec,param)
  prod(tem)  # product
}

#' max likelihood (min negative log-likelihood) for 1-factor copula model 
#'
#' @description
#' min negative log-likelihood for 1-factor copula model 
#'
#' @param nq number of quadrature points
#' @param start starting point (d-vector or m*d vector, e.g. 2*d vector for BB1)
#' @param udata nxd matrix of uniform scores
#' @param dcop name of function for a bivariate copula density
#'       (common for all variables)
#' @param LB lower bound on parameters (scalar or same dimension as start) 
#' @param UB upper bound on parameters (scalar or same dimension as start)
#' @param prlevel printlevel for nlm()
#' @param mxiter maximum number of iteration for nlm()
#'
#' @return nlm object with minimum, estimate, hessian at MLE
#'
#' @examples
#' # See example in r1factor() 
#'
#' @export
#'
ml1factor = function(nq,start,udata,dcop,LB=0,UB=1.e2,prlevel=0,mxiter=100)
{ start = c(start)
  np = length(start)
  gl = gaussLegendre(nq)
  wl = gl$weights
  xl = gl$nodes
  nllkfn1 = function(param)
  { d = ncol(udata)  
    n = nrow(udata)
    nlb = sum(param<=LB | param>=UB)
    if(nlb>0) { return(1.e10) }
    nllk = 0
    if(np==d) par0 = param
    else { par0 = matrix(param,nrow=d,byrow=TRUE) }  # to handle BB1 
    for(i in 1:n)
    { uvec = udata[i,]
      integl = 0; 
      for(iq in 1:nq)
      { fvalgl = d1factcop(xl[iq],uvec,dcop,par0)
        integl = integl+wl[iq]*fvalgl
      }
      nllk = nllk-log(integl);
    }
    return(nllk);
  }
  mle = nlm(nllkfn1,p=start,hessian=T,iterlim=mxiter,print.level=prlevel)
  mle
}  


#' min negative log-likelihood for 1-factor copula model (some parameters can be fixed)
#'
#' @description
#' min negative log-likelihood (nllk) for 1-factor copula model (some parameters can be fixed)
#'
#' @param nq number of quadrature points
#' @param start starting point (d-vector or m*d vector, e.g. 2*d vector for BB1)
#' @param ifixed vector of length(param) of True/False, such that
#'       ifixed[i]=TRUE iff param[i] is fixed at the given value start[j]
#' @param udata nxd matrix of uniform scores
#' @param dcop name of function for a bivariate copula density
#'       (common for all variables)
#' @param LB lower bound on parameters (scalar or same dimension as start) 
#' @param UB upper bound on parameters (scalar or same dimension as start)
#' @param prlevel printlevel for nlm()
#' @param mxiter maximum number of iteration for nlm()
#'
#' @return nlm object with nllk value, estimate, hessian at MLE
#'
#' @export
#'
ml1factor_v2 = function(nq,start,ifixed,udata,dcop,LB=0,UB=1.e2,prlevel=0,mxiter=100)
{ start = c(start)
  np = length(start)
  gl = gaussLegendre(nq)
  wl = gl$weights
  xl = gl$nodes
  nllkfn1 = function(param)
  { d = ncol(udata)  
    n = nrow(udata)
    param0 = start0; param0[!ifixed] = param  # full length parameter for d1factcop
    np0 = length(param0)
    nlb = sum(param0<=LB | param0>=UB)
    if(nlb>0) { return(1.e10) }
    nllk = 0
    if(np==d) par0 = param0
    else { par0 = matrix(param0,nrow=d,byrow=TRUE) }  # to handle BB1 
    for(i in 1:n)
    { uvec = udata[i,]
      integl = 0; 
      for(iq in 1:nq)
      { fvalgl = d1factcop(xl[iq],uvec,dcop,par0)
        integl = integl+wl[iq]*fvalgl
      }
      nllk = nllk-log(integl);
    }
    return(nllk);
  }
  start0 = start # in order to keep the fixed components
  mle = nlm(nllkfn1,p=start[!ifixed],hessian=T,iterlim=mxiter,print.level=prlevel)
  mle
}  

#============================================================

#' negative log-likelihood of 1-factor copula for input to posDefHessMin and posDefHessMinb 
#' 
#' @description
#' negative log-likelihood (nllk) of 1-factor copula for input to posDefHessMin and posDefHessMinb
#'
#' @param  param parameter vector
#' @param  dstruct list with  data set $data, copula name $copname,
#'          $quad is list with quadrature weights and nodes, 
#'          $repar is code for reparametrization (for Gumbel, BB1),
#'          $nu is positive degree of freedom parameter (for t) 
#'       (linking copula is common for all variables).
#' Options for copname are: frank, gumbel, bb1, t.
#' For reflected gumbel or bb1, use something like dstruct$dat = 1-udata
#' @param  iprfn print flag for function and gradient (within Newton-Raphson) iterations)
#'   for BB1, param is 2*d-vector with th1,de1,th2,de2,...
#'
#' @return nllk, grad (gradient), hess (hessian) at MLE
#'
#' @examples
#' cpar_gum = seq(1.9,3.7,0.2)
#' d = length(cpar_gum)
#' n = 300
#' param = c(rbind(cpar_gum,rep(0,d)))  # second par2 is 0 for VineCopula
#' set.seed(111)
#' gum_obj = r1factor(n,d,param,famvec=rep(4,d)) # uses VineCopula
#' udat = gum_obj$udata
#' zdat = qnorm(udat)
#' rmat = cor(zdat)
#' print(round(rmat,3))
#' # run factanal to get loading (rho in normal scale close to Spearman rho)
#' fa1 = factanal(covmat=rmat,factors=1)
#' loadings = c(fa1$loading)
#' # convert loadings to Frank, Gumbel and BB1 parameters 
#' start_frk = frank_rhoS2cpar(loadings)
#' start_gum = gumbel_rhoS2cpar(loadings)
#' tau = bvn_cpar2tau(loadings)
#' start_bb1 = bb1_tau2eqtd(tau)
#' start_bb1 = c(t(start_bb1[,1:2]))
#' gl = gaussLegendre(25)
#' dstrfrk1 = list(copname="frank",data=udat,quad=gl,repar=0, pdf=1)
#' dstrfrk = list(copname="frank",data=udat,quad=gl,repar=0, pdf=0)
#' obj1 = onefactorcop_nllk(start_frk,dstrfrk1) #nllk only
#' obj = onefactorcop_nllk(start_frk,dstrfrk) # nllk, grad, hess
#' print(obj1$fnval)
#' print(obj$grad)
#' #
#' ml_frk = posDefHessMinb(start_frk,onefactorcop_nllk,ifixed=rep(FALSE,d),
#'   dstruct=dstrfrk, LB=rep(-30,d),UB=rep(30,d),iprint=TRUE,eps=1.e-5)
#' dstrgum = list(copname="gumbel",data=udat,quad=gl,repar=0)
#' ml_gum = posDefHessMinb(start_gum,onefactorcop_nllk,ifixed=rep(FALSE,d),
#'   dstruct=dstrgum, LB=rep(-30,d),UB=rep(30,d),iprint=TRUE,eps=1.e-5)
#' dstrgumr = list(copname="gumbel",data=1-udat,quad=gl,repar=0)
#' ml_gumr = posDefHessMinb(start_gum,onefactorcop_nllk,ifixed=rep(FALSE,d),
#'   dstruct=dstrgumr, LB=rep(-30,d),UB=rep(30,d),iprint=TRUE,eps=1.e-5)
#' dstrtnu = list(copname="t",data=udat,quad=gl,repar=0,nu=10)
#' ml_tnu = posDefHessMinb(loadings,onefactorcop_nllk,ifixed=rep(FALSE,d),
#'   dstruct=dstrtnu, LB=rep(-1,d),UB=rep(1,d),iprint=TRUE,eps=1.e-5)
#' dstrbb1 = list(copname="bb1",data=udat,quad=gl,repar=0)
#' ml_bb1 = posDefHessMinb(start_bb1,onefactorcop_nllk,ifixed=rep(FALSE,2*d),
#'   dstruct=dstrbb1, LB=rep(c(0,1),d),UB=rep(20,2*d),iprint=TRUE,eps=1.e-5)
#' dstrbb1r = list(copname="bb1",data=1-udat,quad=gl,repar=0)
#' ml_bb1r = posDefHessMinb(start_bb1,onefactorcop_nllk,ifixed=rep(FALSE,2*d),
#'   dstruct=dstrbb1r, LB=rep(c(0,1),d),UB=rep(20,2*d),iprint=TRUE,eps=1.e-5)
#' cat(ml_frk$fnval, ml_gum$fnval, ml_gumr$fnval, ml_tnu$fnval, ml_bb1$fnval, ml_bb1r$fnval, "\n")
#' # -1342.936 -1560.391 -1198.837 -1454.907 -1560.45 -1549.935
#' cat(ml_frk$iter, ml_gum$iter, ml_gumr$iter, ml_tnu$iter, ml_bb1$iter, ml_bb1r$iter, "\n")
#' # 4 4 5 4 16 6
#'
#' @details
#' linked to Fortran 90 code for speed
#'
#' @export
#'
onefactorcop_nllk = function(param,dstruct,iprfn=FALSE)
{ udata = dstruct$data
  copname = dstruct$copname
  copname = tolower(copname)
  gl = dstruct$quad
  repar = dstruct$repar
  if(copname=="t") nu = dstruct$nu # assumed scalar
  d = ncol(udata);
  n = nrow(udata);
  wl = gl$weights
  xl = gl$nodes
  nq = length(xl)
  npar = length(param)
  if(repar==2) { pr0 = param; param = param^2+rep(c(0,1),d); } 
  if(repar==1) { pr0 = param; param = param^2+1; }
  # check for copula type and call different f90 routines
  if(copname=="frank" | copname=="frk")
  #{ out = .Fortran("frk1fact",
  #    as.integer(npar), as.double(param), as.integer(d), as.integer(n),
  #    as.double(udata), 
  #    as.integer(nq), as.double(wl), as.double(xl), 
  #    nllk=as.double(0.),grad=as.double(rep(0,npar)),
  #    hess=as.double(rep(0,npar*npar))  )
  #}
  #else if(copname=="lfrank" | copname=="lfrk")
  { out = .Fortran("lfrk1fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="gumbel" | copname=="gum")
  { out = .Fortran("gum1fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else if(copname=="t")
  { tl = qt(xl,nu)   # nu is scalar
    tdata = qt(udata,nu) # transform data and nodes to t(nu) scale
    out = .Fortran("t1fact",
      as.integer(npar), as.double(param), as.double(nu),
      as.integer(d), as.integer(n), as.double(tdata), 
      as.integer(nq), as.double(wl), as.double(tl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  # assume param(2xd) has been converted to a column vector
  else if(copname=="bb1")
  { out = .Fortran("bb11fact",
      as.integer(npar), as.double(param), as.integer(d), as.integer(n),
      as.double(udata), 
      as.integer(nq), as.double(wl), as.double(xl), 
      nllk=as.double(0.),grad=as.double(rep(0,npar)),
      hess=as.double(rep(0,npar*npar))  )
  }
  else { message("copname not available"); return (NA); }
  if(iprfn) print(cbind(param,out$grad))
  nllk =out$nllk; hess =matrix(out$hess,npar,npar); grad =out$grad;
  if(repar>0) 
  { tem = diag(2*pr0);  # jacobian
    hess =tem%*%hess%*%tem + 2*diag(grad)  
    grad =2*pr0*grad; # pr0 is transformed parameter
  }
  list(fnval=nllk, grad=grad, hess=hess)
}

#' min negative log-likelihood for 1-factor copula with nlm()
#'
#' @description
#' min negative log-likelihood for 1-factor copula with nlm()
#'
#' @param nq number of quadrature points
#' @param start starting point (d-vector or m*d vector, e.g. 2*d vector for BB1)
#' @param udata nxd matrix of uniform scores
#' @param copname name of copula family such as "gumbel", "frank", "bb1", "t"
#'       (copname common for all variables)
#' @param LB lower bound on parameters (scalar or same dimension as start) 
#' @param UB upper bound on parameters (scalar or same dimension as start)
#' @param ihess flag for hessian option in nlm()
#' @param prlevel printlevel for nlm()
#' @param mxiter max number of iterations for nlm()
#' @param nu degree of freedom parameter if copname ="t"
#'
#' @return MLE as nlm object (estimate, Hessian, SEs, nllk)
#'
#' @examples
#' # See example in r1factor() 
#'
#' @export
#'
ml1factor_f90 = function(nq,start,udata,copname,LB=0,UB=40,ihess=FALSE,prlevel=0,
  mxiter=100,nu=3)
{ # repar hasn't been added to this version
  copname = tolower(copname)
  gl = gaussLegendre(nq)
  wl = gl$weights
  xl = gl$nodes
  nllkfn1 = function(param)
  { d = ncol(udata)  
    n = nrow(udata)
    npar = length(param)
    nlb = sum(param<=LB | param>=UB)
    if(nlb>0) { return(1.e10) }
    if(copname=="frank" | copname=="frk")
    { out = .Fortran("lfrk1fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="gumbel" | copname=="gum")
    { out = .Fortran("gum1fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="t")
    { tl = qt(xl,nu)   # nu is scalar
      tdata = qt(udata,nu) # transform data and nodes to t(nu) scale
      out = .Fortran("t1fact",
        as.integer(npar), as.double(param), as.double(nu),
        as.integer(d), as.integer(n), as.double(tdata), 
        as.integer(nq), as.double(wl), as.double(tl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    # assume param(2xd) has been converted to a column vector
    else if(copname=="bb1")
    { out = .Fortran("bb11fact",
        as.integer(npar), as.double(param), as.integer(d), as.integer(n),
        as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else { message("copname not available"); return (NA); }
    nllk = out$nllk
    if(is.nan(nllk) || is.infinite(nllk) ) { return(1.e10) }
    if(any(is.nan(out$grad)) || any(is.infinite(out$grad)) ) { return(1.e10) }
    attr(nllk,"gradient") = out$grad  # for nlm with analytic gradient
    nllk
  }
  mle = nlm(nllkfn1,p=start,hessian=ihess,iterlim=mxiter,print.level=prlevel,
      check.analyticals=FALSE)
  #if(prlevel>0)
  #{ cat("nllk: \n")
  #  print(mle$minimum)
  #  cat("MLE: \n")
  #  print(mle$estimate)
  #  if(ihess) 
  #  { iposdef = isPosDef(mle$hessian)
  #    if(iposdef)
  #    { cat("SEs: \n")
  #      acov = solve(mle$hessian)
  #      SEs = sqrt(diag(acov))
  #      print(SEs)
  #      mle$SE = SEs
  #    }
  #    else message("Hessian not positive definite\n")
  #  }
  #}
  mle
}  



