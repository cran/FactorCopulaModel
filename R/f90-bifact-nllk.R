# Code written by Pavel Krupskii, function renamed Feb 2025

# structured bi-factor copula model : functions to input to posDefHessMinb
# negative log-likelihoods and derivatives computed in f90

#' negative log-likelihood of bi-factor structured factor copula and derivatives computed in f90
#'  for input to posDefHessMinb
#' 
#' @description
#' negative log-likelihood of bi-factor structured factor copula and derivatives
#' 
#' @param  param parameter vector;
#' parameters for copulas linking U_{ij} and V_0 go first;
#' parameters for copulas linking U_{ij} and V_g (j in group g) given V_0 go next. 
#' For BB1 linking copulas for the global latent, the order is
#'   theta1,..,theta[d],delta1, ...,delta[d];
#' V_0 is the global latent variable that loads on all variables;
#' V_j is a latent variable that loads only for variables in group j
#'  (by group g=1,2,..,mgrp etc)
#' @param  dstruct list with  data set $data, copula name $copname,
#'          $quad is a Gauss-Legendre quadrature object,
#'          $repar is a flag for reparametrization (for Gumbel, BB1),
#'          $nu is a 2-vector with 2 degree of freedom parameters (for t)
#'          $grsize is a vector with group sizes;
#' if dstruct$pdf == 1 the function evaluates nllk only 
#'       (and returns zero gradient and hessian).
#' Options for copname are: frank, gumbel, gumbelfrank, bb1frank, bb1gumbel, t.
#' @param  iprfn indicator for printing of function and gradient (within Newton-Raphson iterations)
#' 
#' @return  nllk, grad, hess (gradient and hessian included)
#' 
#' @examples
#' # donttest
#' grsize = c(4,4,3)
#' d = sum(grsize)
#' n = 500
#' mgrp = length(grsize)
#' par_bi = c(seq(1.4,3.4,0.2),seq(2.0,1.7,-0.1),seq(1.9,1.6,-0.1),rep(1.4,3))
#' set.seed(333)
#' udat_obj = rbifactor(n,grsize,cop=4,par_bi)
#' udat = udat_obj$data
#' summary(udat_obj$v0)
#' summary(udat_obj$vg)
#' zdat = qnorm(udat)
#' rmat = cor(zdat)
#' round(rmat,3)
#' # run bifactor_fa to get bi-factor correlation structure
#' bifa = bifactor_fa(grsize,start=c(rep(0.8,d),rep(0.2,d)),cormat=rmat,n=n,prlevel=0)
#' loading1 = bifa$parmat[,1] # correlations
#' pcor = bifa$parmat[,2] # partial correlations given global latent
#' pcor2 = pcor;  pcor2[pcor<0]=0.05  # for cases with only positive dependence
#' # starting values for different cases
#' # convert loading/pcor to Frank, Gumbel and BB1 parameters etc
#' # Frank for conditional given global latent can allow for conditional negative dependence
#' start_frk1 = frank_rhoS2cpar(loading1)
#' start_frk2 = frank_rhoS2cpar(pcor)
#' start_frk = c(start_frk1,start_frk2)
#' start_gum1 = gumbel_rhoS2cpar(loading1)
#' start_gum2 = gumbel_rhoS2cpar(pcor2)
#' start_gum = c(start_gum1,start_gum2)
#' start_tnu = c(loading1,pcor)
#' tau = bvn_cpar2tau(c(loading1))
#' # order of BB1 parameters has all thetas and then all deltas (different from 1-factor)
#' start_bb1 = bb1_tau2eqtd(tau)
#' start_bb1 = c(start_bb1[,1:2])
#' start_gumfrk = c(start_gum1,start_frk2)
#' start_bb1frk = c(start_bb1,start_frk2)
#' start_bb1gum = c(start_bb1,start_gum2)
#' #
#' gl = gaussLegendre(25)
#' npar = 2*d
#' dstrfrk = list(data=udat,copname="frank",quad=gl,repar=0,grsize=grsize,pdf=0)
#' dstrfrk1 = list(data=udat,copname="frank",quad=gl,repar=0,grsize=grsize,pdf=1)
#' obj1 = bifactorcop_nllk(start_frk,dstrfrk1) # nllk only  
#' obj = bifactorcop_nllk(start_frk,dstrfrk) # nllk, grad, hess 
#' print(obj1$fnval)
#' print(obj$grad)
#' ml_frk = posDefHessMinb(start_frk,bifactorcop_nllk,ifixed=rep(FALSE,npar),
#'   dstrfrk, LB=rep(-20,npar), UB=rep(30,npar), mxiter=30, eps=5.e-5,iprint=TRUE)
#' dstrgum = list(data=udat,copname="gumbel",quad=gl,repar=0,grsize=grsize,pdf=0)
#' ml_gum = posDefHessMinb(start_gum,bifactorcop_nllk,ifixed=rep(FALSE,npar),
#'   dstrgum, LB=rep(1,npar), UB=rep(20,npar), mxiter=30, eps=5.e-5,iprint=TRUE)
#' dstrgumfrk = list(data=udat,copname="gumbelfrank",quad=gl,repar=0,grsize=grsize,pdf=0)
#' ml_gumfrk = posDefHessMinb(start_gumfrk,bifactorcop_nllk,ifixed=rep(FALSE,npar),
#'   dstrgumfrk, LB=c(rep(1,d),rep(-20,d)), UB=rep(25,npar), mxiter=30, eps=5.e-5,iprint=TRUE)
#' dstrtnu = list(data=udat,copname="t",quad=gl,repar=0,grsize=grsize,nu=c(10,20),pdf=0)
#' # slow because of many qt() calculations
#' # numerical issues because data does not have both upper and lower taildep
#' ml_tnu = posDefHessMinb(start_tnu,bifactorcop_nllk,ifixed=rep(FALSE,npar),
#'   dstrtnu, LB=rep(-1,npar), UB=rep(1,npar), mxiter=30, eps=5.e-5,iprint=TRUE)
#' npar3 = 3*d
#' dstrbb1frk = list(data=udat,copname="bb1frank",quad=gl,repar=0,grsize=grsize,pdf=0)
#' ml_bb1frk = posDefHessMinb(start_bb1frk,bifactorcop_nllk,ifixed=rep(FALSE,npar3),
#'   dstrbb1frk, LB=c(rep(0,d),rep(1,d),rep(-20,d)), UB=rep(20,npar3), mxiter=30, eps=5.e-5,iprint=TRUE)
#' dstrbb1gum = list(data=udat,copname="bb1gumbel",quad=gl,repar=0,grsize=grsize,pdf=0)
#' ml_bb1gum = posDefHessMinb(start_bb1gum,bifactorcop_nllk,ifixed=rep(FALSE,npar3),
#'   dstrbb1gum, LB=c(rep(0,d),rep(1,2*d)), UB=rep(20,npar3), mxiter=30, eps=5.e-5,iprint=TRUE)
#' #
#' cat(ml_frk$fnval, ml_gum$fnval, ml_gumfrk$fnval, ml_tnu$fnval, ml_bb1frk$fnval, ml_bb1gum$fnval, "\n")
#' # -2256.602 -2574.16 -2509.274 -6793.963 -2514.124 -2581.214
#' cat(ml_frk$iter, ml_gum$iter, ml_gumfrk$iter, ml_tnu$iter, ml_bb1frk$iter, ml_bb1gum$iter, "\n")
#' # 5 6 5 23 15 12
#' # bi-factor t(10)/t(20) failed because some parameters approached the
#' # upper bound of 1 in which case the numerical integration is inaccurate.
#' 
#' @export 
#' 
bifactorcop_nllk = function(param,dstruct,iprfn=FALSE)
{ udata = dstruct$data
  copname = dstruct$copname
  copname = tolower(copname)
  gl = dstruct$quad
  repar = dstruct$repar
  grsize = dstruct$grsize # vector of group sizes
  mgrp = length(grsize); # total number of groups
  dvar = ncol(udata); # total number of variables
  n = nrow(udata);
  wl = gl$weights
  xl = gl$nodes
  ipdf = dstruct$pdf # ipdf=1 for pdf only, 0 otherwise
  if(copname=="t") # | copname=="tapprox") 
  { nu = dstruct$nu; # assumed a vector of length 2
    tl1 = qt(xl,nu[1]);
    tl2 = qt(xl,nu[2]);
    tl = cbind(tl1,tl2);
  }
  nq = length(xl)
  npar = length(param)
  #if(repar==2) { pr0=param; param=param^2+rep(c(0,1),mgrp); } 
  if(repar==1) { pr0 = param; param = param^2+1; }
  
  if(ipdf==0)  # derivatives
  { # check for copula type and call different f90 routines
    # most options are not available here but are available in version
    # that is companion to "Dependence Modeling with Copulas", by Joe (2014).
    if(copname=="frank" | copname=="frk") 
    { out = .Fortran("strfrk2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="gumbel" | copname=="gum")
    { out = .Fortran("strgum2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="gumfrk" | copname=="gumbelfrank")
    { out = .Fortran("strgumfrk2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else 
    if(copname=="bb1frk" | copname=="bb1frank")
    { out = .Fortran("strbb1frk2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    else if(copname=="bb1gum" | copname=="bb1gumbel")
    { out = .Fortran("strbb1gum2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }    
    else if(copname=="t")
    { out= .Fortran("strt2",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(tl), 
        nllk=as.double(0.),grad=as.double(rep(0,npar)),
        hess=as.double(rep(0,npar*npar))  )
    }
    ## add option of copname=="tapprox") using strt1ipol
    #else if(copname=="tapprox")
    #{ nu1plusone=nu[1]+1 # for conditional
    #  # can replace pp by something better later
    #  pp=c(.0001,.0002,.0005,.001,.002,.005,seq(.01,.99,.01),.995,.998,.999,.9995,.9998,.9999)
    #  nipol=length(pp)
    #  qq=qt(pp,nu1plusone)
    #  pder=pcderiv(qq,pp)  # get derivs for the interpolation
    #  # pass qq,pp,pder to f90 and make call to pchev in fortran
    #  out= .Fortran("strt2ipol",
    #    as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
    #    as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
    #    as.integer(nq), as.double(wl), as.double(tl), 
    #    as.integer(nipol), as.double(qq), as.double(pp), as.double(pder), 
    #    nllk=as.double(0.),grad=as.double(rep(0,npar)),
    #    hess=as.double(rep(0,npar*npar))  )
    #}
    else { message("copname not available"); return (NA); }
    if(iprfn) print(cbind(param,out$grad))
    nllk = out$nllk; hess = matrix(out$hess,npar,npar); grad = out$grad;
    if(repar>0) 
    { tem = diag(2*pr0);  # jacobian
      hess = tem%*%hess%*%tem + 2*diag(grad)  
      grad = 2*pr0*grad; # pr0 is transformed parameter
    }
  }

  if(ipdf==1) # no derivatives
  { # check for copula type and call different f90 routines
    if(copname=="frank" | copname=="frk") 
    { out = .Fortran("strfrk2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.) )
    }
    else if(copname=="gumbel" | copname=="gum")
    { out = .Fortran("strgum2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.) )
    }
    else if(copname=="gumfrk" | copname=="gumbelfrank")
    { out = .Fortran("strgumfrk2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.) )
    }
    else 
    if(copname=="bb1frk" | copname=="bb1frank")
    { out = .Fortran("strbb1frk2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.) )
    }
    else if(copname=="bb1gum" | copname=="bb1gumbel")
    { out = .Fortran("strbb1gum2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(xl), 
        nllk=as.double(0.) )
    }
    else if(copname=="t")
    { out = .Fortran("strt2nllk",
        as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
        as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
        as.integer(nq), as.double(wl), as.double(tl), 
        nllk=as.double(0.) )
    }
    #else if(copname=="tapprox")
    #{ nu1plusone=nu[1]+1 # for conditional
    #  # can replace pp by something better later
    #  pp=c(.0001,.0002,.0005,.001,.002,.005,seq(.01,.99,.01),.995,.998,.999,.9995,.9998,.9999)
    #  nipol=length(pp)
    #  qq=qt(pp,nu1plusone)
    #  pder=pcderiv(qq,pp)  # get derivs for the interpolation
    #  # pass qq,pp,pder to f90 and make call to pchev in fortran
    #  out= .Fortran("strt2ipolnllk",
    #    as.integer(npar), as.double(param), as.integer(mgrp), as.integer(n),
    #    as.integer(dvar), as.double(nu), as.integer(grsize),as.double(udata), 
    #    as.integer(nq), as.double(wl), as.double(tl), 
    #    as.integer(nipol), as.double(qq), as.double(pp), as.double(pder), 
    #    nllk=as.double(0.) )
    #}
    else { message("copname not available"); return (NA); }
    nllk = out$nllk; grad =0; hess =0;
  }
  list(fnval=nllk, grad=grad, hess=hess)
}

