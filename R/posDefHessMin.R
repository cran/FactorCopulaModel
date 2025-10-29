#' Minimization with modified Newton-Raphson iterations,
#' Hessian is modified to be positive definite at each step.
#' Algorithm and code produced by Pavel Krupskii (2013)
#' see PhD thesis Krupskii (2014), UBC and
#' Section 6.2 of # Joe (2014) Dependence Models with Copulas. Chapman&Hall/CRC
#'
#' @description 
#' modified Newton-Raphson minimization with positive Hessian 
#'
#' @param param starting point for minimization 
#' @param objfn function to be minimized with gradient and Hessian
#' @param dstruct  list with  data set and other variables used by objfn
#' @param LB lower bound vector
#' @param UB upper bound vector
#' @param mxiter  max number of iterations
#' @param eps  tolerance for Newton-Raphson iterations
#' @param bdd  bound on difference of 2 consecutive iterations
#'  (useful is starting point is far from solution and func is far from convex)
#' @param iprint  control on amount of printing,
#'          FALSE for no printing of iterations and
#'          TRUE for printing x^{(k)} on each iteration.
#'
#' @return list with
#' fnval = function value at minimum;
#' parmin = param for minimum;
#' invh = inverse Hessian;
#' iconv = 1 if converged, -1 for a boundary point, 0 otherwise;
#' iter = number of iterations.
#'
#' @examples
#' # See examples in onefactorcop_nllk(), bifactorcop_nllk(), nestfactorcop_nllk()
#'
#' @export
#'
posDefHessMin = function(param,objfn,dstruct,LB,UB,mxiter=30,eps=1.e-6,bdd=5,iprint=FALSE)
{ np = length(param)
  iter = 0; mxdif = 1;
  while(mxdif>eps & iter<mxiter)
  { fn = objfn(param,dstruct,iprint);
    fnold = fn$fnval; 
    evd = eigen(fn$hess,F);
    eigval = evd$values;
    eigvec = evd$vectors;
    ind = (eigval<=0);
    #replace negative eigenvalues with a positive number, keeping eigenvectors
    if(sum(ind)>0 & iprint) message(paste("replacing eigenvalues in iteration ", iter))
    eigval[ind] = 1e-1;    
    newhess = eigvec%*%diag(eigval)%*%t(eigvec);
    diff = solve(newhess,fn$grad); 
    mxdif = max(abs(diff))
    iter = iter+1
    param = param-diff
    outb = 0;
    if(any(param<=LB) | any(param>=UB)) outb = 1;
    if(is.infinite(fn$fnval) | is.nan(fn$fnval)) outb = 1;
    if(any(is.infinite(fn$grad)) | any(is.nan(fn$grad))) outb = 1;
    # modified NR to limit size of step to within bounds
    while(outb==1 | mxdif>bdd)
    { mxdif = mxdif/2.;
      outb = 0;
      diff = diff/2; param = param+diff; 
      if(any(param<=LB) | any(param>=UB)) outb = 1;
    }
    #i01 = (param<1e-4  & param>0);
    #i02 = (param>-1e-4 & param<0);
    #param[i01]=1e-4;
    #param[i02]=-1e-4;
    # find a step size to ensure the objective fn decreases
    # later to have option on objfn to not compute derivs
    fnnew = objfn(param,dstruct,F)$fnval; 
    while(fnold<=fnnew) 
    { mxdif = mxdif/2.;
      diff = diff/2; param=param+diff; 
      fnnew = objfn(param,dstruct,F)$fnval;   
    }
    if(iprint) { cat(iter, fnnew, mxdif,"\n") }
  }
  iconv=1;
  if(iter>=mxiter) { message("did not converge\n"); iconv = 0; }
  else { if(max(abs(fn$grad))>eps) iconv = -1 }
  # check for positive definite
  iposdef = isPosDef(fn$hess)
  # inverse Hessian
  invh = solve(fn$hess)
  list(parmin=param,fnval=fn$fnval,invh=invh,iconv=iconv,iposdef=iposdef,iter=iter)
}

 
#' Version with ifixed as argument
#'
#' @description 
#' modified Newton-Raphson minimization with positive Hessian 
#'
#' @param param starting point for minimization 
#' @param objfn function to be minimized with gradient and Hessian
#' @param ifixed vector of length(param) of TRUE/FALSE, such that
#'        ifixed[i]=TRUE iff param[i] is fixed at the given value
#' @param dstruct  list with  data set and other variables used by objfn
#' @param LB lower bound vector
#' @param UB upper bound vector
#' @param mxiter  max number of iterations
#' @param eps  tolerance for Newton-Raphson iterations
#' @param bdd  bound on difference of 2 consecutive iterations
#'  (useful is starting point is far from solution and func is far from convex)
#' @param iprint  control on amount of printing,
#'          FALSE for no printing of iterations and
#'          TRUE for printing x^{(k)} on each iteration.
#'
#' @return list with
#' fnval = function value at minimum;
#' parmin = param for minimum;
#' invh = inverse Hessian; 
#' iconv = 1 if converged, -1 for a boundary point, 0 otherwise;
#' iter = number of iterations.
#'
#' @examples
#' # See examples in onefactorcop_nllk(), bifactorcop_nllk(), nestfactorcop_nllk()
#'
#' @export
#'
posDefHessMinb=function(param,objfn,ifixed,dstruct,LB,UB,mxiter=30,eps=1.e-6,
  bdd=5,iprint=FALSE) 
{ np = length(param)
  iter = 0; mxdif = 1; brk = 0;
  while(mxdif>eps & iter<mxiter)
  { dstruct$pdf = 0;
    fn = objfn(param,dstruct,iprint);
    fnold = fn$fnval;
    # get subset of grad, hess
    sgrad = fn$grad[!ifixed]; shess = fn$hess[!ifixed,!ifixed]
    evd = eigen(shess,F);
    eigval = evd$values;
    eigvec = evd$vectors;
    ind = (eigval<=0);
    #replace negative eigenvalues with a positive number, keeping eigenvectors
    if(sum(ind)>0 & iprint) message(paste("replacing eigenvalues in iteration ", iter))
    eigval[ind] = 1e-1;    
    newhess = eigvec%*%diag(eigval)%*%t(eigvec);
    sdiff = try(solve(newhess,sgrad),T);
    if(is.vector(sdiff) == FALSE) { brk = 1; break; }
    diff = rep(0,np)
    diff[!ifixed] = sdiff
    mxdif = max(abs(diff))
    iter = iter+1
    param = param-diff
    outb = 0;
    if(any(param<=LB) | any(param>=UB)) outb = 1;
    if(is.infinite(fn$fnval) | is.nan(fn$fnval)) outb = 1;
    if(any(is.infinite(fn$grad)) || any(is.nan(fn$grad))) outb = 1;
    # modified NR to limit size of step to within bounds
    while(outb==1 | mxdif>bdd)
    { mxdif = mxdif/2.;
      outb = 0;
      diff = diff/2; param = param+diff; 
      if(any(param<=LB) | any(param>=UB)) outb = 1;
    }
    # find a step size to ensure the objective fn decreases
    LB.ind = param<=LB +0.008
    UB.ind = param>=UB -0.008
    ifixed[LB.ind] = TRUE
    ifixed[UB.ind] = TRUE
    dstruct$pdf = 1; # positive definite
    fnnew = objfn(param,dstruct,F)$fnval;
    brkvr = 0; 
    while(fnold<=fnnew) 
    { mxdif = mxdif/2.; 
      diff = diff/2; param = param+diff; 
      fnnew = objfn(param,dstruct,F)$fnval; 
      brkvr = brkvr+1;
      if(brkvr>10)
      { message("break: step size is very small"); 
        fnnew = fnold-1e-5; mxdif = eps-1e-5; 
      }   
    }
    if(iprint) { cat(iter, fn$fnval, mxdif,"\n") }
  }
  if(brk==0)
  { iconv = 1;
    if(iter>=mxiter) { message("did not converge\n"); iconv = 0; }
    else { if(max(abs(fn$grad))>eps) iconv = -1 }
    # check for positive definite
    iposdef = isPosDef(shess)
    # inverse Hessian
    invh = solve(shess)
  }
  if(brk==1)
  { iconv = 2; #algorithm failed
    param = rep(1e5,np);
    fnnew = 1e5;
    invh = NA;
    isposdef = FALSE;
  }
  list(parmin=param,fnval=fnnew,invh=invh,iconv=iconv,iposdef=iposdef,iter=iter)
}

