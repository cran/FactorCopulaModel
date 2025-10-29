# Front ends to the nllk functions, original version written by Pavel Krupskii

# multivariate t with p-factor and bi-factor correlation structure,
# nested factor as special case of bi-factor
# R interface to tpfactnllk/tbifactnllk nllk+grad 

# bi-factor, nested factor and p-factor

#' negative log-likelihood for the bi-factor Gaussian/t model
#'
#' @description 
#' negative log-likelihood in the bi-factor Gaussian or t model
#'
#' @param rhovec vector of length d*2 for partial correlation representation of loadings,
#'   first d correlations with common factor, then
#'   partial correlations with group factor given common factor
#' @param grsize vector of group sizes for bi-factor model
#' @param tdata nxd data set of scores for t df (e.g., qt(u,df))
#' @param df degree of freedom parameter (positive)
#'
#' @return negative log-likelihood (nllk) of copula for mvt bi-factor model
#'
#' @export
#'
mvtBifact_nllk = function(rhovec,grsize,tdata,df)
{ d = ncol(tdata)
  n = nrow(tdata)
  if(max(abs(rhovec))>0.999) { return(1.e10) }
  mgrp = length(grsize)
  if(df< 300)
  { out = .Fortran("tbifactnllk",
         as.integer(d), as.integer(mgrp), as.integer(grsize),
         as.integer(n), as.double(df),  as.double(rhovec), as.double(tdata),
         nllk=as.double(0), grad=as.double(rep(0,d*2)))
  }
  else # df>=300
  { Robs = cor(tdata)
    out = .Fortran("bifactnllk",
         as.integer(d), as.integer(mgrp), as.double(rhovec), as.integer(grsize),
         as.double(Robs), as.integer(n),
         nllk=as.double(0), grad=as.double(rep(0,d*2)))
    lgdenom = -0.5*sum(tdata^2) -n*d*0.5*log(2*pi)
    out$nllk = out$nllk + lgdenom
  }

  list(nllk=out$nllk,lgrad=out$grad)
}

#' MLE for multivariate normal/t with a bi-factor or nested factor correlation structure
#'
#' @description  
#' MLE for the bi-factor or nested factor structure for multivariate normal/t
#'
#' @param tdata nxd matrix of t-scores or z-scores
#' @param start vector of length 2*d with starting values of partial correlations
#'    values for correlations of observed Z_[gj] and common latent V_0 go first,
#'   then partial correlations of Z_[gj] and V_g given V_0 (j in group g)
#' @param grsize vector of group sizes for bi-factor model
#' @param df degrees of freedom parameter >0
#' @param prlevel print.level for nlm()
#' @param model "bifactor" or "nestfactor"
#' nested-factor is reduced model with fewer parameters
#' @param mxiter maximum number of iterations for nlm()
#'
#' @return nlm object with ($code,$estimate,$gradient,$iterations,$minimum)
#'
#' @examples
#' data(rainstorm)
#' udat = rainstorm$uprecip
#' d = ncol(udat)
#' grsize = rainstorm$grsize
#' df = 10
#' tdata = qt(udat,df)
#' bif = mvtBifact(tdata, c(rep(0.8,d),rep(0.2,d)), grsize, df=df,
#'  prlevel=1, model="bifactor", mxiter=100)
#' #
#' # nested-factor: parameters for group latent linked to global latent
#' # come in the first tree of the 2-truncated vine.
#' nestf = mvtBifact(tdata, c(0.7,0.7,0.7,rep(0.85,d)), grsize, df=df,
#' prlevel=1, model="nestfactor", mxiter=100)
#' # doesn't converge properly, group2 latent matches global latent
#' #
#' st1 = rep(0.7,d)
#' out1t = mvtPfact(tdata,st1,pfact=1,df=df,prlevel=1)
#' st2 = rep(0.7,2*d)
#' out2t = mvtPfact(tdata,st2,pfact=2,df=df,prlevel=1)
#' st3 = c(rep(0.7,grsize[1]),rep(0.1,d),rep(0.7,grsize[2]),rep(0.1,d),rep(0.7,grsize[3]))
#' out3t = mvtPfact(tdata,st3,pfact=3,df=df,prlevel=1)
#' cat(bif$minimum, nestf$minimum, out1t$minimum, out2t$minimum, out3t$minimum,"\n")
#'
#' @details
#' Note the minimum nllk can be the same for different parameter vectors
#'   if some group size values are 1 or 2.
#'
#' @export
#'
mvtBifact = function(tdata,start,grsize,df,prlevel=2,model="bifactor",mxiter=100)
{ d = ncol(tdata);
  n = nrow(tdata);
  #if(full==TRUE & length(start)!=2*d)
  if(model=="bifactor" & length(start)!=2*d)
  { message("start should have length 2*d for bi-factor"); return(0) }
  #if(full==FALSE & length(start)!=d+length(grsize))
  if(model=="nestfactor" & length(start)!=d+length(grsize))
  { message("start should have length d+length(grsize) for nested-factor");
    return(0)
  }
  nllkfn = function(param)
  { if(max(abs(param))>0.999) { return(1.e10) }
    param0 = param;
    #if(full==FALSE) # nested factor model
    if(model=="nestfactor") # nested factor model
    { mgrp = length(grsize);
      ind = 0;
      a2 = param[(mgrp+1):(mgrp+d)]; a1 = rep(0,d);
      for (jg in 1:mgrp)
      { ind1 = ind+1;  ind2 = ind+grsize[jg];
        a1[ind1:ind2] = param[jg];
        ind = ind+grsize[jg];
      }
      rh1 = a1*a2;
      rh2 = a2*sqrt(1-a1^2)/sqrt(1-rh1^2);
      param0 = c(rh1,rh2);
    }
    #out = bifcttlik.vec(tdata,grsize,param0,df);
    out = mvtBifact_nllk(param0,grsize,tdata,df);
    nllk = out$nllk;
    lgrd = out$lgrad;
    #if(full==TRUE) attr(nllk,"gradient")=lgrd;
    if(model=="bifactor") attr(nllk,"gradient") = lgrd;
    #if(full==FALSE)
    if(model=="nestfactor")
    { lgrd0 = rep(0,mgrp+d);
      lgrad1 = lgrd[1:d]; lgrad2 = lgrd[(d+1):(2*d)];
      b1 = (1-(a1*a2)^2)^1.5; b2 = sqrt(1-a1^2);
      tem1 = lgrad1*a2-lgrad2*(1-a2^2)*a1*a2/b1/b2;
      ind0 = 0;
      for (jg in 1:mgrp)
      { ind1 = 1+ind0; ind2 = grsize[jg]+ind0;
        lgrd0[jg] = sum(tem1[ind1:ind2]);
        ind0 = ind0+grsize[jg];
      }
      lgrd0[(mgrp+1):(mgrp+d)] = lgrad1*a1+lgrad2*b2/b1;
      attr(nllk,"gradient") = lgrd0;
    }
    nllk
  }
  # *** should hessian=T be set so that hessian is returned?
  mle = try(nlm(nllkfn,p=start,hessian=FALSE,iterlim=mxiter,print.level=prlevel,
             check.analyticals = FALSE),TRUE)
  if(is.list(mle)==FALSE) { message("mvtBifact() has failed"); return(NULL); }
  bifmle = mle$estimate
  #cat("MLE: \n")
  #if(full==TRUE) bifmle = matrix(bifmle,ncol=2);
  if(model=="bifactor") bifmle = matrix(bifmle,ncol=2);
  #print(bifmle)
  #cat("SEs: \n")
  #if (max(abs(pmle)) < 0.995 & rcond(mle$hessian) > 1e-12)
  #{ acov = solve(mle$hessian); }
  #else {print("NA")}
  #cat("nllk:  \n"); print(mle$minimum)
  mle
}

#======================================================================

#' negative log-likelihood for the p-factor Gaussian/t model
#'
#' @description 
#' negative log-likelihood in the p-factor Gaussian or t model
#'
#' @param rhvec vector of length d*p with partial correlation representation of loadings
#' @param tdata nxd data set of t(df)-scores
#' @param df degree of freedom parameter >0
#'
#' @return negative log-likelihood of copula for mvt p-factor model
#'
#' @export
#'
mvtPfact_nllk = function(rhvec,tdata,df)
{ d = ncol(tdata)
  n = nrow(tdata)
  if(max(abs(rhvec))>0.999) { return(1.e10) }
  p = length(rhvec)/d
  if(df< 300)
  { out = .Fortran("tpfactnllk",
      as.integer(d), as.integer(p),  as.integer(n), as.double(df),
      as.double(rhvec), as.double(tdata),
      nllk=as.double(0), grad=as.double(rep(0,d*p))  )
  }
  else # df>=300
  { Robs = cor(tdata)  
    out = .Fortran("pfactnllk",
      as.integer(d), as.integer(p), as.double(rhvec), as.double(Robs),
      as.integer(n),
      nllk=as.double(0), grad=as.double(rep(0,d*p))  )
    lgdenom = -.5*sum(tdata^2)-n*d*.5*log(2*pi)
    out$nllk = out$nllk + lgdenom
  }
  #nllk=out$nllk; 
  #attr(nllk,"gradient") = out$grad;
  #nllk
  list(nllk=out$nllk,lgrad=out$grad)
}


#' MLE in a MVt model with a p-factor correlation structure
#'
#' @description
#' MLE in a MVt model with a p-factor correlation structure
#'
#' @param tdata nxd matrix of t-scores
#' @param start vector of length p*d with starting values for
#' partial correlation representation of loadings (algberaically independent)
#' @param pfact number of factors (such as 1,2,3,... )
#' @param df degree of freedom parameter >0
#' @param prlevel  print.level for nlm()
#' @param mxiter maximum number of iterations for nlm()
#'
#' @return  nlm object with ($code,$estimate,$gradient,$iterations,$minimum) 
#' Note the minimum nllk can be the same for different parameter vectors
#'   because of invariance of the loading matrix to rotations
#'
#' @examples
#' # See example in mvtBifact()
#'
#' @export
#'
mvtPfact = function(tdata,start,pfact,df,prlevel=0,mxiter=100)
{ d = ncol(tdata);
  n = nrow(tdata);
  if(length(start) != pfact*d) { message("incorrect number of parameters in start"); return(NULL); }
  nllkfn = function(param)
  { if(max(abs(param))>0.999) { return(1.e10) }
    out = mvtPfact_nllk(param,tdata,df);
    nllk = out$nllk; 
    attr(nllk,"gradient")=out$lgrad;
    nllk
  } 

  # *** should hessian=T be set so that hessian is returned?
  mle = try(nlm(nllkfn,p=start,hessian=F,iterlim=mxiter,print.level=prlevel,
       check.analyticals = F),T) 
  if(is.list(mle)==F) { message("mvtPfact() has failed\n"); return(NULL); }
  pmle = mle$estimate
  #cat("MLE: \n")
  pmle = matrix(pmle,ncol=pfact);
  #print(pmle)
  #cat("SEs: \n")
  #if (max(abs(pmle)) < 0.995 & rcond(mle$hessian) > 1e-12) 
  #{ acov=solve(mle$hessian); } 
  #else {print("NA")}
  #cat("nllk:  \n"); cat(mle$minimum,"\n")
  mle
}  

