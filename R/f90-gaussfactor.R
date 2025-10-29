# R interface to pfactor/bifactor nllk+grad 
# functions renamed in Feb 2025, documentation in roxygen2 format
# front end like factanal() 

#' log-likelihood Gaussian p-factor structure correlation matrix
#'
#' @description
#' log-likelihood Gaussian p-factor structure correlation matrix with gradient
#'
#' @param rhvec vector of length d*p with partial corr representation of loadings
#' @param Robs dxd empirical correlation matrix
#' @param nsize sample size
#'
#' @return negative log-likelihood and gradient for Gaussian p-factor model
#'
#' @export
#'
pfactor_nllk = function(rhvec,Robs,nsize)
{ d = nrow(Robs)
  if(max(abs(rhvec))>0.999) { return(1.e10) }
  p = length(rhvec)/d
  out = .Fortran("pfactnllk",
      as.integer(d), as.integer(p), as.double(rhvec), as.double(Robs),
      as.integer(nsize),
      nllk=as.double(0), grad=as.double(rep(0,d*p))  )
  nllk=out$nllk; 
  attr(nllk,"gradient") = out$grad;
  nllk
}

#' log-likelihood Gaussian bi-factor structure correlation matrix
#'
#' @description
#' log-likelihood Gaussian bi-factor structure correlation matrix
#'
#' @param rhvec vector of length d*2 for partial correlation representation of loadings,
#   first d correlations with common factor, then

#   partial correlations with group factor given common factor
#' @param grsize vector of group sizes for bi-factor model
#' @param Robs dxd empirical correlation matrix
#' @param nsize sample size
#'
#' @return negative log-likelihood and gradient for Gaussian bi-factor model
#'
#' @export
#'
bifactor_nllk = function(rhvec,grsize,Robs,nsize)
{ d = nrow(Robs)
  if(max(abs(rhvec))>0.999) { return(1.e10) }
  mgrp = length(grsize)
  out = .Fortran("bifactnllk",
      as.integer(d), as.integer(mgrp), as.double(rhvec), as.integer(grsize),
      as.double(Robs), as.integer(nsize),
      nllk=as.double(0), grad=as.double(rep(0,d*2))  )
  nllk=out$nllk; 
  attr(nllk,"gradient") = out$grad;
  nllk
}


#' Gaussian bi-factor structure correlation matrix
#'
#' @description
#' Gaussian bi-factor structure correlation matrix with quasi-Newton
#'
#' @param grsize vector of group sizes for bi-factor model
#' @param start starting point should have dimension 2*d
#' @param data nsize x d data set to compute the correlation matrix if
#'        correlation matrix (cormat) not given
#' @param cormat dxd empirical correlation matrix
#' @param n sample size
#' @param prlevel print.level for nlm()
#' @param mxiter maximum number of iterations for nlm()
#'
#' @return a list with
#'   $nllk, $parmat= dx2 matrix of correlations and partial correlations
#'
#' @examples
#' data(rainstorm)
#' rmat = rainstorm$cormat
#' n = nrow(rainstorm$zprecip)
#' d = ncol(rmat)
#' grsize = rainstorm$grsize
#' fa1 = pfactor_fa(1,start=rep(0.8,d),cormat=rmat,n=n,prlevel=1)
#' fa2 = pfactor_fa(2,start=rep(0.8,2*d),cormat=rmat,n=n,prlevel=1)
#' st3 = c(rep(0.7,grsize[1]),rep(0.1,d),rep(0.7,grsize[2]),rep(0.1,d),rep(0.7,grsize[3]))
#' fa3 = pfactor_fa(3,start=st3,cormat=rmat,n=n,prlevel=1)
#' names(fa3)
#' loadmat_rotated = fa3$loading%*%fa3$rotmat # more values closer to 0
#' # compare factanal
#' fa3b = factanal(factors=3,covmat=rmat)
#' compare = cbind(loadmat_rotated,fa3b$loadings)
#' print(round(compare,3)) # order of factors is different but interpretation similar
#' #
#' bifa = bifactor_fa(grsize,start=c(rep(0.8,d),rep(0.2,d)),cormat=rmat,n=n,prlevel=1)
#' mgrp = length(grsize)
#' # oblique factor model is much more parsimonious than bi-factor
#' obfa = oblique_fa(grsize,start=rep(0.7,d+mgrp), cormat=rmat, n=n, prlevel=1)
#' #
#' cat(fa1$nllk, fa2$nllk, fa3$nllk, bifa$nllk, obfa$nllk,"\n")
#'
#' @export
#'
bifactor_fa = function(grsize,start,data=1,cormat=NULL,n=100,prlevel=0,mxiter=100)
{ if(is.null(cormat))
  { n = nrow(data); d = ncol(data);
    cormat = cor(data)
  }
  else { d = sum(grsize) }
  if(length(start)!=2*d) { message("start should have length 2*d"); return(0) }
  mle = nlm(bifactor_nllk,p=start,grsize=grsize,Robs=cormat,nsize=n,hessian=F,
    iterlim=mxiter,print.level=prlevel,check.analyticals=F)
  list(nllk=mle$minimum, parmat=matrix(mle$estimate,d,2))
}

# front end like factanal() 

#' Gaussian p-factor structure correlation matrix
#'
#' @description
#' Gaussian p-factor structure correlation matrix with quasi-Newton
#'
#' @param factors p = #factors
#' @param start starting point should have dimension 2*d
#' @param data nsize x d data set to compute the correlation matrix if
#'       correlation matrix (cormat) not given
#' @param cormat dxd empirical correlation matrix
#' @param n sample size
#' @param prlevel print.level for nlm()
#' @param mxiter maximum number of iterations for nlm()
#'
#' @return a list with
#'        $nllk, $rhmat = dxp matrix of partial correlations,
#'        $loading = dxp loading matrix after varimax,
#'        $rotmat = pxp rotation matrix used by varimax
#'
#' @examples
#' # See example in bifactor_fa()
#'
#' @export
#'
pfactor_fa = function(factors,start,data=1,cormat=NULL,n=100,prlevel=0,mxiter=100)
{ if(is.null(cormat))
  { n = nrow(data); d = ncol(data);
    cormat = cor(data)
  }
  else { d = nrow(cormat) }
  p = factors
  if(length(start)!=p*d) 
  { message("start should have length factors*d"); return(0) }
  mle = nlm(pfactor_nllk,p=start,Robs=cormat,nsize=n,hessian=F,
    iterlim=mxiter,print.level=prlevel,check.analyticals=F)
  rhmat = matrix(mle$estimate,d,p)
  if(p>=2)
  { amat = pcor2load(rhmat)
    rotat = varimax(amat)
    loading = as.matrix(rotat$loadings[,1:p])
    rotmat = rotat$rotmat
  }
  else # count #negative loadings and reverse sign if needed
  { nneg = sum(rhmat<0)
    if(nneg>d/2) { rotmat = -1; loading = -rhmat } else { rotmat = 1; loading = rhmat }
  }
  list(nllk=mle$minimum, parmat=rhmat, loading=loading, rotmat=rotmat)
}

#============================================================

