# Tail-weighted dependence measure [zeta(alpha), alpha>0 large].
# Lee D, Joe H and Krupskii P (2018).
# A tail-weighted dependence measure with limit being tail dependence coefficient. 
# J Nonparametric Statistics, 30(2), 262-290.
# http://dx.doi.org/10.1080/10485252.2017.1407414


# Regression limit to get estimates of upper tail dependence parameter.
# Initial version of code written by David Lee , for the method described
# in the above paper. 
# Subsequently there are other methods based on tail expansions of copulas
# to estimate multivariate quantile curves at the same time
# as tail dependence or tail order.


#' Tail dependence parameter estimation
#'
#' @description
#' Tail dependence parameter estimate based on extrapolating zeta(alpha)
#'
#' @param u1 nx1 vector with values (in (0,1) if rank=FALSE)
#' @param u2 nx1 vector with values (in (0,1) if rank=FALSE)
#' @param lowertail TRUE if lower tail-weighted dependence measure, default is FALSE
#' @param eps tolerance (default 0.1) for rate parameter in method 2; non-linear (use method 2) if rate < 1-eps
#' @param semictol tolerance (default 0.1) for exceedance of normal semicorrelation to treat as tail dependent; 
#'  use something like semicol=-0.5 if normal scores plot suggest tail dependence 
#' @param rank TRUE (default) if data matrix needs to be converted to uniform scores in (0,1)
#' @param iprint TRUE for intermediate prints
#'
#' @return upper tail dependence parameter based on extrapolation of zeta(alpha) for large alpha
#'
#' @examples
#' mytest = function(qcond,cpar,n=500,seed=123)
#' { set.seed(seed)
#'   u1 = runif(n)
#'   u2 = qcond(runif(n),u1,cpar)
#'   #convert to uniform scores (marginals are usually not known)
#'   u1 = (rank(u1)-0.5)/n
#'   u2 = (rank(u2)-0.5)/n
#'   alp = c(1,5,10:20)
#'   zetaL = zetaDep(cbind(u1,u2),alp,rank=FALSE,lowertail=FALSE)
#'   zetaU = zetaDep(cbind(u1,u2),alp,rank=FALSE,lowertail=TRUE)
#'   print(cbind(alp,zetaL,zetaU))
#'   utd = tailDep(u1,u2, lowertail=FALSE, eps=0.1, semictol=0.1, rank=FALSE, iprint=TRUE)
#'   ltd = tailDep(u1,u2, lowertail=T, eps=0.1, semictol=0.1, rank=FALSE, iprint=TRUE)
#'   cat(ltd,utd,"\n")
#'   utd = tailDep(u1,u2, lowertail=FALSE, eps=0.1, semictol=0.1, rank=FALSE, iprint=FALSE)
#'   ltd = tailDep(u1,u2, lowertail=T, eps=0.1, semictol=0.1, rank=FALSE, iprint=FALSE)
#'   cat(ltd,utd,"\n")
#'   oldpar = par(no.readonly = TRUE)
#'   on.exit(par(oldpar)) 
#'   par(mfrow=c(2,1))
#'   zetaPlot(cbind(u1,u2),alp,ylim=c(0,1),inverse=FALSE)
#'   zetaPlot(cbind(u1,u2),alp,ylim=c(0,1),inverse=TRUE)
#'   0
#' }
#' mytest(qcondFrank,3)
#' mytest(qcondbvtcop,c(0.6,5))
#'
#' @references
#' Lee D, Joe H, Krupskii P (2018). J Nonparametric Statistics, 30(2), 262-290
#'
#' @export
#'
tailDep = function(u1,u2, lowertail=FALSE, eps=0.1, semictol=0.1, rank=TRUE,
  iprint=FALSE)
{ 
  if(rank)
  { n = length(u1)
    u1 = (rank(u1)-0.5)/n; u2 = (rank(u2)-0.5)/n
  }
  if(lowertail) { u1 = 1-u1; u2 = 1-u2 }  # reflect
  a = 10:20 
  y = rep(0,length(a)) 
  for(j in 1:length(a))
  { y[j] = 2+a[j]*(1-(a[j]/(a[j]+1)-mean(0.5*abs(u1^a[j]-u2^a[j])))^(-1)) }
  # check that above is same as zetaDep(cbind(u1,u2),a[j],rank=FALSE,lowertail)
  #if(iprint) print(cbind(a,y))
  # try different non-linear regressions
  ee1 = nlm(obj1, c(0,1), mat=cbind(a,y), w=a)
  ee2 = nlm(obj2, c(0,1,0.9), mat=cbind(a,y), w=a^0.5)
  ee3 = nlm(obj3, 0.8, mat=cbind(a,y), w=a)
  est1 = ee1$estimate;
  est2 = ee2$estimate;
  est3 = ee3$estimate; # bhat = 2*phat to match publication
  if(iprint)
  { message(paste("method1 coeff ", est1))
    message(paste("method2 coeff ", est2))
    message(paste("method3 coeff ", est3))
  }
  # could add a plot here
  # Decide which estimate/method to use (zeta)
  ainv= 1/a; lmobj= lm(y~ainv); coeff= summary(lmobj)$coefficients
  semic= semiCor(nscore(cbind(u1,u2)))
  semicor= semic[3]; scnorm= bvnSemiCor(semic[1])
  if(iprint & !lowertail) message(paste("ucor, bvnsemic:", semicor, scnorm))
  if(iprint & lowertail) message(paste("lcor, bvnsemic:", semicor, scnorm))
  if(coeff[2,1]<0)
  { #y vs 1/a is decreasing, i.e. increasing to limit
    lambda = 2-2*est3; method= 3 # method 3 in paper
  } 
  else 
  { if(ee2$estimate[3] >= 1-eps | semicor-scnorm > semictol)
    { #can assume linearly increasing, method 1
      lambda = est1[1]; method= 1
    } 
    else 
    { #curved increasing, method 2
      lambda = est2[1]; method = 2
    }
  }
  if(iprint & !lowertail) message(paste("lambdaU: ", method, lambda))
  if(iprint & lowertail) message(paste("lambdaL: ", method, lambda))
  lambda
}

#' @description
#' Regression Method 1 in Lee at al 2018
#'
#' @param b 2-dimensional vector for WLS zetahat ~ b1 + b2/alpha
#' @param mat matrix with alpha values in column 1 and zetas in column 2
#' @param w weight for weighted least squares
#'
#' @return WLS objective
#'
#' not exported
#'
obj1 = function(b,mat,w)
{ if(b[1]<0 | b[1]>1) { return(1e100) } 
  return(mean(((mat[,2]-b[2]/mat[,1]-b[1])^2)/w))
}

# power b[3] is estimated (regression Method 2 in Lee at al 2018)

#' @description
#' Regression method 2 in Lee at al 2018
#'
#' @param b 3-dimensional vector for WLS zetahat ~ b1 + b2/alpha^b3
#' @param mat matrix with alpha values in column 1 and zetas in column 2
#' @param w weight for weighted least squares
#'
#' @return WLS objective
#'
#' not exported
#'
obj2 = function(b,mat,w)
{ if(b[3]<0 | b[3]>1 | b[1]<0) { return(1e100) } 
  return(mean(((mat[,2] - b[2]/(mat[,1])^b[3] - b[1])^2)/w))
}

# increasing measures => only need to estimate 1 parameter

#' @description
#' Regression method 3 in Lee at al 2018
#'
#' @param p scalar for WLS zetahat ~ (2-2*p)+(2*p-4*p*p)/(alpha+1-2*p)
#' @param mat matrix with alpha values in column 1 and zetas in column 2
#' @param w weight for weighted least squares
#'
#' @return WLS objective, b=2*p in Lee at al 2018
#'
#' not exported
#'
obj3 = function(p,mat,w) 
{ if(p<0.5 | p>1) { return(1e100) }
  return(mean(((mat[,2]-(2*p-4*p*p)/(mat[,1]+1-2*p)-(2-2*p))^2)/w))
}


#' Plot zeta(alpha) against alpha
#'
#' @description
#' Plot zeta(alpha) against alpha 
#'
#' @param dat nx2 data matrix with values (u-data in (0,1))
#' @param alpha vector of alpha>0 for zeta measure
#' @param ylim limits for yaxis to pass to plot
#' @param inverse if TRUE, plot zeta against 1/alpha
#'
#' @return nothing is returned, but a plot is produced
#'
#' @export
#'
zetaPlot = function(dat,alpha,ylim=c(0,1),inverse=FALSE)
{ y = rep(NA,length(alpha)); 
  u1 = dat[,1]; u2 = dat[,2]
  a = alpha
  for(j in 1:length(a)) 
  { y[j] = 2+a[j]*(1-(a[j]/(a[j]+1)-mean(.5*abs(u1^a[j]-u2^a[j])))^(-1)) }
  if(!inverse) plot(alpha,y,xlab="alpha",ylab="zeta(alpha)",ylim=ylim)
  if(inverse) plot(1/alpha,y,xlab="1/alpha",ylab="zeta(alpha)",ylim=ylim)
  invisible(0)
}



#======================================================================

#' C_[2|1]^[-1](p|u) for bivariate Student t copula
#'
#' @description
#' bivariate Student t copula conditional quantile
#'
#' @param p 0<p<1, could be a vector
#' @param u 0<u<1, could be a vector 
#' @param cpar  copula parameter: 2-vector with -1<rho<1, df>0
#'
#' @return conditional quantiles of bivariate Student t copula
#'
#' @export
#'
qcondbvtcop = function(p,u,cpar)
{ rho=cpar[1]; df=cpar[2] 
  x2=qt(p,df+1)
  x1=qt(u,df)
  tem=x2*sqrt((df+x1*x1)*(1-rho*rho)/(df+1))+rho*x1
  pt(tem,df)
}

#' C_[2|1]^[-1](p|u) for bivariate Frank copula
#'
#' @description
#' Frank bivariate copula conditional quantile
#'
#' @param p 0<p<1, could be a vector
#' @param u 0<u<1, could be a vector 
#' @param cpar  copula parameter: cpar>0 or cpar<0; cpar=0 input will not work
#'
#' @return conditional quantiles of bivariate Frank copula
#'
#' @details
#'  1-exp(-cpar) becomes 1 in double precision for cpar>37.4;
#'  any argument can be a vector, but all vectors must have same length.
#'  Form of inputs not checked (for readability of code).
#'
#' @export
#'
qcondFrank = function(p,u,cpar)
{ #if(cpar==0) return(p)
  cpar0 = exp(-cpar)
  cpar1 = 1-cpar0
  etem = exp(-cpar*u+log(1./p-1.))   
  tem = 1.-cpar1/(etem+1.);
  v = (-log(tem))/cpar
  isinf = is.infinite(v)
  #print(cbind(v[isinf],tem[isinf],etem[isinf]))
  # v Inf, tem is 0 and etem < 1.e-16
  v[isinf] = (-log(cpar0+etem[isinf]))/cpar
  v
} 

