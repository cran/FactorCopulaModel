# Functions for semi-correlations

#' Semi-correlation for bivariate normal/Gaussian distribution
#'
#' @description
#' semicorrelation assuming bivariate normal/Gaussian copula 
#' 
#' @param rho correlation in (-1,1)
#' 
#' @return Cor(Z1,Z2| Z1>0,Z2>0) when (Z1,Z2)~bivariate standard normal(rho)
#' 
#' @references 
#' Joe (2014), Dependence Modeling with Copulas, Chapman&Hall/CRC; p 71
#' 
#' @export
#' 
bvnSemiCor = function(rho)
{ bp = .25+asin(rho)/(2*pi)
  r1 = sqrt(1-rho^2)
  denom = 2*bp*sqrt(2*pi)
  v10 = (1+rho)/denom
  v20 = 1+rho*r1/(2*bp*pi)
  v11 = rho+r1/(2*bp*pi)
  sc = (v11-v10^2)/(v20-v10^2)
  sc
}

#' Semi-correlations for two variables
#'
#' @description
#' semi-correlations (lower and upper) applied to data after normal scores transform
#' 
#' @param bivdat  nx2 data set
#' @param inscore  TRUE if bivdat has already been converted to normal scores, default FALSE
#' 
#' @return 3-vector with rhoN = correlation of normal scores (vander Waerden correlation) and lower/ upper semi-correlations
#' 
#' @examples
#' # See example in semiCorTable()
#' 
#' @export
#' 
semiCor = function(bivdat,inscore=FALSE)
{ if(!inscore) bivdat=nscore(bivdat)
  ii1 = (bivdat[,1]<0 & bivdat[,2]<0)
  ii2 = (bivdat[,1]>0 & bivdat[,2]>0)
  lcor = cor(bivdat[ii1,1],bivdat[ii1,2])
  ucor = cor(bivdat[ii2,1],bivdat[ii2,2])
  ncor = cor(bivdat[,1],bivdat[,2])
  out = c(ncor,lcor,ucor)
  names(out) = c("ncor","lcor","ucor")
  out
}

# Revised to include variable names in columns
# to do: matrix with upper semicorr in uper triangle of matrix etc

#' Semi-correlation table for a multivariate data set
#'
#' @description
#' Semi-correlation table for several variables
#' 
#' @param mdat nxd multivariate data set with d>=2 columns
#' @param varnames d-vector of (abbreviated) variable names 
#' @param inscore TRUE if mdat has already been converted to normal scores, default FALSE
#' 
#' @return d*(d-1)/2 by 8-column dataframe with columns j1,j2,ncor,lcor,ucor,bvnsemic, varnames[j1] varnames[j2]
#' for 2 variable indices, correlation of normal scores, lower semi-correlation,
#' upper semi-correlation and BVN semi-correlation assuming Gaussian copula
#' Stronger than Gaussian dependence in upper tail if ucor is larger than bvn semicor
#'
#' @examples
#' rmat = toeplitz(c(1,.7,.4))
#' print(rmat)
#' set.seed(1234)
#' zdat = rmvn(n=500,rmat)
#' set.seed(12345)
#' tdat = rmvt(n=500,rmat,nu=5)
#' vnames=c("V1","V2","V3")
#' zsemi = semiCorTable(zdat,vnames)
#' tsemi = semiCorTable(tdat,vnames)
#' cat("trivariate normal\n")
#' print(zsemi)
#' cat("trivariate t(5)\n")
#' print(tsemi)
#' 
#' # data example
#' data(euro07gf)
#' zdat = euro07gf$zscore
#' options(digits=3)
#' tab = semiCorTable(zdat, varnames=colnames(zdat),inscore=TRUE)
#' print(tab)
#'
#' @export
#' 
semiCorTable = function(mdat, varnames, inscore=FALSE)
{ if(!inscore) mdat = nscore(mdat)
  d = ncol(mdat)
  d2 = (d*(d-1))/2
  tab = matrix(0,d2,6)
  ii = 0
  for(j2 in 2:d)
  { for(j1 in 1:(j2-1))
    { scs = semiCor(mdat[,c(j1,j2)])
      bvnsc = bvnSemiCor(scs[1])
      ii = ii+1
      tab[ii,] = c(j1,j2,scs, bvnsc)
    }
  }
  tab = as.data.frame(tab)
  # revise out
  names(tab) = c("j1","j2","ncor","lcor","ucor","bvnsemic")
  tab$var1 = varnames[tab$j1]
  tab$var2 = varnames[tab$j2]
  tab
}


#' Random multivariate normal (standard N(0,1) margins)
#'
#' @description
#' Random multivariate normal (standard N(0,1) margins)
#'
#' @param n simulation sample size
#' @param rmat correlation matrix 
#'
#' @return nxd matrix, where d=nrow(rmat)
#' 
#' @export
#' 
rmvn = function(n,rmat)
{ if(any(diag(rmat)!=1)) return(NA)
  Amat = try(chol(rmat), silent=TRUE)
  if(class(Amat)[1]!="matrix") return(NA)
  d = nrow(rmat)
  z = matrix(rnorm(n*d),n,d)
  z%*%Amat
}

#' Random multivariate t (standard t(nu) margins)
#'
#' @description
#' Random multivariate t (standard t(nu) margins)
#'
#' @param n simulation sample size
#' @param rmat correlation matrix 
#' @param nu degree of freedom parameter 
#'
#' @return nxd matrix, where d=nrow(rmat)
#' 
#' @export
#' 
rmvt = function(n,rmat,nu)
{ if(any(diag(rmat)!=1)) return(NA)
  Amat = try(chol(rmat), silent=TRUE)
  if(class(Amat)[1]!="matrix") return(NA)
  d = nrow(rmat)
  z = matrix(rnorm(n*d),n,d)
  z = z%*%Amat
  w = rchisq(n,nu)
  w = sqrt(w/nu)
  z/w
}


#======================================================================
