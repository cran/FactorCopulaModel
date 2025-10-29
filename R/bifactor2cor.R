# Functions for bi-factor correlation matrices
# and maximum likelihood for multivariate t with these structures
# Code written by Pavel Krupskii, renamed in February 2025

#' Bi-factor partial correlations to correlation matrix 
#'
#' @description
#' Bi-factor partial correlations to correlation matrix, determinant, inverse
#'
#' @param grsize  vector with group sizes: d_1,d_2,...,d_G for G groups
#' @param rh1 vector of length sum(grsize) of correlation with global latent variable, ordered by group index
#' @param rh2 vector of length sum(grsize) of partial correlation with group latent variable given global
#'
#' @return list with Rmat = correlation matrix; det = det(Rmat); 
#'   Rinv = solve(Rmat)
#'
#' @examples
#' grsize = c(5,5,3) 
#' d = sum(grsize)
#' bifpar = c(0.84,0.63,0.58,0.78,0.79, 0.87,0.80,0.74,0.71,0.57, 0.83,0.77,0.80,
#'   0.67,0.58,0.15,0.70,0.47,   0.32,0.27,0.73,0.19,0.12,   0.35,0.23,0.53)
#' bifobj = bifactor2cor(grsize,bifpar[1:d],bifpar[(d+1):(2*d)])
#' rmat = bifobj$Rmat
#' print(det(rmat)-bifobj$det)
#' print(max(abs(solve(rmat)-bifobj$Rinv)))
#' bifobj2 = bifactor2cor_v2(grsize,bifpar[1:d],bifpar[(d+1):(2*d)])
#' rmat2 = bifobj2$Rmat
#' print(det(rmat2)-bifobj2$det)
#' print(max(abs(solve(rmat2)-bifobj2$Rinv)))
#'
#' @export
#'
bifactor2cor = function(grsize,rh1,rh2)
{ mgrp = length(grsize); # mgrp=G=#groups
  d = sum(grsize);
  a2 = rh2*sqrt(1-rh1^2);  # loadings for factor 2
  A = rep(0,mgrp); Aj = rep(0,d);
  eta = rep(0,mgrp); etaa = rep(0,d);
  dzeta1 = rh1^2/(1-rh2^2)/(1-rh1^2); 
  dzeta2 = a2/(1-rh2^2)/(1-rh1^2);
  st = 0;
  fctdiag = matrix(0,nrow=d,ncol=d);

  for(j in 1:mgrp)
  { ind = (st+1):(st+grsize[j]);
    fctdiag[ind,ind] = 1;
    A[j] = 1-grsize[j]+sum(1/(1-rh2[ind]^2));
    Aj[ind] = A[j];
    eta[j] = sum((dzeta2[ind])*(rh1[ind]));  
    etaa[ind] = eta[j]/A[j];     
    st = st+grsize[j];
  }
  fctmat = outer(rh1,rh1)+outer(a2,a2)*fctdiag;
  diag(fctmat) = 1;
  xi = dzeta1/rh1-etaa*dzeta2;
  A0 = 1-sum(eta^2/A)+sum(dzeta1);
  fctdet = A0*prod(A)*prod((1-rh1^2)*(1-rh2^2));
  fctinv = -(1/A0)*outer(xi,xi) - (1/Aj)*outer(dzeta2,dzeta2)*fctdiag;
  diag(fctinv) = 1/(1-rh1^2)/(1-rh2^2) - dzeta2^2/Aj - xi^2/A0;
  list(Rmat=fctmat,det=fctdet,Rinv=fctinv)
}  
    

#' Bi-factor partial correlations to correlation matrix 
#' version 2, using the inverse and determinant of a smaller matrix
#'
#' @description
#' Bi-factor partial correlations to correlation matrix, determinant, inverse
#'
#' @param grsize  vector with group sizes: d_1,d_2,...,d_G for G groups
#' @param rh1 vector of length sum(grsize) of correlation with global latent variable, ordered by group index
#' @param rh2 vector of length sum(grsize) of partial correlation with group latent variable given global
#'
#' @return list with Rmat = correlation matrix; det = det(Rmat); 
#'   Rinv = solve(Rmat)
#'
#' @examples
#' # see examples for bifactor2cor()
#'
#' @export
#'
bifactor2cor_v2 = function(grsize,rh1,rh2)
{ matp = rh1;
  a2 = rh2*sqrt(1-rh1^2);  # loadings for factor 2
  mgrp = length(grsize);
  dvar = sum(grsize);
  st = 0;
  for(jg in 1:mgrp)
  { tem = rep(0,dvar);
    tem.ind = (st+1):(st+grsize[jg]);
    tem[tem.ind] = a2[tem.ind];
    matp = cbind(matp,tem);
    st = st+grsize[jg];
  }
  vrho = 1-rh1^2-a2^2;
  invrho = 1/vrho;
  invrho = diag(invrho);
  prmat = invrho%*%matp;
  nmat = t(matp)%*%prmat+diag(rep(1,mgrp+1));  
  fctmat = matp%*%t(matp) + diag(vrho);
  fctdet = det(nmat)*prod(vrho);
  fctinv = invrho - prmat%*%solve(nmat,t(prmat));
  list(Rmat=fctmat,det=fctdet,Rinv=fctinv)
}


