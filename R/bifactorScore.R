
# version written by HJoe, replacing version of Xinyao Fan

#' Proxies for bi-factor copula model based on Gaussian bi-factor score
#'
#' @description 
#' Proxies (in (0,1)) for bi-factor copula model based on Gaussian bi-factor score
#'
#' @param udata  nxd matrix in (0,1); n is sample size, d is dimension
#' @param start  starting values for fitting the bi-factor Gaussian model
#' @param grsize G-vector with group sizes with G groups
#' @param prlev printlevel in call to nlm
#'
#' @return list with  Aloadmat=estimated loading matrix after N(0,1) transform; 
#'  proxies = nx(G+1) matrix with stage 1 proxies for the latent variables
#' (for global latent in column 1, and then for group latent);
#'  weight = weight matrix based on the correlation matrix of normal score and the loading matrix.
#'
#' @examples
#' # See example in bifactorEstWithProxy()
#'
#' @export
#'
bifactorScore = function(udata, start, grsize, prlev=1)
{ mgrp = length(grsize)
  d = sum(grsize)
  zdata = qnorm(udata)
  Robs = cor(zdata)
  n = nrow(zdata)
  bifact_obj = bifactor_fa(grsize,start,cormat=Robs,n=n,prlevel=prlev)
  bifact_est = bifact_obj$parmat
  pcmat = matrix(0,d,mgrp+1) 
  pcmat[,1] =  bifact_est[,1]
  iend = cumsum(grsize)
  ibeg = iend+1; ibeg= c(1,1+iend[-mgrp])
  for(g in 1:mgrp)
  { pcmat[ibeg[g]:iend[g],g+1] = bifact_est[(ibeg[g]:iend[g]),2] }
  aload = pcor2load(pcmat)
  psi2 = 1- apply(aload^2,1,sum)
  rmat = aload%*%t(aload)
  diag(rmat) = 1
  psi2inv = diag(1/psi2)
  apa = t(aload)%*%psi2inv%*%aload  # A^T %*% Psi2^{-1} %*% A
  values = eigen(apa)$values
  cond = max(values)/min(values)
  if(prlev>0) message(paste("condition number of A^T %*% Psi2^{-1} %*% A  matrix", cond))
  # factor score depends on zdata and aload loading matrix
  wt = solve(rmat,aload) 
  zfs = zdata %*% wt
  ufs = uscore(zfs)
  list(Aloadmat=aload, proxies=ufs, weight=wt)
}

