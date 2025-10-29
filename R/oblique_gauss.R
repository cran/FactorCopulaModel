
# First draft of R functions below were written by Xinyao Fan
# Multivariate Gaussian with structured correlation matrix that is
# oblique factor model

#' oblique factor correlation structure for d variables and m groups 
#'
#' @description
#' For oblique factor correlation structure for d variables and m groups, 
#' convert the vector of parameters theta into
#' the loading matrix and correlation matrix of latent variables
#' The variables are assumed ordered by group.
#'
#' @param theta vector of length d + m*(m-1)/2; d loading parameters followed by
#'      m*(m-1)/2 entries in correlation matrix of latent variables
#'      (lower triangle by row)
#' @param grsize vector of group sizes (variables ordered by group)
#'
#' @return loadings: loading matrix; cor_lat: correlation matrix of latent variables; 
#'  Rmod: correlation matrix based on theta for oblique factor model
#'
#' @examples 
#' theta = c(0.6,0.7,0.8,0.7,0.6,0.5,0.5)
#' oblique_par2load(theta,grsize=c(3,3))
#' #' $loadings
#' #[,1] [,2]
#' #[1,]  0.6  0.0
#' #[2,]  0.7  0.0
#' #[3,]  0.8  0.0
#' #[4,]  0.0  0.7
#' #[5,]  0.0  0.6
#' #[6,]  0.0  0.5
#' # 
#' #$cor_lat
#' #[,1] [,2]
#' #[1,]  1.0  0.5
#' #[2,]  0.5  1.0
#' #
#' #$Rmod
#' #[,1]  [,2] [,3]  [,4] [,5]  [,6]
#' #[1,] 1.00 0.420 0.48 0.210 0.18 0.150
#' #[2,] 0.42 1.000 0.56 0.245 0.21 0.175
#' #[3,] 0.48 0.560 1.00 0.280 0.24 0.200
#' #[4,] 0.21 0.245 0.28 1.000 0.42 0.350
#' #[5,] 0.18 0.210 0.24 0.420 1.00 0.300
#' #[6,] 0.15 0.175 0.20 0.350 0.30 1.000
#'
#' @export
#'
oblique_par2load = function(theta,grsize)
{ d = sum(grsize)
  mgrp = length(grsize)
  loadings = matrix(0,d,mgrp)
  theta_load = theta[1:d]
  theta_cor = theta[(d+1):length(theta)] # theta for correlation  of latents
  cor_lat = corvec2mat(theta_cor)  
  #fill-out the loadings into matrix
  st=0;
  for(j in 1:mgrp)
  { ind = (st+1):(st+grsize[j])
    loadings[ind,j] = theta_load[ind]
    st = st+grsize[j]
  }
  Rmod = loadings%*%cor_lat%*%t(loadings)
  diag(Rmod) = 1
  return(list(loadings=loadings, cor_lat=cor_lat, Rmod=Rmod))
}

# matrix square root of positive definite matrix
msqrt = function(Amat)
{ tem = svd(Amat)
  asqrt = tem$u %*% diag(sqrt(tem$d)) %*% t(tem$v)
  asqrt = 0.5 * (asqrt + t(asqrt))
}


# cor_mat (1,2), (1,3) (2,3), (1,4) etc

#' oblique factor correlation structure for d variables and m groups 
#' include determinant and inverse
#'
#' @description
#' For oblique factor correlation structure for d variables and m groups, 
#' convert the vector of parameters theta into
#' the loading matrix and correlation matrix of latent variables
#' The variables are assumed ordered by group.
#'
#' @param theta vector of length d + m*(m-1)/2; d loading parameters followed by
#'      m*(m-1)/2 entries in correlation matrix of latent variables
#'      (lower triangle by row)
#' @param grsize vector of group sizes (variables ordered by group)
#' @param icheck flag, if TRUE checks are made 
#'
#' @return loadings: loading matrix; cor_lat: correlation matrix of latent variables; 
#'  Rmod: correlation matrix based on theta for oblique factor model
#'
#' @examples 
#' theta = c(0.6,0.7,0.8,0.7,0.6,0.5,0.5)
#' oblique_pp_par2load(theta,grsize=c(3,3))
#' #$loadings
#' #[,1] [,2]
#' #[1,]  0.6  0.0
#' #[2,]  0.7  0.0
#' #[3,]  0.8  0.0
#' #[4,]  0.0  0.7
#' #[5,]  0.0  0.6
#' #[6,]  0.0  0.5
#' #
#' #$cor_lat
#' #[,1] [,2]
#' #[1,]  1.0  0.5
#' #[2,]  0.5  1.0
#' # 
#' #$Rmod
#' #[,1]  [,2] [,3]  [,4] [,5]  [,6]
#' #[1,] 1.00 0.420 0.48 0.210 0.18 0.150
#' #[2,] 0.42 1.000 0.56 0.245 0.21 0.175
#' #[3,] 0.48 0.560 1.00 0.280 0.24 0.200
#' #[4,] 0.21 0.245 0.28 1.000 0.42 0.350
#' #[5,] 0.18 0.210 0.24 0.420 1.00 0.300
#' #[6,] 0.15 0.175 0.20 0.350 0.30 1.000
#'
#' @export
#'
oblique_pp_par2load = function(theta,grsize, icheck=FALSE)
{ d = sum(grsize)
  mgrp = length(grsize)
  loadings = matrix(0,d,mgrp)
  theta_load = theta[1:d]
  theta_cor = theta[(d+1):length(theta)] # theta for correlation  of latents
  cor_lat = corvec2mat(theta_cor)  
  #fill-out the loadings into matrix
  st = 0
  for(j in 1:mgrp)
  { ind = (st+1):(st+grsize[j])
    loadings[ind,j] = theta_load[ind]
    st = st+grsize[j]
  }
  Rmod = loadings%*%cor_lat%*%t(loadings)
  diag(Rmod) = 1
  # get determinant and inverse
  B = loadings %*% msqrt(cor_lat)
  D = diag(1-theta_load^2)
  Dinv = diag(1/(1-theta_load^2))
  tem = diag(rep(1,mgrp)) + t(B)%*% Dinv %*% B
  # size is mgrpxmgrp instead of dxd
  Sdet = prod(1-theta_load^2) * det(tem)
  DB = Dinv%*%B
  Sinv = Dinv - DB %*% solve(tem,t(DB))
  # check
  if(icheck)
  { det2 = det(Rmod)
    message(paste("two dets :", Sdet,det2))
    imat = Rmod%*%Sinv
    tem = max(abs(imat-diag(1,d)))
    message(paste("max abs diff from identity S%*%Sinv :", tem))
  }
  return(list(loadings=loadings, cor_lat=cor_lat, Rmod=Rmod, Rinv=Sinv, Rdet=Sdet))
}

# nllk + gradient for nlm with analytic gradient

#' log-likelihood Gaussian oblique factor structure correlation matrix
#'
#' @description
#' log-likelihood Gaussian oblique factor structure correlation matrix with gradient
#'    for d variables and m groups, 
#'
#' @param theta vector of length d + m*(m-1)/2; d loading parameters followed by
#'      m*(m-1)/2 entries in correlation matrix of latent variables
#'      (lower triangle by row)
#' @param grsize vector of group sizes (variables ordered by group)
#' @param Robs dxd empirical correlation matrix
#' @param nsize sample size if available
#'
#' @return negative log-likelihood and gradient for Gaussian p-factor model
#'
#' @export
#'
oblique_grad_nllk = function(theta, grsize, Robs, nsize=100)
{ if(any(theta>0.999)|any(theta<(-0.999)))
  { # HJ this does not use the positive definite constraint on theta_cor
    # perhaps change to C-vine partial correlation form as input?
    # then convert to correlation matrix:  inputs to oblique_par2load
    #   theta_vec, cor_lat, grsize vector
    return(1e09)
  }
  obj =  oblique_pp_par2load(theta,grsize, icheck=FALSE)
  Rmod = obj$Rmod
  Rinv = obj$Rinv
  corW = obj$cor_lat 
  npar = length(theta); d = sum(grsize)
  cumd = cumsum(grsize)
  mgrp = length(grsize)
  cumd0 = c(0,cumd[-mgrp])+1
  temobs = Rinv %*% Robs 
  grad = rep(0,npar)
  # gradient wrt loadings
  for(j in 1:d)
  { nabla = matrix(0,d,d)
    nabla[j,] = Rmod[j,]/theta[j]
    nabla[,j] = Rmod[,j]/theta[j]
    nabla[j,j] = 0
    temth = Rinv %*% nabla 
    temprod = temth%*%temobs
    grad[j] = (sum(diag(temth)) - sum(diag(temprod))) * 0.5*nsize 
  }
  ip = d
  # gradient wrt latent variable correlations
  for(g2 in 2:mgrp)
  { grp2_ind = (cumd0[g2]:cumd[g2])
    for(g1 in 1:(g2-1))
    { nabla = matrix(0,d,d)
      grp1_ind = (cumd0[g1]:cumd[g1])
      if(corW[g1,g2]!=0)
      { nabla[grp1_ind,grp2_ind] = Rmod[grp1_ind,grp2_ind]/corW[g1,g2] }
      nabla[grp2_ind,grp1_ind] = t(nabla[grp1_ind,grp2_ind]) 
      #print(nabla)
      temth = Rinv %*% nabla 
      temprod = temth%*%temobs
      ip = ip+1
      grad[ip] = (sum(diag(temth)) - sum(diag(temprod))) * 0.5*nsize 
    }
  }
  nllk = 0.5*d*nsize*log(2*pi) + 0.5*nsize*log(obj$Rdet) +
    0.5*nsize*sum(diag(temobs)) 
  attr(nllk,"gradient") = grad
  nllk
}

#' log-likelihood Gaussian oblique factor structure correlation matrix
#'
#' @description
#' negative log-likelihood of the Gaussian oblique factor model 
#'    for d variables and m groups, 
#'
#' @param theta vector of length d + m*(m-1)/2; d loading parameters followed by
#'      m*(m-1)/2 entries in correlation matrix of latent variables
#'      (lower triangle by row)
#' @param grsize vector of group sizes (variables ordered by group)
#' @param Robs dxd (empirical) correlation matrix of normal scores
#' @param nsize sample size used to get Robs if available
#'
#' @return negative log-likelihood value of the oblique Gaussian factor model
#'    with fixed group size at MLE 
#'
#' @examples 
#' rhpar = c(0.81,0.84,0.84, 0.54,0.57,0.49, 0.51,0.54,0.55,0.70, 0.53,0.56,0.53,0.67,0.70)
#' cormat = corvec2mat(rhpar)
#' print(cormat)
#' #    [,1] [,2] [,3] [,4] [,5] [,6]
#' #[1,] 1.00 0.81 0.84 0.54 0.51 0.53
#' #[2,] 0.81 1.00 0.84 0.57 0.54 0.56
#' #[3,] 0.84 0.84 1.00 0.49 0.55 0.53
#' #[4,] 0.54 0.57 0.49 1.00 0.70 0.67
#' #[5,] 0.51 0.54 0.55 0.70 1.00 0.70
#' #[6,] 0.53 0.56 0.53 0.67 0.70 1.00
#' grsize = c(3,3)
#' mgrp = length(grsize)
#' d = sum(grsize)
#' theta = c(rep(0.3,d+mgrp*(mgrp-1)/2))
#' ml_obl = oblique_nllk(theta=theta, grsize, Robs=cormat)
#' print(ml_obl)
#' # 806.7432
#'
#' @export
#'
oblique_nllk = function(theta, grsize, Robs, nsize=100)
{ 
  if(any(theta>0.999)|any(theta<(-0.999)))
  { # HJ this does not use the positive definite constraint on theta_cor
    # perhaps change to C-vine partial correlation form as input?
    # then convert to correlation matrix:  inputs to par2load_oblique
    #   theta_vec, cor_lat, grsize vector
    return(1e09)
  }
  d = nrow(Robs)
  theta_cor = theta[(d+1):length(theta)] # theta for correlation  of latents
  cor_lat = corvec2mat(theta_cor)  
  # check if cor_lat is positive definite
  icheck = isPosDef(cor_lat)
  if(!icheck) { return(1.e10) }
  tem = oblique_par2load(theta,grsize)
  Rmod = tem$Rmod
  #negative log-likelihood for multivariate normal with mu=0, Sigma=Rmod
  nllk = 0.5*d*nsize*log(2*pi) + 0.5*nsize*log(det(Rmod)) +
    0.5*nsize*sum( diag(solve(Rmod,Robs)) )
    #0.5*nsize*sum(diag(solve(Rmod)%*%Robs))
  return(nllk)
}

#' Gaussian oblique factor structure correlation matrix
#'
#' @description
#'  MLE of parameters in the Gaussian oblique factor model
#'     for d variables and m groups, 
#'
#' @param grsize vector of group sizes (variables ordered by group)
#' @param start starting point should have dimension d+m*(m-1)/2
#' @param data nxd data set to compute the correlation matrix if cormat not given
#' @param cormat dxd empirical correlation matrix (of normal scores)
#' @param n sample size, if available
#' @param prlevel print.level for nlm()
#' @param mxiter maximum number of iterations for nlm()
#'
#' @return  list with nllk: negative log-likelihood; 
#' rhovec: the estimated mle;
#' loadings: loading matrix;
#' cor_lat: correlation matrix of the latent variables;
#' Rmod: the correlation matrix with optimized parameters.
#'
#' @examples 
#' # See example in bifact_fa() for a comparison with a data example
#' # Simpler example below
#' rhpar = c(0.81,0.84,0.84, 0.54,0.57,0.49, 0.51,0.54,0.55,0.70,  0.53,0.56,0.53,0.67,0.70)
#' cormat = corvec2mat(rhpar)
#' grsize = c(3,3)
#' mgrp = length(grsize)
#' d = sum(grsize)
#' start = rep(0.7,d+mgrp*(mgrp-1)/2)
#' res = oblique_fa(grsize, start, cormat=cormat, n=100, prlevel=1)
#' # iteration = 18
#' # Parameter:
#' # [1] 0.9005171 0.9064707 0.9270258 0.8202611 0.8480912 0.8243349 0.7039589
#' # Function Value
#' # [1] 622.5533
#' # Gradient:
#' # [1]  1.796252e-05  1.464286e-04  9.424639e-05 -8.344614e-05 -9.640644e-05
#' # [6] -9.799805e-05 -1.705303e-05
#' #
#' # Relative gradient close to zero.
#' # Current iterate is probably solution.
#'  
#' print(res)
#' # $nllk
#' # [1] 622.5533
#' # 
#' # $rhovec
#' # [1] 0.9005171 0.9064707 0.9270258 0.8202611 0.8480912 0.8243349 0.7039589
#' # 
#' # $loadings
#' # [,1]      [,2]
#' # [1,] 0.9005171 0.0000000
#' # [2,] 0.9064707 0.0000000
#' # [3,] 0.9270258 0.0000000
#' # [4,] 0.0000000 0.8202611
#' # [5,] 0.0000000 0.8480912
#' # [6,] 0.0000000 0.8243349
#' # 
#' # $cor_lat
#' # [,1]      [,2]
#' # [1,] 1.0000000 0.7039589
#' # [2,] 0.7039589 1.0000000
#' # 
#' # $Rmod
#' #           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#' # [1,] 1.0000000 0.8162923 0.8348025 0.5199856 0.5376279 0.5225681
#' # [2,] 0.8162923 1.0000000 0.8403217 0.5234234 0.5411823 0.5260230
#' # [3,] 0.8348025 0.8403217 1.0000000 0.5352925 0.5534542 0.5379511
#' # [4,] 0.5199856 0.5234234 0.5352925 1.0000000 0.6956562 0.6761698
#' # [5,] 0.5376279 0.5411823 0.5534542 0.6956562 1.0000000 0.6991112
#' # [6,] 0.5225681 0.5260230 0.5379511 0.6761698 0.6991112 1.0000000
#' 
#' Rmod = res$Rmod  #the average and max abs diff between Rmod and
#' sample correlation matrix
#' print(max(abs(Rmod-cormat))) #0.04657656
#' print(mean(abs(Rmod-cormat))) #0.0120401
#' 
#' @export
#'
oblique_fa = function (grsize, start, data=1, cormat=NULL, 
                       n=100, prlevel=0, mxiter=100)   
{ if (is.null(cormat)) { n = nrow(data); d = ncol(data); cormat = cor(data) }
  else 
  { d = sum(grsize); mgrp = length(grsize); npar = d+mgrp*(mgrp-1)/2 }
  if (length(start) != npar) 
  { message("the dimension of parameter vector is incorrect")
    return(0)
  }
  mle = nlm(oblique_nllk, p=start, grsize=grsize, Robs=cormat, 
            nsize=n, hessian=FALSE, iterlim=mxiter, print.level=prlevel) 
  estimate = mle$estimate
  tem = oblique_par2load(estimate, grsize=grsize)
  list(nllk=mle$minimum, rhovec=mle$estimate, loadings=tem$loadings, 
       cor_lat=tem$cor_lat, Rmod=tem$Rmod)
}

# version with analytic gradient for nlm()

#' Gaussian oblique factor structure correlation matrix
#'
#' @description
#' MLE of parameters in the Gaussian oblique factor model 
#'    for d variables and m groups, 
#'
#' @param grsize vector of group sizes (variables ordered by group)
#' @param start starting vector of length d + m*(m-1)/2; d loading parameters followed by
#'      m*(m-1)/2 entries in correlation matrix of latent variables
#'      (lower triangle by row)
#' @param data n x d data set to compute the correlation matrix if
#'        correlation matrix (cormat) not given
#' @param cormat dxd (empirical) correlation matrix of normal scores
#' @param n sample size, if available
#' @param prlevel print.level for nlm()
#' @param mxiter maximum number of iterations for nlm()
#'
#' @return  list with nllk: negative log-likeilihood;
#' rhovec: the estimated mle;
#' loadings: loading matrix;
#' cor_lat: correlation matrix of the latent variables; 
#' Rmod: the correlation matrix with optimized parameters.
#'
#' @export
#'
oblique_grad_fa = function (grsize, start, data=1, cormat=NULL, 
                            n=100, prlevel=0, mxiter=100)   
{ if (is.null(cormat)) { n = nrow(data); d = ncol(data); cormat = cor(data) }
  else 
  { d = sum(grsize); mgrp = length(grsize); npar = d+mgrp*(mgrp-1)/2 }
  if (length(start) != npar) 
  { message("the dimension of parameter vector is incorrect")
    return(0)
  }
  mle = nlm(oblique_grad_nllk, p=start, grsize=grsize, Robs=cormat, 
            nsize=n, hessian=FALSE, iterlim=mxiter, print.level=prlevel,
            check.analyticals=FALSE) 
  estimate = mle$estimate
  tem = oblique_par2load(estimate, grsize=grsize)
  list(nllk=mle$minimum, rhovec=mle$estimate, loadings=tem$loadings, 
       cor_lat=tem$cor_lat, Rmod=tem$Rmod)
}

