# Currently allowable copula familes using code of VineCopula
# 1 = Gaussian/normal
# 2 = t
# 4 = Gumbel
# 5 = Frank
# 7 = BB1
# 10 = BB8
# 14 = survival Gumbel
# 17 = survival BB1
# 20 = survival BB8
# this covers a range of different tail behaviors, other families could be added later


#' Sequential parameter estimation for bi-factor copula with estimated latent variables
#' using VineCopula::BiCopSelect
#'
#' @description
#' Sequential parameter estimation for bi-factor copula with estimated latent variables
#'
#' @param udata nxd matrix with values in (0,1)
#' @param vglobal n-vector is estimated global latent variables (or test with known values)
#' @param vgroup n*mgrp matrix with estimated group-based latent variables
#' @param grsize G-vector with group sizes with mgrp=G=length(grsize)groups
#' @param famset_global codes for allowable copula families for d global linking copulas
#' @param famset_group codes for allowable copula families for d global linking copulas
#'   VineCopula: current choices to cover a range of tail behavior are:
#'     1 = Gaussian/normal;
#'     2 = t;
#'     4 = Gumbel;
#'     5 = Frank;
#'     7 = BB1;
#'     10 = BB8;
#'     14 = survival Gumbel;
#'     17 = survival BB1;
#'     20 = survival BB8.
#' @param iprint if TRUE print intermediate results
#'
#' @return list with fam = 2*d vector of family codes chosen via BiCopSelect;
#'    dx2 matrix global_par; dx2 matrix group_par;
#'    these contain par,par2 for the selected copula families 
#'    in the 2-truncated vine rooted at the latent variables.
#'
#' @examples
#' # BB1/Frank bi-factor copula
#' set.seed(2024)
#' th1_range = c(0.3,1)
#' th2_range = c(1.1,2.5)
#' th3_range = c(8.5,18.5)
#' grsize = rep(10,3)
#' mgrp = length(grsize)
#' d = sum(grsize)
#' parbi = c(runif(d,th1_range[1],th1_range[2]),  # BB1 theta
#'         runif(d,th2_range[1],th2_range[2]),  # BB1 delta
#'         runif(d,th3_range[1],th3_range[2]))  # Frank
#' n = 500
#' data = rbifactor(n,grsize=grsize,cop=7,parbi)
#' udata = data$data
#' vlat = cbind(data$v0,data$vg)
#' fam_true = c(rep(7,d),rep(5,d))
#' #
#' guess = c(rep(0.7,d),rep(0.5,d)) 
#' bif_obj = bifactorScore(udata, start=guess, grsize, prlev=1)
#' proxy_init = bif_obj$proxies
#' # selection of linking copula families and estimation of their parameters 
#' select1 = bifactorEstWithProxy(udata,proxy_init[,1],proxy_init[,-1], grsize, 
#'   famset_global=c(1,4,5,7,10,17,20), famset_group=c(1,2,4,5))
#' fam1 = select1$fam
#' parglo1 = select1$global_par
#' pargrp1 = select1$group_par
#' print(fam1)  # 7,17,1 for global; 5 for group
#' print(parglo1)
#' print(pargrp1)
#' plot(parbi[(2*d+1):(3*d)],pargrp1[,1])
#' cor(parbi[(2*d+1):(3*d)],pargrp1[,1])
#' #
#' condExpProxy = latentUpdateBifactor(udata=udata, cparvec=c(parglo1,pargrp1),
#'   grsize=grsize, family=fam1, nq=25)
#' proxy_improved = cbind(condExpProxy$v0, condExpProxy$vg)
#' par(mfrow=c(2,2))
#' for(j in 1:(mgrp+1))
#' { plot(proxy_improved[,j],vlat[,j])
#'   print(cor(proxy_improved[,j],vlat[,j]))
#' }
#' rmse_values = sapply(1:ncol(proxy_improved),
#'      function(i) sqrt(mean((proxy_improved[,i] - vlat[,i])^2)))
#' round(rmse_values,3)
#' #
#' # With improved proxies,
#' # selection of linking copula families and estimation of their parameters 
#' select2 = bifactorEstWithProxy(udata,proxy_improved[,1],proxy_improved[,-1],
#'    grsize, famset_global=c(1,4,5,7,10,17,20), famset_group=c(1,2,4,5))
#' parglo2=select2$global_par
#' pargrp2=select2$group_par
#' fam2 = select2$fam  # 7,17 for global; 5 for local
#' cbind(fam_true,fam1,fam2)  
#' plot(parbi[(2*d+1):(3*d)],pargrp2[,1])
#' cor(parbi[(2*d+1):(3*d)],pargrp2[,1]) # higher correlation than pargrp1
#'
#' @references
#' 1. Krupskii P and Joe H (2013).
#' Factor copula models for multivariate data.
#' Journal of Multivariate Analysis, 120, 85-101.
#' 2. Fan X and Joe H (2024).
#' High-dimensional factor copula models with estimation of latent variables
#' Journal of Multivariate Analysis, 201, 105263.
#'
#' @details
#' It is best if variables have been oriented to be positively related
#' to the latent variable
#'
#' @export
#'
bifactorEstWithProxy = function(udata,vglobal,vgroup,grsize, 
  famset_global, famset_group, iprint=FALSE)
{ d = ncol(udata)  # should be equal to sum(grsize)
  mgrp = length(grsize)
  grpindex = rep(1,d)
  cum =  cumsum(grsize)
  for(g in 2:mgrp)
  { ii = (cum[g-1]+1):cum[g]
    grpindex[ii] = g
  }
  tree1par1 = rep(0,d)
  tree2par1 = rep(0,d)
  tree1par2 = rep(0,d)
  tree2par2 = rep(0,d)
  family1 = rep(0,d)
  family2 = rep(0,d)
  for(j in 1:d)
  { tem = BiCopSelect(udata[,j],vglobal, familyset=famset_global, rotations=F) 
    tree1par1[j] = tem$par
    tree1par2[j] = tem$par2
    fam = tem$family
    family1[j] = fam
    if(iprint) cat("tree1", j,tem$par,tem$par2,fam,"\n")
    pseud = BiCopHfunc2(udata[,j],vglobal,family=fam, par=tem$par, par2=tem$par2)
    gr = grpindex[j]
    tem2 = BiCopSelect(pseud, vgroup[,gr], familyset=famset_group, rotations=F)  
    tree2par1[j] = tem2$par
    tree2par2[j] = tem2$par2
    fam = tem2$family
    family2[j] = fam
    if(iprint) cat("tree2", j,tem2$par,tem2$par2,fam,"\n")
  }
  global_par = cbind(tree1par1,tree1par2)
  group_par = cbind(tree2par1,tree2par2)
  fam = c(family1,family2)
  list(fam=fam, global_par=global_par, group_par=group_par)
}

#' Parameter estimation for 1-factor copula with estimated latent variables
#' using VineCopula::BiCopSeelct
#'
#' @description
#' Parameter estimation for 1-factor copula with estimated latent variables
#'
#' @param udata nxd matrix with values in (0,1)
#' @param vlatent vector is estimated latent variables (or test with known values)
#' @param famset 2*d vector of codes for copula families for d global linking copulas
#' and d group-based linking copulas, using those from
#'   VineCopula: current choices to cover a range of tail behavior are:
#'     1 = Gaussian/normal;
#'     2 = t;
#'     4 = Gumbel;
#'     5 = Frank;
#'     7 = BB1;
#'     10 = BB8;
#'     14 = survival Gumbel;
#'     17 = survival BB1;
#'     20 = survival BB8.
#' @param iprint if TRUE print intermediate results
#'
#' @return list with fam = d-vector of family codes chosen via BiCopSelect;
#'    par1 = d-vector; par2 = d-vector of parameters for the selected copula families 
#'    in the 1-truncated vine rooted at the latent variable,
#'
#' @examples
#' # simulate data from 1-factor model with all Frank copulas
#' n = 500
#' d = 40
#' set.seed(20)
#' cpar = runif(d,4.2,18.5)
#' param = c(rbind(cpar,rep(0,d))) #Kendall's tau 0.4 to 0.8
#' data = r1factor(n,d,param,fam=rep(5,d))
#' vlat = data$vlatent # latent variables
#' udata = data$udata
#' proxyMean = uscore(apply(udata,1,mean)) # mean proxy
#' # RMSE of estimated latent variables
#' print(sqrt(mean((proxyMean-vlat)^2)))
#' # first estimation of 1-factor copula parameters
#' # allow for Frank, gaussian, t linking copulas
#' est1 = onefactorEstWithProxy(udata,proxyMean, famset=c(1,2,5))
#' print(est1$fam) # check choices , all 5s (Frank) in this case
#' print(est1$par1)
#' # estimation with only Frank copula as choice
#' est0 = onefactorEstWithProxy(udata,proxyMean, famset=c(5))
#' print(summary(abs(est0$par1-cpar)))  # same as $est1$par1
#' # improved conditional expectation proxies
#' # latentUpdate1factor allows for estimated linking copula with 2-parameters
#' # latentUpdate1factor1 can be used if estimated linking copulas all have par2=0
#' condExpProxy = latentUpdate1factor(c(rbind(est1$par1,est1$par2)),
#'    udata=udata,nq=25,family=rep(5,d))
#' # improved estimation of 1-factor copula parameters
#' est2 = onefactorEstWithProxy(udata,condExpProxy, famset=est1$fam)
#' print(est2$par1)  
#' # simple version of update for 1-parameter linking copulas
#' condExpProxy1 = latentUpdate1factor1(est0$par1,
#'    udata=udata,nq=25,family=rep(5,d))
#' summary(condExpProxy-condExpProxy1)  # 0 because family was chosen as 5 
#' print(summary(abs(est2$par1-cpar)))
#' # rmse of estimated latent variables
#' print(sqrt(mean((condExpProxy-vlat)^2))) 
#' # smaller rmse than initial proxies
#'
#' @references
#' 1. Krupskii P and Joe H (2013).
#' Factor copula models for multivariate data.
#' Journal of Multivariate Analysis, 120, 85-101.
#' 2. Fan X and Joe H (2024).
#' High-dimensional factor copula models with estimation of latent variables
#' Journal of Multivariate Analysis, 201, 105263.
#'
#' @details
#' It is best if variables have been oriented to be positively related
#' to the latent variable
#'
#' @export
#'
onefactorEstWithProxy = function(udata,vlatent, famset, iprint=FALSE)
{ d = ncol(udata)  
  tree1par1 = rep(0,d)
  tree1par2 = rep(0,d)
  family = rep(0,d)
  for(j in 1:d)
  { tem = BiCopSelect(udata[,j],vlatent, familyset=famset, rotations=F) 
    tree1par1[j] = tem$par
    tree1par2[j] = tem$par2
    fam = tem$family
    family[j] = fam
    if(iprint) cat("tree1", j,tem$par,tem$par2,fam,"\n")
  }
  list(fam=family, par1=tree1par1, par2=tree1par2)
}

