
#' simulate from bi-factor copula model
#'
#' @description 
#' simulate from bi-factor copula model  and include corresponding latent variables
#'
#' @param n sample size
#' @param grsize G-vector of group sizes for G groups 
#' @param cop  code for copula families 
#' 1: Gaussian/Gaussian; 2: t/t; 4: Gumbel/Gumbel; 5: Frank/Frank; 7: BB1/Frank;
#' 14: survival Gumbel, 17: survivalBB1 /Frank
#'  if cop = 1, data have standard normal marginals
#'  if cop = 2, data have t marginals
#'  if cop > 2, data have uniform(0,1) marginals
#' @param param vector of parameters (those for the common factor go first)
#' The order in param is the same as in start for mvtbifct(full=T) function.
#' For BB1/Frank: BB1thetas then BB1deltas, then Frank parameters --
#' the order in param is the same as in start for mvtbifct(full=F) function
#'
#' @return list with data: nxd data set with U(0,1) or N(0,1) or t(df) margin; 
#'   v0: n-vector of corresponding global latent variables; and
#'   vg: nxG matrix of corresponding local (group) latent variables.
#'
#' @examples
#' grsize = c(4,3)
#' cop = 4; param4 = c(1.5,1.6,1.7,1.8,1.9,2.0,2.1,  1.1,1.1,1.1,1.1,1.1,1.1,1.1)
#' cop = 14; param14 = c(1.5,1.6,1.7,1.8,1.9,2.0,2.1,  1.1,1.1,1.1,1.1,1.1,1.1,1.1)
#' cop = 5; param5 = c(1.5,1.6,1.7,1.8,1.9,2.0,2.1,  1.1,1.1,1.1,1.1,1.1,1.1,1.1)
#' cop = 1; param1 = c(0.5,0.6,0.7,0.8,0.9,0.4,0.5,  0.1,0.1,0.1,0.1,0.1,0.1,0.1) 
#' cop = 2; param2 = c(0.5,0.6,0.7,0.8,0.9,0.4,0.5,  0.1,0.1,0.1,0.1,0.1,0.1,0.1, 7)
#' cop = 7; param7 = c(0.5,0.6,0.7,0.8,0.9,1.0,1.1, 1.5,1.6,1.7,1.8,1.9,1.4,1.5, 1.1,1.1,1.1,1.1,1.1,1.1,1.1) 
#' cop = 17; param17 = c(0.5,0.6,0.7,0.8,0.9,1.0,1.1, 1.5,1.6,1.7,1.8,1.9,1.4,1.5, 1.1,1.1,1.1,1.1,1.1,1.1,1.1) 
#' set.seed(123)
#' gumdat = rbifactor(n=10, grsize=grsize, cop=4, param=param4)    # U(0,1)
#' gumrdat = rbifactor(n=10, grsize=grsize, cop=14, param=param14) # U(0,1)
#' frkdat = rbifactor(n=10, grsize=grsize, cop=5, param=param5)    # U(0,1)
#' gaudat = rbifactor(n=10, grsize=grsize, cop=1, param=param1)    # N(0,1)
#' bvtdat = rbifactor(n=10, grsize=grsize, cop=2, param=param2)    # t_7
#' bb1frkdat = rbifactor(n=10, grsize=grsize, cop=7, param=param7)    # U(0,1)
#' bb1rfrkdat = rbifactor(n=10, grsize=grsize, cop=17, param=param17) # U(0,1)
#' summary(bb1frkdat$data)
#' summary(bb1frkdat$v0)
#' summary(bb1frkdat$vg)
#'
#' @details
#' The user can modify this code to get other linking copulas.
#'
#' @export
#'
rbifactor = function(n, grsize, cop=5, param) 
{ d = sum(grsize)
  mgrp = length(grsize)
  th1 = param[1:d]  # parameter linking to global latent (Frank, Gumbel, survGumbel)
  th2 = param[(d+1):(2*d)]  # parameter linking to group latent
  #                    or second parameter for global latent (BB1)
  # t/t linking copulas
  if (cop==2) { df0 = param[2*d+1] }  # same dof parameter for all, position 2*d+1
  # BB1 or survival BB1 global, Frank group
  if (cop==7 || cop==17) { th3 = param[(2*d+1):(3*d)] } # Frank parameter linking to group latent
  # Gaussian or t, from partial correlations to loadings
  if (cop==1 || cop==2) { ld2 = th2*sqrt(1-th1^2) }

  zdata = matrix(0, nrow=n, ncol=d)
  # latent is normal for cop=1 or cop=2; 
  if (cop==1 || cop==2) 
  { z0 = rnorm(n)
    z = matrix(rnorm(n*mgrp), ncol=mgrp)
  } 
  else 
  { z0 = runif(n)
    z = matrix(runif(n*mgrp), ncol=mgrp)
  }

  ind = 0
  for (jg in 1:mgrp) 
  { ind1 = ind + 1
    ind2 = ind + grsize[jg]
    ind = ind + grsize[jg]
    for (ij in ind1:ind2) 
    { if (cop==1 || cop==2) 
      { zdata[,ij] = th1[ij]*z0 + z[,jg]*ld2[ij] + sqrt(1-th1[ij]^2-ld2[ij]^2)* rnorm(n)
      } 
      else if (cop==4 || cop==5 || cop==14) # 1-parameter Gumbel, or Frank, or survival Gumbel
      { #bicop1 = BiCop(family=cop, par=th2[ij])
        #q1 = BiCopHinv(runif(n), z[,jg], obj = bicop1)$hinv2
        q1 = VineCopula::BiCopHinv2(runif(n), z[,jg], family=cop, par=th2[ij]) 
        #bicop2 = BiCop(family=cop, par=th1[ij])
        #zdata[,ij] = BiCopHinv(q1, z0, obj = bicop2)$hinv2
        zdata[,ij] = VineCopula::BiCopHinv2(q1, z0, family=cop, par=th1[ij])
      } 
      ## BB1 or survival BB1 for global latent, Frank for group latent
      else if (cop==7 || cop==17) 
      { for (i in 1:n) 
        { #bicop = BiCop(family=5, par=th3[ij])
          #q1 = BiCopHinv(runif(1), z[i,jg], obj=bicop)$hinv2
          q1 = VineCopula::BiCopHinv2(runif(1),z[i,jg], family=5, par=th3[ij]) 
          #bicop1 = BiCop(family=7, par=th1[ij], par2=th2[ij])
          #zdata[i, ij] = BiCopHinv(q1, z0[i], obj = bicop1)$hinv2
          zdata[i,ij] = VineCopula::BiCopHinv2(q1, z0[i], family=cop, par=th1[ij], par2=th2[ij])
        }
      } 
    }
  }

  # maybe can vectorize? if so, need a rchisq vector
  if (cop == 2) 
  { for (i in 1:n) 
    { zdata[i,] = zdata[i,] / sqrt(rchisq(1, df=df0)/df0) }
  }

  return(list(data=zdata, v0=z0, vg=z))
}
