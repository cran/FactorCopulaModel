# Tail-weighted dependence measure [zeta(alpha), alpha>0 large].
# Lee D, Joe H and Krupskii P (2018).
# A tail-weighted dependence measure with limit being tail dependence coefficient. 
# J Nonparametric Statistics, 30(2), 262-290.
# http://dx.doi.org/10.1080/10485252.2017.1407414


#' Upper Tail-weighted dependence measure zeta(C,alpha)
#'
#' @description
#' Upper Tail-weighted dependence measure zeta(C,alpha)
#'
#' @param cpar copula parameter (vector) of pcop
#' @param pcop copula cdf C, assume this takes vectorized form for u,v
#' @param alpha scalar alpha>0 for zeta measure
#' @param zero tolerance to use for zero if integrate is used, e.g., 0.0001
#' @param iGL TRUE to use Gauss-Legendre quadrature
#' @param gl Gauss-Legendre quadratureobject  with nodes/weights if iGL=TRUE
#'
#' @return zeta(C,alpha) ; this is upper tail-weighted is alpha>>1
#'
#' @details
#' zeta(alpha)=2+alpha-alpha/integral
#' integral int_0^1 C( x^{1/alpha}, x^{1/alpha} ) dx.
#' This is a central dependence of measure if alpha =1 and upper tail-weighted is alpha>>1.
#'
#' @references
#' Lee D, Joe H, Krupskii P (2018). J Nonparametric Statistics, 30(2), 262-290
#'
#' @examples
#' # Bivariate margin of 1-factor copula
#' # using conditional cdf via VineCopula::BiCopHfunc2
#' # cpar1 = par,par2 for fam1 for first variable
#' # cpar2 = par,par2 for fam2 for second variable
#' # family codes 1:Gaussian, 2:t, 4:Gumbel, 5:Frank, 7:BB1, 14:survGumbel 17:survBB1
#' p1factbiv = function(u1,u2,cpar1,cpar2,fam1,fam2,nq)
#' { if(length(u1)==1) u1 = rep(u1,nq)
#'   if(length(u2)==1) u2 = rep(u2,nq)
#'   gl = gaussLegendre(nq)
#'   wl = gl$weights
#'   vl = gl$nodes
#'   a1 = BiCopHfunc2(u1,vl,family=fam1,par=cpar1[1], par2=cpar1[2])
#'   a2 = BiCopHfunc2(u2,vl,family=fam2,par=cpar2[1], par2=cpar2[2])
#'   sum(wl*a1*a2)
#' }
#' #
#' # Version of zeta(C) for bivariate margin of 1-factor copula
#' zetaDepC_1factor = function(cpar1,cpar2,fam1,fam2,nq,alpha,zero=0,iGL=FALSE,gl=0)
#' { a1 = 1/alpha
#'   gfn = function(x) 
#'   { nn = length(x)
#'     gval = rep(0,nn)
#'     # need this form because of nesting in p1factbiv()
#'     for(ii in 1:nn) 
#'     { xx = x[ii]^a1
#'       gval[ii] = p1factbiv(xx,xx,cpar1,cpar2,fam1,fam2,nq) 
#'     }
#'     gval
#'   }
#'   if(iGL)
#'   { xq = gl$nodes
#'     wq = gl$weight
#'     tem = sum(wq*gfn(xq))
#'   }
#'   else
#'   { tem = integrate(gfn,zero,1-zero)
#'     tem = tem$value
#'   }
#'   zeta = 2+alpha-alpha/tem
#'   zeta
#' }
#' # Tests for zeta
#' param = matrix(c(0.5,1.5,0.6,1.2),2,2,byrow=TRUE)
#' # Create BB1 copula cdf pbb1 and  BB1 surival copula cdf pbb1r using BiCopCDF
#' pbb1_VC = function(u,v,cpar) { BiCopCDF(u,v,family=7,par=cpar[1],par2=cpar[2]) }
#' pbb1r_VC = function(u,v,cpar) { BiCopCDF(u,v,family=17,par=cpar[1],par2=cpar[2]) }
#' gl21 = gaussLegendre(21)
#' zeta1u_bb1 = zetaDepC(param[1,],pbb1_VC,alpha=10,zero=0.00001,iGL=TRUE,gl=gl21)
#' zeta1l_bb1 = zetaDepC(param[1,],pbb1r_VC,alpha=10,zero=0.00001,iGL=TRUE,gl=gl21)
#' cat(zeta1u_bb1,zeta1l_bb1,"\n")
#' # 0.4504351 0.4776329 
#' zeta2u_bb1 = zetaDepC(param[2,],pbb1_VC,alpha=10,zero=0.00001,iGL=TRUE,gl=gl21)
#' zeta2l_bb1 = zetaDepC(param[2,],pbb1r_VC,alpha=10,zero=0.00001,iGL=TRUE,gl=gl21)
#' cat(zeta2u_bb1,zeta2l_bb1,"\n")
#' # 0.2825654 0.419787 
#' # Bivariate margin of 1-factor copula: linking BB1(param1) and BB1(param2) 
#' # Upper tail
#' zetau_ai = zetaDepC_1factor(param[1,],param[2,],fam1=7,fam2=7,nq=21,alpha=10,
#'    zero=0.00001,iGL=FALSE,gl=0)
#' zetau_gl = zetaDepC_1factor(param[1,],param[2,],fam1=7,fam2=7,nq=21,alpha=10,
#'    zero=0.00001,iGL=TRUE,gl=gl21)
#' cat(zetau_ai,zetau_gl,"\n")
#' # 0.1766584 0.1775621
#' # Lower tail
#' zetal_ai = zetaDepC_1factor(param[1,],param[2,],fam1=17,fam2=17,nq=21,alpha=10,
#'    zero=0.00001,iGL=FALSE,gl=0)
#' zetal_gl = zetaDepC_1factor(param[1,],param[2,],fam1=17,fam2=17,nq=21,alpha=10,
#'    zero=0.00001,iGL=TRUE,gl=gl21)
#' cat(zetal_ai,zetal_gl,"\n")
#' # 0.2664287 0.26737
#' # Ordering is expected based on the individual BB1(param1) and BB1(param2)
#'
#' @export
#'
zetaDepC = function(cpar,pcop,alpha,zero=0,iGL=FALSE,gl=0)
{ a1 = 1/alpha
  gfn = function(x) { pcop(x^a1,x^a1,cpar) }
  if(iGL)
  { xq = gl$nodes
    wq = gl$weight
    tem = sum(wq*gfn(xq))
  }
  else
  { tem = integrate(gfn,zero,1-zero)
    tem = tem$value
  }
  zeta = 2+alpha-alpha/tem
  zeta
}



#' Empirical version of zeta(alpha) tail-weighted dependence measure
#'
#' @description
#' Empirical version of zeta(alpha) tail-weighted dependence measure
#'
#' @param dat nx2 data matrix with values (in (0,1) if rank=FALSE)
#' @param alpha vector of alpha>0 for zeta measure
#' @param rank TRUE (default) if to convert data matrix to uniform scores in (0,1)
#' @param lowertail TRUE if lower tail-weighted dependence measure, default is FALSE
#'
#' @return zeta(C,alpha)
#'
#' @examples
#' data(euro07gf)
#' udat = euro07gf$uscore
#' euro07names = colnames(udat)
#' d = ncol(udat)
#' for(j2 in 2:d)
#' { for(j1 in 1:(j2-1))
#' { zetaU = zetaDep(udat[,c(j1,j2)],alpha=15,rank=FALSE,lowertail=FALSE)
#'   zetaL = zetaDep(udat[,c(j1,j2)],alpha=15,rank=FALSE,lowertail=TRUE)
#'   zeta1 = zetaDep(udat[,c(j1,j2)],alpha=1,rank=FALSE,lowertail=FALSE)
#'   cat(j1,j2,round(zeta1,3),round(zetaL,3),round(zetaU,3),euro07names[j1],euro07names[j2],"\n")
#'   }
#' }
#'
#' @details
#' This is a central dependence of measure if alpha =1 and
#' upper tail-weighted is alpha>>1
#'
#' @references
#' Lee D, Joe H, Krupskii P (2018). J Nonparametric Statistics, 30(2), 262-290
#'
#' @export
#'
zetaDep = function(dat,alpha,rank=TRUE,lowertail=FALSE)
{ if(rank) dat = uscore(dat)
  if(lowertail) dat = 1-dat
  nalp = length(alpha)
  zeta = rep(0,nalp)
  for(i in 1:nalp)
  { alp = alpha[i]
    nu = 0.5*mean(abs(dat[,1]^alp-dat[,2]^alp))
    theta = (alp+alp*(alp+1)*nu) / (alp-(alp+1)*nu)
    zeta[i] = 2-theta
  }
  return(zeta)
}

