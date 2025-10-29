# functions for Spearman rho to copula parameter

#' Gumbel: Blomqvist's beta to copula parameter
#'
#' @description
#' Gumbel: Blomqvist's beta to copula parameter, vectorized
#'
#' @param beta vector of Blomqvist's beta values, 0<beta<1
#'
#' @return vector of Gumbel copula parameters with the given betas
#'
#' @examples
#' b = seq(0.1,0.5,0.1)
#' gumbel_beta2cpar(b)
#'
#' @export
#'
gumbel_beta2cpar = function(beta)
{ ln2 = log(2)
  ln2/log(2-(log(1+beta))/ln2)
}

#' Gumbel: Spearman rho to copula parameter
#'
#' @description
#' Gumbel: Spearman rho to copula parameter
#'
#' @param rho vector of Spearman values, 0<rho<1
#'
#' @return vector of Gumbel copula parameters with the given rho
#'
#' @examples
#' rho = seq(0.1,0.5,0.1)
#' gumbel_rhoS2cpar(rho)
#'
#' @export
#'
gumbel_rhoS2cpar = function(rho)
{
  gumdep = matrix(c(
  1.000000,0.00,0.00000000,
  1.021197,0.02,0.03108610,
  1.043193,0.04,0.06190168,
  1.066042,0.06,0.09244385,
  1.089801,0.08,0.12270912,
  1.114532,0.10,0.15269338,
  1.140303,0.12,0.18239191,
  1.167189,0.14,0.21179939,
  1.195270,0.16,0.24090988,
  1.224635,0.18,0.26971677,
  1.255384,0.20,0.29821286,
  1.287622,0.22,0.32639027,
  1.321470,0.24,0.35424044,
  1.357060,0.26,0.38175415,
  1.394537,0.28,0.40892146,
  1.434065,0.30,0.43573170,
  1.475826,0.32,0.46217347,
  1.520024,0.34,0.48823457,
  1.566890,0.36,0.51390203,
  1.616683,0.38,0.53916205,
  1.669696,0.40,0.56399998,
  1.726265,0.42,0.58840030,
  1.786770,0.44,0.61234656,
  1.851651,0.46,0.63582141,
  1.921413,0.48,0.65880651,
  1.996644,0.50,0.68128252,
  2.078029,0.52,0.70322908,
  2.166371,0.54,0.72462478,
  2.262619,0.56,0.74544712,
  2.367905,0.58,0.76567251,
  2.483585,0.60,0.78527622,
  2.611302,0.62,0.80423239,
  2.753063,0.64,0.82251401,
  2.911346,0.66,0.84009293,
  3.089252,0.68,0.85693986,
  3.290705,0.70,0.87302441,
  3.520753,0.72,0.88831511,
  3.785996,0.74,0.90277949,
  4.095233,0.76,0.91638418,
  4.460465,0.78,0.92909497,
  4.898489,0.80,0.94087701,
  5.433572,0.82,0.95169497,
  6.102112,0.84,0.96151330,
  6.961306,0.86,0.97029646,
  8.106482,0.88,0.97800931,
  9.709230,0.90,0.98461747,
  19.321382,0.95,0.99609174,
  Inf,1.00,1.00000000),48,3,T)
  colnames(gumdep)=c("cpar","beta","rhoS")
  rho2beta = approxfun(gumdep[,3], gumdep[,2])
  betav = rho2beta(rho)
  cpar = gumbel_beta2cpar(betav)
  cpar
}

#' Frank: Spearman rho to copula parameter
#'
#' @description
#' Frank: Spearman rho to copula parameter
#'
#' @param rho vector of Spearman values, -1<rho<1
#'
#' @return vector of Frank copula parameters with the given rho
#'
#' @examples
#' rho = seq(-0.2,0.6,0.1)
#' frank_rhoS2cpar(rho)
#'
#' @export
#'
frank_rhoS2cpar = function(rho)
{
  frkpos = matrix(c(
  0.0000000,0.00,0.00000000,
  0.1600427,0.02,0.02666468,
  0.3203418,0.04,0.05331740,
  0.4811559,0.06,0.07994620,
  0.6427471,0.08,0.10653906,
  0.8053836,0.10,0.13308392,
  0.9693417,0.12,0.15956860,
  1.1349078,0.14,0.18598087,
  1.3023811,0.16,0.21230833,
  1.4720763,0.18,0.23853845,
  1.6443265,0.20,0.26465854,
  1.8194863,0.22,0.29065571,
  1.9979360,0.24,0.31651687,
  2.1800857,0.26,0.34222866,
  2.3663801,0.28,0.36777748,
  2.5573044,0.30,0.39314944,
  2.7533910,0.32,0.41833033,
  2.9552274,0.34,0.44330558,
  3.1634653,0.36,0.46806026,
  3.3788320,0.38,0.49257903,
  3.6021436,0.40,0.51684611,
  3.8343211,0.42,0.54084524,
  4.0764104,0.44,0.56455966,
  4.3296058,0.46,0.58797205,
  4.5952804,0.48,0.61106450,
  4.8750229,0.50,0.63381844,
  5.1706844,0.52,0.65621461,
  5.4844377,0.54,0.67823302,
  5.8188531,0.56,0.69985279,
  6.1769969,0.58,0.72105218,
  6.5625598,0.60,0.74180838,
  6.9800276,0.62,0.76209743,
  7.4349101,0.64,0.78189402,
  7.9340527,0.66,0.80117123,
  8.4860677,0.68,0.81990024,
  9.1019408,0.70,0.83804983,
  9.7959003,0.72,0.85558576,
  10.5866898,0.74,0.87246989,
  11.4994771,0.76,0.88865902,
  12.5687952,0.78,0.90410333,
  13.8432286,0.80,0.91874449,
  15.3931758,0.82,0.93251347,
  17.3243547,0.84,0.94532856,
  19.8027736,0.86,0.95709410,
  23.1045856,0.88,0.96770098,
  27.7258491,0.90,0.97702921,
  Inf,1.00,1.00000000),47,3,T)
  colnames(frkpos)=c("cpar","beta","rhoS")
  frkneg = -frkpos[47:2,]
  frkdep = rbind(frkneg,frkpos)
  rho2beta = approxfun(frkdep[,3], frkdep[,2])
  betav = rho2beta(rho)
  cpar = frank_beta2cpar(betav)
  cpar

}


#' Frank: Blomqvist's beta to copula parameter
#'
#' @description
#' Frank: Blomqvist's beta to copula parameter, vectorized
#'
#' @param beta vector of Blomqvist's beta values, -1<beta<1
#' @param cpar0 starting point for Newton-Raphson iterations
#' @param mxiter maximum number of iterations, default 20
#' @param eps  tolerance for convergence, default 1.e-8
#' @param iprint print flag for iterations, default FALSE
#'
#' @return vector of Frank copula parameters with the given betas
#'
#' @examples
#' b = seq(-0.2,0.5,0.1) 
#' frank_beta2cpar(b)
#' frank_beta2cpar(b,iprint=TRUE)
#' 
#' @details
#' Solve equation to get cpar given Blomqvist's beta, Newton-Raphson iterations;
#'  vectorized input beta is OK, beta=0 fails
#'
#' @export
#'
frank_beta2cpar = function(beta, cpar0=0,mxiter=20,eps=1.e-8,iprint=FALSE)
{ iter = 0
  b = beta
  b[abs(b)<1.e-8] = 0.00000001
  diff = 1.
  cpar = cpar0
  if(cpar0<=0) cpar = 10*b
  ln2 = log(2)
  while(iter<mxiter & max(abs(diff))>eps)
  { tem = exp(cpar/2)
    g = log(1+tem)-ln2-cpar*(1+b)/4
    gp = 0.5/(1+1/tem)-(1+b)/4
    iter = iter+1
    diff = g/gp
    cpar = cpar-diff
    if(iprint) cat(iter," ",cpar," ",diff,"\n")
  }
  if(iter>=mxiter) message("did not converge\n")
  cpar
}

