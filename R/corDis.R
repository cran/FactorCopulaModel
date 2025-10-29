
#' Discrepancy of model-based and observed correlation matrices based on Gaussian log-likelihood
#'
#' @description 
#' Discrepancy of model-based and observed correlation matrices 
#'
#' @param Rmodel model-based correlation matrix
#' @param Rdata empirical correlation matrix (could be observed or polychoric)
#' @param n sample size (if positive integer)
#' @param npar #parameters in the correlation structure
#'
#' @return vector with discrepancy Dfit, and 
#'      also nllk2 (wice negative log-likelihood), BIC, AIC if n and npar are inputted
#'
#' @examples
#' Rmodel = matrix(c(1,.3,.4,.4,.3,1,.5,.6,.4,.5,1,.7,.4,.6,.7,1),4,4)
#' print(Rmodel); print(chol(Rmodel))
#' Rdata = matrix(c(1,.32,.38,.41,.32,1,.53,.61,.38,.53,1,.67,.41,.61,.67,1),4,4)
#' print(corDis(Rmodel,Rdata))
#' print(corDis(Rmodel,Rdata,n=400,npar=3))
#'
#' @export
#'
corDis = function(Rmodel,Rdata,n=0,npar=0)
{ lgdetmod = log(det(Rmodel))
  d = nrow(Rdata)
  tem = sum(diag(solve(Rmodel,Rdata)))
  Dfit = lgdetmod-log(det(Rdata))+tem-d
  if(n==0) return(Dfit)
  nllk2 = n*d*(log(2*pi))+n*lgdetmod + n*tem
  bic = nllk2+log(n)*npar
  aic = nllk2+2*npar
  out = c(Dfit,nllk2,aic,bic)
  names(out) = c("Dfit","nllk2","AIC","BIC")
  out
}

