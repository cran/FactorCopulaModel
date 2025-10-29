#' BB1 tail dependence parameters to copula parameter (theta,delta) 
#'
#' @description
#' BB1 map (lower,upper) tail dependence to copula parameter vector
#'
#' @param taildep tail dependence parameter in mx2 matrix, by row (ltd,utd) in (0,1)^2 
#'
#' @return matrix of copula parameters, by row (theta,delta), theta>0, delta>1
#'
#' @examples
#' cpar = bb1_td2cpar(c(0.4,0.6))
#' print(cpar)
#' #         theta    delta
#' #[1,] 0.3672112 2.060043
#' print(bb1_cpar2td(cpar))
#'
#' @export
#'
bb1_td2cpar = function(taildep)  
{ if(is.vector(taildep)) taildep = matrix(taildep,1,2)
  delta = log(2)/log(2-taildep[,2]); 
  theta = log(2)/(-delta*log(taildep[,1]))
  cbind(theta,delta)
}
