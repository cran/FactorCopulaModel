
#' R interface for Gauss-Legendre quadrature
#'
#' @description 
#' Gauss-Legendre quadrature nodes and weights
#'
#' @param nq number of quadrature points
#'
#' @return structures with
#'
#'   nodes: nq-vector of quadrature nodes,
#'   weights: nq-vector of quadrature weights
#' @examples
#' out = gaussLegendre(15)
#' # same as statmod::gauss.quad.prob(15,dist="uniform") in library(statmod)
#' print(sum(out$weights)) # should be 1
#' print(sum(out$weights*out$nodes)) # should be 0.5  = E(U), U~Uniform(0,1)
#' print(sum(out$weights*out$nodes^2)) # should be 1/3 = E(U^2)
#'
#' @details
#' links to C code translation of jacobi.f in Stroud and Secrest (1966)
#'
#' @export
#'
gaussLegendre = function(nq)
{ # need a shift because C code indexes 1 to nq, not 0 to  nq-1
  if(nq<=0 || nq>70) stop("nq between 1 and 70")
  nq1 = nq+1
  out = .C("gauleg", x1=as.double(0), x2=as.double(1),
         as.integer(nq), xq=as.double(rep(0,nq1)), wq=as.double(rep(0,nq1)) )
  list(nodes=out$xq[-1], weights=out$wq[-1])
}


