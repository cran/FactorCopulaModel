
#' Check if a square symmetric matrix is positive definite
#'
#' @description
#' Check if a square symmetric matrix is positive definite
#'
#' @param amat  symmetric matrix
#'
#' @return TRUE if amat is positive definite, FALSE otherwise
#'
#' @examples
#' a1 = matrix(c(1,.5,.5,1),2,2)
#' a2 = matrix(c(1,1.5,1.5,1),2,2)
#' t1 = try(chol(a1))
#' t2 = try(chol(a2))
#' print(isPosDef(a1))
#' print(isPosDef(a2))
#'
#' @export
#'
isPosDef = function(amat)
{ tt = try(chol(amat), silent=T)
  #ifelse(class(tt)=="matrix",T,F)
  # R v 4.3.3: class of a matrix object is character vector of length 2: "matrix" "array"
  ifelse(class(tt)[1]=="matrix",T,F)
}

