
#' Convert from correlations in vector form to a correlation matrix
#'
#' @description
#' Convert from correlations in vector form to a correlation matrix
#'
#' @param rvec correlations in vector form of length d*(d-1)/2 in the order r12,r13,r23,r14,... r[d-1,d]
#'
#' @return dxd correlation matrix
#'
#' @examples
#' rvec = c(0.3,0.4,0.5,0.4,0.6,0.7)
#' Rmat = corvec2mat(rvec)
#' print(Rmat) # column 1 has 1, 0.3, 0.4, 0.4
#'
#' @export
#'
corvec2mat = function(rvec)
{ dd = length(rvec)
  d = (1+sqrt(1+8*dd))/2
  rmat = matrix(1,d,d)
  k = 0
  for(i in 2:d)
  { for(j in 1:(i-1))
    { k = k+1; rmat[i,j] = rvec[k]; rmat[j,i] = rvec[k] }
  }
  rmat
}

