
#' Rank-based uniform scores transform
#'
#' @description
#' Rank-based uniform scores transform
#'
#' @param data dataframe or matrix, or vector, of reals
#' @param aunif adjustment 'a' for (rank+a)/(n+1+2*a) as scores. n=sample size , default is -0,5
#'
#' @return matrix or vector of uniform scores
#'
#' @export
#' 
uscore = function(data,aunif=-0.5)
{ if(is.vector(data))
  { nr = length(data)
    us = ((1:nr)+aunif)/(nr+1+2*aunif)
    jj = rank(data)
    out = us[jj]
  }
  else
  { nc = ncol(data)
    nr = nrow(data)
    out = matrix(0,nr,nc)
    us = ((1:nr)+aunif)/(nr+1+2*aunif)
    for(j in 1:nc)
    { jj = rank(data[,j])
      tem = us[jj]
      out[,j] = tem
    }
  }
  out
}

