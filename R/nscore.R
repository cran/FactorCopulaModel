#' Rank-based normal scores transform
#'
#' @description
#' Rank-based normal scores transform
#'
#' @param data dataframe or matrix, or vector, of reals
#'
#' @return matrix or vector of normal scores
#'
#' @export
#' 
nscore = function(data)
{ if(is.vector(data))
  { nr = length(data)
    anormal = -0.5  
    # this is a good adjustment 
    # (better than 0 in order than variance transform is closer to 1)
    qn = qnorm(((1:nr)+anormal)/(nr+1+2*anormal))
    jj = rank(data)
    out = qn[jj]
  }
  else
  { nc = ncol(data)
    nr = nrow(data)
    anormal = -0.5
    out = matrix(0,nr,nc)
    qn = qnorm(((1:nr)+anormal)/(nr+1+2*anormal))
    for(j in 1:nc)
    { jj = rank(data[,j])
      tem = qn[jj]
      out[,j] = tem
    }
  }
  out
}

