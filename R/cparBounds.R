# Currently allowable copula familes using code of VineCopula
# 1 = Gussian/normal
# 2 = t
# 4 = Gumbel
# 5 = Frank
# 7 = BB1
# 10 = BB8
# 14 = survival Gumbel
# 17 = survival BB1
# 20 = survival BB8
# this covers a range of different tail behaviors, other families could be added later


#' lower and upper bounds for copula parameters (1-parameter, 2-parameter families) 
#'
#' @description  
#' lower and upper bounds for copula parameters for use in min negative log-likelihood 
#'
#' @param familyvec vector of family codes linking copula families
#'
#' @return lower bound LB1/LB2 and upper bound LB2/UB2 for par1 and par2
#'
#' @examples
#' famvec = c(1,4,5,14,2,7,17,10,20, 4,5,1)
#' out = cparBounds(famvec)
#' print(out)
#' print(out$LB1)
#'
#'@export
#'
cparBounds = function(familyvec)
{ n = length(familyvec)
  bounds = matrix(0,n,4)  # LB1,UB1, LB2,UB2 for par1, par2
  for(i in 1:n)
  { tem = familyvec[i]
    if(tem==2) { bounds[i,] = c(-0.99,0.99,2.01,5000) }
    else if (tem==17 | tem==7) { bounds[i,] = c(0.0001,7,1,7) }
    else if (tem==1) { bounds[i,] = c(-0.99,0.99,0,0) }
    else if (tem==5) { bounds[i,] = c(-99.99,99.99,0,0) }
    else if(tem==10 | tem==20) { bounds[i,] = c(1,8,0.00001,1)  }
    else if(tem==14 | tem==4) { bounds[i,] = c(1,17,0,0) }
  }
  colnames(bounds) = c("LB1","UB1", "LB2","UB2")
  bounds = data.frame(bounds)
  return(bounds)
}

