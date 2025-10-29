# version of method in Brechmann and Joe (2014, CSDA) using igraph
# code mostly written by E C Brechmann in 2012, and then edited by H Joe

# Gaussian correlation matrix with 1-factor and 1-truncation vine for
# residual dependence (conditional dependence given latent variable).

#' correlation matrix for 1-factor plus 1-truncated vine (for residual dependence)
#'
#' @description
#' correlation matrix for 1-factor plus 1-truncated vine (for residual dependence)
#'
#' @param cormat dxd correlation matrix
#' @param loading d-dimensional loading vector (for latent factor), -1<loading[j]<1
#'
#' @return list with R = correlation matrix for structure of 1-factor+Markov tree residual dependence;
#' incl = d*(d-1)/2 binary vector: indicator of edges [1,2], [1,3], [2,3], [1,4], ...[d-1,d] edges in tree with d-1 edges;
#' partcor = conditional correlation matrix given the latent variable.
#'
#' @details
#'  MST algorithm with weights log(1-rho^2), rho's are partial correlations
#'  fiven the latent variable.
#' not exported
#'
residDep = function(cormat,loading)
{ d = dim(cormat)[1]
  # compute conditional/partial correlation matrix
  condcormat = matrix(1,d,d)
  for(i in 2:d) for(j in 1:(i-1)) 
  { condcormat[i,j] = condcormat[j,i] = (cormat[i,j]-loading[i]*loading[j])/
      sqrt(1-loading[i]^2)/sqrt(1-loading[j]^2)
  }
  colnames(condcormat) = rownames(condcormat) = paste("V",1:d,sep="")
  condcormat[condcormat==0] = 1e-16

  # define graph from adjacency matrix (= correlation matrix) and find MST
  g = graph.adjacency(condcormat, mode="lower", weighted=TRUE, diag=FALSE)
  E(g)$name = paste(get.edgelist(g)[,1],get.edgelist(g)[,2],sep=",")
  mst = minimum.spanning.tree(g, weights=log(1-E(g)$weight^2))

  # extract edges
  edges = get.edges(mst,1:(d-1))
  # vector of length d*(d-1)/2 for which edges of the complete graph are included in the MST
  # order of indexing is: [1,2], [1,3], [2,3], [1,4], ...[d-1,d]
  incl = E(g)$name %in% E(mst)$name
  # define order of variables in vine
  node_order = rep(NA,d)
  edge_order = rep(NA,d-1)
  elig = rep(TRUE,d-1)

  # first node (node with few connections to other nodes)
  node_order[1] = order(table(edges))[1]

  # second node (neighbor of first node)
  enum = which(edges %in% node_order[1]) %% (d-1)
  if(enum == 0) enum = d-1
  node_order[2] = edges[enum,][edges[enum,]!=node_order[1]]
  elig[enum] = FALSE
  edge_order[1] = enum

  # remaining nodes (neighbors of already selected nodes)
  for(i in 3:(d-1))
  { poss = rep(FALSE,d-1)
    for(j in 1:(d-1)) poss[j] = sum(edges[j,] %in% node_order[1:(i-1)])>0
    poss[!elig] = FALSE
    enum = min(which(poss))
    node_order[i] = edges[enum,][!(edges[enum,] %in% node_order[1:(i-1)])]
    edge_order[i-1] = enum
    elig[enum] = FALSE
  }
  # last node
  node_order[d] = edges[elig,][!(edges[elig,] %in% node_order[1:(d-1)])]
  edge_order[d-1] = which(elig)

  # sort edges such that first column = j and 
  #  second column = k2(j) in Z_j = phi1*L + phi2*Z_k2(j) + psi*eps
  # equation (2.3) in Brechmann and Joe (2014).
  edges_sort = matrix(NA,d-1,2)
  for(i in 1:(d-1))
  { first = which(edges[edge_order[i],] == node_order[i+1])
    second = which(edges[edge_order[i],] != node_order[i+1])
    edges_sort[i,] = c(edges[edge_order[i],first],edges[edge_order[i],second])
  }
  
  # build R-vine structure matrix with factor = node "d+1"
  #Mat = matrix(NA,d+1,d+1)
  VA = matrix(NA,d+1,d+1)
  #Cor = VA(NA,d+1,d+1)
  PCor = matrix(NA,2,d+1)
  VA[1,] = d+1
  VA[2,2] = node_order[1]
  PCor[1,2] = loading[node_order[1]]
  for(i in 2:d)
  { VA[i+1,i+1] = edges_sort[i-1,1]
    VA[2,i+1] = edges_sort[i-1,2]
    PCor[1,i+1] = loading[node_order[i]]
    PCor[2,i+1] = condcormat[edges_sort[i-1,1],edges_sort[i-1,2]]
  }
  # Vine array has variable 1 for latent variable  and then the original variables (reindexed and incremented by 1)
  # row1 of VA consists 1's (for latent variable)
  # row2 of VA and diagonal of VA indicate edges of Markov tree for residual dependence
  # this tree has edges [j, VA[2,j]; 1] for j=2,...,1+#vars
  #  for partial correlations rho_{j-1,A[2,j]-1;latent} 
  # PCor has loading parameters (cor with latent variable) in PCor[1,2:(d+1)]
  # PCor has partial correlations given latent variable in PCor[2,3:(d+1)]
  #RVM = list(Matrix=VA, PCor=PCor)
  RVobj = list(Varray=VA, PCor=PCor)

  # change labels to 1:(d+1)
  oldVarray = RVobj$Varray
  MM = RVobj$Varray
  oldOrder = diag(MM)
  O = apply(t(1:nrow(MM)),2,"==", MM)
  for(i in 1:nrow(MM)) { MM[O[,oldOrder[i]]] = i }
  RVobj1 = list(Varray=MM, PCor=PCor)

  # compute new correlation matrix from R-vine structure matrix
  new_cormat = RVtrunc2cor(RVobj1)
  # order correlation matrix
  new_cormat = new_cormat[order(oldOrder),order(oldOrder)]
  return(list(R=new_cormat[1:d,1:d],incl=incl,partcor=condcormat))
}

#' compute correlation matrix from 2-truncated R-vine
#'
#' @description
#' compute correlation matrix from 2-truncated R-vine
#'
#' @param RVobj R-vine object with vine array and partial correlation matrix
#' variable 1 is the latent variable;
#' list with $Varray (d+1)x(d+1), $PCor 2x(d+1)
#'
#' @return correlation matrix for 2-truncated vine structure based on 
#' partial correlations in tree 2
#'
#' @details
#' not exported
#'
RVtrunc2cor = function(RVobj)
{ # dimension of matrix is 1 more than number of observed variables
  # variable 1 is the latent variable
  d1 = dim(RVobj$Varray)[1]
  # identify k1 and k2 in Z_j = phi1*Z_k1(j) + phi2*Z_k2(j) + psi*eps
  k1 = RVobj$Varray[1,]
  k2 = RVobj$Varray[2,]
  ksort = cbind(pmin(k1,k2),pmax(k1,k2))

  # given correlations
  rmat = matrix(NA,d1,d1)
  for(j in 2:d1) rmat[k1[j],j] = RVobj$PCor[1,j]
  # compute phi1 and phi2
  phi = matrix(NA,d1,2)
  for(j in 3:d1) phi[j,2] = RVobj$PCor[2,j]*sqrt((1-rmat[k1[j],j]^2)/(1-rmat[ksort[j,1],ksort[j,2]]^2))
  phi[2,1] = rmat[1,2]
  for(j in 3:d1) phi[j,1] = rmat[k1[j],j]-phi[j,2]*rmat[ksort[j,1],ksort[j,2]]

  # compute correlations from partial correlations
  for(j in 3:d1) rmat[k2[j],j] = phi[j,1]*rmat[ksort[j,1],ksort[j,2]]+phi[j,2]
  # compute remaining correlations, maybe this loop can be improved
  for(j in 3:d1)
  { for(i in 1:(j-1))
    { if(is.na(rmat[i,j]))
      { e1min = min(i,k1[j])
        e1max = max(i,k1[j])
        e2min = min(i,k2[j])
        e2max = max(i,k2[j])
        rmat[i,j] = phi[j,1]*rmat[e1min,e1max]+phi[j,2]*rmat[e2min,e2max]
      }
    }
  }
  diag(rmat) = 1
  for(i in 1:(d1-1)) { for(j in (i+1):d1) rmat[j,i]=rmat[i,j] }
  return(rmat)
}

#' Compute correlation matrix according to 1-factor + 1-truncated vine (residual dependence) model 
#'
#' @description
#' Compute correlation matrix according to a 1-factor + 1-truncated vine (for residual dependence) model 
#'
#' @param  cormat dxd correlation matrix
#' @param  start_loading  dx1 loading vector (for latent factor)
#' @param  iter number of iterations for modified EM
#' @param  est "mle" or "mom"
#' @param  plots flag that is TRUE to show plots of EM steps
#' @param  trace flag that is TRUE to print every 100th integer for iter
#'
#' @return components:
#'  loading = final estimate for loading vector; 
#'  R = correlation matrix (from MLE for 1F1T structure); 
#'  Psi = vector of residual variances; 
#'  loadings = matrix where ith row has the ith iteration; 
#'  Rmats = list of correlation matrices; Rmats[[i]] has the ith iteration; 
#'  Rstart = starting value of R based on starting values; 
#'  dists = vector of distance measures as GOF criterion, ith entry for ith iteration; 
#'  incls = iter x d*(d-1)/2 matrix, ith row for ith iteration
#'   columns are whether edge 12, 13, 23, 14, .... (d-1,d) are in residual tree; 
#'  partcor = dxd matrix of partial correlations given latent variable;  
#' loglik = vector of loglik values, ith entry for ith iteration.
#'
#' @examples
#' ##  donttest:
#' library(igraph)
#' data(DJ20142016gf)
#' zdat = dj1416gf$zscore # GARCH-filtered returns that have been transformed to N(0,1)
#' rzmat = cor(zdat)
#' d = ncol(zdat)
#' cat("\n1-factor start for 1f1t\n")
#' fa = factanal(factors=1,covmat=rzmat)
#' start = c(fa$loading)
#' # fitting 1Factor 1Truncated vine residual dependence structure
#' out1f1t = gauss1f1t(rzmat,start_loading=start,iter=20,plots=FALSE)
#' print(out1f1t$loading)
#' cat("\nedges for tree of residual dependence\n")
#' i = 1:d
#' nn = i*(i-1)/2
#' niter = 21  # above iteration bound +1
#' incls = out1f1t$incls[niter,]
#' for(j in 2:d)
#' { for(k in 1:(j-1)) 
#'   { if(incls[nn[j-1]+k]) 
#'     { cat("edge ", k,j," "); cat(out1f1t$partcor[k,j],"\n") } 
#'   }
#' }
#' # extract three columns of this output to use with factor1trvine_nllk
#' # See example in cop1f1t()
#' ## End(donttest)
#'
#' @details
#' A modified EM algorithm is used
#' -- first step of the M-step performed either by MLE or by method of moments
#' -- the second step assumes the moment estimator has been used
#'
#' @references
#' Brechmann EC and Joe H (2014). 
#' Parsimonious parameterization of correlation matrices using truncated
#' vines and factor analysis.
#' Computational Statistics and Data Analysis, 77, 233-251.
#'
#' @export
#'
gauss1f1t = function(cormat, start_loading, iter=10, est="mle", plots=TRUE, trace=TRUE)
{
  d = dim(cormat)[1]
  #loadings = matrix(NA,iter+1,dim(cormat)[1])
  loadings = matrix(NA,iter+1,d)
  Rmats = list()
  dists = rep(NA,iter+1)
  #incls = matrix(NA,iter+1,dim(cormat)[1]*(dim(cormat)[1]-1)/2)
  incls = matrix(NA,iter+1,d*(d-1)/2)
  loglik = rep(NA,iter+1)
  loglik_contr = rep(NA,iter)
  #Psi = rep(NA,dim(cormat)[1])
  Psi = rep(NA,d)
  
  # starting values
  loading = start_loading
  loadings[1,] = loading

  # find initial correlation matrix for starting values
  Rfind = residDep(cormat,start_loading)
  Rstart = Rfind$R
  incls[1,] = Rfind$incl
  Rmat = Rstart
  Rmats[[1]] = Rmat

  # save log likelihood
  loglik[1] = -0.5*log(det(Rmat))-0.5*sum(diag(solve(Rmat,cormat)))
  # save distance to true correlation matrix as goodness-of-fit criterion
  dist0 = log(det(cormat))+dim(cormat)[2]
  dists[1] = -2*loglik[1]-dist0

  for(i in 1:iter)
  { if(trace) { if(i %% 100 == 0) message(paste("iteration",i)) }
    old_loading = loading
    old_Rmat = Rmat
    old_Psi = Psi
    # compute auxiliary quantities
    Rinv = solve(Rmat)
    Sload = Rinv%*%loading
    S_ZF = cormat%*%Sload
    S_FF = as.numeric(t(Sload)%*%cormat%*%Sload + 1 - loading%*%Sload)
    # compute updated loadings
    newton = function(start)
    { invS_ZZ = solve(cormat)
      invR_old = solve(old_Rmat)
      R0load = invR_old%*%old_loading
      func = function(gamma1)
      { g1R0load = gamma1%*%R0load
        Sg1 = invS_ZZ%*%gamma1
        g1Sg1 = gamma1%*%Sg1
        as.numeric(1-g1Sg1)*R0load - as.numeric(S_FF-2*g1R0load+g1Sg1)*Sg1
      }
      grad = function(gamma1)
      { g1R0load = gamma1%*%R0load
        Sg1 = invS_ZZ%*%gamma1
        g1Sg1 = gamma1%*%Sg1
        -2*Sg1 %*% t(R0load) - as.numeric(S_FF-2*g1R0load+g1Sg1)*invS_ZZ - (-2*R0load+2*Sg1) %*% t(Sg1)
      }
      
      newton_iter = 10
      gammas = matrix(NA,newton_iter,d)
      gammas[1,] = start
      for(i in 2:newton_iter) gammas[i,] = gammas[i-1,] - solve(grad(gammas[i-1,]))%*%func(gammas[i-1,])
    
      # diagnostic plot
      #plot(gammas[,1], type="l", ylim=range(gammas))
      #for(i in 2:d) lines(gammas[,i])
      return(gammas[newton_iter,])
    }

    if(est == "mle") { loading = newton(loading) } # ML estimator
    else if(est == "mom")
    { loading = as.vector(S_ZF/S_FF) } # approximative moment estimator
    loadings[i+1,] = loading
    
    # compute idiosyncratic errors
    Psi = diag(cormat - S_ZF%*%t(S_ZF)/S_FF)
    
    # find updated correlation matrix using MST (1-truncated vine)
    Riter = residDep(cormat,loading)
    Rmat = Riter$R
    incls[i+1,] = Riter$incl
    Rmats[[i+1]] = Rmat

    # save log likelihood
    loglik[i+1] = -0.5*log(det(Rmat))-0.5*sum(diag(solve(Rmat,cormat)))
    # save distance
    dists[i+1] = -2*loglik[i+1]-dist0
    # save part that is ignored in optimization of log likelihood
    Rmat_mod = Rmat-loading%*%t(loading)
    loglik_contr[i] = 0.5*(S_FF/(1-loading%*%solve(Rmat)%*%loading)-sum(diag(solve(Rmat_mod,cormat))))
  }
  
  # diagnostic plots
  if(plots == TRUE)
  { cols = rainbow(d)
    plot(loadings[,1], type="l", ylim=range(loadings), xlab="iter", ylab="loading", col=cols[1])
    for(i in 2:d) lines(loadings[,i], col=cols[i])
    plot(dists, type="l", ylim=range(dists,0)*c(1/1.0001,1.0001), xlab="iter", ylab="distance to true correlation matrix")
    abline(h=0, lty=2, col=2)
    true_ll = -0.5*dist0
    plot(loglik, type="l", ylim=range(loglik, true_ll)*c(1/1.0001,1.0001), xlab="iter", ylab="log likelihood")
    abline(h=true_ll, lty=2, col=2)

    dd = d*(d-1)/2
    gg = graph.adjacency(cormat, mode="lower", weighted=TRUE, diag=FALSE)
    E(gg)$name = paste(get.edgelist(gg)[,1],get.edgelist(gg)[,2],sep=",")
    incl_mat = t(incls)*(1:dd)
    plot(incl_mat[1,], ylim=c(1,dd), axes=F, xlab="iter", ylab="", pch=18)
    axis(1)
    axis(2, at=1:dd, labels=E(gg)$name, las=2)
    box()
    for(i in 2:dd) points(incl_mat[i,], pch=18)

    asymp_contr = -(d-1)/2
    plot(loglik_contr, type="l", ylim=range(loglik_contr, asymp_contr)*c(1.0001,1/1.0001), ylab="ignored part of log likelihood")
    abline(h=asymp_contr, lty=2, col=2)
  }
  return(list(loading=loading, R=Rmat,Psi=Psi, loadings=loadings, Rmats=Rmats, 
   Rstart=Rstart,dists=dists,incls=incls,partcor=Riter$partcor,loglik=loglik))
}
