\name{onefactorEstWithProxy}
\alias{onefactorEstWithProxy}
\title{Parameter estimation for 1-factor copula with estimated latent variables using VineCopula::BiCopSeelct 
}
\description{
Parameter estimation for 1-factor copula with estimated latent variables
}

\usage{
onefactorEstWithProxy(udata,vlatent, famset, iprint=FALSE)
}
\arguments{
\item{udata}{nxd matrix with values in (0,1)}
\item{vlatent}{vector is estimated latent variables (or test with known values)}
\item{famset}{2*d vector of codes for copula families for d global linking copulas and d group-based linking copulas, using those from VineCopula: current choices to cover a range of tail behavior are: 1 = Gaussian/normal; 2 = t; 4 = Gumbel; 5 = Frank; 7 = BB1; 10 = BB8; 14 = survival Gumbel; 17 = survival BB1; 20 = survival BB8.}
\item{iprint}{if TRUE print intermediate results}
}
\value{
 list with fam = d-vector of family codes chosen via BiCopSelect;par1 = d-vector; par2 = d-vector of parameters for the selected copula families 
in the 1-truncated vine rooted at the latent variable,
}
\examples{
\dontrun{
# simulate data from 1-factor model with all Frank copulas
n = 500
d = 40
set.seed(20)
cpar = runif(d,4.2,18.5)
param = c(rbind(cpar,rep(0,d))) #Kendall's tau 0.4 to 0.8
data = r1factor(n,d,param,fam=rep(5,d))
vlat = data$vlatent # latent variables
udata = data$udata
proxyMean = uscore(apply(udata,1,mean)) # mean proxy
# RMSE of estimated latent variables
print(sqrt(mean((proxyMean-vlat)^2)))
# first estimation of 1-factor copula parameters
# allow for Frank, gaussian, t linking copulas
est1 = onefactorEstWithProxy(udata,proxyMean, famset=c(1,2,5))
print(est1$fam) # check choices , all 5s (Frank) in this case
print(est1$par1)
# estimation with only Frank copula as choice
est0 = onefactorEstWithProxy(udata,proxyMean, famset=c(5))
print(summary(abs(est0$par1-cpar)))  # same as $est1$par1
# improved conditional expectation proxies
# latentUpdate1factor allows for estimated linking copula with 2-parameters
# latentUpdate1factor1 can be used if estimated linking copulas all have par2=0
condExpProxy = latentUpdate1factor(c(rbind(est1$par1,est1$par2)),
udata=udata,nq=25,family=rep(5,d))
# improved estimation of 1-factor copula parameters
est2 = onefactorEstWithProxy(udata,condExpProxy, famset=est1$fam)
print(est2$par1)  
# simple version of update for 1-parameter linking copulas
condExpProxy1 = latentUpdate1factor1(est0$par1,
udata=udata,nq=25,family=rep(5,d))
summary(condExpProxy-condExpProxy1)  # 0 because family was chosen as 5 
print(summary(abs(est2$par1-cpar)))
# rmse of estimated latent variables
print(sqrt(mean((condExpProxy-vlat)^2))) 
# smaller rmse than initial proxies
}
}
\references{
1. Krupskii P and Joe H (2013).
Factor copula models for multivariate data.
Journal of Multivariate Analysis, 120, 85-101.
2. Fan X and Joe H (2024).
High-dimensional factor copula models with estimation of latent variables
Journal of Multivariate Analysis, 201, 105263.
}
\details{
It is best if variables have been oriented to be positively related
to the latent variable
}
