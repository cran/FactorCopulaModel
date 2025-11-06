\name{mvtBifact}
\alias{mvtBifact}
\title{MLE for multivariate normal/t with a bi-factor or nested factor correlation structure}
\description{
MLE for the bi-factor or nested factor structure for multivariate normal/t
}

\usage{
mvtBifact(tdata,start,grsize,df,prlevel=2,model="bifactor",mxiter=100)
}
\arguments{
\item{tdata}{nxd matrix of t-scores or z-scores}
\item{start}{vector of length 2*d with starting values of partial correlations values for correlations of observed Z_[gj] and common latent V_0 go first, then partial correlations of Z_[gj] and V_g given V_0 (j in group g)}
\item{grsize}{vector of group sizes for bi-factor model}
\item{df}{degrees of freedom parameter >0}
\item{prlevel}{print.level for nlm()}
\item{model}{"bifactor" or "nestfactor" nested-factor is reduced model with fewer parameters}
\item{mxiter}{maximum number of iterations for nlm()}
}
\value{
 nlm object with ($code,$estimate,$gradient,$iterations,$minimum)}
\examples{
\donttest{
data(rainstorm)
udat = rainstorm$uprecip
d = ncol(udat)
grsize = rainstorm$grsize
df = 10
tdata = qt(udat,df)
bif = mvtBifact(tdata, c(rep(0.8,d),rep(0.2,d)), grsize, df=df,
prlevel=1, model="bifactor", mxiter=100)
#
# nested-factor: parameters for group latent linked to global latent
# come in the first tree of the 2-truncated vine.
nestf = mvtBifact(tdata, c(0.7,0.7,0.7,rep(0.85,d)), grsize, df=df,
prlevel=1, model="nestfactor", mxiter=100)
# doesn't converge properly, group2 latent matches global latent
#
st1 = rep(0.7,d)
out1t = mvtPfact(tdata,st1,pfact=1,df=df,prlevel=1)
st2 = rep(0.7,2*d)
out2t = mvtPfact(tdata,st2,pfact=2,df=df,prlevel=1)
st3 = c(rep(0.7,grsize[1]),rep(0.1,d),rep(0.7,grsize[2]),rep(0.1,d),rep(0.7,grsize[3]))
out3t = mvtPfact(tdata,st3,pfact=3,df=df,prlevel=1)
cat(bif$minimum, nestf$minimum, out1t$minimum, out2t$minimum, out3t$minimum,"\n")
}
}
\details{
Note the minimum nllk can be the same for different parameter vectors
if some group size values are 1 or 2.
}

