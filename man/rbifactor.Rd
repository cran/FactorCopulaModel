\name{rbifactor}
\alias{rbifactor}
\title{simulate from bi-factor copula model}
\description{
simulate from bi-factor copula model  and include corresponding latent variables
}

\usage{
rbifactor(n, grsize, cop=5, param) 
}
\arguments{
\item{n}{sample size}
\item{grsize}{G-vector of group sizes for G groups }
\item{cop}{code for copula families  1: Gaussian/Gaussian; 2: t/t; 4: Gumbel/Gumbel; 5: Frank/Frank; 7: BB1/Frank; 14: survival Gumbel, 17: survivalBB1 /Frank if cop = 1, data have standard normal marginals if cop = 2, data have t marginals if cop > 2, data have uniform(0,1) marginals}
\item{param}{vector of parameters (those for the common factor go first) The order in param is the same as in start for mvtbifct(full=T) function. For BB1/Frank: BB1thetas then BB1deltas, then Frank parameters -- the order in param is the same as in start for mvtbifct(full=F) function}
}
\value{
 list with data: nxd data set with U(0,1) or N(0,1) or t(df) margin; v0: n-vector of corresponding global latent variables; and
vg: nxG matrix of corresponding local (group) latent variables.
}
\examples{
grsize = c(4,3)
cop = 4; param4 = c(seq(1.5,2.1,0.1),  rep(1.1,7))
cop = 14; param14 = c(seq(1.5,2.1,0.1),  rep(1.1,7))
cop = 5; param5 = c(seq(1.5,2.1,0.1),  rep(1.1,7))
cop = 1; param1 = c(0.5,0.6,0.7,0.8,0.9,0.4,0.5,  rep(1.1,7)) 
cop = 2; param2 = c(0.5,0.6,0.7,0.8,0.9,0.4,0.5,  rep(1.1,7), 7)
cop = 7; param7 = c(seq(0.5,1.1,0.1), 1.5,1.6,1.7,1.8,1.9,1.4,1.5, rep(1.1,7)) 
cop = 17; param17 = c(seq(0.5,1.1,0.1), 1.5,1.6,1.7,1.8,1.9,1.4,1.5, rep(1.1,7)) 
set.seed(123)
gumdat = rbifactor(n=10, grsize=grsize, cop=4, param=param4)    # U(0,1)
gumrdat = rbifactor(n=10, grsize=grsize, cop=14, param=param14) # U(0,1)
frkdat = rbifactor(n=10, grsize=grsize, cop=5, param=param5)    # U(0,1)
gaudat = rbifactor(n=10, grsize=grsize, cop=1, param=param1)    # N(0,1)
bvtdat = rbifactor(n=10, grsize=grsize, cop=2, param=param2)    # t_7
bb1frkdat = rbifactor(n=10, grsize=grsize, cop=7, param=param7)    # U(0,1)
bb1rfrkdat = rbifactor(n=10, grsize=grsize, cop=17, param=param17) # U(0,1)
summary(bb1frkdat$data)
summary(bb1frkdat$v0)
summary(bb1frkdat$vg)
}
\details{
The user can modify this code to get other linking copulas.
}
