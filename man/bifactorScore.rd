\name{bifactorScore}
\alias{bifactorScore}
\title{Proxies for bi-factor copula model based on Gaussian bi-factor score}
\description{
Proxies (in (0,1)) for bi-factor copula model based on Gaussian bi-factor score
}

\usage{
bifactorScore(udata, start, grsize, prlev=1)
}
\arguments{
\item{udata}{nxd matrix in (0,1); n is sample size, d is dimension}
\item{start}{starting values for fitting the bi-factor Gaussian model}
\item{grsize}{G-vector with group sizes with G groups}
\item{prlev}{printlevel in call to nlm}
}
\value{
 list with  Aloadmat=estimated loading matrix after N(0,1) transform; proxies = nx(G+1) matrix with stage 1 proxies for the latent variables
(for global latent in column 1, and then for group latent);
weight = weight matrix based on the correlation matrix of normal score and the loading matrix.
}
\examples{
#See example in bifactorEstWithProxy()
}
