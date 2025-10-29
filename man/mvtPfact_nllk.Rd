\name{mvtPfact_nllk}
\alias{mvtPfact_nllk}
\title{negative log-likelihood for the p-factor Gaussian/t model}
\description{
negative log-likelihood in the p-factor Gaussian or t model
}

\usage{
mvtPfact_nllk(rhvec,tdata,df)
}
\arguments{
\item{rhvec}{vector of length d*p with partial correlation representation of loadings}
\item{tdata}{nxd data set of t(df)-scores}
\item{df}{degree of freedom parameter >0}
}
\value{
 negative log-likelihood of copula for mvt p-factor model}
