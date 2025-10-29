\name{mvtBifact_nllk}
\alias{mvtBifact_nllk}
\title{negative log-likelihood for the bi-factor Gaussian/t model}
\description{
negative log-likelihood in the bi-factor Gaussian or t model
}

\usage{
mvtBifact_nllk(rhovec,grsize,tdata,df)
}
\arguments{
\item{rhovec}{vector of length d*2 for partial correlation representation of loadings, first d correlations with common factor, then partial correlations with group factor given common factor}
\item{grsize}{vector of group sizes for bi-factor model}
\item{tdata}{nxd data set of scores for t df (e.g., qt(u,df))}
\item{df}{degree of freedom parameter (positive)}
}
\value{
 negative log-likelihood (nllk) of copula for mvt bi-factor model}
