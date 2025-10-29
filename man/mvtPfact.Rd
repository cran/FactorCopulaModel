\name{mvtPfact}
\alias{mvtPfact}
\title{MLE in a MVt model with a p-factor correlation structure}
\description{
MLE in a MVt model with a p-factor correlation structure
}

\usage{
mvtPfact(tdata,start,pfact,df,prlevel=0,mxiter=100)
}
\arguments{
\item{tdata}{nxd matrix of t-scores}
\item{start}{vector of length p*d with starting values for partial correlation representation of loadings (algberaically independent)}
\item{pfact}{number of factors (such as 1,2,3,... )}
\item{df}{degree of freedom parameter >0}
\item{prlevel}{print.level for nlm()}
\item{mxiter}{maximum number of iterations for nlm()}
}
\value{
  nlm object with ($code,$estimate,$gradient,$iterations,$minimum) Note the minimum nllk can be the same for different parameter vectors
because of invariance of the loading matrix to rotations
}
\examples{
# See example in mvtBifact()
}
