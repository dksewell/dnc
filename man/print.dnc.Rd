\name{print.dnc}
\alias{print.dnc}
\title{
Print \command{dnc} object
}
\description{
Print description of a \command{dnc} object.
}
\usage{
\S3method{print}{dnc}(x,printDens=FALSE,...)
}
\arguments{
\item{x}{A \command{dnc} class object.}
\item{printDens}{Logical.  Should the density of the network at each time point be computed and printed?}
\item{...}{additional optional arguments.}
}
\examples{
  data(friendship)
  set.seed(123)
  dncObj <- dnc(friendship,M=4,p=3,method="Gibbs",
                controls=list(nDraws=250,burnin=50,
                              MaxItStg2=25,epsilonStg2=1e-15))
  print(dncObj)
}