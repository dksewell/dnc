\name{BIC.dnc}
\alias{BIC.dnc}
\title{
Compute BIC as in Handcock et al. 2007
}
\description{
The following uses a BIC estimate of \eqn{\pi(Y,\hat{X}|M)} to perform model selection.  Note that this usage is not the typical BIC encountered in simpler contexts.
}
\usage{
\S3method{BIC}{dnc}(object, ...)
}
\arguments{
\item{object}{A \command{dnc} object, a result of running \command{dnc(...)}}
\item{...}{optional additional arguments. None are used.}
}
\details{
Rather than estimating the integrated likelihood \eqn{\pi(Y|G)}, this instead incorporates the MAP estimates of the latent positions and corresponds to \eqn{\pi(Y,\hat{X}|M)}.  The BIC value returned is the following sum:
\deqn{-2 log(\pi(Y|\hat{X},\hat{\theta_1})) + dim(\theta_1)log(\sum y_{ijt})
-2 log(\pi(\hat{X}|\hat{\theta_2})) + dim(\theta_2) \log(nT)}.
See Sewell and Chen (2016) for more details.
}
\value{
A scalar.  Lower values are better.
}
\references{
Handcock, M. S., A.E. Raftery, and J. M. Tantrum (2007).  Model-based clustering for social networks. J.R. Statist. Soc. A, 170, p. 301-354.
}
\examples{
\dontrun{
dncObjs = list()
BICvals = numeric(9)
for(i in 2:10){
  print(i)
  dncObjs[[i]] <- dnc(friendship,M=i,p=3,method="VB",
                      controls=list(nDraws=500,burnin=100,
                                    MaxItStg2=50,epsilonStg2=1e-15))
  BICvals[i-1] <- BIC(dncObjs[[i]])
}
plot(BICvals~c(2:10),type="b",pch=16,
     xlab="Number of communities",ylab="BIC value")
( MBest = which.min(BICvals)+1 )
abline(v=MBest,lty=2,col="blue")
}
}