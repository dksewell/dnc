\name{dnc-package}
\alias{dnc-package}
\docType{package}
\title{
Community detection for dynamic networks.
}
\description{
Dynamic network clustering/community detection using a latent space approach.  Using temporal edge data, network actors are embedded onto a hypersphere and grouped based on direction.
}
\details{
\tabular{ll}{
Package: \tab dnc\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-05-05\cr
License: \tab GPL (>= 2)\cr
}
This package can perform community detection on dynamic (temporal) network data observed at discrete time points.  Communities are assumed fixed, but community membership may change.  The main function is \command{dnc(\ldots)} which can perform variational Bayes estimation or alternatively implement a Gibbs sampler.  A \command{dnc} object is the output, for which there exists the following generic commands: \command{simulate()}, \command{plot()}, \command{print()}, and \command{BIC()}.  Ignorable (MAR, MCAR) missing edge data can be incorporated into the Gibbs sampler.
}
\author{
Daniel K. Sewell

Maintainer: Daniel K. Sewell <daniel-sewell@uiowa.edu>
}
\references{
Sewell, D. K., and Chen, Y. (2016). Latent Space Approaches to Community Detection in Dynamic Networks. Bayesian Analysis. doi: 10.1214/16-BA1000. http://projecteuclid.org/euclid.ba/1461603847
}
\examples{
\dontrun{
VB5 = dnc(Y,M=5,p=3,"VB")
Gibbs5 = dnc(Y,M=5,p=3,"Gibbs")
print(VB5)
print(Gibbs5)
BIC(VB5)
BIC(Gibbs5)
plot(VB5)
plot(Gibbs5,aggregated=FALSE,plotRGL=FALSE)
}
}
