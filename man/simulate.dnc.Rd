\name{simulate.dnc}
\alias{simulate.dnc}
\title{
Simulate DNC
}
\description{
Use an existing \command{dnc} object to simulate new network data.
}
\usage{
\S3method{simulate}{dnc}(object,nsim=1,seed=NULL,...)
}
\arguments{
\item{object}{A \command{dnc} class object}
\item{nsim}{Integer.  Number of simulations to perform.}
\item{seed}{Optional numeric.  If \command{is.null(seed)} is \command{FALSE} then \command{seed} is used in a call to \command{set.seed} before simulating the network data.  The default, \command{NULL}, will not change the random generator state.}
\item{...}{additional optional arguments.}
}
\details{
This function uses the posterior mean of the latent positions and model parameters to simulate new network data sets.  
}
\value{
A n x n x T x \command{nsim} array.
}
\examples{
    data(friendship)
    set.seed(123)
    dncObj <- dnc(friendship,M=4,p=3,method="Gibbs",
                  controls=list(nDraws=250,burnin=50,
                                MaxItStg2=25,epsilonStg2=1e-15))
  newSims = simulate(dncObj,10,123)
}