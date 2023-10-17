\name{dnc}
\alias{dnc}
\title{
Dynamic Network Clustering
}
\description{
Perform dynamic network clustering.  Either variational Bayes or a Gibbs sampler may be implemented.  Setting \command{M=0} performs variational Bayes with no clustering.  Returns posterior parameters (if \command{method="VB"}) or approximate posterior samples (if \command{method="Gibbs"}), as well as the MAP estimates, which may be extracted through \command{dncObj$pm}.
}
\usage{
dnc(Y,M,p=3,method="VB",init=NULL,hyperparms=NULL,Missing=NULL,
    controls=list(MaxIt=500,epsilon=1e-5,MaxItStg2=100,
                  epsilonStg2=1e-15,nDraws=10000,burnin=1000))
}
\arguments{
\item{Y}{Dynamic network data.  This should be in the form of a n x n x T array of 1's and 0's.  Each slice corresponds to a single time point.}
\item{M}{Number of communities (may be zero).}
\item{p}{Dimension of the latent space.}
\item{method}{Method of estimation, either "VB" for variational Bayes, or "Gibbs" for a Gibbs sampler.}
\item{init}{\emph{(Use of this argument is not recommended)} Initial values of the parameters.  A named list containing \command{EOm}, \command{mu}, \command{Sig}, \command{Bi0g}, \command{Bitbar}, \command{Bithk}, \command{Er}, \command{Er2}, \command{ai2}, \command{bi2}, \command{nu}, \command{a3}, \command{b3}, \command{Es}, \command{Es2}, and \command{Gam}. }
\item{hyperparms}{Hyperparameters.  A named list with  \command{cc}, \command{a0Star}, \command{b0Star}, \command{a2Star}, \command{b2Star}, \command{b3Star}, \command{GamStar}.}
\item{Missing}{A matrix whose rows correspond to missing dyads.  \command{Missing} should have three columns: row, column, and time (i.e., the indices for the NA's in \command{Y}).  May be left as \command{NULL} if the missing dyads in \command{Y} are \command{NA}'s.}
\item{controls}{A list of values to control the algorithm.  
\describe{
\item{MaxIt}{The total number of iterations for the VB algorithm.  Ignored if \command{method="Gibbs"} unless \command{M=0}.}
\item{epsilon}{Relative tolerance criteria for evaluating convergence.}
\item{MaxItStg2}{The total number of iterations for the second stage initialization of the VB algorithm/Gibbs sampler.  Ignored if \command{M=0}.}
\item{epsilonStg2}{Relative tolerance criteria for evaluating convergence for the second stage initialization of the VB algorithm/Gibbs sampler.  Ignored if \command{M=0}.}
\item{nDraws}{Total number of post-burn-in samples to be drawn via the Gibbs sampler.  Ignored if \command{method="VB"}.}
\item{burnin}{The number of burn-in samples.  Ignored if \command{method="VB"}.}
} }
}
\details{
This function performs community detection according to the model
\deqn{ logit(P(Y_{ijt} =1)) = \alpha + s_j X_{it}'X_{jt}},
\deqn{\pi(X_{it}|Z_{it}=m) = N(r_i {\bf u}_m,\tau^{-1} I_p)}{\pi(X_{it}|Z_{it}=m) = N(r_{i}*u_{m},\tau_{i}^{-1}I_p).}
While the latent positions, \eqn{X_{it}}'s, live in a p-dim Euclidean space, it is more natural to conceptualize these as living on a (hyper-) sphere with the magnitude of the \eqn{X_{it}}'s as attached attributes that reflect the actors' individual tendency to send and receive edges.

If \command{M=0}, then the prior on \eqn{X_{it}} is given by
\deqn{\pi(X_{i1}) = N(0,\sigma^2 I_p)}
\deqn{\pi(X_{it}|X_{i(t-1)}) = N(X_{i(t-1)},\tau_i^{-1} I_p)}

The variational Bayes approach is typically faster than the Gibbs sampler, but tends to underestimate the spread of the posterior.

Currently, only VB is implemented when \command{M=0} (no clustering), hence \command{method} will be ignored if \command{M=0}.

Ignorable missing data can be estimated within the Gibbs sampler (not using the VB algorithm) by adding the extra step of drawing the missing edges given the latent positions and the model parameters at each iteration.

Using the \command{init} is, in general, strongly discouraged, as this may have a non-negligible negative affect on the performance of the VB or the length of the chain needed to reach convergence.  Unless otherwise specified, both the initialization scheme and the hyperparameters are chosen according to Sewell and Chen (2016).
}
\value{
An object of class \command{dnc}, for which other methods exist (e.g., \command{methods(class="dnc")}).

If \command{method="VB"} and \command{M=0},
\describe{
\item{method}{The estimation algorithm}
\item{Y}{The original data}
\item{mu}{A p x T x n array: Posterior mean of the latent positions}
\item{Sig}{A (Tp) x p x n array: Posterior covariance matrices of the latent positions.  The covariance matrix for \eqn{X_{it}} is \command{dncObj$Sig[(t-1)*p,,i]}}
\item{a0}{Scalar: Posterior shape parameter for \eqn{\sigma^2} in inverse gamma distribution (if \command{M=0}).}
\item{b0}{Scalar: Posterior scale parameter for \eqn{\sigma^2} in inverse gamma distribution (if \command{M=0}).}
\item{ai1}{A n x 1 vector: Posterior mean parameter for the \eqn{r_i}'s in truncated normal distribution (if \command{M>0}).}
\item{bi1}{A n x 1 vector: Posterior variance parameter for the \eqn{r_i}'s in truncated normal distribution (if \command{M>0}).}
\item{Er}{A n x 1 vector: Posterior first moment for the \eqn{r_i}'s (if \command{M>0}).}
\item{Er2}{A n x 1 vector: Posterior second moment for the \eqn{r_i}'s (if \command{M>0}).}
\item{ai2}{A n x 1 vector: Posterior shape parameter for the \eqn{\tau_i}'s in gamma distribution.}
\item{bi2}{A n x 1 vector: Posterior scale parameter for the \eqn{\tau_i}'s in gamma distribution.}
\item{a3}{Scalar: Posterior mean for \eqn{\alpha}.}
\item{b3}{Scalar: Posterior variance for \eqn{\alpha}.}
\item{ai4}{A n x 1 vector: Posterior mean parameter for the \eqn{s_j}'s in truncated normal distribution.}
\item{bi4}{A n x 1 vector: Posterior variance parameter for the \eqn{s_j}'s in truncated normal distribution.}
\item{Es}{A n x 1 vector: Posterior first moment for the \eqn{s_j}'s.}
\item{Es2}{A n x 1 vector: Posterior second moment for the \eqn{s_j}'s.}
\item{nu}{A M x p matrix: Posterior mean directions for the \command{M} clusters/communities, i.e., for the \eqn{u_m}'s (if \command{M>0}).}
\item{kappa}{A M x 1 vector: Posterior concentration parameters for the \command{M} clusters/communities, i.e., for the \eqn{u_m}'s (if \command{M>0}).}
\item{Z}{A n x T matrix: Cluster assignments based on the maximum posterior probabilities, computed marginally at each time point (if \command{M>0}).}
\item{Bi0g}{A n x M matrix: Posterior probabilities of community assignment for each actor at the first observed time point (if \command{M>0}).}
\item{Bithk}{A (MT) x M x n array: Posterior transition probability matrices; \eqn{\pi(Z_{itk}=1|Z_{i(t-1)h}=1,Y)=} \command{dncObj$Bithk[(t-1)*M+h,k,i]}.  Ignore first \command{M} lines (internal use only). (if \command{M>0}).}
\item{Bitbar}{A T x M x n array: Marginal posterior probabilities of community assignments, i.e., \eqn{\pi(Z_{itk}=1|Y)=} \command{dncObj$Bitbar[t,k,i]} (if \command{M>0}).}
\item{Gam}{A (M+1) x M matrix: Posterior concentration parameters for \eqn{\beta_0} (row 1) and for \eqn{\beta_m, m>1} (rows 2 to M+1) in Dirichlet distribution  (if \command{M>0}).}
}

If \command{method="Gibbs"},
\describe{
\item{method}{The estimation algorithm}
\item{Y}{The original data}
\item{X}{A p x T x n x \command{nDraws} array: Posterior samples for the latent positions.}
\item{r}{A n x \command{nDraws} matrix: Posterior samples for the \eqn{r_i}'s.}
\item{tau}{A n x \command{nDraws} matrix: Posterior samples for the \eqn{\tau_i}'s.}
\item{alpha}{A \command{nDraws} x 1 vector: Posterior samples for \eqn{\alpha}.}
\item{s}{A n x \command{nDraws} matrix: Posterior samples for the \eqn{\tau_i}'s.}
\item{u}{A M x p x \command{nDraws} array: Posterior draws for the communities, i.e., the \eqn{u_m}'s.}
\item{Z}{A n x T x \command{nDraws} array: Posterior draws for the community assignments for each actor at each time point.}
\item{beta}{A (M+1) x M x n array: Posterior draws for \eqn{beta_0} (row 1) and \eqn{\beta_m, m>1} (rows 2 to M+1).}
\item{posterior}{A (\command{burnin}+\command{nDraws}) x 1 vector: Posterior values for all iterations of the Gibbs sampler.}
\item{Missing}{A matrix of four columns: The row, column, and time for each missing dyad, as well as the posterior probability that the dyad equals one.}
}

Additionally, each \command{dnc} class object comes with a \command{$pm} value, which is a list of the MAP estimates for \command{alpha}, \command{X}, \command{s}, \command{tau}, \command{r}, \command{u}, \command{Z}, and \command{beta}.
}
\references{
Sewell, D. K., and Chen, Y. (2016). Latent Space Approaches to Community Detection in Dynamic Networks. Bayesian Analysis. doi: 10.1214/16-BA1000. http://projecteuclid.org/euclid.ba/1461603847
}
\examples{
  data(friendship)
  set.seed(123)
  dncObj <- dnc(friendship,M=4,p=3,method="Gibbs",
                controls=list(nDraws=250,burnin=50,
                              MaxItStg2=25,epsilonStg2=1e-15))
  print(dncObj)
  BIC(dncObj)
  par(mar=rep(0,4)+0.05)
  plot(dncObj,plotRGL=FALSE,pch=16,phi=60,lwd=2,cex=1.5)
}