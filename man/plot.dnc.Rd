\name{plot.dnc}
\alias{plot.dnc}
\title{
Plot a \command{dnc} object.
}
\description{
Construct 2D or 3D plots of the latent actor positions using the \command{plot3D} and \command{plot3Drgl} packages.
}
\usage{
\S3method{plot}{dnc}(x,aggregated=TRUE,plotRGL=TRUE,
                    Lines=TRUE,colByComm=TRUE,
                    INDEX=1:min(dim(x$pm$X)[1],3),...)
}
\arguments{
\item{x}{A \command{dnc} class object}
\item{aggregated}{Logical.  Should the time points be aggregated, or should separate plots be made for each time point?}
\item{plotRGL}{Logical.  Should \command{rgl} interactive plot(s) be made?  Ignored if \command{length(INDEX)==2}.}
\item{Lines}{Logical.  Should lines be drawn from the origin to the center of the communities?}
\item{colByComm}{Logical.  Should the actors be colored according to community membership?}
\item{INDEX}{Vector.  Which subset of the \command{p} dimensions should be plotted? \bold{Note: The length of this vector determines whether a 2D or 3D plot should be made.}}
\item{...}{Further arguments to be passed into \command{points()}, \command{surf3D} or \command{points3D} functions.}
}
\details{
This function utlizies functions from the \command{plot3D} and \command{plot3Drgl} packages for 3D functionality.  To see what options you may pass into \command{...}, see documentation for \command{points} in base R, \command{plot3D::surf3D}, or \command{plot3D::points3D}.
}
\examples{
  data(friendship)
  set.seed(123)
  dncObj <- dnc(friendship,M=4,p=3,method="Gibbs",
                controls=list(nDraws=250,burnin=50,
                              MaxItStg2=25,epsilonStg2=1e-15))
  par(mar=rep(0,4)+0.05)
  plot(dncObj,plotRGL=FALSE,pch=16,phi=60,lwd=2,cex=1.5)
}