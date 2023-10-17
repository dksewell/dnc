\name{dncMAP}
\alias{dncMAP}
\title{
MAP estimates
}
\description{
(Re-)Compute the maximum a posteriori (MAP) estimates for a \command{dnc} object.  Primarily intended as an internal function called from \command{dnc(...)}.
}
\usage{
dncMAP(dncObj)
}
\arguments{
\item{dncObj}{A \command{dnc} class object}
}
\value{
A list with the following elements:
\describe{
\item{alpha}{Intercept}
\item{X}{Latent positions}
\item{s}{Receiver scaling effects}
\item{tau}{Latent precision parameters}
\item{sigma2}{Initial latent position variance (if \command{M=0} only)}
\item{r}{Magnitude of mean latent positions (if \command{M>0} only)}
\item{u}{Unit vectors defining communities (if \command{M>0} only)}
\item{Z}{Community assignments (if \command{M>0} only)}
\item{beta}{Initial cluster assignment probabilities (row 1) and transition probability matrix (rows 2 to \command{M}) (if \command{M>0} only)}
}
}
%\examples{
%\dontrun{
%dncObj$pm <- dncMAP(dncObj)
%}
%}