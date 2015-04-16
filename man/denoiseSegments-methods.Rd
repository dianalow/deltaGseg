\name{denoiseSegments-methods}
\docType{methods}
\alias{denoiseSegments-methods}
\alias{denoiseSegments,Trajectories-method}
\alias{denoiseSegments,TransTrajectories-method}
\title{denoiseSegments}
\description{
Performs the function \code{denoiseSegments} on an object of either class "Trajectories" or "TransTrajectories" (classUnion="TrajORTransTraj").
}
\section{Methods}{
\describe{

\item{\code{signature(object = "Trajectories")}}{
Returns an object of class "SegTrajectories"
}
\item{\code{signature(object = "TransTrajectories")}}{
Returns an object of class "SegTrajectories"
}
}}
\keyword{methods}
\references{
None
}
\author{
Diana H.P. Low, Efthimios Motakis
}

\examples{
#data(deltaGseg)
#dt<-denoiseSegments(tt)
}
\keyword{methods}
