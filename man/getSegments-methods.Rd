\name{getSegments-methods}
\docType{methods}
\alias{getSegments-methods}
\alias{getSegments,SegSeriesTrajectories-method}
\alias{getSegments,SegTrajectories-method}
\title{getSegments}
\description{
Returns a matrix from the @smatrix slot (for object of class "SegTrajectories") or @ssmatrix slot (for object of class "SegSeriesTrajectories")
}
\section{Methods}{
\describe{

\item{\code{signature(object = "SegSeriesTrajectories")}}{
Matrix is a result of the function \code{\link{clusterSegments}}
}

\item{\code{signature(object = "SegTrajectories")}}{
Matrix is a result of the function \code{\link{denoiseSegments}}
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
#getSegments(traj1.denoise)
}
\keyword{methods}
