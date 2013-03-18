\name{transformSeries-methods}
\docType{methods}
\alias{transformSeries-methods}
\alias{transformSeries,Trajectories-method}
\title{transformSeries}
\description{
Returns a matrix from the @smatrix slot (for object of class "SegTrajectories") or @ssmatrix slot (for object of class "SegSeriesTrajectories")
}
\section{Methods}{
\describe{

\item{\code{signature(object = "Trajectories")}}{
Returns an object of class "TransTrajectories"
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
#breakpoints<-chooseBreaks(all_breakpoints,3)
#tt<-transformSeries(traj1,breakpoints=breakpoints)
}
\keyword{methods}
