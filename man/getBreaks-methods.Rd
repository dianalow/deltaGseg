\name{getBreaks-methods}
\docType{methods}
\alias{getBreaks-methods}
\alias{getBreaks,TransTrajectories-method}
\title{getBreaks}
\description{
Returns a list from the @breakpoints slot for object of class "TransTrajectories")
}
\section{Methods}{
\describe{

\item{\code{signature(object = "SegSeriesTrajectories")}}{
List of breakpoints per series after transformation of series data using \code{transformSeries}
}

}}
\keyword{methods}
\references{
None
}
\author{
Diana H.P. Low, Efthimios Motakis
}
\seealso{
\code{\link{transformSeries}}
}
\examples{
#data(deltaGseg)
#getBreaks(traj1.tr)
}
\keyword{methods}
