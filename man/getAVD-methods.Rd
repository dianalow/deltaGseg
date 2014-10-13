\name{getAVD-methods}
\docType{methods}
\alias{getAVD-methods}
\alias{getAVD,Trajectories-method}
\alias{getAVD,TransTrajectories-method}
\title{\code{getAVD}}
\description{
Returns numeric vector from the slot @avd in an object of class "Trajectories", or @tavd in an object of class "TransTrajectories."
}
\section{Methods}{
\describe{

\item{\code{signature(object = "Trajectories")}}{
Numeric vector of length equal to series length.
}

\item{\code{signature(object = "TransTrajectories")}}{
Numeric vector of length equal to series length.
}
}}
\seealso{
\code{\link{getAVD}},\code{\link{parseTraj}},\code{\link{transformSeries}}
}
\keyword{methods}
