\name{getBreaks}
\alias{getBreaks}
\title{
getBreaks
}
\description{
Accessor of TransTrajectories object to retrieve computed breakpoints.
}
\usage{
getBreaks(object)
}

\arguments{
  \item{object}{ An object of class 'TransTrajectories'.
}
}
\details{
None.
}
\value{ Returns a list (of numerical breakpoint values) for each trajectory defined in an object of class 'Trajectories'.
}
\references{
None
}
\author{
Diana H.P. Low, Efthimios Motakis
}
\seealso{
\code{\link{splitTraj}},\code{\link{transformSeries}}
}
\examples{
data(deltaGseg)
getBreaks(traj1.tr)
}
