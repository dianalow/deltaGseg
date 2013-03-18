\name{getIntervals}
\alias{getIntervals}
\title{
getIntervals
}
\description{
Helper function to retrieve subpopulations and computes the intervals for each subpopulation after segmentation and clustering.
}
\usage{
getIntervals(object)
}

\arguments{
  \item{object}{ An object of class 'SegSeriesTrajectories'.
}
}
\details{
None.
}
\value{ Returns a list of subpopulations and the intervals.
}
\references{
None
}
\author{
Diana H.P. Low, Efthimios Motakis
}
\seealso{
\code{\link{clusterSegments}}
}
\examples{
data(deltaGseg)
getIntervals(traj1.ss)
}
