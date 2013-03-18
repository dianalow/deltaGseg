\name{getSegments}
\alias{getSegments}
\title{
getSegments
}
\description{
Accessor of SegTrajectories or SegSeriesTrajectories object to retrieve segment matrix.
}
\usage{
getSegments(object)
}
\arguments{
  \item{object}{Object of class "SegTrajectories"" or "SegSeriesTrajectories"
}
}

\value{A matrix with segment information including quantiles.
}
\references{
None.
}
\author{
Diana H.P. Low, Efthimios Motakis
}
\note{
None.
}
\seealso{
\code{\link{denoiseSegments}}, \code{\link{clusterSegments}}
}
\examples{
data(deltaGseg)
segments<-getSegments(traj1.denoise)
}

