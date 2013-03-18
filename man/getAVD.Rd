\name{getAVD}
\alias{getAVD}
\title{
getAVD
}
\description{
Accessor of Trajectories or TransTrajectories object to retrieve adf p-values.
}
\usage{
getAVD(object)
}
\arguments{
  \item{object}{ An object of class 'Trajectories' or 'TransTrajectories'.
}
}

\value{ Returns a vector (of pvalues) for each trajectory defined in the object.
}
\author{
Diana H.P. Low, Efthimios Motakis
}
\note{
None.
}
\seealso{
\code{\link{parseTraj}},\code{\link{transformSeries}}
}
\examples{
data(deltaGseg)
getAVD(traj1)
}