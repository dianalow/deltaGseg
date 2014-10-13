\name{getTraj}
\alias{getTraj}
\title{
getTraj
}
\description{
Accessor of Trajectories or TransTrajectories object to retrieve trajectories.
}
\usage{
getTraj(object)
}

\arguments{
  \item{object}{An object of class "Trajectories" or "TransTrajectories"
}
}
\details{
None.
}
\value{ List of length equal to number of trajectories, each containing matrix with 2 columns. 1: time points, 2: free energy values.
}
\references{
None.
}
\author{
Diana H.P. Low, Efthimios Motakis
}
\examples{
data(deltaGseg)
alltrajectories<-getTraj(traj1)
}
