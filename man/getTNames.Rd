\name{getTNames}
\alias{getTNames}
\title{
getTNames
}
\description{
Accessor of Trajectories and TransTrajectories object to retrieve filename(s) used in trajectory computation.
}
\usage{
getTNames(object)
}

\arguments{
  \item{object}{ An object of class 'Trajectories' or 'TransTrajectories'.
}
}
\details{
None.
}
\value{ Returns character vector of filenames.
}
\references{
None.
}
\author{
Diana H.P. Low, Efthimios Motakis
}

\examples{
data(deltaGseg)
getTNames(traj1)
}
