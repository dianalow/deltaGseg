\name{splitTraj}
\alias{splitTraj}
\title{
Identify Breakpoints in a Trajectory
}
\description{
splitTraj determines the breakpoints to split a given trajectory into a user specified number of segments. This analysis is performed for very long series (more than 20,000 time points) in order to avoid any memory allocation problems in R. Alternatively, it can be used for manual splitting of the series (see \code{transformSeries} with method="override_splitting")
}
\usage{
splitTraj(object, segsplits = rep(5,length(object@filenames)))
}

\arguments{
  \item{object}{ An object of class "Trajectories".
}
  \item{segsplits}{ Numeric vector. The number of breakpoints. The length of segsplits must equal the number of trajectory series in the Trajectories object. Each value specifies the number of splits we want to impose in our long series in order to make it shorter.
}
}
\details{
The output of the function is the estimated points that split the series into smaller sub-series. Typically, the plotted series and the estimated splits are further inspected using plots.
}
\value{ A numeric list of length equal to number of trajectory series, containing the breakpoints for each series.
}
\references{
None
}
\author{
Diana H.P. Low, Efthimios Motakis
}

\examples{
data(deltaGseg)
splitTraj(traj1)
}
