\name{chooseBreaks}
\alias{chooseBreaks}
\title{
Choose a given number of breakpoints in a set of trajectories.
}
\description{
This function automatically chooses a subset of breakpoints from all the estimated breakpoints of function \code{splitTraj}
}
\usage{
chooseBreaks(breakpoints, numbreaks)
}

\arguments{
  \item{breakpoints}{ Numeric list, output from splitTraj.

}
  \item{numbreaks}{ Integer. Number of breakpoints to be returned per trajectory. Breakpoints chosen will be evenly spaced from those defined by splitTraj.
}
}
\details{
None.
}
\value{ Returns a numeric list of breakpoints, one list per trajectory.
}
\author{
Diana H.P. Low, Efthimios Motakis
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{splitTraj}}
}
\examples{
data(deltaGseg)
all_breakpoints<-splitTraj(traj1) #default splits=15 (i.e. 16 segments).
all_breakpoints
chooseBreaks(all_breakpoints,numbreaks=3) 
}

