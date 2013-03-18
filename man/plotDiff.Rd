\name{plotDiff}
\alias{plotDiff}
\title{
plotDiff
}
\description{
Plots (sub)series before and after transformation.
}
\usage{
plotDiff(object,name=NULL)
}

\arguments{
  \item{object}{ An object of class 'TransTrajectories'.
}
  \item{name}{ Character. Name of (sub)series.
}
}
\details{
None.
}
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
data(deltaGseg)
plotDiff(simtraj.tr2,name="1_2")
}
