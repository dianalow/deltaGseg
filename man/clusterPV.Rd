\name{clusterPV}
\alias{clusterPV}
\title{
clusterPV
}
\description{
Wrapper for modified pvclust function.
}
\usage{
clusterPV(object,bootstrap=500)
}

\arguments{
  \item{object}{ An object of class 'SegTrajectories'.
}
  \item{bootstrap}{ Integer. Number of bootstraps to run.
}
}
\details{
This is a wrapper to call the pvclust function that has been modified to suit our deltaGseg computation.
}
\value{ Returns an object of class "pvclust". For use in \code{\link{clusterSegments}} when running the pvclust option.
}
\references{
Shimodaira, H. (2004). Approximately unbiased tests of regions using multistep-multiscale bootstrap resampling, Annals of Statistics 32, 2616-2641.
}
\author{
Diana H.P. Low, Efthimios Motakis
}
\seealso{
\code{\link{clusterSegments}}]
}
\examples{
data(deltaGseg)
clusterPV(traj1.denoise)
}
