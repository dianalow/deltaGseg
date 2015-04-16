\name{parseTraj}
\alias{parseTraj}
\title{
Reads in files containing trajectory data
}
\description{
Reads in files with 2-column, space-separated numerical values containing 1:time points, 2:trajectory(free binding energies).
}
\usage{
parseTraj(path = getwd(), files = NULL,fromfile=TRUE)
}

\arguments{
  \item{path}{Directory containing trajectory files.
}
  \item{files}{Character vector of filenames to read. If not provided, will read all files in given directory and treat them as a set. Can also be used to read in variables if given as list.
}
\item{fromfile}{Logical. If set to FALSE, the files parameter will be used to read in variables.}
}

\details{ 
This is an initialization function for the deltaGseg package. It reads the trajectory files (input) and reports the a short description of the file, the Augmented Dickey-Fuller test p-values for each trajectory in the set and the data plot. The input files should be in tab delimited form with 2 columns: the first column contains the time points 1, 2, ..., T and the second the free binding energies at each time point.
}

\value{ A 'Trajectories' object.
}
\references{
None.
}
\author{
Diana H.P. Low, Efthimios Motakis
}

\seealso{
\code{\link{Trajectories-class}}
}
\examples{
traj1<-parseTraj(path=system.file("extdata",package="deltaGseg"),files=c("D_GBTOT1","D_GBTOT2","D_GBTOT3"))
traj1  #prints summary of traj1 object

# using parseTraj for existing variables ##
subtraj<-getTraj(traj1)[[1]] #extracts first trajectory in the above series
traj2<-parseTraj(files=list(subtraj),fromfile=FALSE)
traj2
}