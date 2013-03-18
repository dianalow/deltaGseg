\name{Trajectories-class}
\Rdversion{1.1}
\docType{class}
\alias{Trajectories-class}
\title{Class \code{"Trajectories"}}
\description{
Objects of this class is a product of the initialization function, \code{parseTraj}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Trajectories", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{path}:}{Object of class \code{"character"} Directory where files were read from. }
    \item{\code{filenames}:}{Object of class \code{"character"} Name of files read.}
    \item{\code{trajlist}:}{Object of class \code{"list"} Trajectories. }
    \item{\code{avd}:}{Object of class \code{"numeric"} adf p-values.
  }
}
}
\section{Accessors}{
In the code snippets below, x is a Trajectories object.
\describe{
  \item{}{\code{getTNames(x)}: Retrieves filenames from slot \code{filenames}. }
  \item{}{\code{getTraj(x)}: Retrieves trajectories list from slot \code{trajlist}.}
  \item{}{\code{getBreaks(x)}: Retrives breakpoints (if any) from slot \code{breakpoints} }
  \item{}{\code{getAVD(x)}: Retrieves adf p-values from slot \code{avd}. }

}
}
\section{Other methods}{
In the code snippets below, x is a Trajectories object.
  \describe{
    \item{}{\code{\link{splitTraj}(x)}: Computes likely breakpoints for the series.}
    \item{}{\code{\link{transformSeries}(x)}: Apply transformation functions for the series if series is not stationary, or to split long series after determining breakpoints with \code{\link{splitTraj}}. }
    \item{}{\code{plot(x,name='all')}: Plots the trajectory series either individually, or combined.  }
    \item{}{\code{show(x)}: Displays summary of object, including inherited classes. This helps in recalling the analysis path taken to produce the current results. }
   }
}

\author{
Diana H.P. Low, Efthimios Motakis
}

\seealso{
\code{\link{parseTraj}}
}
\examples{
showClass("Trajectories")
}
\keyword{classes}
