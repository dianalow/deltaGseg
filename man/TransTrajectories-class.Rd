\name{TransTrajectories-class}
\Rdversion{1.1}
\docType{class}
\alias{TransTrajectories-class}
\title{Class \code{"TransTrajectories"}}
\description{
Object of this class is a product of the \code{transformSeries} function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TransTrajectories", ...)}.
}
\section{Slots}{
  \describe{
      \item{\code{tmethod}:}{Object of class \code{"character"} Transformation method. }
    \item{\code{breakpoints}:}{Object of class \code{"list"} breakpoints, if used. }
    \item{\code{tavd}:}{Object of class \code{"numeric"} adf p-values after transformation.  }
    \item{\code{ttrajlist}:}{Object of class \code{"list"} Transformed trajectories.  }
    \item{\code{difftraj}:}{Object of class \code{"list"} Differentiated trajectories. These trajectories may be produced if \code{\link{transformSeries}} was used with \code{method="differentiation"}. The original trajectory will be kept in \code{ttrajlist} its differentiated version (used only for computation, not presentation) will be stored in this slot. The plotting function \code{\link{plotDiff}} enables the user to compare the original and differentated versions of the subseries.}
    \item{\code{tfilenames}:}{Object of class \code{"character"} Transformed trajectories names. }
    \item{\code{path}:}{Object of class \code{"character"} Inherited from \link{Trajectories-class} }
    \item{\code{filenames}:}{Object of class \code{"character"} Inherited from \link{Trajectories-class} }
    \item{\code{trajlist}:}{Object of class \code{"list"} Trajectories. Inherited from \link{Trajectories-class} }
    \item{\code{avd}:}{Object of class \code{"numeric"} adf p-values. Inherited from \link{Trajectories-class} }

    
  }
}
\section{Extends}{
Class \code{"\linkS4class{Trajectories}"}, directly.
}
\section{Accessors}{
In the code snippets below, x is a TransTrajectories object.
\describe{
  \item{}{\code{getTNames(x)}: Retrieves filenames from slot \code{tfilenames}. }
  \item{}{\code{getTraj(x)}: Retrieves trajectories list from slot \code{ttrajlist}.}
  \item{}{\code{getBreaks(x)}: Retrives breakpoints (if any) from slot \code{breakpoints} }
  \item{}{\code{getAVD(x)}: Retrieves adf p-values from slot \code{tavd}. }

}
}
\section{Other methods}{
In the code snippets below, x is a TransTrajectories object.
  \describe{
    \item{}{\code{\link{denoiseSegments}(x)}: denoising and initial segmentation of trajectory series. }
    \item{}{\code{\link{plotDiff}(x,name='diff_object_name')}: Plots the original and differentiated subseries (one at a time) if \code{method="differentiation"} was used in \code{\link{transformSeries} }}
    \item{}{\code{plot(x)}: Plots the transformed series after \code{\link{transformSeries}}  }
    \item{}{\code{show(x)}: Displays summary of object, including inherited classes. This helps in recalling the analysis path taken to produce the current results. }
   }
}

\author{
Diana H.P. Low, Efthimios Motakis
}

\seealso{
\code{\link{transformSeries}}
}
\examples{
showClass("TransTrajectories")
}
\keyword{classes}
