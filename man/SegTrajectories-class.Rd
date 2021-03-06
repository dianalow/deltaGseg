\name{SegTrajectories-class}
\Rdversion{1.1}
\docType{class}
\alias{SegTrajectories-class}

\title{Class \code{"SegTrajectories"}}
\description{
Objects of this class is a product of the \code{denoiseSegments} function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SegTrajectories", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{sparams}:}{Object of class \code{"character"} Parameters used to run the segmentation algorithm. }
    \item{\code{smatrix}:}{Object of class \code{"matrix"} A matrix containing the preliminary results from segmentation/denoising for each (sub)series generated by \code{\link{transformSeries}}. Each list element contains the following information in matrix form (in columns): the observed data (1st column), the estimated, wavelet denoised data (2nd column), the residuals of the estimation (3rd column), the starting/ending time points of each segment (4th/5th columns), the estimated segment IDs (6th column), the quantiles of the estimated data [minimum, 5\%, 10\%, 15\%, ..., 95\%, maximum] (from columns 7th to 27th) and the series IDs (28th column)  }
    \item{\code{path}:}{Object of class \code{"character"} Inherited from \link{Trajectories-class} }
    \item{\code{filenames}:}{Object of class \code{"character"} Inherited from \link{Trajectories-class} }
    \item{\code{trajlist}:}{Object of class \code{"list"} Trajectories. Inherited from \link{Trajectories-class} }
    \item{\code{avd}:}{Object of class \code{"numeric"} adf p-values. Inherited from \link{Trajectories-class} }
    \item{\code{tmethod}:}{Object of class \code{"character"} Transformation method. Inherited from \link{TransTrajectories-class}  }
    \item{\code{breakpoints}:}{Object of class \code{"list"} breakpoints, if used. }
    \item{\code{tavd}:}{Object of class \code{"numeric"} adf p-values after transformation, if used. Inherited from \link{TransTrajectories-class}  }
    \item{\code{ttrajlist}:}{Object of class \code{"list"} Transformed trajectories. Inherited from \link{TransTrajectories-class} }
    \item{\code{tfilenames}:}{Object of class \code{"character"} Transformed trajectories names. Inherited from \link{TransTrajectories-class} }
  }
}
\section{Extends}{
Class \code{"\linkS4class{TransTrajectories}"}, directly.
Class \code{"\linkS4class{Trajectories}"}, directly.
}

\section{Accessors}{
In the code snippets below, x is a SegTrajectories object.
\describe{
  \item{}{\code{getTNames(x)}: Retrieves filenames from slot \code{filenames} or \code{tfilenames} depending on whether the series has been transformed. }
  \item{}{\code{getTraj(x)}: Retrieves trajectories list from slot \code{trajlist} or \code{ttrajlist} depending on whether the series has been transformed. }
  \item{}{\code{getBreaks(x)}: Retrives breakpoints (if any) from slot \code{breakpoints} }
  \item{}{\code{getAVD(x)}: Retrieves adf p-values from slot \code{avd} or \code{tavd} depending on whether the series has been transformed. }
  \item{}{\code{getSegments(x)}: Retrieves initial segmentation data from slot \code{smatrix} produced by \code{\link{denoiseSegments}} }
}
}
\section{Other methods}{
In the code snippets below, x is a SegTrajectories object.
  \describe{
    \item{}{\code{\link{clusterPV}(x,bootstrap=500)}: Computes p-values to be used with \code{method="pvclust"} in \code{\link{clusterSegments}(x)} }
    \item{}{\code{\link{clusterSegments}(x)}: clustering of segmented trajectories into similar groups. }
    \item{}{\code{plot(x)}: Plots the initial segmentation data after denoising by \code{\link{denoiseSegments}}  }
    \item{}{\code{show(x)}: Displays summary of object, including inherited classes. This helps in recalling the analysis path taken to produce the current results. }
	 }
}

\author{
Diana H.P. Low, Efthimios Motakis
}

\seealso{
\code{\link{denoiseSegments}}
}
\examples{
showClass("SegTrajectories")
}
\keyword{classes}
