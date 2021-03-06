\name{SegSeriesTrajectories-class}
\Rdversion{1.1}
\docType{class}
\alias{SegSeriesTrajectories-class}
\title{Class \code{"SegSeriesTrajectories"}}
\description{
Objects of this class is a product of the \code{clusterSegments} function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SegSeriesTrajectories", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{ssmatrix}:}{Object of class \code{"data.frame"} The output is a data.frame with the following information (in columns): the observed data ("observed"), the estimated, wavelet denoised data ("estimated"), the residuals of the estimation ("residuals"), the estimated subpopulations IDs ("subpopulations"), the series IDs/filenames ("seriesID") }
    \item{\code{ssparams}:}{Object of class \code{"character"} Parameters used to run the clustering algorithm, \code{\link{clusterSegments}}. }
    \item{\code{sparams}:}{Object of class \code{"character"} Parameters used to run the segmentation algorithm, \code{\link{denoiseSegments}}. }
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
    \item{\code{ct}:}{Object of class \code{"numeric"} Grouping information. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{SegTrajectories}"}, directly.
Class \code{"\linkS4class{TransTrajectories}"}, directly.
Class \code{"\linkS4class{Trajectories}"}, directly.
}

\section{Accessors}{
In the code snippets below, x is a SegSeriesTrajectories object.
\describe{
  \item{}{\code{getTNames(x)}: Retrieves filenames from slot \code{filenames} or \code{tfilenames} depending on whether the series has been transformed. }
  \item{}{\code{getTraj(x)}: Retrieves trajectories list from slot \code{trajlist} or \code{ttrajlist} depending on whether the series has been transformed. }
  \item{}{\code{getBreaks(x)}: Retrives breakpoints (if any) from slot \code{breakpoints} }
  \item{}{\code{getAVD(x)}: Retrieves adf p-values from slot \code{avd} or \code{tavd} depending on whether the series has been transformed. }
  \item{}{\code{getSegments(x)}: Retrieves clustered segmentation data from slot \code{smatrix} produced by \code{\link{clusterSegments}} }
}
}
\section{Other methods}{
In the code snippets below, x is a SegSeriesTrajectories object.
  \describe{
    \item{}{\code{\link{diagnosticPlots}(x)}: Generates the diagnostic plots of the wavelet denoising model residuals }
    \item{}{\code{\link{getIntervals}(x)}: Helper function to retrieve subpopulations and computes the intervals for each subpopulation after segmentation and clustering. }
    \item{}{\code{plot(x)}: Plots the final clustered segmentation data after \code{\link{clusterSegments}}  }
    \item{}{\code{show(x)}: Displays summary of object, including inherited classes. This helps in recalling the analysis path taken to produce the current results. }
   }
}
\author{
Diana H.P. Low, Efthimios Motakis
}

\seealso{
\code{\linkS4class{SegTrajectories}}, \code{\link{clusterSegments}} 
}
\examples{
showClass("SegSeriesTrajectories")
}
\keyword{classes}
