\name{denoiseSegments}
\alias{denoiseSegments}
\title{
Wavelet denoising of trajectory series (segments)
}
\description{
This function computes Augmented Dickey-Fuller test for weak-stationarity and carries out segmentation and wavelets denoising.
}
\usage{
denoiseSegments(object,seg_method="BinSeg",maxQ=15,fn=1,factor=0.8,thresh_level=TRUE,minobs=200)
}
\arguments{
  \item{object}{ An object of class "Trajectories" or "TransTrajectories". (the output of functions \code{\link{parseTraj}} or \code{\link{transformSeries}}, respectively).}
  \item{seg_method}{Character. One of "SegNeigh" or "BinSeg". By default it performs the Segment Neighborhood (SegNeigh) method to find multiple changes in mean for data that is assumed to be normally distributed. The value returned is the result of finding the optimal location of up to Q changepoints using the log of the likelihood ratio statistic. Once all changepoint locations have been calculated, the optimal number of changepoints is decided using k*pen as the penalty function where k is the number of changepoints tested (k in (1,Q)). In very large series, the memory demanding "SegNeigh" can be replaced by "BinSeg"" (Binary Segmentation). The segmenation is performed only if the Augmented Dickey-Fuller test P-value is significant at alpha=0.05. If not, an error message appears indicating the need of series splitting or differentiation (see \code{\link{transformSeries}}).
}
  \item{maxQ}{Integer. The maximum number of Q changepoints to be estimated.
}
  \item{fn}{Integer. It specifies the degree of smoothness of the wavelet that you want to use in the decomposition. It takes values from 1 (coarsest/hard smoothing, i.e. Haar's step function) to 10 (finest/soft smoothing). Put simply, the fitted (wavelet denoised/estimated) data of segment q with fn=1 have lower variance than the fitted data of q under fn=10. The former data will resemble a step function over time while the latter will be much closer to the original data.
}
  \item{factor}{Numeric. Between 0.6 and 1 used for re-scaling the denoising threshold.
}
  \item{thresh_level}{Logical. If FALSE then a global threshold is computed on and applied to all scale levels. If TRUE a threshold is computed and applied separately to each scale level (for serious residuals autocorrelation).
  }
  \item{minobs}{Integer. The minimum number of observations a segment should consist of to be accepted
}
}

\value{An object of class "SegTrajectories"

}
\references{
Dickey, D.A. and W.A. Fuller (1979), Distribution of the Estimators for Autoregressive Time Series with a Unit Root, Journal of the American Statistical Association, 74, p. 427-431.

Nason, G.P. (2008) Wavelet methods in Statistics with R. Springer, New York.

Auger, I. E.; Lawrence, C. E. Algorithms for the optimal identification of segment neighborhoods. Bull. Math.Biol. 1989, 51(1), 39-54.

}

\author{
Diana H.P. Low, Efthimios Motakis
}

\seealso{
\code{\link[deltaGseg]{transformSeries}}
}
\examples{
data(deltaGseg)
traj1.denoise<-denoiseSegments(traj1.tr,seg_method="BinSeg",maxQ=15,fn=1,factor=0.8,thresh_level=TRUE,minobs=200)
}
