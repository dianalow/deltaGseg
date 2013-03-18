\name{diagnosticPlots}
\alias{diagnosticPlots}
\title{
diagnosticPlots
}
\description{
This function generates the diagnostic plots of the wavelet denoising model residuals. The assumptions are that the residuals autocorrelation is not significant and that the residuals distribution is approximately normal or, at least, symmetric around 0. We provide plots and test to verify these assumptions (depends on R package fBasics).

}
\usage{
diagnosticPlots(object, norm.test="KS",single.series = FALSE)
}

\arguments{
  \item{object}{An object of class "SegSeriesTrajectories".
}
  \item{norm.test}{Character. One of "KS", "Shapiro", "Agost". Test for residuals normality accepting "KS" (Kolmogorov Smirnov test with Lilliefors correction), "Shapiro"" (Shapiro test for normality) and "Agost"" (D'Agostino test for normality using the skewness and kurtosis of the data ; also gives the skewness and kurtosis p-values for the hypothesis that the respective estimated measures differ from the theoretical values under the normal distribution).
  }
  \item{single.series}{ Logical. If FALSE (default) the residuals of each series are independently analyzed.    
}
}
\details{
The function outputs the standard autocorrelation plots for viewing the residuals autocorrelation, histograms for checking the normality assumptions and the respective P-values to test the normality assumption.
}
\value{A series of plots with printed P-values for the autocorrelation and normality tests.
}
\references{
D'Agostino R.B., Pearson E.S. (1973); Tests for Departure from Normality, Biometrika 60, 613-22.
}
\author{
Diana H.P. Low, Efthimios Motakis
}

\examples{
data(deltaGseg)
diagnosticPlots(traj1.ss,norm.test="KS",single.series=TRUE)
}
