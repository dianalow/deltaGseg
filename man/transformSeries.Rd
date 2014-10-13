\name{transformSeries}
\alias{transformSeries}
\title{
Transforms a Non-Stationary Series into a Weakly-Stationary one
}
\description{
It transforms non-stationary series into weakly stationary (sub)series with three alternative methods (see parameter "methods").
}
\usage{
transformSeries(object, method = "splitting", breakpoints = 1)
}

\arguments{
  \item{object}{An object of class "Trajectories".
}
  \item{method}{Character. One of "differentiation","override_splitting","splitting". See details.
}
  \item{breakpoints}{Integer (for method="splitting") or numeric vector (for method="override_splitting"). See details.
}
}
\details{
(i) method="differentiation": first differences B[t]-B[t-1] are calculated and the first differentiated series is used for further analysis (segmentation and clustering). This option is needed in special cases, only when the series exhibits a trend-like behavior that cannot be removed by splitting. The first differentiations with remove the trend completely (see Appendix in the manual).

(ii) method="splitting": the series are split by automatic data segmentation to subseries. This option divides a non-stationary series to a number of subseries that are weakly stationary;

(iii) method="override_splitting": the series are split into subseries by user-defined cut-offs obtained from the numerical output of the splitTraj function. This option is for long stationary series that cannot be analyzed due to memory limitations. It can be also used for manual splitting when "splitting" option is not satisfactory. Typically, "splitting"" and "override_splitting" generate new data files of subseries.


Determining breakpoints for method
(i) "splitting": an integer specifying the number of split points. This number (a single value applied to all series) denotes the number of subseries that the original series should be divided into.
(ii) "override_splitting": the parameter takes the exact values (time points coordinates) of split points (a list of length equal to the number of series; see manual). One can derive and manually insert these split points after inspecting the output of the splitTraj function (see manual). The user should select a few splits so that the original series is not divided into too many subseries (difficult to process because many new files are generated). Alternatively, function chooseBreaks automatically chooses a subset of breakpoints (not recommended to keep those without inspection)


}
\value{An object of class "TransTrajectories".
}
\references{
Dickey, D.A. and W.A. Fuller (1979). Distribution of the Estimators for Autoregressive Time Series with a Unit Root. Journal of the American Statistical Association 74, 427-431.
}
\author{
Diana H.P. Low, Efthimios Motakis
}

\seealso{
\code{\link{parseTraj}},\code{\link{splitTraj}},\code{\link{chooseBreaks}}
}
\examples{
data(deltaGseg)
trans_series<-transformSeries(traj1,method='splitting',breakpoints=1)
}
