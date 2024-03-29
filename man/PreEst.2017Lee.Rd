% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PreEst.2017Lee.R
\name{PreEst.2017Lee}
\alias{PreEst.2017Lee}
\title{Bayesian Estimation of a Banded Precision Matrix (Lee 2017)}
\usage{
PreEst.2017Lee(X, upperK = floor(ncol(X)/2), logpi = function(k) {
    -k^4
})
}
\arguments{
\item{X}{an \eqn{(n\times p)} data matrix where each row is an observation.}

\item{upperK}{upper bound of bandwidth \eqn{k}.}

\item{logpi}{log of prior distribution for bandwidth \eqn{k}. Default is a function proportional to \eqn{-k^4}.}
}
\value{
a named list containing: \describe{
\item{C}{a \eqn{(p\times p)} MAP estimate for precision matrix.}
}
}
\description{
\code{PreEst.2017Lee} returns a Bayes estimator of the banded precision matrix,
which is defined in subsection 3.3 of Lee and Lee (2017), using the k-BC prior.
The bandwidth is set at the mode of marginal posterior for the bandwidth parameter.
}
\examples{
## generate data from multivariate normal with Identity precision.
pdim = 5
data = matrix(rnorm(100*pdim), ncol=pdim)

## compare different K
out1 <- PreEst.2017Lee(data, upperK=1)
out2 <- PreEst.2017Lee(data, upperK=3)
out3 <- PreEst.2017Lee(data, upperK=5)

## visualize
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,2), pty="s")
image(diag(pdim)[,pdim:1], main="Original Precision")
image(out1$C[,pdim:1],     main="banded2::upperK=1")
image(out2$C[,pdim:1],     main="banded2::upperK=3")
image(out3$C[,pdim:1],     main="banded2::upperK=5")
par(opar)

}
\references{
\insertRef{lee_estimating_2017}{CovTools}
}
