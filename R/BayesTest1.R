#' Bayesian One-Sample Tests for Covariance Matrix
#'
#' Given data, \code{BayesTest1} performs Bayesian version of 1-sample test for Covariance where
#' the null hypothesis is
#' \deqn{H_0 : \Sigma_n = \Sigma_0}
#' where \eqn{\Sigma_n} is the covariance of data model and \eqn{\Sigma_0} is a
#' hypothesized covariance.
#'
#' @param data an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
#' @param method a name of test.
#' @param ... extra arguments to be passed along for each procedure. See below for details.
#' \tabular{lll}{
#' \emph{parameter} \tab \emph{method} \tab \emph{description} \cr
#' \code{a0} \tab \code{"mxPBF"} \tab hyperparameter (see below for details) \cr
#' \code{b0} \tab \code{"mxPBF"} \tab hyperparameter (see below for details) \cr
#' \code{gamma} \tab \code{"mxPBF"} \tab hyperparameter (see below for details)
#' }
#'
#' @return a named list containing one of followings,
#' \tabular{lll}{
#' \emph{element} \tab \emph{method} \tab \emph{description} \cr
#' \code{log.BF.mat} \tab \code{"mxPBF"} \tab a \eqn{(p\times p)} matrix of pairwise log Bayes factors.
#' }
#'
#' @section mxPBF:
#' Let \eqn{X_i} be the \eqn{i}-th column of data matrix. Under the maximum pairwise bayes factor framework, we have following hypothesis,
#' \deqn{H_0: a_{ij}=0~\mathrm{ and }~\tau_{ij}=1 \quad \mathrm{versus. } \quad  H_1: \mathrm{ not }~ H_0.}
#' The model is
#' \deqn{X_i | X_j \sim N_n( a_{ij}X_j, \tau_{ij}^2 I_n )}
#' and the prior is set, under \eqn{H_1},  as
#' \deqn{ a_{ij}|\tau_{ij}^2 \sim N(0, \tau_{ij}^2/(\gamma*||X_j||^2))}
#' \deqn{\tau_{ij}^2 \sim IG(a0, b0).}
#'
#' @examples
#' \dontrun{
#' ## generate data from multivariate normal with trivial covariance.
#' data = matrix(rnorm(100*5), nrow=100)
#'
#' ## run test
#' BayesTest1(data)                 # default setting
#' BayesTest1(data, a0=5.0, b0=5.0) # change hyperparameters
#' }
#'
#' @references
#' \insertRef{lee_maximum_2018}{CovTools}
#'
#' @export
BayesTest1 <- function(data, Sigma0=diag(ncol(data)), method=c("mxPBF"), ...){
  ###########################################################################
  # Preprocessing : Inputs
  # 1. data
  if (!check_datamatrix(data)){
    stop("* BayesTest1 : an input matrix is invalid.")
  }
  n = nrow(data)
  p = ncol(data)
  if ((nrow(data)==1)||(ncol(data)==1)){stop("* BayesTest1 : invalid input matrix X.")}
  # 2. valid Sigma0
  if ((!check_sqmat(Sigma0))||(!check_pd(Sigma0))||(!isSymmetric(Sigma0,tol=sqrt(.Machine$double.eps)))||(nrow(Sigma0)!=p)){
    stop("* BayesTest1 : a given matrix for null hypothess 'Sigma0' is invalid.")
  }
  # 3. method
  method = match.arg(method)
  # 4. extra arguments
  extra.args = (list(...))

  ###########################################################################
  # Preprocessing : Adjust the data for testing I
  scaler = get_invroot(Sigma0)
  X.centered = scale(data, center=TRUE, scale=FALSE)
  X.adjusted = (matrix(X.centered,nrow=n) %*% scaler)

  ###########################################################################
  # Main Computation
  output = switch(method,
                  mxPBF = bayestest1.Lee18(X.adjusted, extra.args))

  ###########################################################################
  return(output)
}

# pack <- "CovTools"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))
