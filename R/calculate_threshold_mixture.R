#' Threshold by Gaussian mixture model crossover
#'
#' @title Gaussian Mixture Model Threshold
#'
#' @description
#' Fits a two-component Gaussian mixture model to \code{marker.data} using
#' \code{\link[mclust]{Mclust}} and finds the analytic crossover point between
#' the two fitted Gaussians via \code{\link[stats]{uniroot}}.  If the model
#' fails to fit or \code{uniroot} fails to converge, falls back to the
#' midpoint of the data range or the mean of component means, respectively.
#'
#' @param marker.data Numeric vector.  Transformed fluorescence values.
#' @param threshold.param Not used by this method; present for interface
#'   consistency with other threshold functions.
#'
#' @return Numeric scalar.  The estimated threshold at the crossover of the
#'   two fitted Gaussian densities.
#'
#' @importFrom mclust Mclust
#'
#' @seealso \code{\link{calculate.threshold.tail}},
#'   \code{\link{calculate.threshold.bimodal}},
#'   \code{\link{calculate.gate.1dsep}}
#'
#' @keywords internal

calculate.threshold.mixture <- function( marker.data, threshold.param ) {
  fit <- mclust::Mclust( marker.data, G = 2, verbose = FALSE )
  if ( is.null( fit ) )
    return( mean( range( marker.data ) ) )
  mu  <- fit$parameters$mean
  sig <- sqrt( fit$parameters$variance$sigmasq )
  pro <- fit$parameters$pro
  # Analytic crossover between two Gaussians
  # Solve: pro[1]*N(x,mu[1],sig[1]) = pro[2]*N(x,mu[2],sig[2])
  # For equal variances this simplifies to the midpoint weighted by means
  # For unequal, solve numerically on the interval between means
  f <- function(x)
    pro[1] * stats::dnorm( x, mu[1], sig[1] ) - pro[2] * stats::dnorm( x, mu[2], sig[2] )
  tryCatch(
    stats::uniroot( f, interval = sort( mu ) )$root,
    error = function(e) mean( mu )
  )
}
