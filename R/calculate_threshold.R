# calculate_threshold.r
#
# Threshold calculation functions used by calculate.gate.1dsep and
# calculate.gate.2dsep.
#
# Each function receives:
#   marker.data      Numeric vector of values along the separation axis.
#   threshold.param  Per-gate override value (character or numeric), or NULL
#                    to use the cytometer default from p$default.tail.probability.
#                    Passed through from the resolved gate parameter list.


#' Threshold by bimodal valley detection
#'
#' @title Bimodal Valley Threshold (Deprecated)
#'
#' @description
#' Finds the valley between two modes in a univariate distribution by
#' converting the data to a one-column \code{\link[flowCore]{flowFrame}} and
#' calling \code{\link[flowStats]{rangeGate}} with \code{alpha = "min"}.  The
#' result is clamped to the observed data range. Deprecated approach.
#'
#' @param marker.data Numeric vector.  Transformed fluorescence or scatter
#'   values along the separation axis.
#' @param threshold.param Not used by this method; present for interface
#'   consistency with other threshold functions.  Pass \code{NULL} or any
#'   value.
#'
#' @return Numeric scalar.  The estimated threshold value, clamped to
#'   \code{[min(marker.data), max(marker.data)]}.
#'
#' @importFrom flowCore flowFrame
#' @importFrom flowStats rangeGate
#'
#' @seealso \code{\link{calculate.threshold.tail}},
#'   \code{\link{calculate.gate.1dsep}}
#'
#' @keywords internal
calculate.threshold.bimodal <- function( marker.data, threshold.param ) {

  flow.data <- as.matrix( marker.data )
  colnames( flow.data ) <- "marker"

  threshold.calc <- rangeGate( flowFrame( flow.data ), "marker",
                               borderQuant = 0, alpha = "min" )

  x.thr <- as.numeric( threshold.calc@min )
  # Clamp to the observed data range
  x.rng <- range( marker.data, na.rm = TRUE )
  x.thr <- max( x.rng[1], min( x.rng[2], x.thr ) )
  x.thr
}

#' Threshold by symmetric tail probability
#'
#' @title Symmetric Tail Threshold
#'
#' @description
#' Estimates the mode of the distribution using the shorth estimator
#' (\code{\link[modeest]{mlv}}, falling back to the median on failure), then
#' constructs a symmetrised distribution around the mode by reflecting the
#' left half of the data.  The threshold is set at the
#' \code{threshold.param} quantile of this symmetric distribution, providing
#' an estimate of the upper tail of the negative population.  A second
#' symmetrisation from the right is computed and the minimum of the two
#' quantiles is returned.
#'
#' @param marker.data Numeric vector.  Must contain at least 10 finite values.
#' @param threshold.param Numeric scalar in \eqn{[0, 1]}.  Quantile of the
#'   symmetrised distribution used as the threshold (e.g. \code{0.99}).
#'
#' @return Numeric scalar.  The estimated threshold.
#'
#' @importFrom modeest mlv
#'
#' @seealso \code{\link{calculate.threshold.bimodal}},
#'   \code{\link{calculate.threshold.mixture}},
#'   \code{\link{calculate.gate.1dsep}}
#'
#' @keywords internal
calculate.threshold.tail <- function( marker.data, threshold.param ) {
  x <- marker.data

  if ( length( x ) == 0L || all( is.na( x ) ) )
    stop( "calculate.threshold.tail: marker.data is empty or all NA" )

  x <- x[ is.finite( x ) ]
  if ( length( x ) < 10L )
    stop( sprintf( "calculate.threshold.tail: too few finite values (%d)", length( x ) ) )

  x.mode <- tryCatch(
    as.numeric( mlv( x, method = "shorth" ) ),
    error = function(e) {
      warning( sprintf( "calculate.threshold.tail: mlv() failed (%s); falling back to median", conditionMessage(e) ) )
      stats::median( x )
    }
  )

  if ( length( x.mode ) != 1L || !is.finite( x.mode ) ) {
    warning( "calculate.threshold.tail: mode is non-scalar or non-finite; falling back to median" )
    x.mode <- stats::median( x )
  }

  x.symm.low  <- x[ x < x.mode ]
  x.symm.high <- 2 * x.mode - x.symm.low
  x.symm      <- c( x.symm.low, x.symm.high )

  x.symm.rev.high <- x[ x > x.mode ]
  x.symm.rev.low  <- 2 * x.mode - x.symm.rev.high
  x.symm.rev      <- c( x.symm.rev.low, x.symm.rev.high )

  if ( length( x.symm ) == 0L || length( x.symm.rev ) == 0L )
    stop( sprintf(
      "calculate.threshold.tail: symmetrised data empty (mode=%.3f at %s of range [%.3f, %.3f])",
      x.mode,
      if ( x.mode <= min(x) ) "bottom" else if ( x.mode >= max(x) ) "top" else "interior",
      min(x), max(x)
    ))

  if ( !is.numeric( threshold.param ) || length( threshold.param ) != 1L ||
       threshold.param < 0 || threshold.param > 1 )
    stop( sprintf( "calculate.threshold.tail: threshold.param must be in [0,1], got: %s",
                   deparse( threshold.param ) ) )

  x.thr     <- as.numeric( stats::quantile( x.symm,     threshold.param ) )
  x.thr.rev <- as.numeric( stats::quantile( x.symm.rev, threshold.param ) )

  result <- min( x.thr, x.thr.rev )

  if ( !is.finite( result ) )
    stop( sprintf(
      "calculate.threshold.tail: non-finite threshold (x.thr=%.3f, x.thr.rev=%.3f)",
      x.thr, x.thr.rev
    ))

  result
}
