#' Threshold at the nearest kernel density minimum to the data centre
#'
#' @title KDE Local Minimum Threshold
#'
#' @description
#' Estimates the density of \code{marker.data} with \code{\link[stats]{density}}
#' and locates all interior local minima.  Among these, the minimum closest to
#' the centre of the data range is returned as the threshold.  Falls back to
#' the midpoint of the data range if no local minimum is found (e.g.
#' monotone or unimodal distribution).
#'
#' @param marker.data Numeric vector.  Values to threshold.
#' @param threshold.param Numeric scalar or \code{NULL}.  Bandwidth adjustment
#'   factor passed to \code{\link[stats]{density}} as \code{adjust}.  If
#'   \code{NULL}, defaults to \code{1}.
#'
#' @return Numeric scalar.  The x-coordinate of the selected local density
#'   minimum, or the midpoint of \code{range(marker.data)} if no minimum
#'   exists.
#'
#' @seealso \code{\link{calculate.threshold.tail}},
#'   \code{\link{calculate.threshold.bimodal}},
#'   \code{\link{calculate.gate.1dsep}}
#'
#' @keywords internal

calculate.threshold.kde.min <- function( marker.data, threshold.param ) {
  d    <- stats::density( marker.data, adjust = threshold.param %||% 1 )
  # Find local minima in the interior
  dy   <- diff( d$y )
  mins <- which( dy[-1] < 0 & dy[-length(dy)] >= 0 ) + 1L
  if ( length( mins ) == 0L )
    return( mean( range( marker.data ) ) )  # fallback
  # Pick the minimum closest to the centre of the data
  mid  <- mean( range( marker.data ) )
  d$x[ mins[ which.min( abs( d$x[ mins ] - mid ) ) ] ]
}
