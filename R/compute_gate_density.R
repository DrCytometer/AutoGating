
# ==============================================================================
# compute.gate.density
#
# Trim outliers from gate.data for kde2d stability, compute a 2-D kernel
# density, interpolate back to the full data, and return a plotting data frame.
#
# @param gate.data  n×2 numeric matrix (parent-population events).
# @param p          Resolved parameter list (density.quantile.x/y,
#                   density.trimneg.x/y).
# @param ap         Autogating parameter object (density.grid.n).
# @return data.frame with columns x, y, z (z = interpolated density value).
# ==============================================================================

#' Compute 2-D kernel density for gate scatter plotting
#'
#' @title 2-D Kernel Density Estimation for Gate Plots
#'
#' @description
#' Trims outlying events to improve kernel density stability, estimates a 2-D
#' kernel density over the trimmed events using
#' \code{AutoSpectralRcpp::fast_kde2d_cpp}, then interpolates density values
#' back to all parent-population events via \code{\link[fields]{interp.surface}}.
#' Returns a data frame suitable for density-coloured scatter plotting with
#' \pkg{ggplot2}.
#'
#' @param gate.data Numeric matrix (n × 2).  Parent-population events
#'   (columns = x channel, y channel).
#' @param p Named list.  Resolved gate parameters; uses
#'   \code{density.quantile.x}, \code{density.quantile.y},
#'   \code{density.trimneg.x}, and \code{density.trimneg.y}.
#' @param ap Autogating parameter object; uses \code{density.grid.n} and
#'   \code{bandwidth.factor}.
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{\code{x}}{Numeric. x-axis values for all parent events.}
#'     \item{\code{y}}{Numeric. y-axis values for all parent events.}
#'     \item{\code{z}}{Numeric. Interpolated kernel density value for each
#'       event (used for colour mapping).}
#'   }
#'   Returns a zero-row data frame if \code{gate.data} contains no events.
#'
#' @importFrom MASS bandwidth.nrd kde2d
#' @importFrom fields interp.surface
#'
#' @seealso \code{\link{build.gate.plot}}, \code{\link{gate.plot}}
#'
#' @keywords internal
compute.gate.density <- function( gate.data, p, ap ) {

  if ( nrow( gate.data ) == 0L )
    return( data.frame( x = numeric(0), y = numeric(0), z = numeric(0) ) )

  # ---- trim for kde2d -------------------------------------------------------

  x <- gate.data[ , 1L ]
  y <- gate.data[ , 2L ]

  xmin <- min( x );  xmax <- max( x )
  ymin <- min( y );  ymax <- max( y )

  dx.lo <- stats::quantile( x, p$density.quantile.x )
  dx.hi <- stats::quantile( x, 1 - p$density.quantile.x )
  dy.lo <- stats::quantile( y, p$density.quantile.y )
  dy.hi <- stats::quantile( y, 1 - p$density.quantile.y )

  eps <- .Machine$double.eps^0.5

  if ( dx.lo <= xmin ) dx.lo <- ( 1 + eps ) * xmin
  if ( p$density.trimneg.x && dx.lo < 0 ) dx.lo <- 0
  if ( dx.hi >= xmax ) dx.hi <- ( 1 - eps ) * xmax

  if ( dy.lo <= ymin ) dy.lo <- ( 1 + eps ) * ymin
  if ( p$density.trimneg.y && dy.lo < 0 ) dy.lo <- 0
  if ( dy.hi >= ymax ) dy.hi <- ( 1 - eps ) * ymax

  keep <- x >= dx.lo & x <= dx.hi & y >= dy.lo & y <= dy.hi
  dens.data <- gate.data[ keep, , drop = FALSE ]

  # ---- interpolate density to full data ------------------------------------
# this should be changed to match AutoSpectral, I think
  if ( nrow( dens.data ) > 1L ) {
    bw <- apply( dens.data, 2, MASS::bandwidth.nrd )
    bw[ bw == 0 ] <- 0.1

    if ( requireNamespace("AutoSpectralRcpp", quietly = TRUE) &&
         "fast_kde2d_cpp" %in% ls(getNamespace("AutoSpectralRcpp")) ) {
      dens <- AutoSpectralRcpp::fast_kde2d_cpp(
        dens.data[, 1 ],
        dens.data[, 2 ],
        n = ap$density.grid.n,
        h = ap$bandwidth.factor * bw * 0.1,
        x_limits = range( dens.data[, 1 ] ),
        y_limits = range( dens.data[, 2 ] )
      )
    } else {
      dens <- MASS::kde2d(
        dens.data[, 1 ],
        dens.data[, 2 ],
        n = ap$density.grid.n,
        h = ap$bandwidth.factor * bw )
    }

    z <- fields::interp.surface(dens, gate.data)

  } else {
    z <- rep( 1 / nrow( gate.data ), nrow( gate.data ) )
  }

  return( data.frame( x = x, y = y, z = z ) )
}







