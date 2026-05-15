# gate_singlet.r
#
# Fit a robust linear model to the FSC-A / FSC-H relationship and return a
# polygonGate that captures the singlet population.
#

#' Fit a singlet-discrimination polygon gate using robust linear regression
#'
#' @title Singlet Gate via Robust Linear Model
#'
#' @description
#' Identifies singlet events by modelling the linear relationship between
#' scatter area and scatter height using \code{\link[MASS]{rlm}}.  A grid of
#' representative points is sampled along the area axis, local spread
#' (residual width) is estimated at each grid point, and the spread envelope
#' is modulated by a cubic spline to produce a tapered polygon gate.  The
#' result is a \code{\link[flowCore]{polygonGate}} object ready for use with
#' the \pkg{flowCore} framework.
#'
#' @param x Numeric matrix.  Must contain at least columns named by
#'   \code{area} and \code{height} (and \code{sidescatter} if provided).
#' @param area Character scalar.  Name of the scatter-area channel column.
#'   Default \code{"FSC-A"}.
#' @param height Character scalar.  Name of the scatter-height channel column.
#'   Default \code{"FSC-H"}.
#' @param sidescatter Character scalar or \code{NULL}.  Optional side-scatter
#'   channel to include as an additional covariate in the regression.
#'   Default \code{NULL}.
#' @param filterId Character scalar.  Identifier for the resulting
#'   \code{polygonGate}.  Default \code{"singlet"}.
#' @param maxit Integer.  Maximum IRLS iterations passed to
#'   \code{\link[MASS]{rlm}}.  Default \code{100}.
#' @param gate.grid.n Integer.  Number of evenly spaced grid points along the
#'   area axis used to estimate the spread envelope.  Default \code{20}.
#' @param gate.skip.low Numeric in \eqn{[0, 1]}.  Fraction of the area range
#'   to exclude at the low end before fitting (removes near-zero debris).
#'   Default \code{0.015}.
#' @param gate.skip.high Numeric in \eqn{[0, 1]}.  Fraction of the area range
#'   to exclude at the high end before fitting (removes saturated events).
#'   Default \code{0.15}.
#' @param spread.span Numeric.  Half-width (in area units, as a multiple of
#'   the grid spacing) of the local window used to estimate residual spread at
#'   each grid point.  Default \code{0.5}.
#' @param spread.factor Numeric.  Multiplier applied to the local spread to
#'   set gate width.  Default \code{10}.
#' @param spread.spline.x Numeric vector.  x-knots (in \eqn{[0, 1]},
#'   normalised grid position) for the spline envelope.
#' @param spread.spline.y Numeric vector.  y-values of the spline envelope
#'   at each knot (spread scale factors).
#' @param ... Additional arguments passed to \code{\link[MASS]{rlm}}.
#'
#' @return A \code{\link[flowCore]{polygonGate}} object whose
#'   \code{boundaries} slot contains the singlet polygon in the
#'   area-vs-height plane.
#'
#' @importFrom MASS rlm
#' @importFrom flowCore polygonGate
#'
#' @seealso \code{\link{calculate.boundary.singlets}},
#'   \code{\link[flowCore]{polygonGate}}
#'
#' @export
gate.singlet <- function( x,
    area        = "FSC-A",
    height      = "FSC-H",
    sidescatter = NULL,
    filterId    = "singlet",
    maxit           = 100,
    gate.grid.n     = 20,
    gate.skip.low   = 0.01,
    gate.skip.high  = 0.15,
    spread.span     = 0.5,
    spread.factor   = 10,
    spread.spline.x = c( 0, 0.2, 0.4, 1 ),
    spread.spline.y = c( 1, 0.5, 1, 0 ),
    ... )
{
  stopifnot( is.matrix( x ) )
  stopifnot( is.character( area )   && length( area )   == 1L )
  stopifnot( is.character( height ) && length( height ) == 1L )
  if ( !is.null( sidescatter ) )
    stopifnot( is.character( sidescatter ) )

    # ---- filter out extrema --------------------------------------------------

    x <- data.frame( x[
        x[ , area ] > gate.skip.low  * max( x[ , area ] ) +
            ( 1 - gate.skip.low  ) * min( x[ , area ] ) &
        x[ , area ] < gate.skip.high * min( x[ , area ] ) +
            ( 1 - gate.skip.high ) * max( x[ , area ] ) &
        x[ , height ] > gate.skip.low  * max( x[ , height ] ) +
            ( 1 - gate.skip.low  ) * min( x[ , height ] ) &
        x[ , height ] < gate.skip.high * min( x[ , height ] ) +
            ( 1 - gate.skip.high ) * max( x[ , height ] ),
        c( area, height, sidescatter ) ] )

    channel_names <- c( area, height )
    area        <- make.names( area )
    height      <- make.names( height )

    rlm_formula <- paste( height, area, sep = " ~ " )
    if ( !is.null( sidescatter ) ) {
        sidescatter <- make.names( sidescatter )
        ssc_ratio   <- paste0( "I(", sidescatter, " / ", area, ")" )
        rlm_formula <- paste( rlm_formula, sidescatter, ssc_ratio, sep = " + " )
    }
    rlm_formula <- stats::as.formula( rlm_formula )

    rlm_fit <- MASS::rlm( rlm_formula, data = x, maxit = maxit, ... )
    if ( !rlm_fit$converged )
        warning( "The IRLS algorithm in rlm() did not converge." )

    # ---- build grid of representative points ---------------------------------

    x.delta     <- ( max( x[[ area ]] ) - min( x[[ area ]] ) ) / ( gate.grid.n - 1 )
    x.grid      <- seq( min( x[[ area ]] ), max( x[[ area ]] ), by = x.delta )
    x.point.idx <- sapply( x.grid, function( xg ) which.min( abs( x[[ area ]] - xg ) ) )
    x.predict   <- x[ x.point.idx, ]
    predictions <- stats::predict( rlm_fit, x.predict, interval = "none" )

    # ---- estimate spread at each grid point ----------------------------------

    x.spread <- sapply( seq_len( gate.grid.n ), function( idx ) {
        x.range.idx <- which( abs( x[[ area ]] - x.grid[ idx ] ) < x.delta * spread.span )
        x.pred      <- x[ x.range.idx, ]
        pred        <- stats::predict( rlm_fit, x.pred, interval = "none" )

        pos.bol  <- x[[ height ]][ x.range.idx ] > pred
        pos.idx  <- x.range.idx[ pos.bol ]
        pred.pos <- pred[ pos.bol ]

        stats::quantile( x[[ height ]][ pos.idx ] - pred.pos, 0.5 )
    } )

    # Any NAs (from empty grid cells) replaced with minimum observed spread
    x.spread[ is.na( x.spread ) ] <- min( x.spread, na.rm = TRUE )

    # Scale spread by the spline envelope
    x.spread <- x.spread *
      stats::splinefun( spread.spline.x, spread.spline.y, method = "natural" )(
            ( seq_len( gate.grid.n ) - 1L ) / ( gate.grid.n - 1L ) )

    # ---- construct polygon boundary ------------------------------------------

    gate_vertices <- rbind(
        t( sapply( seq_len( gate.grid.n ), function( idx )
            c( x.predict[[ area ]][ idx ],
               predictions[ idx ] - spread.factor * x.spread[ idx ] ) ) ),
        c( x.predict[[ area ]][ gate.grid.n ], predictions[ gate.grid.n ] ),
        t( sapply( rev( seq_len( gate.grid.n ) ), function( idx )
            c( x.predict[[ area ]][ idx ],
               predictions[ idx ] + spread.factor * x.spread[ idx ] ) ) )
    )

    colnames( gate_vertices ) <- channel_names
    flowCore::polygonGate( gate_vertices, filterId = filterId )
}
