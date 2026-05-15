# gate_fs_singlet.R
#
# Singlet gate via flowStats::gate_singlet().
#
# This is an alternative to gate.singlet() (which uses a bespoke robust linear
# model).  gate.fs.singlet() is a thin wrapper around the flowStats
# implementation, normalising input/output so that it is a drop-in replacement
# from the perspective of calculate.boundary.fs.singlet().
#

#' Singlet gate using \code{flowStats::gate_singlet}
#'
#' @title Singlet Gate via flowStats
#'
#' @description
#' Identifies singlet events by calling \code{\link[flowStats]{gate_singlet}},
#' which fits a robust linear model to the area-vs-height scatter relationship
#' and constructs an elliptical polygon gate around the singlet population.
#' This function accepts the same \code{matrix} input and returns the same
#' \code{\link[flowCore]{polygonGate}} output as \code{\link{gate.singlet}},
#' making it a drop-in replacement within the AutoGating pipeline when the
#' algorithm is set to \code{"fs.singlet"}.
#'
#' @param x Numeric matrix.  Must contain at least the columns named by
#'   \code{area} and \code{height}.
#' @param area Character scalar.  Name of the scatter-area channel column.
#'   Default \code{"FSC-A"}.
#' @param height Character scalar.  Name of the scatter-height channel column.
#'   Default \code{"FSC-H"}.
#' @param filterId Character scalar.  Identifier for the resulting
#'   \code{polygonGate}.  Default \code{"singlet"}.
#' @param wider Numeric scalar passed to \code{\link[flowStats]{gate_singlet}}
#'   as the \code{wider} argument.  Controls the width of the gate envelope;
#'   larger values admit more doublets but reduce false exclusions.
#'   Default \code{0.1}.
#' @param prediction.level Numeric scalar in \eqn{(0, 1)} passed to
#'   \code{\link[flowStats]{gate_singlet}}.  The coverage probability used to
#'   set the gate boundaries (analogous to a prediction interval).
#'   Default \code{0.99}.
#' @param subsample Integer or \code{NULL}.  Number of events to subsample
#'   before calling \code{flowStats::gate_singlet} to limit computation time
#'   on large files.  The gate is fitted on the subsample but the returned
#'   \code{polygonGate} is applied to all events downstream.  Set to
#'   \code{NULL} to use all events.  Default \code{50000L}.
#' @param ... Additional arguments passed to
#'   \code{\link[flowStats]{gate_singlet}}.
#'
#' @return A \code{\link[flowCore]{polygonGate}} object in the area-vs-height
#'   plane, directly usable with \code{\link[flowCore]{filter}} or the
#'   AutoGating point-in-polygon assignment via \code{\link{apply.gate}}.
#'
#' @importFrom flowStats gate_singlet
#' @importFrom flowCore flowFrame polygonGate
#'
#' @seealso \code{\link{gate.singlet}}, \code{\link[flowStats]{gate_singlet}}
#'
#' @export
gate.fs.singlet <- function(
    x,
    area             = "FSC-A",
    height           = "FSC-H",
    filterId         = "singlet",
    wider            = 0.1,
    prediction.level = 0.99,
    subsample        = 50000L,
    ...
) {
  stopifnot( is.matrix( x ) )
  stopifnot( is.character( area )   && length( area )   == 1L )
  stopifnot( is.character( height ) && length( height ) == 1L )
  stopifnot( area   %in% colnames( x ) )
  stopifnot( height %in% colnames( x ) )

  # ---- optional subsampling ------------------------------------------------
  # flowStats::gate_singlet can be slow on large FCS files; fitting on a
  # representative subsample is usually sufficient because the singlet ridge
  # is highly stable.

  if ( !is.null( subsample ) && nrow( x ) > subsample ) {
    set.seed( 42L )
    fit.idx <- sample.int( nrow( x ), subsample )
    x.fit   <- x[ fit.idx, , drop = FALSE ]
  } else {
    x.fit <- x
  }

  # ---- build a minimal flowFrame for the two required channels -------------
  # flowStats::gate_singlet requires a flowFrame; we create one containing
  # only the area and height columns to keep memory usage minimal.

  mat <- x.fit[ , c( area, height ), drop = FALSE ]

  # flowFrame requires a numeric matrix with channel metadata; build the
  # minimal AnnotatedDataFrame that flowCore expects.
  param.df <- data.frame(
    name        = c( area, height ),
    desc        = c( area, height ),
    range       = c( max( mat[ , 1L ] ), max( mat[ , 2L ] ) ),
    minRange    = c( min( mat[ , 1L ] ), min( mat[ , 2L ] ) ),
    maxRange    = c( max( mat[ , 1L ] ), max( mat[ , 2L ] ) ),
    row.names   = c( area, height ),
    stringsAsFactors = FALSE
  )
  param.adf <- Biobase::AnnotatedDataFrame( param.df )

  ff <- flowCore::flowFrame(
    exprs      = mat,
    parameters = param.adf
  )

  # ---- call flowStats::gate_singlet ----------------------------------------

  fs.gate <- flowStats::gate_singlet(
    x                = ff,
    area             = area,
    height           = height,
    filterId         = filterId,
    wider            = wider,
    prediction.level = prediction.level,
    ...
  )

  # flowStats::gate_singlet returns a polygonGate directly — return as-is.
  # The boundaries slot uses the original channel names, which is consistent
  # with what gate.singlet() produces and what apply.gate() expects.
  fs.gate
}
