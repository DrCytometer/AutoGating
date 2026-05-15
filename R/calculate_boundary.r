# calculate_boundary.r
#
# Boundary calculation functions for "free"-type gates.
#
# calculate.boundary.all()       — bounding box of all data (gate 0)
# calculate.boundary.singlets()  — singlet polygon via gate.singlet()
#
# Both receive (marker.data, p) where p is the resolved parameter list from
# resolve.gate.param().  make.rect.boundary() is defined in analyze_aux.r.

#' Calculate a bounding-box boundary encompassing all events
#'
#' @title All-Events Bounding Box Gate Boundary
#'
#' @description
#' Constructs a slightly padded rectangular boundary that encompasses the full
#' range of the two-column event data.  Used as gate 0 (the root gate) to
#' include every event in the dataset.
#'
#' @param marker.data Numeric matrix (n × 2).  Full event data.
#' @param p Named list.  Resolved gate parameters (currently unused but
#'   required for interface consistency).
#' @param ap Autogating parameter object (currently unused but required for
#'   interface consistency).
#'
#' @return A 4 × 2 numeric matrix representing the closed rectangle, as
#'   produced by \code{\link{make.rect.boundary}}.
#'
#' @seealso \code{\link{calculate.gate.free}}, \code{\link{make.rect.boundary}}
#'
#' @keywords internal
calculate.boundary.all <- function( marker.data, p, ap ) {

  xmin <- min( marker.data[ , 1L ] )
  xmax <- max( marker.data[ , 1L ] )
  ymin <- min( marker.data[ , 2L ] )
  ymax <- max( marker.data[ , 2L ] )

  make.rect.boundary(
    1.01 * xmin - 0.01 * xmax,
    -0.01 * xmin + 1.01 * xmax,
    1.01 * ymin - 0.01 * ymax,
    -0.01 * ymin + 1.01 * ymax
  )
}

#' Calculate a singlet-discrimination gate boundary
#'
#' @title Singlet Gate Boundary
#'
#' @description
#' Fits a robust linear model to the area-vs-height scatter relationship using
#' \code{\link{gate.singlet}} and returns the resulting polygon boundary
#' matrix.  Scatter channel names are taken from
#' \code{ap$singlet.scatter.parameter}.
#'
#' @param marker.data Numeric matrix (n × 2).  Event data with columns
#'   corresponding to the area and height scatter channels specified in
#'   \code{ap$singlet.scatter.parameter}.
#' @param p Named list.  Resolved gate parameters (passed through to
#'   \code{\link{gate.singlet}} for any future extensibility).
#' @param ap Autogating parameter object; uses
#'   \code{ap$singlet.scatter.parameter} to identify the area
#'   (\code{[1]}) and height (\code{[2]}) channels.
#'
#' @return A numeric matrix of polygon vertices (rows = vertices, 2 columns)
#'   with dimension names removed.
#'
#' @seealso \code{\link{gate.singlet}}, \code{\link{calculate.gate.free}}
#'
#' @keywords internal
calculate.boundary.singlets <- function( marker.data, p, ap ) {

  scatter <- colnames( marker.data )
  stopifnot( length( unique( scatter ) ) == 2 ) # must be separate parameters
  # is it FSC or SSC?
  fsc.active <- length( grep( "FSC", scatter, ignore.case = TRUE ) ) == 2L
  ssc.active <- length( grep( "SSC", scatter, ignore.case = TRUE ) ) == 2L
  stopifnot( fsc.active | ssc.active ) # must be the same parameter's Area vs H/W
  stopifnot( any( grepl( "-A", scatter, ignore.case = TRUE ) ) ) # one must be area

  # could reverse order here if Area isn't first

  res <- gate.singlet(
    marker.data,
    area           = scatter[ 1L ],
    height         = scatter[ 2L ],
    maxit          = 100,
    gate.grid.n    = 20,
    gate.skip.low  = 0.06,
    gate.skip.high = 0.10,
    spread.span    = 1.0,
    spread.factor  = 7,
    spread.spline.x = c( 0.00, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 0.90, 1.00 ),
    spread.spline.y = c( 2.00, 2.10, 1.90, 1.30, 0.90, 0.45, 0.30, 0.25, 0.25 )
  )

  boundary <- res@boundaries
  dimnames( boundary ) <- NULL
  boundary
}

# ---- Addition to calculate_boundary.r ---------------------------------------
#
# Paste this block immediately after calculate.boundary.singlets() in
# calculate_boundary.r.
# ------------------------------------------------------------------------------

#' Calculate a singlet-discrimination gate boundary using flowStats
#'
#' @title Singlet Gate Boundary via flowStats
#'
#' @description
#' Fits a singlet gate using \code{\link{gate.fs.singlet}} (which wraps
#' \code{\link[flowStats]{gate_singlet}}) and returns the polygon boundary
#' matrix.  This function is the \code{"fs.singlet"} counterpart to
#' \code{\link{calculate.boundary.singlets}} (\code{"singlets"} algorithm) and
#' is discovered automatically by \code{\link{calculate.gate.free}} via
#' \code{get("calculate.boundary.fs.singlet")}.
#'
#' Scatter channel names are taken from \code{ap$singlet.scatter.parameter},
#' exactly as in \code{calculate.boundary.singlets}.  The optional tuning
#' parameters \code{wider}, \code{prediction.level}, and \code{subsample} may
#' be supplied in the resolved gate parameter list \code{p}; if absent,
#' defaults from \code{\link{gate.fs.singlet}} are used.
#'
#' @param marker.data Numeric matrix (n × 2).  Event data with columns
#'   corresponding to the area and height scatter channels specified in
#'   \code{ap$singlet.scatter.parameter}.
#' @param p Named list.  Resolved gate parameters.  Recognised optional fields:
#'   \describe{
#'     \item{\code{fs.wider}}{Passed as \code{wider} to
#'       \code{\link{gate.fs.singlet}}.  Default \code{0.1}.}
#'     \item{\code{fs.prediction.level}}{Passed as \code{prediction.level}.
#'       Default \code{0.99}.}
#'     \item{\code{fs.subsample}}{Passed as \code{subsample}.  Default
#'       \code{50000L}.}
#'   }
#' @param ap Autogating parameter object; uses
#'   \code{ap$singlet.scatter.parameter} to identify the area (\code{[1]}) and
#'   height (\code{[2]}) channels.
#'
#' @return A numeric matrix of polygon vertices (rows = vertices, 2 columns)
#'   with dimension names removed, ready for storage in \code{fcs.gates}.
#'
#' @seealso \code{\link{gate.fs.singlet}}, \code{\link{calculate.gate.free}},
#'   \code{\link{calculate.boundary.singlets}}
#'
#' @keywords internal
calculate.boundary.fs.singlets <- function( marker.data, p, ap ) {

  scatter <- colnames( marker.data )
  stopifnot( length( unique( scatter ) ) == 2 ) # must be separate parameters
  # is it FSC or SSC?
  fsc.active <- length( grep( "FSC", scatter, ignore.case = TRUE ) ) == 2L
  ssc.active <- length( grep( "SSC", scatter, ignore.case = TRUE ) ) == 2L
  stopifnot( fsc.active | ssc.active ) # must be the same parameter's Area vs H/W
  stopifnot( any( grepl( "-A", scatter, ignore.case = TRUE ) ) ) # one must be area

  # Pull optional tuning overrides from the resolved parameter list, falling
  # back to gate.fs.singlet() defaults when absent.
  wider            <- if ( !is.null( p$fs.wider ) )            p$fs.wider            else 0.1
  prediction.level <- if ( !is.null( p$fs.prediction.level ) ) p$fs.prediction.level else 0.99
  subsample        <- if ( !is.null( p$fs.subsample ) )        p$fs.subsample        else 50000L

  res <- gate.fs.singlet(
    marker.data,
    area             = scatter[ 1L ],
    height           = scatter[ 2L ],
    wider            = wider,
    prediction.level = prediction.level,
    subsample        = subsample
  )

  boundary <- res@boundaries
  dimnames( boundary ) <- NULL
  boundary
}

