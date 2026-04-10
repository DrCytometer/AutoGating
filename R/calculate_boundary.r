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

  scatter <- ap$singlet.scatter.parameter   # e.g. c( "FSC-A", "FSC-H" )

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
