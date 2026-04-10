
# ==============================================================================
# do.gate.generic  (orchestrator)
# ==============================================================================

#' Orchestrate boundary calculation, event assignment, and plotting for one gate
#'
#' @title Generic Gate Orchestrator
#'
#' @description
#' The central workhorse of the gating pipeline.  For a single gate defined by
#' \code{gd}, it:
#' \enumerate{
#'   \item Resolves the full parameter set via \code{\link{resolve.gate.param}}.
#'   \item Retrieves parent-population event indices from \code{fcs.gates} (or
#'     uses all events for gate number 0).
#'   \item Optionally calls the appropriate \code{calculate.gate.<type>()}
#'     function to compute and store polygon boundary matrices.
#'   \item Calls \code{\link{apply.gate}} to run point-in-polygon assignment.
#'   \item Calls \code{\link{gate.plot}} to generate and save a diagnostic PNG
#'     (skipped for gate number 0).
#' }
#'
#' @param fcs.gates Named list.  Accumulator for gate results.
#' @param fcs.gate.data Numeric matrix with exactly two columns corresponding
#'   to \code{gd$gate.marker.x} and \code{gd$gate.marker.y}.
#' @param gd Named list.  Gate definition.  Required fields:
#'   \code{gate.number}, \code{gate.name}, \code{gate.type},
#'   \code{gate.marker.x}, \code{gate.marker.y}, \code{parent.gate},
#'   \code{parent.popul}, \code{populations}, \code{algorithm},
#'   \code{gate.param}.
#' @param ap Autogating parameter object returned by
#'   \code{\link{get.autogating.param}}.
#' @param calculate.gate Logical.  If \code{TRUE}, compute gate boundaries;
#'   if \code{FALSE}, reuse boundaries already stored in \code{fcs.gates}.
#' @param figure.dir Character.  Output directory for diagnostic PNG figures.
#'
#' @return Updated \code{fcs.gates} list with boundary, event indices,
#'   fraction, and label populated for each non-ignored population of this
#'   gate.
#'
#' @seealso \code{\link{apply.gate}}, \code{\link{gate.plot}},
#'   \code{\link{resolve.gate.param}}, \code{\link{calculate.gate.free}},
#'   \code{\link{calculate.gate.1dsep}}, \code{\link{calculate.gate.2dsep}}
#'
#' @keywords internal
do.gate.generic <- function( fcs.gates, fcs.gate.data, gd, ap,
                             calculate.gate, figure.dir ) {

  # ---- resolve full parameter set ------------------------------------------

  p <- resolve.gate.param( gd, ap )

  gate.name    <- gd$gate.name
  gate.number  <- gd$gate.number
  gate.marker.x <- gd$gate.marker.x
  gate.marker.y <- gd$gate.marker.y
  popul.names   <- names( gd$populations )

  stopifnot( colnames( fcs.gate.data ) == c( gate.marker.x, gate.marker.y ) )

  # ---- get parent event indices --------------------------------------------
  parent.gate  <- gd$parent.gate
  parent.popul <- gd$parent.popul

  if ( gate.number > 0L ) {

    stopifnot( length( parent.gate ) == length( parent.popul ) )

    gate.parent.idx <- sort(
      unlist(
        lapply( seq_along( parent.gate ),
                function( i )
                  fcs.gates[[ parent.gate[ i ] ]][[ parent.popul[ i ] ]]$idx
        )
      )
    )

    stopifnot( anyDuplicated( gate.parent.idx ) == 0L )
  } else {
    gate.parent.idx <- seq_len( nrow( fcs.gate.data ) )
  }

  gate.data <- fcs.gate.data[ gate.parent.idx, , drop = FALSE ]

  # ---- calculate gate boundary if required ---------------------------------

  if ( calculate.gate ) {
    calc.fn <- get( paste0( "calculate.gate.", gd$gate.type ) )
    message( sprintf( "  -> %s()", paste0( "calculate.gate.", gd$gate.type ) ) )
    gate.boundary <- calc.fn(
      gate.data, popul.names, gd$algorithm, p,
      gate.number, gate.name, fcs.gates, figure.dir, ap
    )

    fcs.gates[[ gate.name ]] <- vector( "list", length( gate.boundary ) )
    names( fcs.gates[[ gate.name ]] ) <- names( gate.boundary )

    for ( pop in names( gate.boundary ) )
      fcs.gates[[ gate.name ]][[ pop ]] <- list(
        parent.gate  = parent.gate,
        parent.popul = parent.popul,
        boundary     = gate.boundary[[ pop ]]
      )
  }

  # ---- apply gate boundaries (PIP) -----------------------------------------

  fcs.gates <- apply.gate( fcs.gates, gate.data, gate.parent.idx,
                           gate.name, gd$populations, ap )

  # ---- plot ----------------------------------------------------------------

  if ( gate.number > 0L )
    gate.plot( fcs.gates, gate.data, gate.name, gate.number,
               gate.marker.x, gate.marker.y, gd$populations, p,
               figure.dir, ap )

  return(fcs.gates)
}
