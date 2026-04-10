
# ==============================================================================
# do.population.gates
# ==============================================================================

#' Run all population gates in sequence
#'
#' @title Execute Population Gating Pipeline
#'
#' @description
#' Iterates over a list of population gate definitions and calls
#' \code{\link{do.gate.generic}} for each one, passing the appropriate
#' two-column data slice from the full expression matrix.  Progress is logged
#' to the console via \code{message()}.
#'
#' @param fcs.gates Named list.  Accumulator for gate results; updated in
#'   place across iterations and returned on completion.
#' @param fcs.gates.data Numeric matrix.  Full expression matrix (rows =
#'   events, columns = channels/antigens) for all samples in the current
#'   batch.
#' @param gate.defs List of gate definition lists, each following the
#'   standard gate definition schema (fields: \code{gate.number},
#'   \code{gate.name}, \code{gate.type}, \code{gate.marker.x},
#'   \code{gate.marker.y}, \code{parent.gate}, \code{parent.popul},
#'   \code{populations}, \code{algorithm}, \code{gate.param},
#'   \code{stats.parent}).
#' @param ap Autogating parameter object returned by
#'   \code{\link{get.autogating.param}}.
#' @param calculate.gates Logical.  If \code{TRUE}, gate boundaries are
#'   (re-)calculated from data; if \code{FALSE}, previously stored boundaries
#'   in \code{fcs.gates} are used directly and only event assignment is
#'   repeated.
#' @param figure.dir Character.  Output directory for diagnostic PNG figures.
#'
#' @return Updated \code{fcs.gates} list with all population gates applied.
#'
#' @seealso \code{\link{do.gate.generic}}, \code{\link{do.activation.gates}}
#'
#' @export
do.population.gates <- function( fcs.gates, fcs.gates.data, gate.defs, ap,
                                 calculate.gates, figure.dir ) {

  for ( gd in gate.defs ) {
    gate.marker <- c( gd$gate.marker.x, gd$gate.marker.y )

    message(
      sprintf( "%0*d - %s", ap$gate.number.width, gd$gate.number, gd$gate.name )
    )

    fcs.gates <- do.gate.generic(
      fcs.gates,
      fcs.gates.data[ , gate.marker ],
      gd,
      ap,
      calculate.gate = calculate.gates,
      figure.dir     = figure.dir
    )
  }

  return(fcs.gates)
}
