
# ==============================================================================
# apply.gate
#
# Run point.in.polygon for every non-ignored population and store idx /
# fraction / label back into fcs.gates.
# ==============================================================================

#' Assign events to gate populations using point-in-polygon testing
#'
#' @title Apply Gate Boundaries to Events
#'
#' @description
#' For each non-ignored population within a gate, runs
#' \code{\link[sp]{point.in.polygon}} against the stored boundary matrix and
#' writes back the event indices, population fraction, and display label into
#' \code{fcs.gates}.
#'
#' @param fcs.gates Named list.  Gate results accumulator; must already
#'   contain a \code{fcs.gates[[gate.name]][[pop.name]]$boundary} matrix for
#'   each population being evaluated.
#' @param gate.data Numeric matrix.  Two-column event data for the parent
#'   population (rows = events, columns = x and y channels).
#' @param gate.parent.idx Integer vector.  Row indices (into the full event
#'   matrix) corresponding to the rows of \code{gate.data}.
#' @param gate.name Character.  Key into \code{fcs.gates} identifying this
#'   gate.
#' @param populations Named list of population metadata lists.  Each element
#'   must contain at least \code{label} (character display label) and the
#'   population name must match a key already present in
#'   \code{fcs.gates[[gate.name]]}.  Populations whose name equals
#'   \code{ap$ignore} are skipped.
#' @param ap Autogating parameter object; uses \code{ap$ignore} sentinel.
#'
#' @return Updated \code{fcs.gates} list.  For each non-ignored population,
#'   the sub-list \code{fcs.gates[[gate.name]][[pop.name]]} gains three new
#'   elements:
#'   \describe{
#'     \item{\code{label}}{Character.  Display label copied from
#'       \code{populations}.}
#'     \item{\code{idx}}{Integer vector.  Absolute event row indices passing
#'       the gate.}
#'     \item{\code{fraction}}{Numeric scalar in \eqn{[0, 1]}.  Proportion of
#'       parent events inside the boundary.}
#'   }
#'
#' @importFrom sp point.in.polygon
#'
#' @seealso \code{\link{do.gate.generic}}
#'
#' @keywords internal
apply.gate <- function( fcs.gates, gate.data, gate.parent.idx, gate.name,
                        populations, ap ) {

  for ( pop.name in names( populations ) ) {

    if ( pop.name == ap$ignore )
      next

    pop.info     <- populations[[ pop.name ]]
    pop.boundary <- fcs.gates[[ gate.name ]][[ pop.name ]]$boundary

    pip <- point.in.polygon(
      gate.data[ , 1 ], gate.data[ , 2 ],
      pop.boundary[ , 1 ], pop.boundary[ , 2 ]
    )

    fcs.gates[[ gate.name ]][[ pop.name ]]$label    <- pop.info$label
    fcs.gates[[ gate.name ]][[ pop.name ]]$idx      <-
      gate.parent.idx[ pip == 1L ]
    fcs.gates[[ gate.name ]][[ pop.name ]]$fraction <- mean( pip == 1L )
  }

  return(fcs.gates)
}
