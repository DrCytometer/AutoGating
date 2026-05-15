# ==============================================================================
# do.activation.gates
# ==============================================================================

#' Run all activation marker gates in sequence
#'
#' @title Execute Activation Gating Pipeline
#'
#' @description
#' Iterates over a list of activation gate definitions.  For each definition,
#' iterates over each activation marker and constructs a synthetic gate
#' definition compatible with \code{\link{do.gate.generic}}, then calls that
#' function to calculate boundaries and assign events.  Display-only rows
#' (those with \code{gate.type == ""} or an empty \code{markers} vector) are
#' silently skipped.  Per-marker algorithm overrides are respected via
#' \code{marker.algorithm.override}.
#'
#' @param fcs.gates Named list.  Accumulator for gate results; updated in
#'   place and returned on completion.
#' @param fcs.gates.data Numeric matrix.  Full expression matrix (rows =
#'   events, columns = channels/antigens).
#' @param act.defs List of activation gate definition lists.  Each element
#'   must contain: \code{gate.type}, \code{markers} (character vector),
#'   \code{parent.gate}, \code{parent.popul}, \code{popul.label},
#'   \code{label.pos}, \code{gate.marker.x}, \code{gate.marker.y},
#'   \code{algorithm}, \code{gate.param}, \code{popul.number}, and optionally
#'   \code{marker.algorithm.override} (named list of marker → algorithm
#'   strings).
#' @param ap Autogating parameter object returned by
#'   \code{\link{get.autogating.param}}.
#' @param calculate.gates Logical.  Whether to (re-)calculate gate boundaries
#'   from data.
#' @param figure.dir Character.  Output directory for diagnostic PNG figures.
#' @param plot.gates Logical.  Whether to produce figures for the gates.
#'
#' @return Updated \code{fcs.gates} list with all activation gates applied.
#'
#' @seealso \code{\link{do.gate.generic}}, \code{\link{do.population.gates}}
#'
#' @export
do.activation.gates <- function(
    fcs.gates,
    fcs.gates.data,
    act.defs,
    ap,
    calculate.gates,
    figure.dir,
    plot.gates = TRUE
  ) {

  for ( ad in act.defs ) {

    # Skip display-only rows
    if ( ad$gate.type == "" || length( ad$markers ) == 0L )
      next

    parent.gate  <- ad$parent.gate
    parent.popul <- ad$parent.popul
    popul.label  <- ad$popul.label
    label.pos    <- ad$label.pos

    for ( marker in ad$markers ) {

      gate.name <- sprintf( "%s.%s.%s", parent.gate, parent.popul, marker )

      # Determine algorithm for this marker
      alg <- if ( marker %in% names( ad$marker.algorithm.override ) )
        ad$marker.algorithm.override[[ marker ]]
      else
        ad$algorithm

      # Build an activation gate definition compatible with do.gate.generic.
      # "*" (ap$ignore) is a wildcard meaning "use the current marker on this
      # axis"; a non-"*" value is the fixed population axis marker (e.g. "CD4").
      act.gd <- list(
        gate.number   = ad$popul.number,
        gate.name     = gate.name,
        gate.type     = ad$gate.type,
        gate.marker.x = if ( ad$gate.marker.x == ap$ignore ) marker
        else ad$gate.marker.x,
        gate.marker.y = if ( ad$gate.marker.y == ap$ignore ) marker
        else ad$gate.marker.y,
        parent.gate   = parent.gate,
        parent.popul  = parent.popul,
        populations = list(
          `*`       = list( label = "*", label.pos = 0L ),
          activated = list(
            label     = sprintf( "%s+ %s", marker, popul.label ),
            label.pos = label.pos
          )
        ),
        algorithm     = alg,
        gate.param    = ad$gate.param,
        stats.parent  = list( gate = ap$ignore, popul = ap$ignore )
      )

      # Columns to pass to do.gate.generic: the resolved x and y markers
      # (both wildcards now substituted), preserving order x then y so that
      # colnames(fcs.gate.data) == c(gate.marker.x, gate.marker.y) holds.
      gate.marker <- unique( c( act.gd$gate.marker.x, act.gd$gate.marker.y ) )

      message( sprintf( "%0*d - %s",
                        ap$gate.number.width, ad$popul.number, gate.name ) )

      fcs.gates <- do.gate.generic(
        fcs.gates,
        fcs.gates.data[ , gate.marker ],
        act.gd,
        ap,
        calculate.gate = calculate.gates,
        figure.dir     = figure.dir,
        plot.gates     = plot.gates
      )
    }
  }

  return(fcs.gates)
}
