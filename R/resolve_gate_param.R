# ==============================================================================
# resolve.gate.param
#
# Merge cytometer defaults with per-gate overrides.  Returns a flat list
# containing every parameter key the pipeline needs.
# ==============================================================================

#' Resolve the full parameter set for a single gate
#'
#' @title Resolve Gate Parameters
#'
#' @description
#' Merges cytometer-level gate defaults (\code{ap$gate.defaults}) with the
#' per-gate overrides stored in \code{gd$gate.param}.  The gate-level values
#' take precedence.  The resulting flat list is passed as \code{p} to all
#' downstream gate calculation and plotting functions.
#'
#' @param gd Named list.  A single gate definition as found in the population-
#'   or activation-gate definition list.  Must contain a \code{gate.param}
#'   element (may be an empty list).
#' @param ap Autogating parameter object returned by
#'   \code{\link{get.autogating.param}}.
#'
#' @return A named list containing every parameter key in
#'   \code{ap$gate.defaults}, with any matching keys in \code{gd$gate.param}
#'   overriding the defaults.
#'
#' @seealso \code{\link{do.gate.generic}}, \code{\link{merge.gate.param}}
#'
#' @keywords internal
resolve.gate.param <- function( gd, ap ) {
  p <- ap$gate.defaults

  for ( key in names( gd$gate.param ) )
    p[[ key ]] <- gd$gate.param[[ key ]]

  return(p)
}
