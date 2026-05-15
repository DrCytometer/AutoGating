# get_autogating_param_aurora.r
#
#' Instrument-specific parameter overrides for the Cytek Aurora
#'
#' @title Cytek Aurora Autogating Parameters
#'
#' @description
#' Applies Cytek Aurora-specific values to the base autogating parameter list:
#' scatter/fluorescence data ranges, figure axis limits, default
#' \code{flowjo_biexp} transformation parameters, and gate algorithm defaults.
#' Called automatically by \code{\link{get.autogating.param}} when
#' \code{cytometer = "aurora"}.
#'
#' @param ap Named list.  The base parameter list returned by
#'   \code{get.autogating.param.minimal()}.
#'
#' @return Updated \code{ap} list with Aurora-specific values applied.
#'
#' @seealso \code{\link{get.autogating.param}}
#'
#' @keywords internal
get.autogating.param.aurora <- function( ap ) {

    ap$cytometer <- "Aurora"

    # ---- raw data ranges (scatter channels) -----------------------------------
    ap$default.scatter.param <- c( "FSC-A", "SSC-A" )
    ap$time.param <- c( "Time" )
    ap$scatter.data.min <- 0
    ap$scatter.data.max <- 4194304

    # ---- raw data ranges (fluorescence channels) ------------------------------

    ap$expr.data.min <- -111
    ap$expr.data.max <- 4194304

    # ---- figure axis limits --------------------------------------------------

    ap$figure.limits <- list(
        scatter.min = 0,
        scatter.max = 4200000,
        expr.min    = -5000,
        expr.max    = 4200000
    )

    # ---- default transformation parameters ------------------------------------
    # All fluorescence channels use these unless overridden by a per-experiment
    # transformation CSV.

    ap$transformation.defaults <- list(
        length    = 256,
        max.range = 4194304,
        pos       = 5.62,
        neg       = 0,
        width     = -1000
    )

    # ---- gate algorithm defaults ----------------------------------------------
    # Override only what differs from the universal minimal defaults.

    ap$gate.defaults$density.trimneg.x <- TRUE
    ap$gate.defaults$density.trimneg.y <- FALSE

    ap$gate.defaults$break.x        <- c( -1, -0.1, 0, 0.1, 1, 10, 100 )
    ap$gate.defaults$break.x.factor <- 1e4
    ap$gate.defaults$break.y        <- c( -1, -0.1, 0, 0.1, 1, 10, 100 )
    ap$gate.defaults$break.y.factor <- 1e4

    # Scatter-channel plot breaks (FSC / SSC axes)
    ap$gate.defaults$break.scatter        <- c( 0, 1000, 2000, 3000, 4000, 4500 )
    ap$gate.defaults$break.scatter.factor <- 1000
    ap$gate.defaults$break.fluor          <- c( -1e3, 0, 1e3, 1e4, 1e5, 1e6 )

    ap$gate.defaults$label.size <- 4

    ap$gate.defaults$algorithm.1dsep         <- "tail"
    ap$gate.defaults$algorithm.2dsep         <- c( "mixture", "mixture" )
    ap$gate.defaults$default.tail.probability <- 0.99

    ap
}
