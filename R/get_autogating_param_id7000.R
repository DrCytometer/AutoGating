# get_autogating_param_id7000.r
#
#' Instrument-specific parameter overrides for the Sony ID7000
#'
#' @title Sony ID7000 Autogating Parameters
#'
#' @description
#' Applies Sony ID7000-specific values to the base autogating parameter list:
#' scatter/fluorescence data ranges, figure axis limits, default
#' \code{flowjo_biexp} transformation parameters, and gate algorithm defaults.
#' Called automatically by \code{\link{get.autogating.param}} when
#' \code{cytometer = "id7000"}.
#'
#' @param ap Named list.  The base parameter list returned by
#'   \code{get.autogating.param.minimal()}.
#'
#' @return Updated \code{ap} list with ID7000-specific values applied.
#'
#' @seealso \code{\link{get.autogating.param}}
#'
#' @keywords internal

get.autogating.param.id7000 <- function( ap ) {

    ap$cytometer <- "ID7000"

    # ---- raw data ranges (scatter channels) -----------------------------------
    ap$default.scatter.param <- c( "FSC-A", "SSC-A" )
    ap$time.param <- c( "TIME" )
    ap$scatter.data.min <- 0
    ap$scatter.data.max <- 1048576

    # ---- raw data ranges (fluorescence channels) ------------------------------

    ap$expr.data.min <- -111
    ap$expr.data.max <- 1000000

    # ---- figure axis limits --------------------------------------------------

    ap$figure.limits <- list(
        scatter.min = 0,
        scatter.max = 1050000,
        expr.min    = -2000,
        expr.max    = 1050000
    )

    # ---- default transformation parameters ------------------------------------
    # All fluorescence channels use these unless overridden by a per-experiment
    # transformation CSV.  The only channel that differs in practice for st1 is
    # BYG790-A (CCR7): neg = 1, width = -20, pos = 4.418539922.

    ap$transformation.defaults <- list(
        length    = 256,
        max.range = 1000000,
        pos       = 5,
        neg       = 0,
        width     = -250
    )

    # ---- gate algorithm defaults ----------------------------------------------
    # Override only what differs from the universal minimal defaults.

    ap$gate.defaults$density.trimneg.x <- TRUE
    ap$gate.defaults$density.trimneg.y <- FALSE

    ap$gate.defaults$break.x        <- c( -1, 0, 1, 10, 100 )
    ap$gate.defaults$break.x.factor <- 1000
    ap$gate.defaults$break.y        <- c( -1, 0, 1, 10, 100 )
    ap$gate.defaults$break.y.factor <- 1000

    # Scatter-channel plot breaks (FSC / SSC axes)
    ap$gate.defaults$break.scatter        <- c( 0, 200, 400, 600, 800, 1000 )
    ap$gate.defaults$break.scatter.factor <- 1000
    ap$gate.defaults$break.fluor          <- c( -1e3, 0, 1e3, 1e4, 1e5, 1e6 )

    ap$gate.defaults$label.size <- 4

    ap$gate.defaults$algorithm.1dsep         <- "tail"
    ap$gate.defaults$algorithm.2dsep         <- c( "bimodal", "bimodal" )
    ap$gate.defaults$default.tail.probability <- 0.99

    return(ap)
}
