# get_autogating_param_aurora.r
#
#' Instrument-specific parameter overrides for the Cytek Aurora
#'
#' @title Cytek Aurora Autogating Parameters
#'
#' @description
#' Applies Cytek Aurora-specific values to the base autogating parameter list:
#' figure axis limits and Gating-ML transformation parameters.
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

    # ---- figure axis limits --------------------------------------------------

    ap$figure.limits <- list(
        scatter.min = 0,
        scatter.max = 4200000,
        expr.min    = -5000,
        expr.max    = 4200000
    )

    # ---- default transformation parameters ------------------------------------
    # Gating-ML logicle for fluorescence (W from FlowJo width basis -1000),
    # linear for scatter; both map the top of the range to 1.  A
    # per-experiment transformation CSV overrides individual channels.

    ap$transform.defaults <- list(
        fluorescence = list( id = 1L, scale.t = 4194304,
                             logicle.w = log10( 1000 ) / 2,
                             logicle.m = 5.62, logicle.a = 0 ),
        scatter      = list( id = 0L, scale.t = 4194304, linear.a = 0 )
    )

    ap
}
