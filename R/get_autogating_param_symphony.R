# get_autogating_param_symphony.r
#
#' Instrument-specific parameter overrides for the BD FACSymphony
#'
#' @title BD FACSymphony Autogating Parameters
#'
#' @description
#' Applies BD FACSymphony-specific values to the base autogating parameter
#' list: the singlet scatter pair, figure axis limits and Gating-ML
#' transformation parameters.
#' Called automatically by \code{\link{get.autogating.param}} when
#' \code{cytometer = "symphony"}.
#'
#' @param ap Named list.  The base parameter list returned by
#'   \code{get.autogating.param.minimal()}.
#'
#' @return Updated \code{ap} list with Symphony-specific values applied.
#'
#' @seealso \code{\link{get.autogating.param}}
#'
#' @keywords internal

get.autogating.param.symphony <- function( ap ) {

    ap$cytometer <- "Symphony"

    # ---- singlet scatter pair ------------------------------------------------
    ap$singlet.scatter.parameter <- c( "FSC-A", "FSC-H" )

    # ---- figure axis limits --------------------------------------------------

    ap$figure.limits <- list(
        scatter.min = 0,
        scatter.max = 265000,
        expr.min    = -2000,
        expr.max    = 265000
    )

    # ---- default transformation parameters ------------------------------------
    # Gating-ML logicle for fluorescence (W from FlowJo width basis -100),
    # linear for scatter; both map the top of the range to 1.  A
    # per-experiment transformation CSV overrides individual channels.

    ap$transform.defaults <- list(
        fluorescence = list( id = 1L, scale.t = 262144,
                             logicle.w = log10( 100 ) / 2,
                             logicle.m = 4.418539922, logicle.a = 0 ),
        scatter      = list( id = 0L, scale.t = 262144, linear.a = 0 )
    )

    ap
}
