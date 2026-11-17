# get_autogating_param_id7000.r
#
#' Instrument-specific parameter overrides for the Sony ID7000
#'
#' @title Sony ID7000 Autogating Parameters
#'
#' @description
#' Applies Sony ID7000-specific values to the base autogating parameter list:
#' figure axis limits and Gating-ML transformation parameters.
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

    # ---- figure axis limits --------------------------------------------------

    ap$figure.limits <- list(
        scatter.min = 0,
        scatter.max = 1050000,
        expr.min    = -2000,
        expr.max    = 1050000
    )

    # ---- default transformation parameters ------------------------------------
    # Gating-ML logicle for fluorescence (W from FlowJo width basis -250),
    # linear for scatter; both map the top of the range to 1.  A
    # per-experiment transformation CSV overrides individual channels.

    ap$transform.defaults <- list(
        fluorescence = list( id = 1L, scale.t = 1000000,
                             logicle.w = log10( 250 ) / 2,
                             logicle.m = 5, logicle.a = 0 ),
        scatter      = list( id = 0L, scale.t = 1048576, linear.a = 0 )
    )

    ap
}
