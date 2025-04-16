# get_autogating_param_symphony.r

#' @title get.autogating.param.symphony
#'
#' @description
#' Get cytometer-specific parameters for FACSymphony files in order to perform
#' autogating.
#'
#' @param autogating.param The initial general parameter list provided by calling
#' get.autogating.param.minimal.
#' @return description
#' @export


get.autogating.param.symphony <- function( autogating.param ) {

  # add cytometer-specific parameters

  autogating.param$cytometer <- "Symphony"

  autogating.param$scatter.data.min.x <- 0

  autogating.param$scatter.data.max.x <- 262144

  autogating.param$scatter.data.min.y <- 0

  autogating.param$scatter.data.max.y <- 262144

  autogating.param$expr.data.min <- -111

  autogating.param$expr.data.max <- 262144

  autogating.param$default.scatter.parameter <- c( "FSC-A", "SSC-A" )

  autogating.param$default.time.parameter <- "Time"

  autogating.param$default.transformation.param <- list(
    length = 256,
    max.range = 262144,
    pos = 4.42,
    neg = 0,
    width = -20
  )

  autogating.param$data.step <- 5e4

  autogating.param$fcs.parent.regexp <- sprintf( "^%s\\.[0-9]+$",
                                                 autogating.param$fcs.parent.string )
  autogating.param$fcs.parent.regexp.sub <- sprintf( "^%s\\.([0-9]+)$",
                                                     autogating.param$fcs.parent.string )

  return( autogating.param )
}
