# get_autogating_param_aurora.r

#' @title get.autogating.param.aurora
#'
#' @description
#' Get cytometer-specific parameters for Aurora files in order to perform autogating.
#'
#' @param autogating.param The initial general parameter list provided by calling
#' get.autogating.param.minimal.
#' @return autogating.param Returns updated autogating parameters
#' @export
#'


get.autogating.param.aurora <- function( autogating.param ) {

  # add cytometer-specific parameters

  autogating.param$cytometer <- "Aurora"

  autogating.param$scatter.data.min.x <- 0

  autogating.param$scatter.data.max.x <- 4194304

  autogating.param$scatter.data.min.y <- 0

  autogating.param$scatter.data.max.y <- 4194304

  autogating.param$expr.data.min <- -111

  autogating.param$expr.data.max <- 4194304

  autogating.param$default.scatter.parameter <- c( "FSC-A", "SSC-A" )

  autogating.param$default.time.parameter <- "Time"

  autogating.param$default.transformation.param <- list(
    length = 256,
    max.range = 4194304,
    pos = 5.62,
    neg = 0,
    width = -1000
  )

  autogating.param$data.step <- 5e5

  autogating.param$fcs.parent.regexp <- sprintf( "^%s\\.[0-9]+$",
                                                 autogating.param$fcs.parent.string )
  autogating.param$fcs.parent.regexp.sub <- sprintf( "^%s\\.([0-9]+)$",
                                                     autogating.param$fcs.parent.string )

  return( autogating.param )
}
