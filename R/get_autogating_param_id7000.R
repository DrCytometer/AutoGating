# get_autogating_param_id7000.r

#' @title get.autogating.param.id7000
#'
#' @description
#' Get cytometer-specific parameters for ID7000 files in order to perform autogating.
#'
#' @param autogating.param The initial general parameter list provided by calling
#' get.autogating.param.minimal.
#' @return description
#' @export


get.autogating.param.id7000 <- function( autogating.param ) {

  # add cytometer-specific parameters

  autogating.param$cytometer <- "ID7000"

  autogating.param$scatter.data.min.x <- 0

  autogating.param$scatter.data.max.x <- 1048576

  autogating.param$scatter.data.min.y <- 0

  autogating.param$scatter.data.max.y <- 1048576

  autogating.param$expr.data.min <- -111

  autogating.param$expr.data.max <- 1000000

  autogating.param$default.scatter.parameter <- c( "FSC-A", "SSC-A" )

  autogating.param$default.time.parameter <- "TIME"

  autogating.param$default.transformation.param <- list(
    length = 256,
    max.range = 1000000,
    pos = 5,
    neg = 0,
    width = -250
  )

  autogating.param$data.step <- 1e5

  autogating.param$fcs.parent.regexp <- sprintf( "^%s\\.[0-9]+$",
                                                 autogating.param$fcs.parent.string )
  autogating.param$fcs.parent.regexp.sub <- sprintf( "^%s\\.([0-9]+)$",
                                                 autogating.param$fcs.parent.string )

  return( autogating.param )
}
