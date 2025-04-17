# get_autogating_param.r

get.autogating.param <- function( cytometer, fcs.experiment )
{
  autogating.param <- get.autogating.param.minimal()

  # cytometer-specific parameters

  get.param.function <- get0( sprintf( "get.autogating.param.%s", cytometer ) )

  if ( is.null( get.param.function ) ) {
    cat( "\033[31m", "Unsupported cytometer", "\033[0m\n" )
    stop()
  }

  autogating.param <- get.param.function( autogating.param )

  directories <- create.directories( autogating.param, fcs.experiment )

  autogating.param$fcs.figure.popgate.dir <- directories$pop.dir
  autogating.param$fcs.figure.actgate.dir <- directories$act.dir

  return( autogating.param )

}

