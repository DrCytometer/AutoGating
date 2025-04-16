# get_autogating_param.r

get.autogating.param <- function( cytometer )
{
  autogating.param <- get.autogating.param.minimal()

  # cytometer-specific parameters

  get.param.function <- get0( sprintf( "get.autogating.param.%s", cytometer ) )

  if ( is.null( get.param.function ) ) {
    cat( "\033[31m", "Unsupported cytometer", "\033[0m\n" )
    stop()
  }

  autogating.param <- get.param.function( autogating.param )

  autogating.param

}

