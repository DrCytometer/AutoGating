# calculate_gate_free.r

calculate.gate.free <- function( gate.data, popul.name, gate.algorithm,
                                 gate.param, gate.number, gate.name,
                                 fcs.population.gates, figure.dir, agp )
{
  stopifnot( length( gate.algorithm ) == 1 && length( popul.name ) == 1 )

  gate.boundary.list <- list()

  if ( gate.algorithm[ 1 ] == "replicate" )
  {
    origin.gate <- fcs.population.gates[[ gate.param$origin.gate ]]

    stopifnot( length( origin.gate ) == 1 )

    gate.boundary.list[[ popul.name[ 1 ] ]] <- origin.gate[[ 1 ]]$boundary
  }
  else
  {
    # get algorithm parameters

    quantile.x <- as.numeric( gate.param$quantile.x )
    quantile.y <- as.numeric( gate.param$quantile.y )

    boundary.param <- gate.param$boundary.param

    # trim data

    gate.data.trim <- gate.data[
      gate.data[ , 1 ] >= quantile( gate.data[ , 1 ], quantile.x ) &
        gate.data[ , 1 ] <= quantile( gate.data[ , 1 ], 1 - quantile.x ) &
        gate.data[ , 2 ] >= quantile( gate.data[ , 2 ], quantile.y ) &
        gate.data[ , 2 ] <= quantile( gate.data[ , 2 ], 1 - quantile.y ), ]

    # calculate boundary

    calculate.boundary.function <- get(
      paste0( "calculate.boundary.", gate.algorithm ) )

    gate.boundary <- calculate.boundary.function( gate.data.trim,
                                                  boundary.param )

    gate.boundary.list[[ popul.name[ 1 ] ]] <- gate.boundary
  }

  gate.boundary.list
}
