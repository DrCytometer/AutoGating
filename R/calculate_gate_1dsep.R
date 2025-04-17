# calculate_gate_1dsep.r

calculate.gate.1dsep <- function( gate.data, popul.name, gate.algorithm,
                                  gate.param, gate.number, gate.name, fcs.population.gates, figure.dir )
{
  stopifnot( length( gate.algorithm ) == 1 )

  gate.boundary.list <- list()

  if ( gate.algorithm[ 1 ] == "replicate" )
  {
    origin.gate <- fcs.population.gates[[ gate.param$origin.gate ]]

    stopifnot( length( popul.name ) == length( origin.gate ) )

    for ( pop.idx in 1 : length( popul.name ) )
    {
      pop.name <- popul.name[ pop.idx ]
      if ( pop.name != fcs.gate.parameter.ignore )
        gate.boundary.list[[ pop.name ]] <-
          origin.gate[[ pop.idx ]]$boundary
    }
  }
  else
  {
    stopifnot( length( popul.name ) == 2 )

    # get algorithm parameters

    quantile.x <- as.numeric( gate.param$quantile.x )
    quantile.y <- as.numeric( gate.param$quantile.y )

    limit.xmin <- as.numeric( gate.param$limit.xmin )
    limit.xmax <- as.numeric( gate.param$limit.xmax )

    limit.ymin <- as.numeric( gate.param$limit.ymin )
    limit.ymax <- as.numeric( gate.param$limit.ymax )

    sep.axis <- as.numeric( gate.param$sep.axis )

    region.min <- as.numeric( gate.param$region.min )
    region.max <- as.numeric( gate.param$region.max )

    threshold.param <- gate.param$threshold.param

    # trim data

    gate.data.trim <- gate.data[
      gate.data[ , 1 ] >= quantile( gate.data[ , 1 ], quantile.x ) &
        gate.data[ , 1 ] <= quantile( gate.data[ , 1 ], 1 - quantile.x ) &
        gate.data[ , 2 ] >= quantile( gate.data[ , 2 ], quantile.y ) &
        gate.data[ , 2 ] <= quantile( gate.data[ , 2 ], 1 - quantile.y ), ]

    # work always on x axis

    if ( sep.axis == 2 )
      gate.data.trim <- gate.data.trim[ , c( 2, 1 ) ]

    # get data limits

    gate.data.trim.xmin <- min( gate.data.trim[ , 1 ] )
    gate.data.trim.xmax <- max( gate.data.trim[ , 1 ] )

    gate.data.trim.ymin <- min( gate.data.trim[ , 2 ] )
    gate.data.trim.ymax <- max( gate.data.trim[ , 2 ] )

    # extract region for gate calculation

    gate.data.calc <- gate.data.trim[
      gate.data.trim[ , 1 ] >= ( 1 - region.min ) * gate.data.trim.xmin +
        region.min * gate.data.trim.xmax &
        gate.data.trim[ , 1 ] <= ( 1 - region.max ) * gate.data.trim.xmin +
        region.max * gate.data.trim.xmax, ]

    # calculate threshold

    calculate.threshold.function.x <- get(
      paste0( "calculate.threshold.", gate.algorithm ) )

    gate.xthr <- calculate.threshold.function.x( gate.data.calc[ , 1 ],
                                                 threshold.param )

    # plot density and threshold

    png( filename = file.path( figure.dir, sprintf(
      "%0*d - %s - calculation.png", fcs.gate.number.width,
      gate.number, gate.name ) ),
      width = 1024, height = 768 )
    par( mar = c( 5, 5, 4, 1 ) )
    plot( density( gate.data.trim[ , 1 ] ),
          xlab = colnames( gate.data.trim )[ 1 ],
          main = sprintf( "%s  -  1dsep %s", gate.name, gate.algorithm ),
          lwd = 3, cex.lab = 2, cex.axis = 2, cex.main = 2.5, font.main = 1 )
    abline( v = range( gate.data.calc[ , 1 ] ), lwd = 2, lty = 2 )
    abline( v = gate.xthr, lwd = 2, col = "blue3" )
    dev.off()

    # calculate boundaries

    if ( popul.name[ 1 ] != fcs.gate.parameter.ignore )
    {
      # negative population

      gate.xmin <- ( 1 - limit.xmin ) * gate.data.trim.xmin +
        limit.xmin * gate.data.trim.xmax
      gate.xmax <- gate.xthr
      gate.ymin <- ( 1 - limit.ymin ) * gate.data.trim.ymin +
        limit.ymin * gate.data.trim.ymax
      gate.ymax <- ( 1 - limit.ymax ) * gate.data.trim.ymin +
        limit.ymax * gate.data.trim.ymax

      gate.boundary <- rbind(
        c( gate.xmin, gate.ymin ),
        c( gate.xmax, gate.ymin ),
        c( gate.xmax, gate.ymax ),
        c( gate.xmin, gate.ymax ) )

      if ( sep.axis == 2 )
        gate.boundary <- gate.boundary[ , c( 2, 1 ) ]

      gate.boundary.list[[ popul.name[ 1 ] ]] <- gate.boundary
    }

    if ( popul.name[ 2 ] != fcs.gate.parameter.ignore )
    {
      # positive population

      gate.xmin <- gate.xthr
      gate.xmax <- ( 1 - limit.xmax ) * gate.data.trim.xmin +
        limit.xmax * gate.data.trim.xmax
      gate.ymin <- ( 1 - limit.ymin ) * gate.data.trim.ymin +
        limit.ymin * gate.data.trim.ymax
      gate.ymax <- ( 1 - limit.ymax ) * gate.data.trim.ymin +
        limit.ymax * gate.data.trim.ymax

      gate.boundary <- rbind(
        c( gate.xmin, gate.ymin ),
        c( gate.xmax, gate.ymin ),
        c( gate.xmax, gate.ymax ),
        c( gate.xmin, gate.ymax ) )

      if ( sep.axis == 2 )
        gate.boundary <- gate.boundary[ , c( 2, 1 ) ]

      gate.boundary.list[[ popul.name[ 2 ] ]] <- gate.boundary
    }
  }

  gate.boundary.list
}
