# calculate_gate_2dsep.r


calculate.gate.2dsep <- function( gate.data, popul.name, gate.algorithm,
    gate.param, gate.number, gate.name, fcs.population.gates, figure.dir, agp )
{
    stopifnot( length( gate.algorithm ) == 2 ||
        ( length( gate.algorithm ) == 1 &&
            gate.algorithm[ 1 ] == "replicate" ) )

    join.boundary <- gate.param$join

    gate.boundary.list <- list()

    if ( gate.algorithm[ 1 ] == "replicate" )
    {
        origin.gate <- fcs.population.gates[[ gate.param$origin.gate ]]

        stopifnot( is.null( join.boundary ) &&
            length( popul.name ) == length( origin.gate ) )

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
        stopifnot( ( is.null( join.boundary ) && length( popul.name ) == 4 ) ||
            ( ! is.null( join.boundary ) && length( popul.name ) == 1 ) )

        # get algorithm parameters

        quantile.x <- as.numeric( gate.param$quantile.x )
        quantile.y <- as.numeric( gate.param$quantile.y )

        limit.xmin <- as.numeric( gate.param$limit.xmin )
        limit.xmax <- as.numeric( gate.param$limit.xmax )

        limit.ymin <- as.numeric( gate.param$limit.ymin )
        limit.ymax <- as.numeric( gate.param$limit.ymax )

        region.xmin <- as.numeric( gate.param$region.xmin )
        region.xmax <- as.numeric( gate.param$region.xmax )

        region.ymin <- as.numeric( gate.param$region.ymin )
        region.ymax <- as.numeric( gate.param$region.ymax )

        threshold.param <- gate.param$threshold.param

        # trim data

        gate.data.trim <- gate.data[
            gate.data[ , 1 ] >= quantile( gate.data[ , 1 ], quantile.x ) &
            gate.data[ , 1 ] <= quantile( gate.data[ , 1 ], 1 - quantile.x ) &
            gate.data[ , 2 ] >= quantile( gate.data[ , 2 ], quantile.y ) &
            gate.data[ , 2 ] <= quantile( gate.data[ , 2 ], 1 - quantile.y ), ]

        # get data limits

        gate.data.trim.xmin <- min( gate.data.trim[ , 1 ] )
        gate.data.trim.xmax <- max( gate.data.trim[ , 1 ] )

        gate.data.trim.ymin <- min( gate.data.trim[ , 2 ] )
        gate.data.trim.ymax <- max( gate.data.trim[ , 2 ] )

        # extract region for gate calculation

        gate.data.calc <- gate.data.trim[
            gate.data.trim[ , 1 ] >=
                ( 1 - region.xmin ) * gate.data.trim.xmin +
                    region.xmin * gate.data.trim.xmax &
            gate.data.trim[ , 1 ] <=
                ( 1 - region.xmax ) * gate.data.trim.xmin +
                    region.xmax * gate.data.trim.xmax &
            gate.data.trim[ , 2 ] >=
                ( 1 - region.ymin ) * gate.data.trim.ymin +
                    region.ymin * gate.data.trim.ymax &
            gate.data.trim[ , 2 ] <=
                ( 1 - region.ymax ) * gate.data.trim.ymin +
                    region.ymax * gate.data.trim.ymax, ]

        # calculate thresholds

        calculate.threshold.function.x <- get(
            paste0( "calculate.threshold.", gate.algorithm[[ 1 ]] ) )

        gate.xthr <- calculate.threshold.function.x( gate.data.calc[ , 1 ],
            threshold.param )

        calculate.threshold.function.y <- get(
            paste0( "calculate.threshold.", gate.algorithm[[ 2 ]] ) )

        gate.ythr <- calculate.threshold.function.y( gate.data.calc[ , 2 ],
            threshold.param, agp )

        # plot densities and thresholds

        png( filename = file.path( figure.dir, sprintf(
                "%0*d - %s - calculation.png", agp$fcs.gate.number.width,
                gate.number, gate.name ) ),
            width = 1024, height = 1536 )
        par( mfrow = c( 2, 1 ), mar = c( 5, 5, 1, 1 ), oma = c( 0, 0, 3, 0 ) )
        plot( density( gate.data.trim[ , 1 ] ),
            xlab = colnames( gate.data.trim )[ 1 ], main = "", lwd = 3,
            cex.lab = 2, cex.axis = 2 )
        abline( v = range( gate.data.calc[ , 1 ] ), lwd = 2, lty = 2 )
        abline( v = gate.xthr, lwd = 2, col = "blue3" )
        plot( density( gate.data.trim[ , 2 ] ),
            xlab = colnames( gate.data.trim )[ 2 ], main = "", lwd = 3,
            cex.lab = 2, cex.axis = 2 )
        abline( v = range( gate.data.calc[ , 2 ] ), lwd = 2, lty = 2 )
        abline( v = gate.ythr, lwd = 2, col = "blue3" )
        mtext( sprintf( "%s  -  2dsep %s", gate.name,
                paste0( gate.algorithm, collapse = ":" ) ),
            outer = TRUE, cex = 2.5 )
        dev.off()

        # calculate boundaries

        if ( is.null( join.boundary ) )
        {
            if ( popul.name[ 1 ] != agp$fcs.gate.parameter.ignore )
            {
                # double negative population

                gate.xmin <- ( 1 - limit.xmin ) * gate.data.trim.xmin +
                    limit.xmin * gate.data.trim.xmax
                gate.xmax <- gate.xthr
                gate.ymin <- ( 1 - limit.ymin ) * gate.data.trim.ymin +
                    limit.ymin * gate.data.trim.ymax
                gate.ymax <- gate.ythr

                gate.boundary <- rbind(
                    c( gate.xmin, gate.ymin ),
                    c( gate.xmax, gate.ymin ),
                    c( gate.xmax, gate.ymax ),
                    c( gate.xmin, gate.ymax ) )

                gate.boundary.list[[ popul.name[ 1 ] ]] <- gate.boundary
            }

            if ( popul.name[ 2 ] != agp$fcs.gate.parameter.ignore )
            {
                # x-positive population

                gate.xmin <- gate.xthr
                gate.xmax <- ( 1 - limit.xmax ) * gate.data.trim.xmin +
                    limit.xmax * gate.data.trim.xmax
                gate.ymin <- ( 1 - limit.ymin ) * gate.data.trim.ymin +
                    limit.ymin * gate.data.trim.ymax
                gate.ymax <- gate.ythr

                gate.boundary <- rbind(
                    c( gate.xmin, gate.ymin ),
                    c( gate.xmax, gate.ymin ),
                    c( gate.xmax, gate.ymax ),
                    c( gate.xmin, gate.ymax ) )

                gate.boundary.list[[ popul.name[ 2 ] ]] <- gate.boundary
            }

            if ( popul.name[ 3 ] != agp$fcs.gate.parameter.ignore )
            {
                # y-positive population

                gate.xmin <- ( 1 - limit.xmin ) * gate.data.trim.xmin +
                    limit.xmin * gate.data.trim.xmax
                gate.xmax <- gate.xthr
                gate.ymin <- gate.ythr
                gate.ymax <- ( 1 - limit.ymax ) * gate.data.trim.ymin +
                    limit.ymax * gate.data.trim.ymax

                gate.boundary <- rbind(
                    c( gate.xmin, gate.ymin ),
                    c( gate.xmax, gate.ymin ),
                    c( gate.xmax, gate.ymax ),
                    c( gate.xmin, gate.ymax ) )

                gate.boundary.list[[ popul.name[ 3 ] ]] <- gate.boundary
            }

            if ( popul.name[ 4 ] != agp$fcs.gate.parameter.ignore )
            {
                # double positive population

                gate.xmin <- gate.xthr
                gate.xmax <- ( 1 - limit.xmax ) * gate.data.trim.xmin +
                    limit.xmax * gate.data.trim.xmax
                gate.ymin <- gate.ythr
                gate.ymax <- ( 1 - limit.ymax ) * gate.data.trim.ymin +
                    limit.ymax * gate.data.trim.ymax

                gate.boundary <- rbind(
                    c( gate.xmin, gate.ymin ),
                    c( gate.xmax, gate.ymin ),
                    c( gate.xmax, gate.ymax ),
                    c( gate.xmin, gate.ymax ) )

                gate.boundary.list[[ popul.name[ 4 ] ]] <- gate.boundary
            }
        }
        else
        {
            # join boundaries

            if ( join.boundary == "1:2:3" )
            {
                gate.xmin <- ( 1 - limit.xmin ) * gate.data.trim.xmin +
                    limit.xmin * gate.data.trim.xmax
                gate.xmax <- ( 1 - limit.xmax ) * gate.data.trim.xmin +
                    limit.xmax * gate.data.trim.xmax
                gate.ymin <- ( 1 - limit.ymin ) * gate.data.trim.ymin +
                    limit.ymin * gate.data.trim.ymax
                gate.ymax <- ( 1 - limit.ymax ) * gate.data.trim.ymin +
                    limit.ymax * gate.data.trim.ymax

                gate.boundary <- rbind(
                    c( gate.xmin, gate.ymin ),
                    c( gate.xmax, gate.ymin ),
                    c( gate.xmax, gate.ythr ),
                    c( gate.xthr, gate.ythr ),
                    c( gate.xthr, gate.ymax ),
                    c( gate.xmin, gate.ymax ) )
            }
            else
                stopifnot( FALSE )

            gate.boundary.list[[ popul.name[ 1 ] ]] <- gate.boundary
        }
    }

    gate.boundary.list
}

