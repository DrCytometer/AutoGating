# calcalate_population_statistics.r

calculate.population.statistics <- function( fcs.gates.sample,
    fcs.population.gates.definition )
{
    samples <- names( fcs.gates.sample )

    population.stats <- data.frame( Sample = samples )

    for ( popul.gate.def in fcs.population.gates.definition )
    {
        gate.name <- popul.gate.def$gate.name

        popul.name <- strsplit( popul.gate.def$popul.name,
            fcs.gate.parameter.delimiter )[[ 1 ]]

        stats.parent.gate <- strsplit( popul.gate.def$stats.parent.gate,
            fcs.gate.parameter.delimiter )[[ 1 ]]

        stats.parent.popul <- strsplit( popul.gate.def$stats.parent.popul,
            fcs.gate.parameter.delimiter )[[ 1 ]]

        if ( length( stats.parent.gate ) == 0 &&
             length( stats.parent.popul ) == 0 )
        {
            stats.parent.gate <- paste0( fcs.parent.string, ".1" )
            stats.parent.popul <- paste0( fcs.parent.string, ".1" )
        }

        # calculate population fractions for statistics

        if ( stats.parent.gate != fcs.gate.parameter.ignore &&
             stats.parent.popul != fcs.gate.parameter.ignore )
        {
            stopifnot( length( stats.parent.gate ) ==
                    length( stats.parent.popul ) )

            for ( sp.idx in 1 : length( stats.parent.gate ) )
            {
                stats.parent.ga <- stats.parent.gate[ sp.idx ]
                stats.parent.pop <- stats.parent.popul[ sp.idx ]

                stopifnot( length( stats.parent.ga ) == 1 &&
                    length( stats.parent.pop ) == 1 )

                # obtain actual parent gate and population

                if ( grepl( fcs.parent.regexp, stats.parent.ga ) &&
                        stats.parent.ga == stats.parent.pop )
                {
                    stats.parent.n <- as.numeric(
                        sub( fcs.parent.regexp.sub, "\\1", stats.parent.ga ) )

                    stats.gate.curr <- gate.name

                    while ( stats.parent.n > 0 )
                    {
                        stats.parent.ga <- fcs.gates.sample[[
                            1 ]][[ stats.gate.curr ]][[ 1 ]]$parent.gate
                        stats.parent.pop <- fcs.gates.sample[[
                            1 ]][[ stats.gate.curr ]][[ 1 ]]$parent.popul

                        stopifnot( length( stats.parent.ga ) == 1 &&
                            length( stats.parent.pop ) == 1 )

                        stats.gate.curr <- stats.parent.ga

                        stats.parent.n <- stats.parent.n - 1
                    }
                }

                stats.parent.label <- fcs.gates.sample[[ 1 ]][[
                    stats.parent.ga ]][[ stats.parent.pop ]]$label

                for ( pn.idx in 1 : length( popul.name ) )
                {
                    pop.name <- popul.name[ pn.idx ]

                    if ( pop.name != fcs.gate.parameter.ignore )
                    {
                        stats.pop.label <- fcs.gates.sample[[ 1 ]][[
                            gate.name ]][[ pop.name ]]$label

                        stats.fraction.label <- paste0( stats.pop.label,
                            fcs.stats.separator, stats.parent.label )

                        for ( samp in samples )
                        {
                            population.idx <- fcs.gates.sample[[
                                samp ]][[ gate.name ]][[ pop.name ]]$idx
                            parent.idx <- fcs.gates.sample[[
                                samp ]][[ stats.parent.ga ]][[
                                    stats.parent.pop ]]$idx

                            stopifnot( population.idx %in% parent.idx )

                            fraction.in.parent <- mean( parent.idx %in%
                                population.idx )

                            population.stats[
                                population.stats$Sample == samp,
                                stats.fraction.label ] <- fraction.in.parent
                        }
                    }
                }
            }
        }
    }

    population.stats
}

