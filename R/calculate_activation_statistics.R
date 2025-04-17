# calculate_activation_statistics.r

calculate.activation.statistics <- function( fcs.gates.sample,
                                             fcs.activation.gates.definition, fcs.expr.data )
{
  samples <- names( fcs.gates.sample )

  activation.stats <- data.frame( Sample = samples )

  for ( activ.gate.def in fcs.activation.gates.definition )
  {
    marker.activation <- fcs.marker.activation[
      ! is.na( activ.gate.def[ fcs.marker.activation ] ) &
        activ.gate.def[ fcs.marker.activation ] != "" ]

    if ( length( marker.activation ) > 0 )
    {
      parent.gate <- activ.gate.def$parent.gate
      parent.popul <- activ.gate.def$parent.popul

      for ( marker.act in marker.activation )
      {
        gate.name <- sprintf( "%s.%s.%s", parent.gate, parent.popul,
                              marker.act )

        stats.pop.label <- fcs.gates.sample[[ 1 ]][[
          parent.gate ]][[ parent.popul ]]$label

        stats.fraction.label <- sprintf( "%s+ fraction %s %s",
                                         marker.act, fcs.stats.separator, stats.pop.label )

        stats.median.label <- sprintf( "%s+ median %s %s",
                                       marker.act, fcs.stats.separator, stats.pop.label )

        stats.iqr.label <- sprintf( "%s+ iqr %s %s",
                                    marker.act, fcs.stats.separator, stats.pop.label )

        for ( samp in samples )
        {
          population.idx <- fcs.gates.sample[[ samp ]][[
            gate.name ]][[ fcs.activation.population.label ]]$idx

          parent.idx <- fcs.gates.sample[[ samp ]][[ parent.gate ]][[
            parent.popul ]]$idx

          stopifnot( population.idx %in% parent.idx )

          marker.act.fraction <-
            mean( parent.idx %in% population.idx )

          marker.act.expr.data <- fcs.transform.inv[[ marker.act ]](
            fcs.expr.data[ population.idx, marker.act ] )

          marker.act.median <- median( marker.act.expr.data )
          marker.act.iqr <- IQR( marker.act.expr.data )

          activation.stats[ activation.stats$Sample == samp,
                            c( stats.fraction.label, stats.median.label,
                               stats.iqr.label ) ] <-
            c( marker.act.fraction, marker.act.median,
               marker.act.iqr )
        }
      }
    }
  }

  activation.stats
}
