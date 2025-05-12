# autogate.r

#' @title autogate

autogate <- function( param.dir, data.dir, transform.filename,
                      fcs.experiment, fcs.panel, agp,
                      fcs.gates.all,
                      compensate = FALSE, compensation.filename,
                      setup = FALSE ) {

  # set up parallel processing

  if( agp$parallel ){
    plan( multisession, workers = agp$worker.process.n )
    options( future.globals.maxSize = agp$max.memory.n )
    lapply.function <- future_lapply
  } else {
    lapply.function <- lapply.sequential
  }

  fcs.sample <- fcs.experiment$sample.filename
  sample.n <- fcs.experiment$sample.n

  # read in fcs files

  flow.data <- read.fcs.data( param.dir, data.dir, transform.filename,
                              fcs.experiment, fcs.panel, agp,
                              compensate, compensation.filename,
                              setup )

  # apply population and activation gates to each individual sample

  fcs.gates.sample <- lapply.function( fcs.sample, function( samp ) {
    fcs.gates.samp <- fcs.gates.all

    fcs.gates.samp.data <- fcs.expr.data[ fcs.event.sample == samp, ]

    fcs.gates.samp <- do.population.gates( fcs.gates.samp, fcs.gates.samp.data,
                                           fcs.population.gates.definition, FALSE,
                                           fcs.figure.popgate.dir[ samp ] )

    fcs.gates.samp <- do.activation.gates( fcs.gates.samp, fcs.gates.samp.data,
                                           fcs.activation.gates.definition, FALSE,
                                           fcs.figure.actgate.dir[ samp ] )

    fcs.gates.samp
  } )

  names( fcs.gates.sample ) <- fcs.sample

  cat( str( fcs.gates.sample, max.level = 1 ) )


  # calculate population statistics

  population.stats <- calculate.population.statistics( fcs.gates.sample,
                                                       fcs.population.gates.definition )

  activation.stats <- calculate.activation.statistics( fcs.gates.sample,
                                                       fcs.activation.gates.definition,
                                                       fcs.expr.data )

  cat( str( population.stats ) )
  cat( str( activation.stats ) )

  write.table( population.stats,
               file = file.path( agp$fcs.statistics.dir, fcs.population.statistics.filename ),
               sep = ",", quote = FALSE, row.names = FALSE )

  write.table( activation.stats,
               file = file.path( agp$fcs.statistics.dir, fcs.activation.statistics.filename ),
               sep = ",", quote = FALSE, row.names = FALSE )


  # calculate threshold for activation gates on all samples together

  fcs.activation.thresholds <- do.activation.thresholds(
    fcs.expr.data[ fcs.event.sample %in% fcs.gate.definition.sample, ],
    fcs.population.gates, fcs.figure.dir[ "all" ], fcs.activation.thresholds,
    fcs.activation.thresholds.param )


  # calculate activation gates on all samples together

  fcs.activation.gates.calculate.name <- as.vector( unlist(
    sapply( names( fcs.activation.gates ), function( act.marker ) {
      act.popul <- names( fcs.activation.gates[[ act.marker ]] )
      sprintf( "%s.%s", rep( act.marker, length( act.popul ) ), act.popul )
    } ) ) )

  fcs.activation.gates.calculate <- rep( TRUE,
                                         length.out = length( fcs.activation.gates.calculate.name ) )
  names( fcs.activation.gates.calculate ) <- fcs.activation.gates.calculate.name

  fcs.activation.gates <- do.activation.gates(
    fcs.expr.data[ fcs.event.sample %in% fcs.gate.definition.sample, ],
    fcs.population.gates, fcs.figure.dir[ "all" ], fcs.activation.thresholds,
    fcs.activation.gates, fcs.activation.gates.calculate,
    fcs.activation.gates.param )


  # apply activation gates on each individual sample

  fcs.activation.gates.sample <- list()
  for ( samp in fcs.sample )
    fcs.activation.gates.sample[[ samp ]] <- fcs.activation.gates

  fcs.activation.gates.calculate.sample <- list()
  for ( samp in fcs.sample ) {
    fcs.activation.gates.calculate.sample[[ samp ]] <-
      rep( FALSE, length.out = length( fcs.activation.gates.calculate ) )
    names( fcs.activation.gates.calculate.sample[[ samp ]] ) <-
      names( fcs.activation.gates.calculate )
  }

  for ( samp in fcs.sample )
    fcs.activation.gates.sample[[ samp ]] <- do.activation.gates(
      fcs.expr.data[ fcs.event.sample == samp, ],
      fcs.population.gates.sample[[ samp ]],
      fcs.figure.dir[ samp ],
      fcs.activation.thresholds, fcs.activation.gates.sample[[ samp ]],
      fcs.activation.gates.calculate.sample[[ samp ]],
      fcs.activation.gates.param )


  # calculate population statistics

  fcs.population.statistics.sample <- list()

  for ( samp in fcs.sample )
    fcs.population.statistics.sample[[ samp ]] <- do.population.statistics(
      fcs.expr.data[ fcs.event.sample == samp, ],
      fcs.population.gates.sample[[ samp ]],
      fcs.activation.gates.sample[[ samp ]],
      fcs.activation.gates.param,
      fcs.population.statistics.param,
      fcs.statistics.dir )

  fcs.batch.statistics <- do.batch.statistics( fcs.population.statistics.sample,
                                               fcs.batch.statistics.param,
                                               fcs.statistics.dir )

}
