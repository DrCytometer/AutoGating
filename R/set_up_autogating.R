# set_up_autogating.r

set.up.autogating <- function( param.dir, data.dir, transform.filename,
                               fcs.experiment, fcs.panel, agp,
                               fcs.pop.gates.def.filename, fcs.act.gates.def.filename,
                               compensate = FALSE, compensation.filename,
                               setup = TRUE ) {
  # read transform

  fcs.transform.both <- read.transformation( file.path( param.dir,
                                                        transform.filename ),
                                             fcs.panel )

  fcs.transform <- fcs.transform.both[[ 1 ]]

  # read in fcs files

  fcs.data <- read.fcs.data( param.dir, data.dir, fcs.transform.both,
                              fcs.experiment, fcs.panel, agp,
                              compensate, compensation.filename,
                              setup )

  # more fixes for lazy evaluation preventing reading/writing of flowSet

  name.update.idx <- c( fcs.panel$lineage.idx, fcs.panel$activation.idx )
  print( names( fcs.transform ))
  names( fcs.transform )[ 1:length( name.update.idx ) ] <- fcs.panel$panel$antigen[ name.update.idx ]
  print( names( fcs.transform ))

  # read definition of population gates

  cat( "\033[34m", "Reading gate definitions...", "\033[0m\n" )

  fcs.pop.gates.def.table <- read.csv( file.path( param.dir,
                                            fcs.pop.gates.def.filename ),
                                 stringsAsFactors = FALSE,
                                 check.names = FALSE )

  fcs.pop.gates.def <- lapply( 1 : nrow( fcs.pop.gates.def.table ),
    function( pgd.idx ) {
      pop.gate.def <- unclass( fcs.pop.gates.def.table[
        pgd.idx, ] )
      attr( pop.gate.def, "row.names" ) <- NULL
      pop.gate.def
    } )

  # read definition of activation gates

  fcs.act.gates.def.table <- read.csv( file.path( param.dir,
                                            fcs.act.gates.def.filename ),
                                 stringsAsFactors = FALSE,
                                 check.names = FALSE,
                                 row.names = 1 )

  fcs.act.gates.def <- lapply( 1 : nrow( fcs.act.gates.def.table ),
    function( agd.idx ) {
      act.gate.def <- unclass( fcs.act.gates.def.table[
        agd.idx, ] )
      attr( act.gate.def, "row.names" ) <- NULL
      act.gate.def
    } )

  # calculate population and activation gates on joined selected samples

  cat( "\033[34m", "Calculating optimal gate positions...", "\033[0m\n" )

  fcs.gates.all <- list()

  fcs.gates.all <- do.population.gates( fcs.gates.all, fcs.data$expr.data,
                                        fcs.pop.gates.def, TRUE,
                                        agp$fcs.figure.popgate.dir[ "all" ],
                                        agp, fcs.transform, fcs.panel )

  fcs.gates.all <- do.activation.gates.def( fcs.gates.all, fcs.data$expr.data,
                                        fcs.act.gates.def, TRUE,
                                        agp$fcs.figure.actgate.dir[ "all" ],
                                        agp, fcs.transform, fcs.panel )

  str( fcs.gates.all, max.level = 1 )

  return( fcs.gates.all )

}
