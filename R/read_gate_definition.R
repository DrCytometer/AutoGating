# read_gate_definition.r

#' @title read.gate.definition
#'

read.gate.definition <- function( path, pop.gate.def.filename,
                                  act.gate.def.filename, agp ) {

  # read definition of population gates
  fcs.population.gates.definition.table <- read.csv( file.path( new.param.dir,
                                                                fcs.population.gates.definition.filename ), stringsAsFactors = FALSE )

  fcs.population.gates.definition <- lapply(
    1 : nrow( fcs.population.gates.definition.table ),
    function( pgd.idx ) {
      pop.gate.def <- unclass( fcs.population.gates.definition.table[
        pgd.idx, ] )
      attr( pop.gate.def, "row.names" ) <- NULL
      pop.gate.def
    } )

  # read definition of activation gates

  fcs.activation.gates.definition.table <- read.csv( file.path( param.dir,
                                                                fcs.activation.gates.definition.filename ), stringsAsFactors = FALSE )

  # correct names for activation markers changed with table reading

  fcs.marker.activation.idx <- which(
    names( fcs.activation.gates.definition.table ) %in%
      make.names( fcs.marker ) )

  fcs.marker.activation.raw <- names( fcs.activation.gates.definition.table )[
    fcs.marker.activation.idx ]

  fcs.marker.activation <- fcs.marker[ match( fcs.marker.activation.raw,
                                              make.names( fcs.marker ) ) ]

  fcs.marker.activation.raw
  fcs.marker.activation

  names( fcs.activation.gates.definition.table )[ fcs.marker.activation.idx ] <-
    fcs.marker.activation

  fcs.activation.gates.definition <- lapply(
    1 : nrow( fcs.activation.gates.definition.table ),
    function( agd.idx ) {
      act.gate.def <- unclass( fcs.activation.gates.definition.table[
        agd.idx, ] )
      attr( act.gate.def, "row.names" ) <- NULL
      act.gate.def
    } )

  return( list(
    pop.gates = fcs.population.gates.definition,
    act.gates = fcs.activation.gates.definition ) )
}
