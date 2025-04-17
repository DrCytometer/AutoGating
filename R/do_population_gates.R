# do.population.gates.r

#' @title do.population.gates
#'
#' @description
#' A function to calculate and execute gating for population (cell type) gates
#' in the panel.
#'
#' @param fcs.gates The gates to be generated.
#' @param fcs.gates.data Data to be used to generate the gates.
#' @param fcs.popul.gates.definition Definitions of the gates, obtained from the
#' gate parameter table.
#' @param calculate.gates Logical?
#' @param figure.popgate.dir The directory where the figures will be generated.
#' @return fcs.gates
#' @export
#'

do.population.gates <- function( fcs.gates, fcs.gates.data,
                                 fcs.popul.gates.definition, calculate.gates, figure.popgate.dir )
{
  for ( popul.gate.def in fcs.popul.gates.definition )
  {
    gate.number <- popul.gate.def$gate.number
    gate.name <- popul.gate.def$gate.name
    gate.marker <- c( popul.gate.def$gate.marker.x,
                      popul.gate.def$gate.marker.y )

    cat( sprintf( "%0*d - %s\n", fcs.gate.number.width, gate.number,
                  gate.name ) )

    fcs.gates <- do.gate.generic( fcs.gates,
                                  fcs.gates.data[ , gate.marker ], popul.gate.def, calculate.gates,
                                  figure.popgate.dir )

    str( fcs.gates[[ gate.name ]] )
    cat( "\n" )
  }

  fcs.gates
}
