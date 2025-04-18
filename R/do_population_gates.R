# do.population.gates.r

#' @title do.population.gates
#'
#' @description
#' A function to calculate and execute gating for population (cell type) gates
#' in the panel.
#'
#' @param fcs.gates The gates to be generated.
#' @param fcs.gates.data Data to be used to generate the gates.
#' @param fcs.pop.gates.def Definitions of the gates, obtained from the
#' gate parameter table.
#' @param calculate.gates Logical?
#' @param figure.popgate.dir The directory where the figures will be generated.
#' @param agp The AutoGating parameter list.
#' @param fcs.transform The biexponential transformation list.
#' @param fcs.panel Description of the cytometry panel. Produced by read.panel.design().
#' @return fcs.gates
#' @export
#'

do.population.gates <- function( fcs.gates, fcs.gates.data, fcs.pop.gates.def,
                                 calculate.gates, figure.popgate.dir, agp,
                                 fcs.transform, fcs.panel )
{
  for ( popul.gate.def in fcs.pop.gates.def )
  {
    gate.number <- popul.gate.def$gate.number

    gate.name <- popul.gate.def$gate.name

    gate.marker <- c( popul.gate.def$gate.marker.x,
                      popul.gate.def$gate.marker.y )

    cat( sprintf( "%0*d - %s\n", agp$fcs.gate.number.width, gate.number,
                  gate.name ) )

    fcs.gates <- do.gate.generic( fcs.gates,
                                  fcs.gates.data[ , gate.marker ], popul.gate.def,
                                  calculate.gates, figure.popgate.dir, agp,
                                  fcs.transform, fcs.panel )

    print( str( fcs.gates[[ gate.name ]] ) )
    cat( "\n" )
  }

  fcs.gates
}
