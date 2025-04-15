# do_activation_gates_main.r

#' @title do.activation.gates.main
#'
#' @description
#' A function to calculate and execute gating for activation markers in the panel.
#'
#' @param fcs.gates The gates for flow cytometry data.
#' @param fcs.gates.data Data used for generating? gates.
#' @param fcs.activ.gates.definition Definitions of the activation marker gates.
#' @param calculate.gates Logical?
#' @param figure.actgate.dir Directory for output?
#' @return fcs.gates
#' @export
#'

do.activation.gates.main <- function( fcs.gates, fcs.gates.data,
                                 fcs.activ.gates.definition, calculate.gates,
                                 figure.actgate.dir )
{
  for ( activ.gate.def in fcs.activ.gates.definition )
  {
    marker.activation <- fcs.marker.activation[
      ! is.na( activ.gate.def[ fcs.marker.activation ] ) &
        activ.gate.def[ fcs.marker.activation ] != "" ]

    if ( length( marker.activation ) > 0 )
    {
      popul.number <- activ.gate.def$popul.number

      parent.gate <- activ.gate.def$parent.gate
      parent.popul <- activ.gate.def$parent.popul

      popul.label.base <- activ.gate.def$popul.label
      popul.label.pos <- activ.gate.def$popul.label.pos

      gate.algorithm.base <- activ.gate.def$gate.algorithm

      stopifnot( xor(
        activ.gate.def$gate.marker.x == fcs.gate.parameter.ignore,
        activ.gate.def$gate.marker.y == fcs.gate.parameter.ignore
      ) )

      marker.act.axis <- ifelse(
        activ.gate.def$gate.marker.x == fcs.gate.parameter.ignore,
        1, 2 )

      activ.gate.def$gate.number <- popul.number

      activ.gate.def$popul.label.pos <- paste( fcs.gate.parameter.ignore,
                                               popul.label.pos, sep = fcs.gate.parameter.delimiter )

      for ( marker.act in marker.activation )
      {
        gate.name <- sprintf( "%s.%s.%s", parent.gate, parent.popul,
                              marker.act )

        popul.label <- sprintf( "%s+ %s", marker.act,
                                popul.label.base )

        activ.gate.def$gate.name <- gate.name

        activ.gate.def$popul.name <- paste( fcs.gate.parameter.ignore,
                                            fcs.activation.population.label,
                                            sep = fcs.gate.parameter.delimiter )

        activ.gate.def$popul.label <- paste( fcs.gate.parameter.ignore,
                                             popul.label,
                                             sep = fcs.gate.parameter.delimiter )

        if ( marker.act.axis == 1 )
          activ.gate.def$gate.marker.x <- marker.act
        else
          activ.gate.def$gate.marker.y <- marker.act

        activ.gate.def$gate.algorithm <- ifelse(
          activ.gate.def[[ marker.act ]] ==
            fcs.activation.marker.selected,
          gate.algorithm.base, activ.gate.def[[ marker.act ]] )

        gate.marker <- c( activ.gate.def$gate.marker.x,
                          activ.gate.def$gate.marker.y )

        cat( sprintf( "%0*d - %s\n", fcs.gate.number.width,
                      popul.number, gate.name ) )

        fcs.gates <- do.gate.generic( fcs.gates,
                                      fcs.gates.data[ , gate.marker ], activ.gate.def,
                                      calculate.gates, figure.actgate.dir )

        str( fcs.gates[[ gate.name ]] )
        cat( "\n" )
      }
    }
  }

  fcs.gates
}
