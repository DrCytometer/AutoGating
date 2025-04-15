# do.activation.gates.r

#' @title do.activation.gates
#'
#' @description
#' A function to calculate and execute gating for activation markers in the panel.
#'
#' @param fcs.expr.data The marker expression data from the fcs files.
#' @param fcs.population.gates List of population-defining (cell type) gates.
#' @param fcs.figure.dir Directory where the figures will be generated.
#' @param fcs.activation.thresholds The thresholds used to define the positive/negative
#' cutoffs for the activation marker gates.
#' @param fcs.activation.gates.calculate Logical?
#' @param fcs.activation.gates.param Parameters for calculating the activation gates.
#' @return description
#' @export
#'
#'

do.activation.gates <- function( fcs.expr.data, fcs.population.gates,
                                 fcs.figure.dir, fcs.activation.thresholds, fcs.activation.gates,
                                 fcs.activation.gates.calculate, fcs.activation.gates.param )
{
  for ( fcs.activation.marker in names( fcs.activation.gates ) )
    for ( fcs.activation.population in
          names( fcs.activation.gates[[ fcs.activation.marker ]] ) )
      fcs.activation.gates <- do.act.gate( fcs.expr.data,
                                           fcs.population.gates, fcs.figure.dir,
                                           fcs.activation.thresholds[[ fcs.activation.marker ]]$threshold,
                                           fcs.activation.marker, fcs.activation.population,
                                           fcs.activation.gates,
                                           fcs.activation.gates.calculate[ sprintf( "%s.%s",
                                                                                    fcs.activation.marker, fcs.activation.population ) ],
                                           fcs.activation.gates.param[[ fcs.activation.marker ]][[
                                             fcs.activation.population ]] )

  fcs.activation.gates
}
