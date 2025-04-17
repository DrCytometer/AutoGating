# do_activation_thresholds.r

#' @title do.activation.thresholds
#'
#' @description
#' A function to calculate the thresholds for gating activation markers in the panel.
#'
#' @param fcs.expr.data The marker expression data from the fcs files.
#' @param fcs.population.gates List of population-defining (cell type) gates.
#' @param fcs.activation.thresholds The thresholds used to define the positive/negative
#' cutoffs for the activation marker gates.
#' @param fcs.activation.thresholds.param Parameters for calculating the
#' activation gate thresholds.
#' @return description
#' @export
#'

do.activation.thresholds <- function( fcs.expr.data,
                                      fcs.population.gates, fcs.figure.dir, fcs.activation.thresholds,
                                      fcs.activation.thresholds.param )
{
  for ( activation.marker in names( fcs.activation.thresholds ) )
    if ( fcs.activation.thresholds[[ activation.marker ]] == 0 ) {
      # general threshold function
      fcs.activation.thresholds <- do.act.threshold(
        fcs.expr.data, fcs.population.gates, fcs.figure.dir,
        fcs.activation.thresholds, activation.marker,
        fcs.activation.thresholds.param[[ activation.marker ]] )
    }
  else {
    # special threshold function
    threshold.function <-
      get( paste0( "do.act.threshold.", activation.marker ) )
    fcs.activation.thresholds <- threshold.function( fcs.expr.data,
                                                     fcs.population.gates, fcs.figure.dir, fcs.activation.thresholds )

  }

  fcs.activation.thresholds
}
