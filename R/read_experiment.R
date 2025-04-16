# read_experiment.r

#' @title read.experiment
#'
#' @description
#' The function reads in the user's definition of their experiment from the csv
#' file and outputs a list with the information structured for downstream use.
#'
#' @param path The file path to the location of the fcs.experiment.filename file.
#' @param fcs.experiment.filename The csv file defining the experiment.
#' @return fcs.experiment
#' \describe{
#'  \item{experiment}{Definition of the experiment.}
#'  \item{sample}{Character vector of the fcs samples.}
#'  \item{sample.n}{Numeric. Number of samples.}
#'  \item{sample.filename}{Character vector of the fcs file names.}
#'  \item{sample.gate.calculation}{Samples to be used in calculating the gates.}
#'  \item{sample.gate.calc.filename}{Filenames of samples to be use for gate calculation.}
#'  \item{sample.gate.calculation.n}{Numeric. Number of samples for gate calculation.}
#' }
#' @export
#'
#'

read.experiment <- function( path, fcs.experiment.filename ) {

  # read definition of experiment

  fcs.experiment <- read.csv( file.path( path, fcs.experiment.filename ),
                              stringsAsFactors = FALSE )

  cat( "\033[34m", paste( "Here is the definition of your experiment:" ),
       "\033[0m\n" )
  print( fcs.experiment )

  stopifnot( anyDuplicated( fcs.experiment$filename ) == 0 )
  stopifnot( anyDuplicated( fcs.experiment$sample ) == 0 )

  fcs.sample <- fcs.experiment$sample
  fcs.sample.n <- length( fcs.sample )

  fcs.sample.filename <- fcs.experiment$filename

  fcs.sample.gate.calculation <- fcs.sample[
    which( fcs.experiment$gate.calculation ) ]
  fcs.sample.gate.calc.filename <- fcs.sample.filename[
    which( fcs.experiment$gate.calculation ) ]
  fcs.sample.gate.calculation.n <- length( fcs.sample.gate.calculation )

  cat( "\033[34m", paste( "Your experiment contains these", fcs.sample.n,
                          "samples." ), "\033[0m\n" )
  cat( fcs.sample, "\n"  )

  cat( "\033[34m", paste( "To calculate the gate definitions, these",
                          fcs.sample.gate.calculation.n,
                          "samples will be used." ), "\033[0m\n" )
  cat( fcs.sample.gate.calculation, "\n" )

  fcs.experiment <- list(
    experiment = fcs.experiment,
    sample = fcs.sample,
    sample.n = fcs.sample.n,
    sample.filename = fcs.sample.filename,
    sample.gate.calculation = fcs.sample.gate.calculation,
    sample.gate.calc.filename = fcs.sample.gate.calc.filename,
    sample.gate.calculation.n = fcs.sample.gate.calculation.n
  )

  return( fcs.experiment )
}
