# create_directories.r

#' @title create.directories
#'
#' @description
#' Creates the output folders for the pipeline.
#'
#' @param agp The AutoGating parameter list.
#' @export
#'
#'


# move this so it gets called automatically within something else

create.directories <- function( agp, fcs.experiment ) {

  fcs.figure.popgate.all.dir <- paste0( agp$fcs.figure.dir.basename, agp$pop.gate.all.label )
  fcs.figure.popgate.sample.dir <- sprintf( agp$pop.gate.label,
                                            agp$fcs.figure.dir.basename, fcs.experiment$sample )

  fcs.figure.popgate.dir <- c( fcs.figure.popgate.all.dir,
                               fcs.figure.popgate.sample.dir )
  names( fcs.figure.popgate.dir ) <- c( "all", fcs.experiment$sample )

  fcs.figure.actgate.all.dir <- paste0( agp$fcs.figure.dir.basename, agp$act.gate.all.label )
  fcs.figure.actgate.sample.dir <- sprintf( agp$act.gate.label,
                                            agp$fcs.figure.dir.basename, fcs.experiment$sample )

  fcs.figure.actgate.dir <- c( fcs.figure.actgate.all.dir,
                               fcs.figure.actgate.sample.dir )
  names( fcs.figure.actgate.dir ) <- c( "all", fcs.experiment$sample )

  for ( fig.dir in c( fcs.figure.popgate.dir, fcs.figure.actgate.dir ) ){
    #full.path <- file.path( output.dir, fig.dir )
    if ( !file.exists( fig.dir ) ) {
      dir.create( fig.dir, recursive = TRUE )
    }
  }


  # create statistics directory

  if ( ! file.exists( agp$fcs.statistics.dir ) ){
    #full.path <- file.path( output.dir, fcs.statistics.dir )
    dir.create( agp$fcs.statistics.dir, recursive = TRUE )
  }

}
