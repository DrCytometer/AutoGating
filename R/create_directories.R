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

# move these to agp
# "_popgate_all"
# "%s_popgate_%s"
# "all"
# "_actgate_all"
#  "%s_actgate_%s"

# move this so it gets called automatically within something else

create.directories <- function( agp ) {

  fcs.figure.popgate.all.dir <- paste0( agp$fcs.figure.dir.basename, "_popgate_all" )
  fcs.figure.popgate.sample.dir <- sprintf( "%s_popgate_%s",
                                            agp$fcs.figure.dir.basename, fcs.sample )

  fcs.figure.popgate.dir <- c( fcs.figure.popgate.all.dir,
                               fcs.figure.popgate.sample.dir )
  names( fcs.figure.popgate.dir ) <- c( "all", fcs.sample )

  fcs.figure.actgate.all.dir <- paste0( fcs.figure.dir.basename, "_actgate_all" )
  fcs.figure.actgate.sample.dir <- sprintf( "%s_actgate_%s",
                                            fcs.figure.dir.basename, fcs.sample )

  fcs.figure.actgate.dir <- c( fcs.figure.actgate.all.dir,
                               fcs.figure.actgate.sample.dir )
  names( fcs.figure.actgate.dir ) <- c( "all", fcs.sample )

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
