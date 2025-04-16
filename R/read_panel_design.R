# read_panel_design.r

#' @title read.panel.design
#'
#' @description
#' Reads in the panel design csv file.
#'
#' @param param.dir File path to the fcs.panel.filename
#' @param fcs.panel.filename CSV file describing the panel structure.
#' @return fcs.panel
#' @export
#'
#'


read.panel.design <- function( param.dir, fcs.panel.filename ) {

  # read panel design

  fcs.panel <- read.csv( file.path( param.dir, fcs.panel.filename ),
                         stringsAsFactors = FALSE )

  fcs.panel$type <- factor( fcs.panel$type, levels = c( "scatter", "time",
                                                        "lineage", "activation" ) )

  # define markers

  fcs.panel.scatter.idx <- which( fcs.panel$type == "scatter" )
  fcs.panel.time.idx <- which( fcs.panel$type == "time" )
  fcs.panel.lineage.idx <- which( fcs.panel$type == "lineage" )
  fcs.panel.activation.idx <- which( fcs.panel$type == "activation" )

  fcs.panel.marker.idx <- c( fcs.panel.scatter.idx, fcs.panel.time.idx,
                             fcs.panel.lineage.idx, fcs.panel.activation.idx )

  stopifnot( sort( fcs.panel.marker.idx ) == 1 : length( fcs.panel$dye ) )

  fcs.marker.scatter <- fcs.panel$antigen[ fcs.panel.scatter.idx ]
  fcs.marker.time <- fcs.panel$antigen[ fcs.panel.time.idx ]
  fcs.marker.lineage <- fcs.panel$antigen[ fcs.panel.lineage.idx ]
  fcs.marker.activation <- fcs.panel$antigen[ fcs.panel.activation.idx ]

  fcs.marker.scatter.label <- fcs.marker.scatter
  fcs.marker.time.label <- fcs.marker.time
  fcs.marker.lineage.label <- paste0(
    fcs.panel$dye[ fcs.panel.lineage.idx ], " - ",
    fcs.panel$antigen[ fcs.panel.lineage.idx ] )
  fcs.marker.activation.label <- paste0(
    fcs.panel$dye[ fcs.panel.activation.idx ], " - ",
    fcs.panel$antigen[ fcs.panel.activation.idx ] )

  fcs.marker <- c( fcs.marker.scatter, fcs.marker.time, fcs.marker.lineage,
                   fcs.marker.activation )

  stopifnot( fcs.marker == fcs.panel$antigen[ fcs.panel.marker.idx ] )

  fcs.marker.label <- c( fcs.marker.scatter.label, fcs.marker.time.label,
                         fcs.marker.lineage.label, fcs.marker.activation.label )

  fcs.marker.label.figure <- fcs.marker.label
  names( fcs.marker.label.figure ) <- fcs.marker

  fcs.marker.label.table <- fcs.marker
  names( fcs.marker.label.table ) <- fcs.marker

  cat("\033[34m", "Markers", "\033[0m\n")
  cat( paste( fcs.marker, "\n" ) )
  cat("\033[34m", "Markers & fluorophores", "\033[0m\n")
  cat( paste( fcs.marker.label, "\n" ) )

  fcs.panel <- list(
    panel = fcs.panel,
    scatter = fcs.marker.scatter,
    scatter.idx = fcs.panel.scatter.idx,
    scatter.label = fcs.marker.scatter.label,
    time = fcs.marker.time,
    time.idx = fcs.panel.time.idx,
    time.label = fcs.marker.time.label,
    activation = fcs.marker.activation,
    activation.idx = fcs.panel.activation.idx,
    activation.label = fcs.marker.activation.label,
    lineage = fcs.marker.lineage,
    lineage.idx = fcs.panel.lineage.idx,
    lineage.label = fcs.marker.lineage.label,
    marker = fcs.marker,
    marker.idx = fcs.panel.marker.idx,
    marker.label = fcs.marker.label,
    marker.label.figure = fcs.marker.label.figure,
    marker.label.table = fcs.marker.label.table
  )

  return( fcs.panel )
}
