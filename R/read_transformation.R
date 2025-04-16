# read_transformation.r

#' @title read.transformation
#'
#' @description
#' Read in a csv file defining the biexponential transformation parameters.
#'
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom flowCore linearTransform
#'
#' @param transformation.filepath Path (including filename) to the location of
#' the file defining the transformation.
#' @param fcs.panel Description of the cytometry panel. Produced by read.panel.design().
#'
#' @return The transformation and inverse transformation lists for the channels.
#'
#' @export
#'

read.transformation <- function( transformation.filepath, fcs.panel ) {
  fcs.transform.param <- read.csv( transformation.filepath,
                                   stringsAsFactors = FALSE )

  fcs.transform.dye <- fcs.panel$dye[ c( fcs.panel$lineage.idx,
                                         fcs.panel$activation.idx ) ]

  fcs.transform.antigen <- fcs.panel$antigen[ c( fcs.panel$lineage.idx,
                                                 fcs.panel$activation.idx ) ]

  stopifnot( sort( fcs.transform.param$dye ) == sort( fcs.transform.dye ) )

  fcs.transform <- lapply( 1 : nrow( fcs.transform.param ),
                           function( ftp.idx ) flowWorkspace::flowjo_biexp(
                             channelRange = fcs.transform.param$length[ ftp.idx ],
                             maxValue = fcs.transform.param$max.range[ ftp.idx ],
                             pos = fcs.transform.param$pos[ ftp.idx ],
                             neg = fcs.transform.param$neg[ ftp.idx ],
                             widthBasis = fcs.transform.param$width[ ftp.idx ],
                             inverse = FALSE
                           ) )

  fcs.transform.inv <- lapply( 1 : nrow( fcs.transform.param ),
                               function( ftp.idx ) flowWorkspace::flowjo_biexp(
                                 channelRange = fcs.transform.param$length[ ftp.idx ],
                                 maxValue = fcs.transform.param$max.range[ ftp.idx ],
                                 pos = fcs.transform.param$pos[ ftp.idx ],
                                 neg = fcs.transform.param$neg[ ftp.idx ],
                                 widthBasis = fcs.transform.param$width[ ftp.idx ],
                                 inverse = TRUE
                               ) )

  fcs.transform.reorder <- match( fcs.transform.dye, fcs.transform.param$dye )

  fcs.transform <- fcs.transform[ fcs.transform.reorder ]
  fcs.transform.inv <- fcs.transform.inv[ fcs.transform.reorder ]

  names( fcs.transform ) <- fcs.transform.antigen
  names( fcs.transform.inv ) <- fcs.transform.antigen

  for ( fms in fcs.panel$scatter )
  {
    fcs.transform[[ fms ]] <- flowCore::linearTransform()
    fcs.transform.inv[[ fms ]] <- flowCore::linearTransform()
  }

  return( list( fcs.transform, fcs.transform.inv ) )
}
