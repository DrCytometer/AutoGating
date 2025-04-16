# read_compensation_csv.r

read.compensation.csv <- function( compensation.filepath, fcs.panel )
{
  fcs.compensation <- read.csv( compensation.filepath, stringsAsFactors = FALSE,
                                check.names = FALSE, row.names = 1 )

  fcs.compensation.dye <- fcs.panel$panel$dye[ c( fcs.panel$lineage.idx,
                                                  fcs.panel$activation.idx ) ]
  fcs.compensation.antigen <- fcs.panel$panel$antigen[ c( fcs.panel$lineage.idx,
                                                          fcs.panel$activation.idx ) ]

  stopifnot( sort( colnames( fcs.compensation ) ) == sort( fcs.compensation.dye ) )

  fcs.compensation.reorder <- match( fcs.compensation.dye, colnames( fcs.compensation ) )

  fcs.compensation <- fcs.compensation[ fcs.compensation.reorder,
                                        fcs.compensation.reorder ]

  return( fcs.compensation )
}
