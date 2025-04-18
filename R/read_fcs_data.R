# read_fcs_data.r

#' @title read.fcs.data
#'
#' @description
#' Read in the fcs files for autogating.
#'
#' @importFrom flowCore read.flowSet compensate compensation transformList exprs
#' @importFrom flowCore read.FCS fsApply
#'
#' @param param.dir File path to the transformation and compensation files.
#' @param data.dir File path to the fcs files.
#' @param fcs.transform.both The biexponential and inverse transformation
#' to be applied.
#' @param fcs.experiment Description of the cytometry experiment. Produced by
#' read.experiment().
#' @param fcs.panel Description of the cytometry panel. Produced by
#' read.panel.design().
#' @param agp AutoGating parameter list. Produced by get.autogating.param().
#' @param compensate Logical. If TRUE, compensation will be applied and a csv
#' file must be provided to compensation.filename.
#' @param compensation.filename A csv file providing the compensation matrix to
#' be applied. Required if compensate is set to TRUE.
#' @param setup Logical. If TRUE, uses only the files selected for defining
#' (calculating) the gates. If FALSE, processes all files.
#'
#' @return fcs.data The expression (time, scatter and marker) data from the fcs files.
#'
#' @export
#'


read.fcs.data <- function( param.dir, data.dir, fcs.transform.both,
                           fcs.experiment, fcs.panel, agp,
                           compensate = FALSE, compensation.filename,
                           setup = TRUE ) {

  # select samples to be used

  if ( setup ) {
    fcs.sample <- fcs.experiment$sample.gate.calc.filename
    sample.name <- fcs.experiment$sample.gate.calculation
    sample.n <- fcs.experiment$sample.gate.calculation.n
  } else {
    fcs.sample <- fcs.experiment$sample.filename
    sample.name <- fcs.experiment$sample
    sample.n <- fcs.experiment$sample.n
  }

  # read compensation and transformation

  if ( compensate )
    fcs.compensation <- read.compensation.csv( file.path( param.dir,
                                                          compensation.filename ),
                                               fcs.panel )

  fcs.transform <- fcs.transform.both[[ 1 ]]
  fcs.transform.inv <- fcs.transform.both[[ 2 ]]

  # read fcs files and apply compensation and transformation

  cat( "\033[34m", "Reading fcs files...", "\033[0m\n" )

  # suppressWarnings to shut down truncate_max_range warnings

  fcs.flow.set <- suppressWarnings( read.flowSet( files = fcs.sample,
                                path = data.dir, transformation = FALSE,
                                truncate_max_range = TRUE ) )


  fs.colnames <- colnames( exprs( fcs.flow.set[[ 1 ]] ) )

  stopifnot( fcs.panel$panel$dye %in% fs.colnames )

  fcs.flow.set <- fcs.flow.set[ , fcs.panel$panel$dye ]

  # compensate
  if ( compensate )
    fcs.flow.set <- flowCore::compensate( fcs.flow.set,
                                          flowCore::compensation( fcs.compensation ) )

  # transform
  transformed.list <- flowCore::fsApply( fcs.flow.set, function( ff ) {
    flowCore::transform( ff, flowCore::transformList( names( fcs.transform ), fcs.transform ) )
  }, simplify = FALSE)

  fcs.flow.set <- as( transformed.list, "flowSet" )

  # extract expression data

  fcs.sample.event.number.max <- 0

  for ( flow.idx in 1 : sample.n )
  {
    flow.sample.data <- fcs.flow.set[[ flow.idx ]]

    fcs.sample.event.number <- nrow( exprs( flow.sample.data ) )

    if ( fcs.sample.event.number > fcs.sample.event.number.max )
      fcs.sample.event.number.max <- fcs.sample.event.number
  }

  fcs.event.number.width <- floor( log10( fcs.sample.event.number.max ) ) + 1

  fcs.event.regexp <- sprintf( "\\.[0-9]{%d}$", fcs.event.number.width )

  cat( "\033[34m", "Extracting marker expression data...", "\033[0m\n" )

  fcs.expr.data <- lapply( 1 : sample.n, function( flow.idx ) {
    expr.data <- exprs( fcs.flow.set[[ flow.idx ]] )

    expr.data <- expr.data[ , fcs.panel$panel$dye ]
    colnames( expr.data ) <- fcs.panel$panel$antigen

    fcs.the.sample <- sample.name[ flow.idx ]

    fcs.the.event <- sprintf( "%s.%0*d", fcs.the.sample,
                              fcs.event.number.width, 1 : nrow( expr.data ) )
    rownames( expr.data ) <- fcs.the.event

    stopifnot( sort( colnames( expr.data ) ) == sort( fcs.panel$panel$antigen ) )

    expr.data <- expr.data[ , fcs.panel$marker ]

    attr( expr.data, "ranges" ) <- NULL
    attr( colnames( expr.data ), "names" ) <- NULL

    expr.data
  } )

  fcs.expr.data <- do.call( rbind, fcs.expr.data )

  fcs.expr.data[ , agp$default.time.parameter ] <- fcs.expr.data[ , agp$default.time.parameter ] / 100

  rm( fcs.flow.set )

  # set events

  fcs.event <- rownames( fcs.expr.data )

  fcs.event.n <- length( fcs.event )

  fcs.event.sample <- sub( fcs.event.regexp, "", fcs.event )
  fcs.event.sample <- factor( fcs.event.sample, levels = sample.name )

  fcs.data <- list(
    expr.data = fcs.expr.data,
    event = fcs.event,
    event.n = fcs.event.n,
    event.sample = fcs.event.sample
  )

  cat( "\033[34m", "Flow data extracted.", "\033[0m\n" )

  return( fcs.data )
}
