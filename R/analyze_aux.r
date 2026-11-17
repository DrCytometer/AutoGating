# analyze_aux.r
#
# Utility functions for the FCS autogating pipeline.
#
# Contents:
#   seeded.call()             - reproducible subsamples
#   read.transformation()     - build Gating-ML transform lists
#   load.gate.definitions()   - load an R list or RDS gate definition file


#' Call a function with a fixed random seed
#'
#' @title Seeded Call
#'
#' @description
#' Sets the seed, calls \code{f} and restores the caller's random number
#' stream, so that subsamples are reproducible without changing the random
#' numbers drawn elsewhere in the session.
#'
#' @param seed Integer seed.
#' @param f Function of no arguments.
#'
#' @return The value of \code{f()}.
#'
#' @keywords internal
seeded.call <- function( seed, f ) {
    env <- globalenv()
    had.seed <- exists( ".Random.seed", envir = env, inherits = FALSE )
    if ( had.seed ) old.seed <- get( ".Random.seed", envir = env, inherits = FALSE )
    on.exit( {
        if ( had.seed ) assign( ".Random.seed", old.seed, envir = env )
        else if ( exists( ".Random.seed", envir = env, inherits = FALSE ) )
            rm( ".Random.seed", envir = env )
    } )
    set.seed( seed )
    f()
}


#' Build forward and inverse transform lists for all panel channels
#'
#' @title Build Channel Transformations
#'
#' @description
#' Builds the Gating-ML 2.0 transform (see \code{\link{make.transform}}) for
#' every panel channel: logicle for fluorescence channels, linear for
#' scatter channels and the identity for the time channel, starting from the
#' cytometer defaults in \code{ap$transform.defaults}. A per-experiment CSV
#' can override fluorescence or scatter channels. Two CSV layouts are
#' accepted, keyed by the \code{dye} column:
#' \itemize{
#'   \item Gating-ML parameters: any of \code{id}, \code{scale.t},
#'     \code{linear.a}, \code{logicle.w}, \code{logicle.m},
#'     \code{logicle.a}, \code{log.m}.
#'   \item FlowJo biex settings (\code{max.range}, \code{pos}, \code{neg},
#'     \code{width}; \code{length} is ignored), converted with
#'     \code{\link{flowjo.biex.to.logicle}}.
#' }
#'
#' @param transformation.filepath Character or \code{NULL}.  Path to the
#'   override CSV.  \code{NULL} uses cytometer defaults throughout.
#' @param fcs.panel \code{data.frame}.  Panel definition with columns
#'   \code{dye} and \code{antigen}.
#' @param fcs.panel.scatter.idx Integer vector.  Scatter-channel row indices.
#' @param fcs.panel.lineage.idx Integer vector.  Lineage-channel row indices.
#' @param fcs.panel.activation.idx Integer vector.  Activation-channel row
#'   indices.
#' @param ap Autogating parameter object; uses \code{ap$transform.defaults}.
#' @param fcs.panel.time.idx Integer vector.  Time-channel row indices
#'   (identity transform).  Default none.
#'
#' @return A list:
#'   \describe{
#'     \item{\code{[[1]]}}{Named list of forward functions (raw to display),
#'       keyed by antigen name.}
#'     \item{\code{[[2]]}}{Named list of inverse functions, keyed by antigen.}
#'     \item{\code{[[3]]}}{Named list of transform parameter lists, keyed by
#'       antigen.}
#'   }
#'
#' @seealso \code{\link{get.autogating.param}}, \code{\link{make.transform}}
#'
#' @keywords internal
read.transformation <- function( transformation.filepath, fcs.panel,
    fcs.panel.scatter.idx, fcs.panel.lineage.idx, fcs.panel.activation.idx,
    ap, fcs.panel.time.idx = integer( 0 ) ) {

    fluor.idx <- c( fcs.panel.lineage.idx, fcs.panel.activation.idx )

    params <- list()
    for ( i in fluor.idx )
        params[[ fcs.panel$antigen[ i ] ]] <- ap$transform.defaults$fluorescence
    for ( i in fcs.panel.scatter.idx )
        params[[ fcs.panel$antigen[ i ] ]] <- ap$transform.defaults$scatter
    for ( i in fcs.panel.time.idx )
        params[ fcs.panel$antigen[ i ] ] <- list( NULL )

    if ( !is.null( transformation.filepath ) ) {

        ovr <- utils::read.csv( transformation.filepath, stringsAsFactors = FALSE )
        stopifnot( "dye" %in% colnames( ovr ) )

        gml.cols    <- c( "id", "scale.t", "linear.a", "logicle.w",
                          "logicle.m", "logicle.a", "log.m" )
        flowjo.cols <- c( "max.range", "pos", "neg", "width" )
        is.flowjo   <- all( flowjo.cols %in% colnames( ovr ) )
        if ( !is.flowjo && !any( gml.cols %in% colnames( ovr ) ) )
            stop( sprintf(
                "Transformation CSV '%s' has neither Gating-ML columns (%s) nor FlowJo biex columns (%s).",
                transformation.filepath, paste( gml.cols, collapse = ", " ),
                paste( flowjo.cols, collapse = ", " ) ) )

        unknown <- setdiff( ovr$dye, fcs.panel$dye )
        if ( length( unknown ) > 0L )
            message( sprintf( "Transformation CSV: ignoring dyes not in the panel: %s",
                paste( unknown, collapse = ", " ) ) )

        for ( r in which( ovr$dye %in% fcs.panel$dye[ c( fluor.idx, fcs.panel.scatter.idx ) ] ) ) {
            ant <- fcs.panel$antigen[ match( ovr$dye[ r ], fcs.panel$dye ) ]
            if ( is.flowjo ) {
                params[[ ant ]] <- flowjo.biex.to.logicle( ovr$max.range[ r ],
                    ovr$pos[ r ], ovr$neg[ r ], ovr$width[ r ] )
            } else {
                p <- params[[ ant ]]
                for ( col in intersect( gml.cols, colnames( ovr ) ) )
                    if ( !is.na( ovr[[ col ]][ r ] ) )
                        p[[ col ]] <- ovr[[ col ]][ r ]
                params[[ ant ]] <- p
            }
        }
    }

    fns <- lapply( params, make.transform )

    list(
        lapply( fns, `[[`, "forward" ),
        lapply( fns, `[[`, "inverse" ),
        params
    )
}


#' Load a gate definition list from an R file or RDS file
#'
#' @title Load Gate Definitions from File
#'
#' @description
#' Loads a gate definition list from either an R source file or an RDS file.
#' For R source files, the file is sourced into a fresh environment and the
#' single variable it defines is returned.  If the file defines more than one
#' variable, an error is raised.  If \code{filepath} is \code{NULL}, an empty
#' list is returned.
#'
#' @param filepath Character scalar or \code{NULL}.  Path to the gate
#'   definition file.  Extension determines loading method: \code{".rds"}
#'   uses \code{\link[base]{readRDS}}; any other extension uses
#'   \code{\link[base]{source}}.
#'
#' @return A list containing the gate definitions (population or activation),
#'   or an empty list if \code{filepath} is \code{NULL}.
#'
#' @seealso \code{\link{get.autogating.param}},
#'   \code{\link{convert.gate.definitions}}, \code{\link{define.gates}}
#'
#' @keywords internal
load.gate.definitions <- function( filepath ) {
    if ( is.null( filepath ) ) return( list() )

    ext <- tolower( tools::file_ext( filepath ) )

    if ( ext == "rds" ) {
        readRDS( filepath )
    } else {
        env <- new.env( parent = baseenv() )
        source( filepath, local = env )
        objs <- ls( env )
        if ( length( objs ) != 1L )
            stop( sprintf(
                "Gate definition file '%s' must assign exactly one variable.",
                filepath ) )
        env[[ objs[ 1 ] ]]
    }
}
