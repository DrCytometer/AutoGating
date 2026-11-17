# run_autogating.R
#
# Applies a trained gating model to every sample and writes the population
# counts, the marker medians of each population and the gating status of
# each gate as CSV files.


#' Apply a gating model to every sample
#'
#' @title Run Gating on Every Sample
#'
#' @description
#' Applies the model returned by \code{\link{define.gates}} (or read with
#' \code{\link{read.agmodel}}) to each sample on its own, as Honeychrome
#' Automated Gating does: each gate uses its application mode from
#' \code{model$application} ("fixed" applies the stored boundary,
#' "recalculate" calculates it on the sample and flags a move larger than
#' the drift limit).  Stored boundaries are first converted to the data's
#' axis transforms (\code{\link{remap.model}}).  With
#' \code{ap$worker.n > 1}, samples are gated in parallel through
#' \code{\link[AutoSpectral]{create.parallel.lapply}}.
#'
#' With \code{ap$figure.samples}, each gate is drawn on each sample in
#' \code{ap$figure.dir/<sample>} (\code{\link{write.gate.figures}}).
#'
#' Three CSV files are written to \code{ap$statistics.dir}:
#' \describe{
#'   \item{\code{ap$population.statistics.filename}}{One row per sample and
#'     population (see \code{\link{population.table}}), with a
#'     \code{sample} column.}
#'   \item{\code{ap$median.statistics.filename}}{Median of each lineage and
#'     activation channel in each population: \code{median} on the display
#'     scale and \code{median.raw} in raw units.}
#'   \item{\code{ap$status.filename}}{Per sample and gate: source of the
#'     boundary, drift, flag and message.}
#' }
#'
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}.
#' @param model Gating model with trained boundaries.
#' @param ap Autogating parameter object returned by
#'   \code{\link{get.autogating.param}}.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{results}}{Gating result of each sample (see
#'       \code{\link{run.gating}}); row indices refer to that sample's rows
#'       of \code{fcs.setup$fcs.expr.data}, in order.}
#'     \item{\code{population.table}, \code{medians}, \code{status}}{The
#'       tables written to disk.}
#'     \item{\code{model}}{The model as applied (boundaries in the data's
#'       transforms).}
#'   }
#'
#' @seealso \code{\link{setup.autogating}}, \code{\link{define.gates}},
#'   \code{\link{run.gating}}, \code{\link{population.table}},
#'   \code{\link{population.medians}}
#'
#' @importFrom AutoSpectral create.parallel.lapply
#'
#' @export
run.autogating <- function( fcs.setup, model, ap ) {

    transforms <- fcs.setup$ap$transform.params
    transforms <- transforms[ lengths( transforms ) > 0L ]
    remapped <- remap.model( model, transforms )
    model <- remapped$model
    for ( msg in remapped$messages )
        message( msg )

    gate.defs <- model$gate.definitions
    data <- fcs.setup$fcs.expr.data
    missing <- missing.gate.channels( gate.defs, colnames( data ) )
    if ( length( missing ) )
        stop( sprintf( "Gate channels not in the data: %s", paste( missing, collapse = ", " ) ) )

    app <- model$application
    reference <- model$trained.boundaries
    rows.by.sample <- split( seq_len( nrow( data ) ), fcs.setup$fcs.event.sample )
    marker.channels <- c( fcs.setup$fcs.marker.lineage, fcs.setup$fcs.marker.activation )
    inverse <- fcs.setup$ap$transform.inv
    figures <- isTRUE( ap$figure.samples )
    figure.dir <- ap$figure.dir %||% "figures"


    # ==========================================================================
    # 1. Gate each sample
    # ==========================================================================

    gate.sample <- function( s ) {
        d <- data[ rows.by.sample[[ s ]], , drop = FALSE ]
        res <- run.gating( gate.defs, d, reference = reference,
                           modes = app$mode.per.gate, default.mode = app$default.mode,
                           drift.limit = app$drift.limit, drift.action = app$drift.action,
                           transforms = transforms )
        pt <- population.table( gate.defs, res, nrow( d ) )
        if ( nrow( pt ) )
            pt <- cbind( sample = s, pt, stringsAsFactors = FALSE )
        med <- population.medians( gate.defs, res, d, marker.channels )
        if ( nrow( med ) ) {
            med$median.raw <- mapply( function( ch, v ) inverse[[ ch ]]( v ),
                                      med$channel, med$median, USE.NAMES = FALSE )
            med <- cbind( sample = s, med, stringsAsFactors = FALSE )
        }
        st <- do.call( rbind, lapply( names( res$status ), function( g ) {
            x <- res$status[[ g ]]
            data.frame( sample = s, gate = g, source = x$source, drift = x$drift,
                        flag = x$flag, message = x$message, stringsAsFactors = FALSE )
        } ) )
        if ( figures )
            write.gate.figures( gate.defs, res, d, fcs.setup, ap,
                                file.path( figure.dir, gsub( "[^A-Za-z0-9._+-]+", "_", s ) ), title = s )
        list( result = res, populations = pt, medians = med, status = st )
    }

    samples <- fcs.setup$fcs.sample
    worker.n <- as.integer( ap$worker.n %||% 1L )
    if ( worker.n > 1L ) {
        par <- create.parallel.lapply( list( worker.process.n = worker.n, bird.seed = ap$seed ),
                                       exports = character( 0 ), parallel = TRUE,
                                       threads = worker.n, export.env = environment(),
                                       allow.mclapply.mac = TRUE )
        if ( !is.null( par$cleanup ) )
            on.exit( par$cleanup(), add = TRUE )
        out <- par$lapply( samples, gate.sample )
    } else {
        out <- lapply( samples, gate.sample )
    }
    failed <- vapply( out, function( o ) inherits( o, "try-error" ) || !is.list( o ), logical( 1 ) )
    if ( any( failed ) )
        stop( sprintf( "Gating failed for sample(s): %s", paste( samples[ failed ], collapse = ", " ) ) )
    names( out ) <- samples


    # ==========================================================================
    # 2. Tables and files
    # ==========================================================================

    stack <- function( key ) {
        tbl <- do.call( rbind, lapply( out, `[[`, key ) )
        if ( is.null( tbl ) ) return( data.frame() )
        rownames( tbl ) <- NULL
        tbl
    }
    population.tbl <- stack( "populations" )
    medians        <- stack( "medians" )
    status         <- stack( "status" )

    flagged <- if ( nrow( status ) ) status[ !is.na( status$flag ), , drop = FALSE ] else status
    if ( nrow( flagged ) )
        message( sprintf( "%d gate result(s) flagged (%s); see %s.", nrow( flagged ),
                          paste( sprintf( "%s %d", names( table( flagged$flag ) ),
                                          as.integer( table( flagged$flag ) ) ), collapse = ", " ),
                          ap$status.filename ) )

    dir.create( ap$statistics.dir, recursive = TRUE, showWarnings = FALSE )
    utils::write.csv( population.tbl, file.path( ap$statistics.dir, ap$population.statistics.filename ),
                      row.names = FALSE )
    utils::write.csv( medians, file.path( ap$statistics.dir, ap$median.statistics.filename ),
                      row.names = FALSE )
    utils::write.csv( status, file.path( ap$statistics.dir, ap$status.filename ),
                      row.names = FALSE )
    message( sprintf( "Statistics written to %s.", ap$statistics.dir ) )
    if ( figures )
        message( sprintf( "Gate figures of each sample written to %s.", figure.dir ) )

    list(
        results          = lapply( out, `[[`, "result" ),
        population.table = population.tbl,
        medians          = medians,
        status           = status,
        model            = model
    )
}
