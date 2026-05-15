# setup_autogating.R
#
# Reads all experiment inputs, applies compensation and transformation, and
# returns the combined expression matrix together with the sample metadata
# structures needed by the rest of the pipeline.


#' Read, compensate, and transform FCS data for an experiment
#'
#' @title Set Up FCS Data for Autogating
#'
#' @description
#' Orchestrates the input stage of the AutoGating pipeline:
#' \enumerate{
#'   \item Sets the random seed from \code{ap$seed}.
#'   \item Reads and validates the experiment manifest CSV.
#'   \item Reads and validates the panel CSV; drops the \code{AutoF}
#'     autofluorescence channel if present.
#'   \item Builds forward and inverse biexponential transform lists via
#'     \code{\link{read.transformation}}.
#'   \item Reads all FCS files into a \code{flowSet}, subsets and renames
#'     channels to antigen names, then applies spectral compensation (if a
#'     compensation matrix is supplied) and biexponential transformation.
#'   \item Concatenates all per-sample expression matrices into a single
#'     event-labelled numeric matrix.
#' }
#'
#' The function returns a named list (\code{fcs.setup}) that bundles the
#' expression matrix with all derived metadata vectors and panel indices
#' required by \code{\link{define.gates}} and \code{\link{run.autogating}}.
#' Transforms are also written back into \code{ap} (via the \code{ap}
#' argument, passed by reference using \code{<<-} would be inappropriate;
#' instead the caller should update \code{ap} from
#' \code{fcs.setup$ap}).
#'
#' @param ap Autogating parameter object returned by
#'   \code{\link{get.autogating.param}}.  Must contain valid values for
#'   \code{ap$experiment.filename}, \code{ap$panel.filename},
#'   \code{ap$data.dir}, \code{ap$compensation.filename},
#'   \code{ap$transformation.filename}, and \code{ap$seed}.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{fcs.expr.data}}{Numeric matrix, all events × all markers.
#'       Row names are \code{"<sample>.<event_number>"} strings.  Column order
#'       matches \code{fcs.marker} (scatter, time, lineage, activation).}
#'     \item{\code{fcs.sample}}{Character vector of unique sample identifiers
#'       (one per FCS file), in manifest order.}
#'     \item{\code{fcs.sample.n}}{Integer.  Total number of samples.}
#'     \item{\code{fcs.sample.filename}}{Character vector of FCS filenames
#'       (relative to \code{ap$data.dir}), in manifest order.}
#'     \item{\code{fcs.sample.gate.calc}}{Character vector of sample names
#'       for which \code{gate.calculation == TRUE} in the manifest.}
#'     \item{\code{fcs.event}}{Character vector of all event identifiers
#'       (row names of \code{fcs.expr.data}).}
#'     \item{\code{fcs.event.n}}{Integer.  Total number of events.}
#'     \item{\code{fcs.event.sample}}{Factor.  Per-event sample membership,
#'       levels in manifest order.}
#'     \item{\code{fcs.event.regexp}}{Character scalar.  Regular expression
#'       matching the \code{".<event_number>"} suffix so event identifiers can
#'       be stripped back to sample names.}
#'     \item{\code{fcs.panel}}{data.frame.  Validated panel table (columns
#'       \code{dye}, \code{antigen}, \code{type}) with \code{type} as a factor
#'       and \code{AutoF} rows removed.}
#'     \item{\code{fcs.panel.scatter.idx}}{Integer vector.  Row indices of
#'       scatter channels in \code{fcs.panel}.}
#'     \item{\code{fcs.panel.lineage.idx}}{Integer vector.  Row indices of
#'       lineage channels.}
#'     \item{\code{fcs.panel.activation.idx}}{Integer vector.  Row indices of
#'       activation channels.}
#'     \item{\code{fcs.marker}}{Character vector.  All antigen names in
#'       pipeline column order (scatter, time, lineage, activation).}
#'     \item{\code{fcs.marker.scatter}}{Character vector.  Scatter antigen
#'       names.}
#'     \item{\code{fcs.marker.lineage}}{Character vector.  Lineage antigen
#'       names.}
#'     \item{\code{fcs.marker.activation}}{Character vector.  Activation
#'       antigen names.}
#'     \item{\code{fcs.marker.label.figure}}{Named character vector.  Display
#'       labels for axis titles (format: \code{"<dye> - <antigen>"} for
#'       fluorescence channels), named by antigen.}
#'     \item{\code{fcs.marker.label.table}}{Named character vector.  Antigen
#'       names as column labels for statistics tables, named by antigen.}
#'     \item{\code{ap}}{The input \code{ap} list, augmented with
#'       \code{ap$transform} and \code{ap$transform.inv} (the forward and
#'       inverse transform function lists keyed by antigen name).}
#'   }
#'
#' @seealso \code{\link{get.autogating.param}}, \code{\link{define.gates}},
#'   \code{\link{run.autogating}}, \code{\link{read.transformation}},
#'   \code{\link{read.compensation.csv}}
#'
#' @importFrom flowCore read.flowSet compensate transform transformList exprs
#'   compensation
#'
#' @export
setup.autogating <- function( ap ) {

    set.seed( ap$seed )

    # ==========================================================================
    # 1. Read and validate experiment manifest
    # ==========================================================================

    fcs.experiment <- utils::read.csv( ap$experiment.filename,
                                       stringsAsFactors = FALSE )

    stopifnot( anyDuplicated( fcs.experiment$filename ) == 0 )
    stopifnot( anyDuplicated( fcs.experiment$sample )   == 0 )

    # gate.calculation: treat missing / empty as FALSE (apply but do not
    # re-calculate boundaries for unlabelled samples)
    fcs.experiment$gate.calculation[
        is.na( fcs.experiment$gate.calculation ) |
        fcs.experiment$gate.calculation == "" ] <- FALSE
    fcs.experiment$gate.calculation <-
        as.logical( fcs.experiment$gate.calculation )

    fcs.sample           <- fcs.experiment$sample
    fcs.sample.n         <- length( fcs.sample )
    fcs.sample.filename  <- fcs.experiment$filename
    fcs.sample.gate.calc <- fcs.sample[ fcs.experiment$gate.calculation ]

    message( sprintf( "%d samples total, %d used for gate calculation.",
                      fcs.sample.n, length( fcs.sample.gate.calc ) ) )


    # ==========================================================================
    # 2. Read and validate panel
    # ==========================================================================

    fcs.panel <- utils::read.csv( ap$panel.filename, stringsAsFactors = FALSE )

    # Remove the spectral autofluorescence channel — it is not gated
    fcs.panel <- fcs.panel[ fcs.panel$antigen != "AutoF", ]

    fcs.panel$type <- factor( fcs.panel$type,
        levels = c( "scatter", "time", "lineage", "activation" ) )

    fcs.panel.scatter.idx    <- which( fcs.panel$type == "scatter" )
    fcs.panel.time.idx       <- which( fcs.panel$type == "time" )
    fcs.panel.lineage.idx    <- which( fcs.panel$type == "lineage" )
    fcs.panel.activation.idx <- which( fcs.panel$type == "activation" )

    # Every panel row must belong to exactly one type
    stopifnot(
        sort( c( fcs.panel.scatter.idx, fcs.panel.time.idx,
                 fcs.panel.lineage.idx, fcs.panel.activation.idx ) ) ==
        seq_len( nrow( fcs.panel ) )
    )

    fcs.marker.scatter    <- fcs.panel$antigen[ fcs.panel.scatter.idx ]
    fcs.marker.time       <- fcs.panel$antigen[ fcs.panel.time.idx ]
    fcs.marker.lineage    <- fcs.panel$antigen[ fcs.panel.lineage.idx ]
    fcs.marker.activation <- fcs.panel$antigen[ fcs.panel.activation.idx ]

    # Pipeline column order: scatter, time, lineage, activation
    fcs.marker <- c( fcs.marker.scatter, fcs.marker.time,
                     fcs.marker.lineage, fcs.marker.activation )

    # Figure axis labels: fluorescence channels prefixed with dye name
    fcs.marker.label <- c(
        fcs.marker.scatter,
        fcs.marker.time,
        paste0( fcs.panel$dye[ fcs.panel.lineage.idx ],
                " - ", fcs.marker.lineage ),
        paste0( fcs.panel$dye[ fcs.panel.activation.idx ],
                " - ", fcs.marker.activation )
    )

    fcs.marker.label.figure          <- fcs.marker.label
    names( fcs.marker.label.figure ) <- fcs.marker

    fcs.marker.label.table          <- fcs.marker
    names( fcs.marker.label.table ) <- fcs.marker


    # ==========================================================================
    # 3. Build transformation (and optionally read compensation matrix)
    # ==========================================================================

    if ( !is.null( ap$compensation.filename ) ) {
        fcs.compensation <- read.compensation.csv(
            ap$compensation.filename, fcs.panel,
            fcs.panel.lineage.idx, fcs.panel.activation.idx )
    }

    fcs.transform.both <- read.transformation(
        ap$transformation.filename, fcs.panel,
        fcs.panel.scatter.idx, fcs.panel.lineage.idx,
        fcs.panel.activation.idx, ap )

    fcs.transform     <- fcs.transform.both[[ 1 ]]
    fcs.transform.inv <- fcs.transform.both[[ 2 ]]

    # Attach transforms to ap so downstream functions (e.g.
    # calculate.activation.statistics) can back-transform expression values
    ap$transform     <- fcs.transform
    ap$transform.inv <- fcs.transform.inv


    # ==========================================================================
    # 4. Read FCS files; apply compensation and transformation
    # ==========================================================================

    fcs.flow.set <- suppressWarnings(
        flowCore::read.flowSet(
            files               = fcs.sample.filename,
            path                = ap$data.dir,
            transformation      = NULL,
            truncate_max_range  = TRUE
        )
    )

    stopifnot( fcs.panel$dye %in% colnames( fcs.flow.set ) )

    fcs.flow.set <- fcs.flow.set[ , fcs.panel$dye ]
    colnames( fcs.flow.set ) <- fcs.panel$antigen

    if ( !is.null( ap$compensation.filename ) ) {
        fcs.flow.set <- flowCore::compensate(
            fcs.flow.set,
            flowCore::compensation( fcs.compensation )
        )
    }

    fcs.flow.set <- flowCore::transform(
        fcs.flow.set,
        flowCore::transformList( names( fcs.transform ), fcs.transform )
    )


    # ==========================================================================
    # 5. Extract and label expression matrix
    # ==========================================================================

    # Determine zero-padded event-ID width from the largest sample so that
    # row names sort lexicographically
    fcs.sample.event.n.max <- max( sapply(
        seq_len( fcs.sample.n ),
        function( i ) nrow( flowCore::exprs( fcs.flow.set[[ i ]] ) )
    ) )
    fcs.event.number.width <- floor( log10( fcs.sample.event.n.max ) ) + 1L
    fcs.event.regexp       <- sprintf( "\\.[0-9]{%d}$",
                                       fcs.event.number.width )

    fcs.expr.data <- do.call( rbind, lapply(
        seq_len( fcs.sample.n ),
        function( i ) {
            ed   <- flowCore::exprs( fcs.flow.set[[ i ]] )
            samp <- fcs.sample[ i ]

            rownames( ed ) <- sprintf(
                "%s.%0*d", samp, fcs.event.number.width,
                seq_len( nrow( ed ) )
            )

            stopifnot( sort( colnames( ed ) ) == sort( fcs.panel$antigen ) )

            ed <- ed[ , fcs.marker ]
            attr( ed, "ranges" )            <- NULL
            attr( colnames( ed ), "names" ) <- NULL
            ed
        }
    ) )

    # Rescale time channel from raw units to seconds (÷ 100)
    fcs.expr.data[ , "Time" ] <- fcs.expr.data[ , "Time" ] / 100

    fcs.event        <- rownames( fcs.expr.data )
    fcs.event.n      <- length( fcs.event )
    fcs.event.sample <- factor(
        sub( fcs.event.regexp, "", fcs.event ),
        levels = fcs.sample
    )

    message( sprintf( "%d total events across %d samples.",
                      fcs.event.n, fcs.sample.n ) )


    # ==========================================================================
    # Return
    # ==========================================================================

    list(
        fcs.expr.data            = fcs.expr.data,
        fcs.sample               = fcs.sample,
        fcs.sample.n             = fcs.sample.n,
        fcs.sample.filename      = fcs.sample.filename,
        fcs.sample.gate.calc     = fcs.sample.gate.calc,
        fcs.event                = fcs.event,
        fcs.event.n              = fcs.event.n,
        fcs.event.sample         = fcs.event.sample,
        fcs.event.regexp         = fcs.event.regexp,
        fcs.panel                = fcs.panel,
        fcs.panel.scatter.idx    = fcs.panel.scatter.idx,
        fcs.panel.lineage.idx    = fcs.panel.lineage.idx,
        fcs.panel.activation.idx = fcs.panel.activation.idx,
        fcs.marker               = fcs.marker,
        fcs.marker.scatter       = fcs.marker.scatter,
        fcs.marker.lineage       = fcs.marker.lineage,
        fcs.marker.activation    = fcs.marker.activation,
        fcs.marker.label.figure  = fcs.marker.label.figure,
        fcs.marker.label.table   = fcs.marker.label.table,
        ap                       = ap
    )
}
