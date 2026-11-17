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
#'   \item Reads and validates the experiment manifest CSV.
#'   \item Reads and validates the panel CSV; drops the \code{AutoF}
#'     autofluorescence channel if present.
#'   \item Builds the Gating-ML transforms for every panel channel via
#'     \code{\link{read.transformation}} (logicle for fluorescence, linear
#'     for scatter, identity for time).
#'   \item Reads each FCS file with \code{\link[AutoSpectral]{readFCS}},
#'     keeping only the panel channels (raw values, no truncation at
#'     \code{$PnR}), applies compensation if a spillover matrix is supplied,
#'     then transforms every channel onto the 0 to 1 display scale.
#'   \item Stacks the per-sample matrices into one matrix, with a factor
#'     giving each event's sample.
#' }
#'
#' Compensation follows the FCS \code{$SPILLOVER} convention: the CSV holds the spillover
#' matrix (rows = fluorochromes, columns = detectors, diagonal 1), and the
#' compensated values are \eqn{X S^{-1}} on the lineage and activation
#' channels.
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
#'   \code{ap$data.dir}, \code{ap$compensation.filename}
#'   and \code{ap$transformation.filename}.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{fcs.expr.data}}{Numeric matrix, all events x all markers,
#'       on the display scale, without row names.  Column order matches
#'       \code{fcs.marker} (scatter, time, lineage, activation).  Events are
#'       in sample order, then file order.}
#'     \item{\code{fcs.sample}}{Character vector of unique sample identifiers
#'       (one per FCS file), in manifest order.}
#'     \item{\code{fcs.sample.n}}{Integer.  Total number of samples.}
#'     \item{\code{fcs.sample.filename}}{Character vector of FCS filenames
#'       (relative to \code{ap$data.dir}), in manifest order.}
#'     \item{\code{fcs.sample.gate.calc}}{Character vector of sample names
#'       for which \code{gate.calculation == TRUE} in the manifest.}
#'     \item{\code{fcs.event.n}}{Integer.  Total number of events.}
#'     \item{\code{fcs.event.sample}}{Factor.  Per-event sample membership
#'       (one entry per row of \code{fcs.expr.data}), levels in manifest
#'       order.}
#'     \item{\code{fcs.panel}}{data.frame.  Validated panel table (columns
#'       \code{dye}, \code{antigen}, \code{type}) with \code{type} as a factor
#'       and \code{AutoF} rows removed.}
#'     \item{\code{fcs.panel.scatter.idx}}{Integer vector.  Row indices of
#'       scatter channels in \code{fcs.panel}.}
#'     \item{\code{fcs.panel.lineage.idx}}{Integer vector.  Row indices of
#'       lineage channels.}
#'     \item{\code{fcs.panel.activation.idx}}{Integer vector.  Row indices of
#'       activation channels.}
#'     \item{\code{fcs.panel.time.idx}}{Integer vector.  Row indices of the
#'       time channel (raw values, not rescaled).}
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
#'       \code{ap$transform} and \code{ap$transform.inv} (forward and
#'       inverse transform functions keyed by antigen name) and
#'       \code{ap$transform.params} (the transform parameters, keyed by
#'       antigen).}
#'   }
#'
#' @seealso \code{\link{get.autogating.param}}, \code{\link{define.gates}},
#'   \code{\link{run.autogating}}, \code{\link{read.transformation}},
#'   \code{\link{read.compensation.csv}}
#'
#' @importFrom AutoSpectral readFCS
#'
#' @export
setup.autogating <- function( ap ) {

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

    # Remove the spectral autofluorescence channel - it is not gated
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
    # 3. Build transformations; read compensation matrix
    # ==========================================================================

    fcs.compensation <- NULL
    if ( !is.null( ap$compensation.filename ) ) {
        fcs.compensation <- read.compensation.csv(
            ap$compensation.filename, fcs.panel,
            fcs.panel.lineage.idx, fcs.panel.activation.idx )
        fcs.compensation.inv <- solve( fcs.compensation )
    }

    fcs.transform.all <- read.transformation(
        ap$transformation.filename, fcs.panel,
        fcs.panel.scatter.idx, fcs.panel.lineage.idx,
        fcs.panel.activation.idx, ap, fcs.panel.time.idx )

    # Attach transforms to ap so downstream functions can convert between
    # raw and display values
    ap$transform        <- fcs.transform.all[[ 1 ]]
    ap$transform.inv    <- fcs.transform.all[[ 2 ]]
    ap$transform.params <- fcs.transform.all[[ 3 ]]


    # ==========================================================================
    # 4. Read FCS files; compensate and transform
    # ==========================================================================

    fcs.expr.list <- lapply( seq_len( fcs.sample.n ), function( i ) {

        ed <- readFCS( file.path( ap$data.dir, fcs.sample.filename[ i ] ),
                       columns = fcs.panel$dye )
        colnames( ed ) <- fcs.panel$antigen

        if ( !is.null( fcs.compensation ) ) {
            comp.cols <- colnames( fcs.compensation )
            ed[ , comp.cols ] <- ed[ , comp.cols, drop = FALSE ] %*%
                fcs.compensation.inv
        }

        for ( ant in fcs.panel$antigen )
            ed[ , ant ] <- ap$transform[[ ant ]]( ed[ , ant ] )

        ed[ , fcs.marker, drop = FALSE ]
    } )


    # ==========================================================================
    # 5. Stack per-sample matrices
    # ==========================================================================

    fcs.sample.event.n <- vapply( fcs.expr.list, nrow, integer( 1L ) )

    fcs.expr.data <- do.call( rbind, fcs.expr.list )
    rm( fcs.expr.list )
    dimnames( fcs.expr.data ) <- list( NULL, fcs.marker )

    fcs.event.n      <- nrow( fcs.expr.data )
    fcs.event.sample <- factor(
        rep( fcs.sample, times = fcs.sample.event.n ),
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
        fcs.event.n              = fcs.event.n,
        fcs.event.sample         = fcs.event.sample,
        fcs.panel                = fcs.panel,
        fcs.panel.scatter.idx    = fcs.panel.scatter.idx,
        fcs.panel.lineage.idx    = fcs.panel.lineage.idx,
        fcs.panel.activation.idx = fcs.panel.activation.idx,
        fcs.panel.time.idx       = fcs.panel.time.idx,
        fcs.marker               = fcs.marker,
        fcs.marker.scatter       = fcs.marker.scatter,
        fcs.marker.lineage       = fcs.marker.lineage,
        fcs.marker.activation    = fcs.marker.activation,
        fcs.marker.label.figure  = fcs.marker.label.figure,
        fcs.marker.label.table   = fcs.marker.label.table,
        ap                       = ap
    )
}
