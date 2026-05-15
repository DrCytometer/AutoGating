# run_autogating.R
#
# Applies pre-calculated gate boundaries to every individual sample, computes
# population and activation statistics, writes output CSVs, and optionally
# runs differential abundance analysis.


#' Apply gates per sample and export statistics
#'
#' @title Run Per-Sample Gating and Export Statistics
#'
#' @description
#' Executes the per-sample application stage of the AutoGating pipeline:
#' \enumerate{
#'   \item For each sample, subsets the expression matrix and calls
#'     \code{\link{do.population.gates}} and \code{\link{do.activation.gates}}
#'     with \code{calculate.gates = FALSE}, reusing boundaries stored in
#'     \code{fcs.gates.all}.  Diagnostic figures are written per sample.
#'   \item Calls \code{\link{calculate.population.statistics}} and
#'     \code{\link{calculate.activation.statistics}} to produce wide-format
#'     data frames of per-sample statistics.
#'   \item Writes the statistics tables as comma-separated files to
#'     \code{ap$statistics.dir}.
#'   \item If \code{experiment.da.filename} is not \code{NULL}, runs
#'     \code{\link{calculate.differential.abundance}} and returns its output.
#' }
#'
#' Sample-level gating runs sequentially by default.  Parallel execution is
#' not implemented inside this function; if parallelism is needed, callers
#' may replace the internal \code{lapply} with \code{parallel::mclapply}
#' after setting \code{ap$worker.n}.
#'
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}.
#'   Required elements:
#'   \describe{
#'     \item{\code{fcs.expr.data}}{Numeric matrix. Full expression matrix for
#'       all events (rows) across all samples.}
#'     \item{\code{fcs.sample}}{Character vector of sample identifiers.}
#'     \item{\code{fcs.event.sample}}{Factor of per-event sample membership.}
#'   }
#' @param fcs.defined Named list returned by \code{\link{define.gates}}.
#'   Required elements:
#'   \describe{
#'     \item{\code{fcs.gates.all}}{Nested gate results list carrying the
#'       shared boundary polygons.}
#'     \item{\code{fcs.population.gates}}{List of population gate
#'       definitions.}
#'     \item{\code{fcs.activation.gates}}{List of activation gate
#'       definitions.}
#'     \item{\code{fcs.figure.popgate.dir}}{Named character vector of
#'       per-sample population-gate figure directory paths.}
#'     \item{\code{fcs.figure.actgate.dir}}{Named character vector of
#'       per-sample activation-gate figure directory paths.}
#'   }
#' @param ap Autogating parameter object (as returned by
#'   \code{\link{get.autogating.param}} and augmented by
#'   \code{\link{setup.autogating}}).  Required fields:
#'   \code{ap$statistics.dir}, \code{ap$population.statistics.filename},
#'   \code{ap$activation.statistics.filename}, \code{ap$transform.inv}.
#' @param experiment.da.filename Character scalar or \code{NULL}.  Path to a
#'   supplementary experiment CSV with \code{Group} (and optionally
#'   \code{Block} and \code{Covariate.*}) columns for differential abundance
#'   analysis.  When \code{NULL} (default), the differential abundance step
#'   is skipped.
#' @param plot.gates Logical.  If \code{TRUE} (default), per-sample
#'   diagnostic scatter plots are written to the sample figure directories.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{fcs.gates.sample}}{Named list of per-sample gate results
#'       (same structure as \code{fcs.gates.all} but one element per sample).}
#'     \item{\code{population.stats}}{data.frame of population fraction
#'       statistics, one row per sample.}
#'     \item{\code{activation.stats}}{data.frame of activation statistics
#'       (fraction positive, median, IQR), one row per sample.}
#'     \item{\code{differential.abundance}}{The object returned by
#'       \code{\link{calculate.differential.abundance}}, or \code{NULL} if
#'       \code{experiment.da.filename} was not supplied.}
#'   }
#'
#' @seealso \code{\link{setup.autogating}}, \code{\link{define.gates}},
#'   \code{\link{do.population.gates}}, \code{\link{do.activation.gates}},
#'   \code{\link{calculate.population.statistics}},
#'   \code{\link{calculate.activation.statistics}},
#'   \code{\link{calculate.differential.abundance}}
#'
#' @export
run.autogating <- function(
    fcs.setup,
    fcs.defined,
    ap,
    experiment.da.filename = NULL,
    plot.gates             = TRUE
) {

    fcs.expr.data          <- fcs.setup$fcs.expr.data
    fcs.sample             <- fcs.setup$fcs.sample
    fcs.event.sample       <- fcs.setup$fcs.event.sample

    fcs.gates.all          <- fcs.defined$fcs.gates.all
    fcs.population.gates   <- fcs.defined$fcs.population.gates
    fcs.activation.gates   <- fcs.defined$fcs.activation.gates
    fcs.figure.popgate.dir <- fcs.defined$fcs.figure.popgate.dir
    fcs.figure.actgate.dir <- fcs.defined$fcs.figure.actgate.dir


    # ==========================================================================
    # 1. Apply gates to each individual sample
    # ==========================================================================

    fcs.gates.sample <- lapply( fcs.sample, function( samp ) {

        samp.data  <- fcs.expr.data[ fcs.event.sample == samp, ]

        # Start from the shared boundaries; boundaries are not re-calculated
        samp.gates <- fcs.gates.all

        samp.gates <- do.population.gates(
            samp.gates,
            samp.data,
            fcs.population.gates,
            ap,
            calculate.gates = FALSE,
            figure.dir      = fcs.figure.popgate.dir[ samp ],
            plot.gates      = plot.gates
        )

        samp.gates <- do.activation.gates(
            samp.gates,
            samp.data,
            fcs.activation.gates,
            ap,
            calculate.gates = FALSE,
            figure.dir      = fcs.figure.actgate.dir[ samp ],
            plot.gates      = plot.gates
        )

        samp.gates
    } )

    names( fcs.gates.sample ) <- fcs.sample


    # ==========================================================================
    # 2. Calculate statistics
    # ==========================================================================

    population.stats <- calculate.population.statistics(
        fcs.gates.sample, fcs.population.gates, ap )

    activation.stats <- calculate.activation.statistics(
        fcs.gates.sample, fcs.activation.gates, fcs.expr.data, ap )


    # ==========================================================================
    # 3. Write statistics to disk
    # ==========================================================================

    utils::write.table(
        population.stats,
        file      = file.path( ap$statistics.dir,
                               ap$population.statistics.filename ),
        sep       = ",",
        quote     = FALSE,
        row.names = FALSE
    )

    utils::write.table(
        activation.stats,
        file      = file.path( ap$statistics.dir,
                               ap$activation.statistics.filename ),
        sep       = ",",
        quote     = FALSE,
        row.names = FALSE
    )

    message( sprintf( "Statistics written to %s.", ap$statistics.dir ) )


    # ==========================================================================
    # 4. Optional differential abundance analysis
    # ==========================================================================

    da.result <- NULL

    if ( !is.null( experiment.da.filename ) ) {
        da.result <- calculate.differential.abundance(
            population.stats,
            activation.stats,
            ap,
            experiment.filename = experiment.da.filename
        )
        message( "Differential abundance analysis complete." )
    }


    # ==========================================================================
    # Return
    # ==========================================================================

    list(
        fcs.gates.sample      = fcs.gates.sample,
        population.stats      = population.stats,
        activation.stats      = activation.stats,
        differential.abundance = da.result
    )
}
