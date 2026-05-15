# define_gates.R
#
# Loads gate definition files, creates output directories, and calculates
# gate boundaries on the pooled gate-calculation samples.  Returns the
# fcs.gates.all accumulator ready for per-sample application by
# run.autogating().


#' Load gate definitions and calculate boundaries on pooled samples
#'
#' @title Define and Calculate Gate Boundaries
#'
#' @description
#' Performs the gate definition and boundary-calculation stage of the
#' AutoGating pipeline:
#' \enumerate{
#'   \item Loads population and activation gate definition lists from the
#'     files referenced in \code{ap} via \code{\link{load.gate.definitions}}.
#'   \item Creates the figure and statistics output directories on disk.
#'   \item Pools the expression data from all samples that have
#'     \code{gate.calculation == TRUE} in the experiment manifest.
#'   \item Calls \code{\link{do.population.gates}} and
#'     \code{\link{do.activation.gates}} on the pooled data with
#'     \code{calculate.gates = TRUE} to compute shared gate boundaries.
#' }
#'
#' The returned object (\code{fcs.gates.all}) is a nested named list whose
#' structure is documented in \code{CLAUDE.md}.  It carries both the computed
#' boundary polygons and the event-index assignments for the pooled data.
#' Pass it, together with the gate definition lists, to
#' \code{\link{run.autogating}} for per-sample application.
#'
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}.
#'   Required elements:
#'   \describe{
#'     \item{\code{fcs.expr.data}}{Numeric matrix. Full expression matrix for
#'       all events.}
#'     \item{\code{fcs.sample}}{Character vector of all sample identifiers.}
#'     \item{\code{fcs.sample.gate.calc}}{Character vector of sample
#'       identifiers used for boundary calculation.}
#'     \item{\code{fcs.event.sample}}{Factor of per-event sample membership.}
#'   }
#' @param ap Autogating parameter object returned by
#'   \code{\link{get.autogating.param}} (and augmented by
#'   \code{\link{setup.autogating}} with \code{ap$transform} and
#'   \code{ap$transform.inv}).  Required fields: \code{ap$figure.dir.basename},
#'   \code{ap$statistics.dir}, \code{ap$population.gates.filename},
#'   \code{ap$activation.gates.filename}.
#' @param plot.gates Logical.  If \code{TRUE} (default), diagnostic density
#'   scatter plots are saved to the \code{"all"} figure directories.  Set to
#'   \code{FALSE} to suppress plot generation (e.g. in automated testing).
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{fcs.gates.all}}{Nested gate results list (see
#'       \code{CLAUDE.md}) for the pooled gate-calculation data.  Contains
#'       \code{boundary}, \code{idx}, \code{fraction}, \code{label},
#'       \code{parent.gate}, and \code{parent.popul} for every non-ignored
#'       population in every gate.}
#'     \item{\code{fcs.population.gates}}{List of population gate definitions
#'       loaded from \code{ap$population.gates.filename}.}
#'     \item{\code{fcs.activation.gates}}{List of activation gate definitions
#'       loaded from \code{ap$activation.gates.filename}.}
#'     \item{\code{fcs.figure.popgate.dir}}{Named character vector of
#'       population-gate figure directory paths, one entry per sample plus an
#'       \code{"all"} entry for the pooled figures.}
#'     \item{\code{fcs.figure.actgate.dir}}{Named character vector of
#'       activation-gate figure directory paths, same structure as above.}
#'   }
#'
#' @seealso \code{\link{setup.autogating}}, \code{\link{run.autogating}},
#'   \code{\link{load.gate.definitions}}, \code{\link{do.population.gates}},
#'   \code{\link{do.activation.gates}}
#'
#' @export
define.gates <- function( fcs.setup, ap, plot.gates = TRUE ) {

    fcs.expr.data        <- fcs.setup$fcs.expr.data
    fcs.sample           <- fcs.setup$fcs.sample
    fcs.sample.gate.calc <- fcs.setup$fcs.sample.gate.calc
    fcs.event.sample     <- fcs.setup$fcs.event.sample


    # ==========================================================================
    # 1. Load gate definitions
    # ==========================================================================

    fcs.population.gates <- load.gate.definitions( ap$population.gates.filename )
    fcs.activation.gates <- load.gate.definitions( ap$activation.gates.filename )


    # ==========================================================================
    # 2. Create output directories
    # ==========================================================================

    fcs.figure.popgate.dir <- c(
        all = paste0( ap$figure.dir.basename, "_popgate_all" ),
        stats::setNames(
            sprintf( "%s_popgate_%s", ap$figure.dir.basename, fcs.sample ),
            fcs.sample
        )
    )

    fcs.figure.actgate.dir <- c(
        all = paste0( ap$figure.dir.basename, "_actgate_all" ),
        stats::setNames(
            sprintf( "%s_actgate_%s", ap$figure.dir.basename, fcs.sample ),
            fcs.sample
        )
    )

    for ( d in c( fcs.figure.popgate.dir, fcs.figure.actgate.dir,
                  ap$statistics.dir ) ) {
        if ( !file.exists( d ) )
            dir.create( d, recursive = TRUE )
    }


    # ==========================================================================
    # 3. Calculate boundaries on pooled gate-calculation samples
    # ==========================================================================

    fcs.gates.all.data <- fcs.expr.data[
        fcs.event.sample %in% fcs.sample.gate.calc, ]

    fcs.gates.all <- list()

    fcs.gates.all <- do.population.gates(
        fcs.gates.all,
        fcs.gates.all.data,
        fcs.population.gates,
        ap,
        calculate.gates = TRUE,
        figure.dir      = fcs.figure.popgate.dir[ "all" ],
        plot.gates      = plot.gates
    )

    fcs.gates.all <- do.activation.gates(
        fcs.gates.all,
        fcs.gates.all.data,
        fcs.activation.gates,
        ap,
        calculate.gates = TRUE,
        figure.dir      = fcs.figure.actgate.dir[ "all" ],
        plot.gates      = plot.gates
    )


    # ==========================================================================
    # Return
    # ==========================================================================

    list(
        fcs.gates.all          = fcs.gates.all,
        fcs.population.gates   = fcs.population.gates,
        fcs.activation.gates   = fcs.activation.gates,
        fcs.figure.popgate.dir = fcs.figure.popgate.dir,
        fcs.figure.actgate.dir = fcs.figure.actgate.dir
    )
}
