# get_autogating_param.r
#
# Public entry point: get.autogating.param()
# Internal helper:    get.autogating.param.minimal()


#' Build the complete autogating parameter list for a given cytometer
#'
#' @title Build Autogating Parameter List
#'
#' @description
#' Public entry point for the autogating parameter system.  Calls
#' \code{get.autogating.param.minimal()} to establish universal defaults, then
#' dispatches to the instrument-specific override function
#' (e.g. \code{get.autogating.param.symphony()}) identified by the
#' \code{cytometer} argument, and returns the completed list.
#'
#' @param cytometer Character string identifying the cytometer (e.g.
#'   \code{"symphony"}, \code{"aurora"}, \code{"id7000"}).  Must match a
#'   loaded \code{get.autogating.param.<cytometer>()} function
#'   (matched case-insensitively via \code{get0()}).
#' @param data.dir Character.  Path to the directory containing FCS files.
#' @param panel.filename Character.  Path to the panel CSV file with columns
#'   \code{dye}, \code{antigen}, and \code{type}.
#' @param experiment.filename Character.  Path to the experiment CSV file with
#'   columns \code{filename}, \code{sample}, and \code{gate.calculation}.
#' @param compensation.filename Character or \code{NULL}.  Path to a square
#'   compensation matrix CSV.  \code{NULL} disables spectral compensation.
#' @param transformation.filename Character or \code{NULL}.  Path to a
#'   per-experiment transformation override CSV keyed by \code{dye}, with
#'   either Gating-ML columns (\code{id}, \code{scale.t}, \code{linear.a},
#'   \code{logicle.w}, \code{logicle.m}, \code{logicle.a}, \code{log.m}) or
#'   FlowJo biex columns (\code{max.range}, \code{pos}, \code{neg},
#'   \code{width}); see \code{\link{read.transformation}}.  \code{NULL}
#'   uses the cytometer defaults in \code{ap$transform.defaults}.
#' @param population.gates.filename Character or \code{NULL}.  Path to an R
#'   source file or RDS file containing the population gate definition list.
#'   \code{NULL} uses an empty list (no population gating).
#' @param activation.gates.filename Character or \code{NULL}.  Path to an R
#'   source file or RDS file containing the activation gate definition list.
#'   \code{NULL} uses an empty list (no activation gating).
#' @param model.filename Character or \code{NULL}.  Path to a gating model
#'   (\code{.agmodel}, e.g. saved by Honeychrome).  When set, it is used
#'   instead of the gate definition files; see \code{\link{define.gates}}.
#' @param events.per.sample Integer.  Events drawn from each training
#'   sample for gate calculation (0 = all).  Default \code{50000L}, as in
#'   Honeychrome.
#' @param figure.dir Character.  Directory of the gate figures.  Default
#'   \code{"figures"}.
#' @param figure.training Logical.  Draw each gate on the pooled training
#'   events in \code{\link{define.gates}} (\code{figure.dir/training}).
#'   Default \code{TRUE}.
#' @param figure.samples Logical.  Draw each gate on each sample in
#'   \code{\link{run.autogating}} (\code{figure.dir/<sample>}); one PNG
#'   per gate and sample.  Default \code{FALSE}.
#' @param statistics.dir Character.  Directory where output CSV statistics are
#'   written.  Default \code{"statistics"}.
#' @param population.statistics.filename Character.  Filename (within
#'   \code{statistics.dir}) for population statistics output.
#'   Default \code{"population_statistics.csv"}.
#' @param median.statistics.filename Character.  Filename (within
#'   \code{statistics.dir}) for the median of each channel in each
#'   population.  Default \code{"population_medians.csv"}.
#' @param status.filename Character.  Filename (within
#'   \code{statistics.dir}) for the per-sample, per-gate gating status.
#'   Default \code{"gating_status.csv"}.
#' @param seed Integer.  Random seed of the training subsample.
#'   Default \code{20180101L}.
#' @param worker.n Integer.  Number of samples gated in parallel (via
#'   \code{AutoSpectral::create.parallel.lapply}); 1 gates them one after
#'   another.  Default \code{1L}.
#'
#' @return A named list (the \emph{autogating parameter object}, \code{ap})
#'   containing all pipeline configuration: file paths, the singlet scatter
#'   pair, transformation defaults, figure settings and output file names.
#'
#' @seealso \code{\link{get.autogating.param.symphony}},
#'   \code{\link{get.autogating.param.aurora}},
#'   \code{\link{get.autogating.param.id7000}},
#'   \code{\link{define.gates}}, \code{\link{run.autogating}}
#'
#' @export

get.autogating.param <- function(
    cytometer,
    data.dir,
    panel.filename,
    experiment.filename,
    compensation.filename          = NULL,
    transformation.filename        = NULL,
    population.gates.filename      = NULL,
    activation.gates.filename      = NULL,
    model.filename                 = NULL,
    events.per.sample              = 50000L,
    figure.dir                     = "figures",
    figure.training                = TRUE,
    figure.samples                 = FALSE,
    statistics.dir                 = "statistics",
    population.statistics.filename = "population_statistics.csv",
    median.statistics.filename     = "population_medians.csv",
    status.filename                = "gating_status.csv",
    seed                           = 20180101L,
    worker.n                       = 1L
) {
    ap <- get.autogating.param.minimal()

    # ---- file paths -----------------------------------------------------------

    ap$data.dir                        <- data.dir
    ap$panel.filename                  <- panel.filename
    ap$experiment.filename             <- experiment.filename
    ap$compensation.filename           <- compensation.filename
    ap$transformation.filename         <- transformation.filename
    ap$population.gates.filename       <- population.gates.filename
    ap$activation.gates.filename       <- activation.gates.filename
    ap$model.filename                  <- model.filename
    ap$events.per.sample               <- as.integer( events.per.sample )
    ap$figure.dir                      <- figure.dir
    ap$figure.training                 <- isTRUE( figure.training )
    ap$figure.samples                  <- isTRUE( figure.samples )
    ap$statistics.dir                  <- statistics.dir
    ap$population.statistics.filename  <- population.statistics.filename
    ap$median.statistics.filename      <- median.statistics.filename
    ap$status.filename                 <- status.filename

    # ---- execution ------------------------------------------------------------

    ap$seed     <- seed
    ap$worker.n <- as.integer( worker.n )

    # ---- cytometer-specific overrides -----------------------------------------

    fn <- get0( sprintf( "get.autogating.param.%s", tolower( cytometer ) ) )
    if ( is.null( fn ) )
        stop( sprintf( "No parameter function found for cytometer '%s'.", cytometer ) )

    fn( ap )
}


# ---------------------------------------------------------------------------
# Internal: universal defaults
# ---------------------------------------------------------------------------

#' Build the universal minimal autogating parameter defaults
#'
#' @title Universal Minimal Autogating Defaults
#'
#' @description
#' Internal constructor that returns a named list of cytometer-agnostic
#' defaults: the singlet scatter pair, gate-numbering width, figure
#' settings, transformation defaults and figure axis limits.  Every
#' cytometer-specific parameter function receives this list as its starting
#' point and overrides the fields that differ for that instrument.
#'
#' @return A named list with the following top-level elements (among others):
#'   \describe{
#'     \item{\code{singlet.scatter.parameter}}{Character vector
#'       \code{c(area, height)} of forward scatter, used to recognise the
#'       singlet gate of an imported FlowJo workspace.  \code{NULL} until
#'       overridden.}
#'     \item{\code{transform.defaults}}{Named list with elements
#'       \code{fluorescence} and \code{scatter}, each a Gating-ML transform
#'       parameter list (see \code{\link{make.transform}}).}
#'     \item{\code{figure.limits}}{Raw-unit axis limits of gate figures:
#'       \code{scatter.min}/\code{scatter.max} for channels of panel type
#'       "scatter", \code{expr.min}/\code{expr.max} for fluorescence.}
#'     \item{\code{figure.*}}{Size, resolution, point and label sizes of
#'       gate figures (see \code{\link{gate.figure}}).}
#'   }
#'
#' @keywords internal
get.autogating.param.minimal <- function() {

    list(
        # ---- singlet scatter pair (set per cytometer) --------------------------
        # c( area, height ) of forward scatter, used to recognise the singlet
        # gate of an imported FlowJo workspace.

        singlet.scatter.parameter = NULL,

        # ---- gate numbering ---------------------------------------------------

        gate.number.width = 3L,       # zero-padded width for figure filenames

        # ---- figures ----------------------------------------------------------

        figure.max.events   = 100000L,  # events drawn per gate figure
        figure.width        = 6,        # inches
        figure.height       = 6,
        figure.resolution   = 150,      # dots per inch
        figure.point.size   = 1,        # scattermore point size (pixels)
        figure.label.size   = 3.5,      # population labels (mm)
        figure.label.alpha  = 0.7,      # label background opacity
        figure.text.size    = 11,       # axis text and titles (pt)
        density.grid.n      = 128L,     # density estimate grid per axis
        bandwidth.factor    = 1,        # multiplies the normal-reference bandwidth

        # ---- transformation defaults (overridden per cytometer) ---------------
        # Gating-ML transform parameters (see transform.R) for fluorescence
        # and scatter channels, used unless a per-experiment transformation
        # CSV overrides a channel.

        transform.defaults = list(
            fluorescence = list( id = 1L, scale.t = 262144,
                                 logicle.w = 1, logicle.m = 4.418539922,
                                 logicle.a = 0 ),
            scatter      = list( id = 0L, scale.t = 262144, linear.a = 0 )
        ),

        # ---- figure axis limits (overridden per cytometer) --------------------

        figure.limits = list(
            scatter.min = 0,
            scatter.max = 265000,
            expr.min    = -2000,
            expr.max    = 265000
        )
    )
}
