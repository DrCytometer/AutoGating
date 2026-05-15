# get_autogating_param.r
#
# Public entry point: get.autogating.param()
# Internal helpers:   get.autogating.param.minimal()
#                     merge.gate.param()


#' Build the complete autogating parameter list for a given cytometer
#'
#' @title Build Autogating Parameter List
#'
#' @description
#' Public entry point for the autogating parameter system.  Calls
#' \code{get.autogating.param.minimal()} to establish universal defaults, then
#' dispatches to the instrument-specific override function
#' (e.g. \code{get.autogating.param.symphony()}) identified by the
#' \code{cytometer} argument.  After the cytometer-specific function returns,
#' derived regular-expression fields are appended and the completed list is
#' returned.
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
#'   per-experiment transformation override CSV whose columns
#'   (\code{dye}, \code{length}, \code{max.range}, \code{pos}, \code{neg},
#'   \code{width}) override the cytometer defaults in
#'   \code{ap$transformation.defaults}.  \code{NULL} uses cytometer defaults
#'   for every channel.
#' @param population.gates.filename Character or \code{NULL}.  Path to an R
#'   source file or RDS file containing the population gate definition list.
#'   \code{NULL} uses an empty list (no population gating).
#' @param activation.gates.filename Character or \code{NULL}.  Path to an R
#'   source file or RDS file containing the activation gate definition list.
#'   \code{NULL} uses an empty list (no activation gating).
#' @param figure.dir.basename Character.  Base string used to construct
#'   per-gate figure directory names.  Default \code{"figure"}.
#' @param statistics.dir Character.  Directory where output CSV statistics are
#'   written.  Default \code{"statistics"}.
#' @param population.statistics.filename Character.  Filename (within
#'   \code{statistics.dir}) for population statistics output.
#'   Default \code{"population_statistics.csv"}.
#' @param activation.statistics.filename Character.  Filename (within
#'   \code{statistics.dir}) for activation statistics output.
#'   Default \code{"activation_statistics.csv"}.
#' @param seed Integer.  Random seed used for reproducible downsampling.
#'   Default \code{20180101L}.
#' @param worker.n Integer.  Number of parallel workers passed to
#'   \code{mclapply}.  Default \code{1L}.
#'
#' @return A named list (the \emph{autogating parameter object}, \code{ap})
#'   containing all pipeline configuration: file paths, cytometer-specific
#'   ranges, gate algorithm defaults, transformation defaults, figure
#'   aesthetics, sentinel values, and derived regular expressions.
#'
#' @seealso \code{\link{get.autogating.param.symphony}},
#'   \code{\link{get.autogating.param.aurora}},
#'   \code{\link{get.autogating.param.id7000}},
#'   \code{\link{do.population.gates}}, \code{\link{do.activation.gates}}
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
    figure.dir.basename            = "figure",
    statistics.dir                 = "statistics",
    population.statistics.filename = "population_statistics.csv",
    activation.statistics.filename = "activation_statistics.csv",
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
    ap$figure.dir.basename             <- figure.dir.basename
    ap$statistics.dir                  <- statistics.dir
    ap$population.statistics.filename  <- population.statistics.filename
    ap$activation.statistics.filename  <- activation.statistics.filename

    # ---- execution ------------------------------------------------------------

    ap$seed     <- seed
    ap$worker.n <- as.integer( worker.n )

    # ---- cytometer-specific overrides -----------------------------------------

    fn <- get0( sprintf( "get.autogating.param.%s", tolower( cytometer ) ) )
    if ( is.null( fn ) )
        stop( sprintf( "No parameter function found for cytometer '%s'.", cytometer ) )

    ap <- fn( ap )

    # ---- derived regexp (requires fcs.parent.string to be set) ----------------

    ap$parent.regexp     <- sprintf( "^%s\\.[0-9]+$",    ap$parent.string )
    ap$parent.regexp.sub <- sprintf( "^%s\\.([0-9]+)$",  ap$parent.string )

    ap
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
#' defaults covering sentinel values, gate-numbering conventions, figure
#' aesthetics, gate algorithm defaults, transformation defaults, and figure
#' axis limits.  Every cytometer-specific parameter function receives this
#' list as its starting point and overrides the fields that differ for that
#' instrument.
#'
#' @return A named list with the following top-level elements (among others):
#'   \describe{
#'     \item{\code{singlet.scatter.parameter}}{Character vector
#'       \code{c(area, height)} identifying the scatter channels used for
#'       singlet discrimination.  \code{NULL} until overridden.}
#'     \item{\code{leukocyte.scatter.parameter}}{Character vector
#'       \code{c(forward, side)} for the leukocyte gate.  \code{NULL} until
#'       overridden.}
#'     \item{\code{ignore}}{Sentinel string (\code{"*"}) marking a population
#'       slot or gate field that should be skipped.}
#'     \item{\code{gate.defaults}}{Named list of per-gate algorithm defaults
#'       (quantile trimming, region fractions, axis break vectors, threshold
#'       algorithm names, and the default tail probability).}
#'     \item{\code{transformation.defaults}}{Named list of
#'       \code{flowjo_biexp} parameter defaults (\code{length},
#'       \code{max.range}, \code{pos}, \code{neg}, \code{width}).}
#'     \item{\code{figure.limits}}{Named list of scatter and expression axis
#'       limits for plot generation.}
#'   }
#'
#' @keywords internal
get.autogating.param.minimal <- function() {

    list(
      # ---- scatter channel pairs --------------------------------------------
      # Must be overridden by each cytometer-specific function.
      # singlet.scatter.parameter:   c( area.channel, height.channel )
      #   used by gate.singlet() — area vs height on the same detector
      # leukocyte.scatter.parameter: c( forward.channel, side.channel )
      #   used by the leukocyte gate — two orthogonal scatter detectors

      singlet.scatter.parameter   = NULL,
      leukocyte.scatter.parameter = NULL,

        # ---- sentinel / delimiter values used throughout the pipeline ---------

        ignore          = "*",        # population slot to skip
        parent.string   = "parent",   # prefix for relative parent references
        stats.separator = " :in: ",   # used when building statistic column names

        # ---- gate numbering ---------------------------------------------------

        gate.number.width = 3L,       # zero-padded width for figure filenames

        # ---- activation gating ------------------------------------------------

        activation.marker.selected   = "X",        # value meaning "use default alg"
        activation.population.label  = "activated",

        # ---- figure defaults --------------------------------------------------

        density.grid.n          = 100L,
        density.barheight       = 30,
        bandwidth.factor        = 1,
        palette.n               = 1000L,
        palette.base.n          = 1000000L,
        figure.label.alpha      = 0.7,
        figure.label.width.factor  = 0.010,
        figure.label.height.factor = 0.026,
        figure.width            = 7.5,
        figure.height           = 7.0,
        figure.margin           = c( up = 5, right = 2, down = 0, left = 10 ),
        figure.point.size       = 0.1,
        figure.axis.label.size  = 12,
        figure.axis.title.size  = 12,

        # ---- gate algorithm defaults (overridden per cytometer / per gate) ----

        gate.defaults = list(
            quantile.x         = 0,
            quantile.y         = 0.0001,
            sep.axis           = 2,
            limit.xmin         = -0.01,
            limit.xmax         =  1.01,
            limit.ymin         = -0.01,
            limit.ymax         =  1.01,
            region.xmin        = 0.01,
            region.xmax        = 0.99,
            region.ymin        = 0.01,
            region.ymax        = 0.99,
            region.min         = 0.01,
            region.max         = 0.99,
            density.quantile.x = 0,
            density.quantile.y = 0.0001,
            density.trimneg.x  = TRUE,
            density.trimneg.y  = FALSE,
            break.x            = c( -0.1, 0, 1, 2, 3, 4 ),
            break.x.factor     = 1e6,
            break.y            = c( -0.1, 0, 1, 2, 3, 4 ),
            break.y.factor     = 1e6,
            label.size         = 4,
            algorithm.1dsep    = "tail",
            algorithm.2dsep    = c( "bimodal", "bimodal" ),
            default.tail.probability = 0.99
        ),

        # ---- transformation defaults (overridden per cytometer) ---------------
        # These are applied to all fluorescence channels unless a per-experiment
        # transformation CSV provides channel-level overrides.

        transformation.defaults = list(
            length    = 256,
            max.range = 262144,
            pos       = 4.418539922,
            neg       = 0,
            width     = -100
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


# ---------------------------------------------------------------------------
# Helper: merge a gate's own param list with cytometer gate.defaults.
# Returns a flat list with every key present; gate-level values win.
# ---------------------------------------------------------------------------

#' Merge per-gate parameter overrides with cytometer defaults
#'
#' @title Merge Gate Parameters with Cytometer Defaults
#'
#' @description
#' Takes a gate definition's own \code{gate.param} list and applies each
#' named element on top of the cytometer-level \code{ap$gate.defaults},
#' returning a single flat list that contains every parameter key the gating
#' pipeline requires.  Gate-level values take precedence over cytometer
#' defaults.
#'
#' @param gate.param Named list of per-gate parameter overrides (may be
#'   empty).
#' @param ap Autogating parameter object (the \code{ap} list returned by
#'   \code{\link{get.autogating.param}}).
#'
#' @return A named list equivalent to \code{ap$gate.defaults} with any
#'   entries in \code{gate.param} overwritten.
#'
#' @seealso \code{\link{resolve.gate.param}}
#'
#' @keywords internal
merge.gate.param <- function( gate.param, ap ) {
    defaults <- ap$gate.defaults
    for ( key in names( gate.param ) )
        defaults[[ key ]] <- gate.param[[ key ]]
    defaults
}
