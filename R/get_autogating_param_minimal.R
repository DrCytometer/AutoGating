# get.autogating.param.minimal.r

#' @title get.autogating.param.minimal
#'
#' @description
#' Get a list of parameter essential for performing autogating.
#' To be called via get.autogating.param().
#' @return The AutoGating parameter list.
#' @export
#'
#' @examples
#' agp <- get.autogating.param( cytometer = "Aurora" )
#'


get.autogating.param.minimal <- function() {

  list(

    # inputs

    param.dir = "./Parameter",
    fcs.data.dir = "./Data",

    fcs.compensation.filename = "compensation.csv",
    fcs.transformation.filename = "transformation.csv",
    fcs.experiment.filename = "experiment.csv",
    fcs.panel.filename = "panel.csv",

    fcs.population.gates.definition.filename = "population_gates.csv",
    fcs.activation.gates.definition.filename = "activation_gates.csv",

    # output directories

    fcs.figure.dir.basename = "figure",
    fcs.figure.popgate.dir = NULL,
    fcs.figure.actgate.dir = NULL,
    fcs.statistics.dir = "statistics",
    fcs.population.statistics.filename = "population_statistics.csv",
    fcs.activation.statistics.filename = "activation_statistics.csv",

    # labels
    pop.gate.label = "%s_popgate_%s",
    pop.gate.all.label = "_popgate_all",
    act.gate.label = "%s_actgate_%s",
    act.gate.all.label = "_actgate_all",

    # figure parameters
    # set configuration parameters for graphics

    # fcs.histogram.bin.n = 256,
    fcs.bandwidth.factor = 20,
    fcs.density.grid.n = 200,
    fcs.density.barheight = 30,

    fcs.palette.n = 1000,
    fcs.palette.base.n = 1000000,

    fcs.figure.label.size = 4.5,
    fcs.figure.label.alpha = 0.7,
    fcs.figure.label.width.factor = 0.010,
    fcs.figure.label.height.factor = 0.026,

    fcs.figure.width = 7.5,
    fcs.figure.height = 7,
    fcs.figure.margin.up = 5,
    fcs.figure.margin.ri = 2,
    fcs.figure.margin.do = 0,
    fcs.figure.margin.le = 10,
    fcs.figure.point.size = 0.1,
    fcs.figure.axis.label.size = 12,
    fcs.figure.axis.title.size = 12,

    # gating parameters
    # set configuration parameters for gates

    fcs.default.density.trimneg.x = FALSE,
    fcs.default.density.trimneg.y = FALSE,

    fcs.default.tail.probability = 0.99,

    fcs.gate.parameter.delimiter = ";",
    fcs.gate.parameter.delimiter.ex = ":",
    fcs.gate.parameter.ignore = "*",
    fcs.gate.parameter.operator = "=",

    fcs.gate.number.width = 3,

    fcs.activation.marker.selected = "X",
    fcs.activation.population.label = "activated",

    fcs.parent.string = "parent",
    fcs.parent.regexp = NULL,
    fcs.parent.regexp.sub = NULL,

    fcs.stats.separator = " :in: ",

    # parallel processing parameters

    parallel = TRUE,
    fcs.seed.base = 42,
    worker.process.n = parallelly::availableCores() - 1,
    max.memory.n = 4 * 1024^3,

    # cytometer parameters

    cytometer = NULL,
    scatter.data.min.x = NULL,
    scatter.data.max.x = NULL,
    scatter.data.min.y = NULL,
    scatter.data.max.y = NULL,
    expr.data.min = NULL,
    expr.data.max = NULL,
    default.scatter.parameter = NULL,
    default.time.parameter = "Time",
    default.transformation.param = NULL,
    data.step = NULL
  )

}
