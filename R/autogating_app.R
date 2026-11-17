# autogating_app.R
#
# Launches the AutoGating app (inst/shiny/app.R), an interactive editor for
# gating models: load an experiment, adjust the display transforms, build
# or edit the gate hierarchy on the pooled training events and see each
# gate as it is calculated, import gates from an .agmodel file, a FlowJo
# workspace or the earlier R gate lists, train, save the model for
# Honeychrome and run it on every sample.
#
# The app is evaluated in an environment whose parent is the package
# namespace, so it calls the same functions as the pipeline.  shiny, bslib,
# DT and shinyFiles (file and folder browsers) are needed only for the app
# (Suggests).


#' Open the AutoGating app
#'
#' @title AutoGating App
#'
#' @description
#' Opens the app in a browser.  Start it empty and load an experiment from
#' its files, or pass an experiment already read with
#' \code{\link{setup.autogating}} (and optionally a model) to edit it.
#' Gates are calculated on the pooled training events
#' (\code{ap$events.per.sample} events from each sample marked
#' \code{gate.calculation}), as \code{\link{define.gates}} does.
#'
#' The app has four tabs:
#' \describe{
#'   \item{Data}{Load an experiment (cytometer, FCS folder, panel,
#'     experiment, compensation and transformation files) and a gating
#'     model: empty, an \code{.agmodel} file, a FlowJo workspace (optionally
#'     converted to calculated gates, \code{\link{convert.imported.gates}})
#'     or the earlier R gate lists.}
#'   \item{Transforms}{Change the logicle, log or linear transform of a
#'     channel, previewed on the training events or any population;
#'     stored boundaries follow the change.  Save the transforms as a CSV
#'     (\code{\link{write.transformation.csv}}).}
#'   \item{Gates}{The gate hierarchy with each population's share of its
#'     parent; the figure of the selected gate; an editor for its type,
#'     parent, channels, populations, algorithm and parameters.  Clicking
#'     the figure sets a manual threshold.  Changes are recalculated at
#'     once and can be undone.}
#'   \item{Save and run}{Save the model (\code{\link{write.agmodel}}),
#'     write the training figures, run the model on every sample
#'     (\code{\link{run.autogating}}) and return to R.}
#' }
#'
#' @param fcs.setup Optional experiment returned by
#'   \code{\link{setup.autogating}}.
#' @param ap Optional autogating parameter object; with \code{fcs.setup}
#'   it sets the training events, figure and output settings.
#' @param model Optional gating model; without one, the files named in
#'   \code{ap} are used (as in \code{\link{define.gates}}), else the model
#'   starts empty.
#' @param launch.browser Passed to \code{\link[shiny]{runApp}}.
#'
#' @return When closed with "Finish", a list with the \code{model} and the
#'   \code{fcs.setup} (whose transforms may have changed); otherwise
#'   \code{NULL}.
#'
#' @export
autogating.app <- function( fcs.setup = NULL, ap = NULL, model = NULL,
                            launch.browser = interactive() ) {
    needed <- c( "shiny", "bslib", "DT", "shinyFiles" )
    missing <- needed[ !vapply( needed, requireNamespace, logical( 1 ), quietly = TRUE ) ]
    if ( length( missing ) )
        stop( sprintf( "The app needs the package(s) %s: install.packages(c(%s)).",
                       paste( missing, collapse = ", " ),
                       paste( sprintf( "\"%s\"", missing ), collapse = ", " ) ) )
    file <- system.file( "shiny", "app.R", package = "AutoGating" )
    if ( !nzchar( file ) )
        stop( "The app file inst/shiny/app.R is not installed." )
    env <- new.env( parent = environment( autogating.app ) )
    env$app.start <- list( fcs.setup = fcs.setup, ap = ap, model = model )
    app <- eval( parse( file, keep.source = FALSE ), envir = env )
    shiny::runApp( app, launch.browser = launch.browser )
}
