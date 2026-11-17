# define_gates.R
#
# Builds the gating model of an experiment: gate definitions from a
# Honeychrome .agmodel file or from the earlier R gate files (converted by
# convert.gate.definitions), with boundaries trained on an equal-sized
# subsample of each training sample, as Honeychrome's Gating Model Builder
# does.


#' Pooled training events
#'
#' @title Training Events
#'
#' @description
#' Draws at most \code{events.per.sample} events from each training sample
#' (all of a smaller sample), so that every sample weighs the same in the
#' pooled data.  The caller's random stream is restored afterwards.
#'
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}.
#' @param samples Training sample names.
#' @param events.per.sample Events per sample; 0 keeps all.
#' @param seed Seed of the subsample.
#'
#' @return List with \code{rows} (row indices of
#'   \code{fcs.setup$fcs.expr.data}) and \code{per.sample} (per sample:
#'   \code{n_events_used}, \code{n_events_file}, \code{n_events_kept}).
#'
#' @keywords internal
training.events <- function( fcs.setup, samples, events.per.sample, seed ) {
    rows.by.sample <- split( seq_len( fcs.setup$fcs.event.n ), fcs.setup$fcs.event.sample )
    draw <- function() {
        rows <- list()
        per.sample <- list()
        for ( s in samples ) {
            r <- rows.by.sample[[ s ]]
            n <- length( r )
            if ( events.per.sample > 0L && n > events.per.sample )
                r <- sort( r[ sample.int( n, events.per.sample ) ] )
            rows[[ s ]] <- r
            per.sample[[ s ]] <- list( n_events_used = length( r ), n_events_file = n,
                                       n_events_kept = n )
        }
        list( rows = rows, per.sample = per.sample )
    }
    drawn <- seeded.call( seed, draw )
    rows <- drawn$rows
    per.sample <- drawn$per.sample
    list( rows = unlist( rows, use.names = FALSE ), per.sample = per.sample )
}


#' Gate channels missing from the data
#'
#' @title Missing Gate Channels
#'
#' @param gate.defs List of gate definitions.
#' @param channels Data column names.
#'
#' @return Character vector of the channels gates need that are not in
#'   \code{channels}.
#'
#' @keywords internal
missing.gate.channels <- function( gate.defs, channels ) {
    by.name <- gates.by.name( gate.defs )
    used <- unlist( lapply( gate.defs, function( g ) gate.channels( effective.gate.def( g, by.name ) ) ) )
    setdiff( unique( used ), channels )
}


#' Gating model from the files named in the parameters
#'
#' @title Model from Files
#'
#' @description
#' The \code{.agmodel} file \code{ap$model.filename}, else the R gate files
#' \code{ap$population.gates.filename} and \code{ap$activation.gates.filename}
#' converted by \code{\link{convert.gate.definitions}} (its messages are
#' shown).
#'
#' @param ap Autogating parameter object.
#' @param panel Panel \code{data.frame}, for the channel names of an
#'   \code{.agmodel}.
#'
#' @return Gating model; stops when neither is set.
#'
#' @keywords internal
model.from.files <- function( ap, panel ) {
    if ( !is.null( ap$model.filename ) )
        return( read.agmodel( ap$model.filename, panel ) )
    up <- convert.gate.definitions(
        load.gate.definitions( ap$population.gates.filename ),
        load.gate.definitions( ap$activation.gates.filename ) )
    for ( msg in up$messages )
        message( msg )
    if ( length( up$gate.definitions ) == 0L )
        stop( "No gate definitions: set ap$model.filename or ap$population.gates.filename." )
    new.agmodel( up$gate.definitions )
}


#' Fit a gating model to an experiment's data
#'
#' @title Prepare Model
#'
#' @description
#' Converts stored boundaries to the data's axis transforms
#' (\code{\link{remap.model}}), sets the model's transforms and channel map
#' to the data's, and gives every gate an application mode.  Stops when a
#' gate needs a channel the data do not have.
#'
#' @param model Gating model.
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}.
#'
#' @return List with \code{model} and \code{messages}.
#'
#' @keywords internal
prepare.model <- function( model, fcs.setup ) {
    transforms <- fcs.setup$ap$transform.params
    transforms <- transforms[ lengths( transforms ) > 0L ]
    remapped <- remap.model( model, transforms )
    model <- remapped$model
    model$transforms <- transforms
    for ( ch in setdiff( names( transforms ), names( model$channel.map ) ) )
        model$channel.map[[ ch ]] <- ch
    gate.defs <- model$gate.definitions
    missing <- missing.gate.channels( gate.defs, colnames( fcs.setup$fcs.expr.data ) )
    if ( length( missing ) )
        stop( sprintf( "Gate channels not in the data: %s", paste( missing, collapse = ", " ) ) )
    gate.names <- vapply( gate.defs, function( g ) g$gate.name, character( 1 ) )
    for ( g in setdiff( gate.names, names( model$application$mode.per.gate ) ) )
        model$application$mode.per.gate[[ g ]] <- model$application$default.mode %||% "fixed"
    list( model = model, messages = remapped$messages )
}


#' Training samples and events of an experiment
#'
#' @title Training Set
#'
#' @description
#' The samples marked \code{gate.calculation} in the experiment file (every
#' sample when none is marked) and \code{ap$events.per.sample} events drawn
#' from each (\code{\link{training.events}}, seed \code{ap$seed}).
#'
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}.
#' @param ap Autogating parameter object.
#'
#' @return List with samples, events.per.sample, seed, rows, per.sample,
#'   data (the pooled events) and messages.
#'
#' @keywords internal
training.set <- function( fcs.setup, ap ) {
    messages <- character( 0 )
    samples <- fcs.setup$fcs.sample.gate.calc
    if ( length( samples ) == 0L ) {
        samples <- fcs.setup$fcs.sample
        messages <- "No samples are marked for gate calculation; training on every sample."
    }
    events.per.sample <- as.integer( ap$events.per.sample %||% 50000L )
    tr <- training.events( fcs.setup, samples, events.per.sample, ap$seed )
    list( samples = samples, events.per.sample = events.per.sample, seed = ap$seed,
          rows = tr$rows, per.sample = tr$per.sample,
          data = fcs.setup$fcs.expr.data[ tr$rows, , drop = FALSE ], messages = messages )
}


#' Store a training run in a model
#'
#' @title Record Training
#'
#' @param model Gating model.
#' @param result Gating result on the training events
#'   (\code{\link{train.boundaries}}).
#' @param ts Training set (\code{\link{training.set}}).
#'
#' @return The model with the trained boundaries and its training summary
#'   (samples, events per sample, seed, status of each gate).
#'
#' @keywords internal
record.training <- function( model, result, ts ) {
    model$trained.boundaries <- result$boundaries
    model$training.summary <- list(
        training_sample_names = as.list( ts$samples ),
        n_events_per_sample   = lapply( ts$per.sample, `[[`, "n_events_used" ),
        events_cap_per_sample = ts$events.per.sample,
        seed                  = ts$seed,
        n_events_pooled       = length( ts$rows ),
        per_sample            = ts$per.sample,
        gate_status           = lapply( result$status, function( s )
            list( source = s$source, drift = s$drift, flag = s$flag, message = s$message ) ),
        complete              = TRUE )
    model
}


#' Build and train the gating model of an experiment
#'
#' @title Define and Train Gates
#'
#' @description
#' Builds the gating model (see \code{\link{read.agmodel}}) and trains its
#' boundaries:
#' \enumerate{
#'   \item The gate definitions come from \code{model}, else from the
#'     \code{.agmodel} file \code{ap$model.filename}, else from the R gate
#'     files \code{ap$population.gates.filename} and
#'     \code{ap$activation.gates.filename}, converted by
#'     \code{\link{convert.gate.definitions}}.
#'   \item Stored boundaries are converted to the data's axis transforms
#'     (\code{\link{remap.model}}).
#'   \item Gates without a stored boundary (every gate when
#'     \code{retrain = TRUE}) are calculated on the pooled training events:
#'     \code{ap$events.per.sample} events drawn from each sample marked
#'     \code{gate.calculation} in the experiment file (every sample when none
#'     is marked).  Gates with a stored boundary are applied as they are.
#' }
#' The model's transforms are set to the data's, and its training summary
#' records the samples and events used.  With \code{ap$figure.training},
#' each gate is drawn on the pooled training events in
#' \code{ap$figure.dir/training} (\code{\link{write.gate.figures}}).  Save
#' the model for Honeychrome with \code{\link{write.agmodel}}, or edit it in
#' \code{\link{autogating.app}}.
#'
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}.
#' @param ap Autogating parameter object returned by
#'   \code{\link{get.autogating.param}}.
#' @param model Optional gating model to train or reuse.
#' @param retrain Logical.  Recalculate every gate, including those with a
#'   stored boundary.  Default \code{FALSE}.
#'
#' @return Gating model with trained boundaries.
#'
#' @seealso \code{\link{run.autogating}}, \code{\link{write.agmodel}},
#'   \code{\link{convert.gate.definitions}}, \code{\link{train.boundaries}}
#'
#' @export
define.gates <- function( fcs.setup, ap, model = NULL, retrain = FALSE ) {

    # ==========================================================================
    # 1. Gate definitions
    # ==========================================================================

    if ( is.null( model ) )
        model <- model.from.files( ap, fcs.setup$fcs.panel )

    prepared <- prepare.model( model, fcs.setup )
    model <- prepared$model
    for ( msg in prepared$messages )
        message( msg )
    gate.defs <- model$gate.definitions
    gate.names <- vapply( gate.defs, function( g ) g$gate.name, character( 1 ) )


    # ==========================================================================
    # 2. Train gates without a stored boundary
    # ==========================================================================

    done <- if ( retrain ) list() else
        model$trained.boundaries[ intersect( names( model$trained.boundaries ), gate.names ) ]
    to.train <- setdiff( gate.names, names( done ) )
    if ( length( to.train ) == 0L ) {
        message( "Every gate has a stored boundary; nothing to train." )
        return( model )
    }

    ts <- training.set( fcs.setup, ap )
    for ( msg in ts$messages )
        message( msg )
    message( sprintf( "Training %d gate(s) on %d events from %d sample(s).",
                      length( to.train ), length( ts$rows ), length( ts$samples ) ) )
    result <- train.boundaries( gate.defs, ts$data, done = done, transforms = model$transforms )
    for ( g in names( result$status ) ) {
        s <- result$status[[ g ]]
        if ( isTRUE( s$flag %in% c( "failed", "skipped" ) ) )
            message( sprintf( "Gate '%s' not trained: %s", g, s$message ) )
    }
    model <- record.training( model, result, ts )

    if ( isTRUE( ap$figure.training ) ) {
        dir <- file.path( ap$figure.dir %||% "figures", "training" )
        files <- write.gate.figures( gate.defs, result, ts$data, fcs.setup, ap, dir,
                                     title = "Training" )
        message( sprintf( "%d gate figure(s) written to %s.", length( files ), dir ) )
    }
    model
}
