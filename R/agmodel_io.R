# agmodel_io.R
#
# Reads and writes Honeychrome Automated Gating models (.agmodel files,
# JSON, format 2), so that gating models move between AutoGating and
# Honeychrome.
#
# In R a model is a named list mirroring the file:
#   format.version      2
#   metadata            as stored (JSON key names)
#   gate.definitions    list of gate definitions (see gate_model.R)
#   trained.boundaries  gate -> population -> boundary entry (see
#                       gate_calculation.R); polygons are n x 2 matrices and
#                       missing thresholds NA
#   channel.map         channel -> antigen
#   transforms          channel -> transform parameters (see transform.R);
#                       the transforms the boundaries are in
#   training.summary    as stored (JSON key names)
#   af, qc              as stored (JSON key names)
#   application         mode.per.gate (gate -> "fixed" or "recalculate"),
#                       default.mode, drift.limit, drift.action
#
# Key names are snake_case in the file and use dots in R.  Names that are
# data (gate, population, channel and sample names) and all string values
# are kept exactly.  The metadata, training_summary, af and qc blocks and
# each gate's import record (source) are not used by AutoGating and are
# kept as read.
#
# Channels: Honeychrome names channels by the FCS $PnN keyword (the panel
# dye column), while AutoGating names data columns by antigen.  Given a
# panel, reading renames model channels to antigens and writing renames
# them back to dyes.


agmodel.format.version <- 2L


#' File key to R key
#' @title JSON Key to R Key
#' @param k Character vector of keys.
#' @return Keys with underscores replaced by dots.
#' @keywords internal
json.to.r.key <- function( k ) gsub( "_", ".", k, fixed = TRUE )


#' R key to file key
#' @title R Key to JSON Key
#' @param k Character vector of keys.
#' @return Keys with dots replaced by underscores.
#' @keywords internal
r.to.json.key <- function( k ) gsub( ".", "_", k, fixed = TRUE )


#' Rename the keys of nested lists
#' @title Rename Keys
#' @param x Any value.
#' @param fn Function applied to the names of every named list in \code{x}.
#' @return \code{x} with renamed keys (\code{NULL} elements kept).
#' @keywords internal
rename.keys <- function( x, fn ) {
    if ( !is.list( x ) )
        return( x )
    out <- lapply( x, rename.keys, fn = fn )
    if ( !is.null( names( x ) ) )
        names( out ) <- fn( names( x ) )
    out
}


#' A list written as a JSON object
#' @title JSON Object
#' @param x Named list, empty list or \code{NULL}.
#' @return \code{x}, or an empty named list (written as \code{\{\}}).
#' @keywords internal
json.object <- function( x ) {
    if ( length( x ) == 0L ) stats::setNames( list(), character( 0 ) ) else x
}


#' Boundary entry read from a model file
#' @title Boundary Entry from JSON
#' @param e Entry as parsed from JSON.
#' @return Entry with R keys, the polygon as a matrix and missing
#'   thresholds \code{NA}.
#' @keywords internal
entry.from.json <- function( e ) {
    e <- rename.keys( e, json.to.r.key )
    b <- e[[ "boundary" ]]
    if ( is.list( b ) )
        e[[ "boundary" ]] <- if ( length( b ) )
            matrix( as.numeric( unlist( b ) ), ncol = 2L, byrow = TRUE ) else NULL
    for ( key in intersect( c( "threshold.x", "threshold.y" ), names( e ) ) )
        e[ key ] <- list( if ( is.null( e[[ key ]] ) ) NA_real_ else as.numeric( e[[ key ]] ) )
    if ( is.list( e[[ "range" ]] ) )
        e[[ "range" ]] <- as.numeric( unlist( e[[ "range" ]] ) )
    if ( is.list( e[[ "shape.channels" ]] ) )
        e[[ "shape.channels" ]] <- as.character( unlist( e[[ "shape.channels" ]] ) )
    e
}


#' Boundary entry written to a model file
#' @title Boundary Entry to JSON
#' @param e Boundary entry.
#' @return Entry with file keys, ready for \code{jsonlite::toJSON}.
#' @keywords internal
entry.to.json <- function( e ) {
    if ( is.matrix( e[[ "boundary" ]] ) )
        e[[ "boundary" ]] <- unname( e[[ "boundary" ]] )
    if ( !is.null( e[[ "range" ]] ) )
        e[[ "range" ]] <- as.list( as.numeric( e[[ "range" ]] ) )
    if ( !is.null( e[[ "shape.channels" ]] ) )
        e[[ "shape.channels" ]] <- as.list( e[[ "shape.channels" ]] )
    rename.keys( e, r.to.json.key )
}


#' Gate definition read from a model file
#' @title Gate Definition from JSON
#' @param g Gate definition as parsed from JSON.
#' @return Gate definition with R keys.
#' @keywords internal
gate.from.json <- function( g ) {
    out <- list()
    for ( k in names( g ) ) {
        v <- g[[ k ]]
        out[ json.to.r.key( k ) ] <- list( switch( k,
            populations       = json.object( lapply( v, rename.keys, fn = json.to.r.key ) ),
            stats_parent      = json.object( v ),
            template          = ,
            template_original = json.object( lapply( v, entry.from.json ) ),
            source            = v,
            rename.keys( v, json.to.r.key ) ) )
    }
    out
}


#' Gate definition written to a model file
#' @title Gate Definition to JSON
#' @param g Gate definition.
#' @return Gate definition with file keys.
#' @keywords internal
gate.to.json <- function( g ) {
    out <- list()
    for ( k in names( g ) ) {
        v <- g[[ k ]]
        out[ r.to.json.key( k ) ] <- list( switch( k,
            populations       = json.object( lapply( v, rename.keys, fn = r.to.json.key ) ),
            stats.parent      = json.object( v ),
            gate.param        = json.object( rename.keys( v, r.to.json.key ) ),
            template          = ,
            template.original = json.object( lapply( v, entry.to.json ) ),
            source            = v,
            rename.keys( v, r.to.json.key ) ) )
    }
    for ( key in c( "populations", "gate_param", "stats_parent" ) )
        if ( is.null( out[[ key ]] ) )
            out[[ key ]] <- json.object( NULL )
    out
}


#' Channels a model refers to
#' @title Model Channels
#' @param model Gating model.
#' @return Character vector of channel names.
#' @keywords internal
model.channels <- function( model ) {
    shapes <- function( entries ) unlist( lapply( entries, function( e ) e[[ "shape.channels" ]] ) )
    ch <- c( unlist( lapply( model$gate.definitions, function( g )
                 c( referenced.channels( g ), shapes( g[[ "template.original" ]] ) ) ) ),
             unlist( lapply( model$trained.boundaries, shapes ) ),
             names( model$transforms ), names( model$channel.map ) )
    unique( ch[ vapply( ch, has.value, logical( 1 ) ) ] )
}


#' Rename the channels of a model
#' @title Rename Model Channels
#' @param model Gating model.
#' @param to Named character vector, old channel name -> new.
#' @return Model with renamed gate axes, shape channels, transforms and
#'   channel map.
#' @keywords internal
rename.model.channels <- function( model, to ) {
    rn <- function( ch ) {
        hit <- ch %in% names( to )
        ch[ hit ] <- to[ ch[ hit ] ]
        unname( ch )
    }
    rn.entries <- function( entries ) lapply( entries, function( e ) {
        if ( !is.null( e[[ "shape.channels" ]] ) )
            e[[ "shape.channels" ]] <- rn( e[[ "shape.channels" ]] )
        e
    } )
    model$gate.definitions <- lapply( model$gate.definitions, function( g ) {
        for ( key in c( "gate.marker.x", "gate.marker.y" ) )
            if ( has.value( g[[ key ]] ) )
                g[[ key ]] <- rn( g[[ key ]] )
        for ( key in c( "template", "template.original" ) )
            if ( length( g[[ key ]] ) )
                g[[ key ]] <- rn.entries( g[[ key ]] )
        g
    } )
    model$trained.boundaries <- lapply( model$trained.boundaries, rn.entries )
    if ( length( model$transforms ) )
        names( model$transforms ) <- rn( names( model$transforms ) )
    if ( length( model$channel.map ) )
        names( model$channel.map ) <- rn( names( model$channel.map ) )
    model
}


#' Assemble a gating model
#'
#' @title New Gating Model
#'
#' @description
#' Builds a format-2 gating model (see \code{\link{read.agmodel}}) with
#' every gate applied with its stored boundary ("fixed") by default.
#'
#' @param gate.definitions List of gate definitions.
#' @param trained.boundaries Named list gate -> population -> boundary
#'   entry.
#' @param transforms Named list channel -> transform parameters of the
#'   boundaries.
#' @param channel.map Named list channel -> antigen; by default each
#'   channel maps to itself.
#' @param training.summary Named list (JSON key names) describing the
#'   training data.
#' @param metadata Named list (JSON key names).
#' @param application Named list overriding the default application
#'   settings (\code{mode.per.gate}, \code{default.mode},
#'   \code{drift.limit}, \code{drift.action}).
#'
#' @return Gating model.
#'
#' @export
new.agmodel <- function( gate.definitions, trained.boundaries = list(),
                         transforms = list(), channel.map = NULL,
                         training.summary = list(), metadata = list(),
                         application = list() ) {
    gate.names <- vapply( gate.definitions, function( g ) g$gate.name, character( 1 ) )
    app <- list( mode.per.gate = stats::setNames( as.list( rep( "fixed", length( gate.names ) ) ),
                                                  gate.names ),
                 default.mode = "fixed", drift.limit = ag.defaults$drift.limit,
                 drift.action = "flag" )
    app[ names( application ) ] <- application
    transforms <- transforms[ lengths( transforms ) > 0L ]
    if ( is.null( channel.map ) )
        channel.map <- stats::setNames( as.list( names( transforms ) ), names( transforms ) )
    meta <- list( created = format( Sys.time(), "%Y-%m-%dT%H:%M:%S" ),
                  autogating_version = tryCatch( as.character( utils::packageVersion( "AutoGating" ) ),
                                                 error = function( e ) "" ),
                  n_gates = length( gate.definitions ),
                  n_populations = sum( vapply( gate.definitions,
                      function( g ) length( g$populations ), integer( 1 ) ) ) )
    meta[ names( metadata ) ] <- metadata
    model <- list(
        format.version     = agmodel.format.version,
        metadata           = meta,
        gate.definitions   = unname( gate.definitions ),
        trained.boundaries = trained.boundaries,
        channel.map        = channel.map,
        transforms         = transforms,
        training.summary   = training.summary,
        af                 = list( sample_af_profiles = json.object( NULL ) ),
        qc                 = list( enabled = FALSE, settings = json.object( NULL ),
                                   fingerprints = json.object( NULL ) ),
        application        = app
    )
    normalize.agmodel( model, agmodel.format.version )$model
}


#' Validate and complete a model read from a file
#' @title Normalise Model
#' @description
#' Fills the defaults Honeychrome fills (empty populations, parameters and
#' statistics parents, a "root" parent, per-gate application modes),
#' repairs statistics parents ending in "/", checks the hierarchy for
#' cycles and names untrained gates.
#' @param model Gating model.
#' @param version Format version of the file.
#' @return List with \code{model} and \code{messages}.
#' @keywords internal
normalize.agmodel <- function( model, version ) {
    messages <- character( 0 )
    repaired <- character( 0 )
    for ( i in seq_along( model$gate.definitions ) ) {
        g <- model$gate.definitions[[ i ]]
        if ( !has.value( g$gate.name ) )
            stop( "A gate definition has no gate_name." )
        for ( key in c( "populations", "gate.param", "stats.parent" ) )
            if ( is.null( g[[ key ]] ) )
                g[[ key ]] <- json.object( NULL )
        if ( !has.value( g$parent.gate ) ) {
            g[ "parent.gate" ] <- list( NULL )
            g$parent.popul <- "root"
        }
        for ( p in names( g$stats.parent ) ) {
            v <- g$stats.parent[[ p ]]
            if ( is.character( v ) && length( v ) == 1L && endsWith( v, "/" ) ) {
                g$stats.parent[[ p ]] <- sub( "/+$", "", v )
                repaired <- union( repaired, g$gate.name )
            }
        }
        model$gate.definitions[[ i ]] <- g
    }
    if ( length( repaired ) )
        messages <- c( messages, paste0(
            "Statistics parents ending in '/' named no population; they now name the gate itself: ",
            paste( repaired, collapse = ", " ) ) )
    ordered.gate.defs( model$gate.definitions )
    if ( version == 1L ) {
        messages <- c( messages, paste(
            "This model was saved before transform settings were recorded; gates are",
            "applied as stored, without checking the data's axis transforms." ) )
        model$transforms <- list()
    }
    app <- model$application
    app$default.mode <- app$default.mode %||% "fixed"
    app$drift.limit  <- app$drift.limit %||% ag.defaults$drift.limit
    app$drift.action <- app$drift.action %||% "flag"
    modes <- app$mode.per.gate %||% list()
    for ( g in model$gate.definitions )
        if ( !isTRUE( modes[[ g$gate.name ]] %in% c( "fixed", "recalculate" ) ) )
            modes[[ g$gate.name ]] <- app$default.mode
    app$mode.per.gate <- modes
    model$application <- app
    untrained <- setdiff( vapply( model$gate.definitions, function( g ) g$gate.name, character( 1 ) ),
                          names( model$trained.boundaries ) )
    if ( length( untrained ) )
        messages <- c( messages, paste0( length( untrained ), " untrained gate(s): ",
            paste( utils::head( untrained, 8L ), collapse = ", " ),
            if ( length( untrained ) > 8L ) sprintf( " and %d more", length( untrained ) - 8L ) ) )
    list( model = model, messages = messages )
}


#' Read a Honeychrome gating model
#'
#' @title Read .agmodel File
#'
#' @description
#' Reads a gating model saved by Honeychrome Automated Gating or by
#' \code{\link{write.agmodel}} (format 1 or 2) into the R structure
#' described in the file header of \code{agmodel_io.R}.  With a panel, the
#' model's channels are renamed to the panel's antigens: a channel matching
#' a panel dye takes its antigen, otherwise a channel whose
#' \code{channel_map} antigen is in the panel takes that antigen.
#' Information about the model (untrained gates, repairs, channels missing
#' from the panel) is reported with \code{message()}.
#'
#' @param path Path to the \code{.agmodel} file.
#' @param panel Optional panel \code{data.frame} with \code{dye} and
#'   \code{antigen} columns (e.g. \code{fcs.setup$fcs.panel}).
#'
#' @return Gating model (named list).
#'
#' @seealso \code{\link{write.agmodel}}, \code{\link{remap.model}}
#'
#' @export
read.agmodel <- function( path, panel = NULL ) {
    m <- jsonlite::read_json( path, simplifyVector = FALSE )
    version <- as.integer( m$format_version %||% NA )
    if ( !isTRUE( version %in% 1:2 ) )
        stop( sprintf( "Unsupported model format_version %s in '%s'.",
                       format( m$format_version ), path ) )
    if ( !is.list( m$gate_definitions ) )
        stop( sprintf( "Model '%s' has no gate_definitions list.", path ) )

    app <- m$application %||% list()
    model <- list(
        format.version     = agmodel.format.version,
        metadata           = json.object( m$metadata ),
        gate.definitions   = lapply( m$gate_definitions, gate.from.json ),
        trained.boundaries = lapply( json.object( m$trained_boundaries ),
                                     function( entries ) lapply( entries, entry.from.json ) ),
        channel.map        = json.object( m$channel_map ),
        transforms         = lapply( json.object( m$transforms ), function( p ) {
            p <- rename.keys( p, json.to.r.key )
            if ( is.list( p[[ "limits" ]] ) )
                p[[ "limits" ]] <- as.numeric( unlist( p[[ "limits" ]] ) )
            p
        } ),
        training.summary   = json.object( m$training_summary ),
        af                 = m$af %||% list( sample_af_profiles = json.object( NULL ) ),
        qc                 = m$qc %||% list( enabled = FALSE, settings = json.object( NULL ),
                                             fingerprints = json.object( NULL ) ),
        application        = list( mode.per.gate = json.object( app$mode_per_gate ),
                                   default.mode  = app$default_mode,
                                   drift.limit   = app$drift_limit,
                                   drift.action  = app$drift_action )
    )
    res <- normalize.agmodel( model, version )
    model <- res$model
    messages <- res$messages

    if ( !is.null( panel ) ) {
        chans <- model.channels( model )
        to <- vapply( chans, function( ch ) {
            if ( ch %in% panel$dye )
                return( panel$antigen[ match( ch, panel$dye ) ] )
            ant <- model$channel.map[[ ch ]]
            if ( has.value( ant ) && ant %in% panel$antigen )
                return( ant )
            ch
        }, character( 1 ) )
        model <- rename.model.channels( model, to )
        gate.chans <- unique( unlist( lapply( model$gate.definitions, referenced.channels ) ) )
        missing <- setdiff( gate.chans, panel$antigen )
        if ( length( missing ) )
            messages <- c( messages, paste( "Gate channels not in the panel:",
                                            paste( missing, collapse = ", " ) ) )
    }
    for ( msg in messages )
        message( msg )
    model
}


#' Write a gating model for Honeychrome
#'
#' @title Write .agmodel File
#'
#' @description
#' Writes a gating model as a format-2 \code{.agmodel} file that Honeychrome
#' Automated Gating opens.  With a panel, channels named by antigen are
#' written under their dye names (the FCS \code{$PnN} names Honeychrome
#' uses), and \code{channel_map} maps each dye to its antigen.  Numbers are
#' written with 17 significant digits, so they read back exactly.  The file
#' is written beside \code{path} and then moved into place.
#'
#' @param model Gating model (see \code{\link{read.agmodel}},
#'   \code{\link{new.agmodel}}).
#' @param path Output path.
#' @param panel Optional panel \code{data.frame} with \code{dye} and
#'   \code{antigen} columns.
#'
#' @return \code{path}, invisibly.
#'
#' @seealso \code{\link{read.agmodel}}
#'
#' @export
write.agmodel <- function( model, path, panel = NULL ) {
    if ( !is.null( panel ) ) {
        chans <- model.channels( model )
        to <- ifelse( chans %in% panel$antigen, panel$dye[ match( chans, panel$antigen ) ], chans )
        model <- rename.model.channels( model, stats::setNames( to, chans ) )
        for ( i in which( !( to %in% names( model$channel.map ) ) ) )
            model$channel.map[[ to[ i ] ]] <- chans[ i ]
    }
    app <- model$application
    transforms <- model$transforms[ lengths( model$transforms ) > 0L ]
    out <- list(
        format_version     = agmodel.format.version,
        metadata           = json.object( model$metadata ),
        gate_definitions   = unname( lapply( model$gate.definitions, gate.to.json ) ),
        trained_boundaries = json.object( lapply( model$trained.boundaries,
            function( entries ) json.object( lapply( entries, entry.to.json ) ) ) ),
        channel_map        = json.object( model$channel.map ),
        transforms         = json.object( lapply( transforms, function( p ) {
            if ( !is.null( p[[ "limits" ]] ) )
                p[[ "limits" ]] <- as.list( as.numeric( p[[ "limits" ]] ) )
            rename.keys( p, r.to.json.key )
        } ) ),
        training_summary   = json.object( model$training.summary ),
        af                 = model$af %||% list( sample_af_profiles = json.object( NULL ) ),
        qc                 = model$qc %||% list( enabled = FALSE, settings = json.object( NULL ),
                                                 fingerprints = json.object( NULL ) ),
        application        = list( mode_per_gate = json.object( app$mode.per.gate ),
                                   default_mode  = app$default.mode %||% "fixed",
                                   drift_limit   = app$drift.limit %||% ag.defaults$drift.limit,
                                   drift_action  = app$drift.action %||% "flag" )
    )
    txt <- jsonlite::toJSON( out, auto_unbox = TRUE, null = "null", na = "null",
                             digits = I( 17 ), pretty = TRUE )
    tmp <- paste0( path, ".tmp" )
    writeLines( txt, tmp, useBytes = TRUE )
    if ( !file.rename( tmp, path ) )
        stop( sprintf( "Could not write '%s'.", path ) )
    invisible( path )
}
