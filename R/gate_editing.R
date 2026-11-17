# gate_editing.R
#
# Editing a gating model gate by gate, as the gate editor of the AutoGating
# app does: a gate definition is shown as a flat set of form fields and
# built back from them; adding, renaming or removing a gate keeps every
# reference to it in the rest of the model consistent.
#
# Form fields (gate.form.values / gate.from.form):
#   name, type, parent ("gate/pop", "" for all events), x, y, sep.axis
#   algorithm             1dsep, singlets and free gates
#   algorithm.x/.y        2dsep axes
#   populations           named character: region -> population name.
#                         1dsep "below", "above"; 2dsep the four quadrants
#                         (quadrants with the same name form one population,
#                         an empty name leaves the quadrant out); singlets
#                         and free gates "inside"
#   axis.lo/.hi, tail.fraction, threshold.factor, fixed.threshold
#                         per axis "x" and "y" (1dsep: "x" only); NA leaves
#                         the engine default
#   replicate             per axis: "gate:<name>" copies that gate's
#                         threshold on the same channel; "population:<key>"
#                         calculates it on that population with
#                         replicate.algorithm
#   width.sd              singlets band half-width (robust SDs)
#   debris.algorithm, edge.density   cells gate
#   origin                origin gate of a "replicate" gate
#   mode                  "fixed" or "recalculate" when applied to samples
#
# Parameters of a 1dsep gate are stored without an axis suffix
# (axis.lo, tail.fraction), those of a 2dsep gate with one (axis.lo.x), as
# the engine reads them.


#' Population keys of a gate hierarchy
#'
#' @title Population Keys
#'
#' @param gate.defs List of gate definitions.
#'
#' @return Character vector of \code{"gate/pop"} keys in hierarchy order.
#'
#' @keywords internal
population.keys <- function( gate.defs ) {
    by.name <- gates.by.name( gate.defs )
    keys <- character( 0 )
    for ( g in ordered.gate.defs( gate.defs ) ) {
        eff <- effective.gate.def( g, by.name )
        keys <- c( keys, population.key( g$gate.name, names( eff$populations ) ) )
    }
    keys
}


#' Regions a gate type names populations by
#'
#' @title Form Regions
#'
#' @param gate.type Character gate type.
#'
#' @return Character vector of regions (empty for replicate gates).
#'
#' @keywords internal
form.regions <- function( gate.type ) {
    switch( gate.type,
            `1dsep` = c( "below", "above" ),
            `2dsep` = ag.quadrant.regions,
            singlets = "inside",
            free = "inside",
            character( 0 ) )
}


#' Form fields of a gate definition
#'
#' @title Gate Form Values
#'
#' @description
#' The gate as the flat fields described in the file header of
#' \code{gate_editing.R}.
#'
#' @param gate.def Gate definition, or \code{NULL} for a new 1dsep gate.
#' @param mode Application mode of the gate.
#'
#' @return Named list of form fields.
#'
#' @keywords internal
gate.form.values <- function( gate.def = NULL, mode = "fixed" ) {
    if ( is.null( gate.def ) )
        gate.def <- list( gate.name = "", gate.type = "1dsep", algorithm = "mixture",
                          populations = default.populations.for.type( "1dsep" ) )
    type <- gate.def$gate.type %||% "1dsep"
    gp <- gate.def$gate.param %||% list()
    two <- identical( type, "2dsep" )
    num <- function( key ) {
        v <- suppressWarnings( as.numeric( gp[[ key ]] ) )
        if ( length( v ) == 1L ) v else NA_real_
    }
    suffix <- function( axis ) if ( two ) paste0( ".", axis ) else ""
    regions <- population.regions( gate.def )
    names.by.region <- stats::setNames( rep( "", length( form.regions( type ) ) ), form.regions( type ) )
    for ( p in names( regions ) ) {
        r <- regions[[ p ]]
        parts <- if ( two ) region.quadrants( r ) else r
        names.by.region[ intersect( parts, names( names.by.region ) ) ] <- p
    }
    v <- list(
        name = gate.def$gate.name %||% "", type = type,
        parent = if ( has.value( gate.def$parent.gate ) )
            population.key( gate.def$parent.gate, gate.def$parent.popul %||% "" ) else "",
        x = gate.def$gate.marker.x %||% "", y = gate.def$gate.marker.y %||% "",
        sep.axis = if ( identical( gate.def$sep.axis, "y" ) ) "y" else "x",
        algorithm = gate.def$algorithm %||% "",
        algorithm.x = gp$algorithm.x %||% gate.def$algorithm %||% "",
        algorithm.y = gp$algorithm.y %||% gate.def$algorithm %||% "",
        populations = names.by.region,
        width.sd = num( "width.sd" ), debris.algorithm = gp$debris.algorithm %||% "",
        edge.density = num( "edge.density" ),
        origin = gate.def$origin.gate %||% "", mode = mode )
    for ( axis in c( "x", "y" ) ) {
        s <- suffix( axis )
        for ( key in c( "axis.lo", "axis.hi", "tail.fraction", "threshold.factor", "fixed.threshold" ) ) {
            stored <- if ( key == "fixed.threshold" ) paste0( key, ".", axis ) else paste0( key, s )
            v[[ paste0( key, ".", axis ) ]] <- if ( !two && axis == "y" ) NA_real_ else num( stored )
        }
        src <- gp[[ paste0( "replicate.", axis ) ]]
        v[[ paste0( "replicate.", axis ) ]] <-
            if ( is.list( src ) && has.value( src$gate ) ) paste0( "gate:", src$gate )
            else if ( is.list( src ) && has.value( src$population ) ) paste0( "population:", src$population )
            else ""
        v[[ paste0( "replicate.algorithm.", axis ) ]] <-
            if ( is.list( src ) && has.value( src$algorithm ) ) src$algorithm else "mixture"
    }
    v
}


#' Gate definition from form fields
#'
#' @title Gate from Form
#'
#' @description
#' Builds a gate definition from the fields described in the file header of
#' \code{gate_editing.R}.  Fields of \code{old} the form does not show
#' (number, statistics parents, drawn shapes, other parameters, population
#' flags) are kept.
#'
#' @param v Named list of form fields.
#' @param old Gate definition being edited, or \code{NULL}.
#'
#' @return Gate definition.
#'
#' @keywords internal
gate.from.form <- function( v, old = NULL ) {
    type <- v$type
    if ( !( type %in% ag.gate.types ) )
        stop( sprintf( "Unknown gate type '%s'.", type ) )
    name <- trimws( v$name %||% "" )
    if ( !nzchar( name ) || grepl( "/", name, fixed = TRUE ) )
        stop( "A gate needs a name without '/'." )
    g <- old %||% list()
    g$gate.name <- name
    g$gate.type <- type
    value <- function( x ) if ( has.value( x ) ) x else NULL
    g$gate.marker.x <- value( v$x )
    g$gate.marker.y <- value( v$y )
    g$sep.axis <- if ( identical( type, "1dsep" ) && has.value( v$y ) ) v$sep.axis %||% "x" else NULL
    parent <- if ( has.value( v$parent ) ) parse.population.key( v$parent ) else NULL
    g$parent.gate <- parent[ 1L ]
    g$parent.popul <- if ( is.null( parent ) ) "root" else parent[ 2L ]
    g$origin.gate <- if ( type == "replicate" ) value( v$origin ) else NULL
    if ( type == "replicate" && is.null( g$origin.gate ) )
        stop( "A replicate gate needs an origin gate." )

    gp <- g$gate.param %||% list()
    if ( !identical( old$gate.type, type ) ) gp <- list()
    two <- type == "2dsep"
    set <- function( key, x ) {
        x <- suppressWarnings( as.numeric( x ) )
        gp[[ key ]] <<- if ( length( x ) == 1L && is.finite( x ) ) x else NULL
    }
    if ( type %in% ag.threshold.types ) {
        for ( axis in if ( two ) c( "x", "y" ) else "x" ) {
            s <- if ( two ) paste0( ".", axis ) else ""
            for ( key in c( "axis.lo", "axis.hi", "tail.fraction", "threshold.factor" ) )
                set( paste0( key, s ), v[[ paste0( key, ".", axis ) ]] )
            set( paste0( "fixed.threshold.", axis ), v[[ paste0( "fixed.threshold.", axis ) ]] )
            algo <- if ( two ) v[[ paste0( "algorithm.", axis ) ]] else v$algorithm
            src <- v[[ paste0( "replicate.", axis ) ]] %||% ""
            gp[[ paste0( "replicate.", axis ) ]] <- if ( identical( algo, "replicate" ) ) {
                if ( startsWith( src, "gate:" ) ) list( gate = substring( src, 6L ) )
                else if ( startsWith( src, "population:" ) )
                    list( population = substring( src, 12L ),
                          algorithm = v[[ paste0( "replicate.algorithm.", axis ) ]] %||% "mixture" )
                else stop( sprintf( "Axis %s replicates a threshold: choose its source.", axis ) )
            } else NULL
        }
        if ( two ) {
            gp$algorithm.x <- v$algorithm.x
            gp$algorithm.y <- v$algorithm.y
            g$algorithm <- v$algorithm.x
        } else {
            g$algorithm <- v$algorithm
        }
    } else if ( type == "singlets" ) {
        g$algorithm <- v$algorithm %||% "singlets"
        set( "width.sd", v$width.sd )
    } else if ( type == "free" ) {
        g$algorithm <- v$algorithm %||% "cells"
        gp$debris.algorithm <- value( v$debris.algorithm )
        set( "edge.density", v$edge.density )
    } else {
        g$algorithm <- "replicate"
    }
    allowed <- ag.algorithms.by.type[[ type ]]
    for ( a in c( g$algorithm, if ( two ) c( gp$algorithm.x, gp$algorithm.y ) ) )
        if ( !( a %in% allowed ) )
            stop( sprintf( "Algorithm '%s' is not available for %s gates.", a, type ) )
    g$gate.param <- gp

    regions <- form.regions( type )
    if ( length( regions ) ) {
        pop.names <- trimws( v$populations[ regions ] )
        pop.names[ is.na( pop.names ) ] <- ""
        names( pop.names ) <- regions
        if ( any( grepl( "/", pop.names, fixed = TRUE ) ) )
            stop( "Population names cannot contain '/'." )
        used <- unique( pop.names[ nzchar( pop.names ) ] )
        if ( length( used ) == 0L )
            stop( "Name at least one population." )
        if ( type == "1dsep" && length( used ) < sum( nzchar( pop.names ) ) )
            stop( "The two sides of a 1dsep gate need different names." )
        pops <- list()
        for ( p in used ) {
            region <- if ( two ) join.quadrants( regions[ pop.names == p ] ) else regions[ pop.names == p ]
            keep <- old$populations[[ p ]] %||% list()
            keep$label <- p
            keep$region <- region
            pops[[ p ]] <- keep
        }
        g$populations <- pops
    } else {
        g$populations <- NULL
    }
    g
}


#' Add or replace a gate in a model
#'
#' @title Put Gate
#'
#' @description
#' Adds \code{gate.def}, or replaces the gate \code{old.name}.  A renamed
#' gate or population is renamed wherever it is referred to: children's
#' parents, replicate origins and sources, statistics parents, the
#' application mode and the trained boundary.  Populations are matched to
#' the old ones by region.  The gate's trained boundary is kept when only
#' names changed, and dropped otherwise (so that it is calculated again).
#'
#' @param model Gating model.
#' @param gate.def Gate definition.
#' @param old.name Name of the gate replaced, or \code{NULL} to add.
#' @param mode Application mode ("fixed" or "recalculate"), or \code{NULL}
#'   to keep the current one.
#'
#' @return The model.
#'
#' @keywords internal
put.gate <- function( model, gate.def, old.name = NULL, mode = NULL ) {
    defs <- model$gate.definitions
    nm <- vapply( defs, function( g ) g$gate.name, character( 1 ) )
    new.name <- gate.def$gate.name
    i <- if ( is.null( old.name ) ) NA_integer_ else match( old.name, nm )
    if ( !is.null( old.name ) && is.na( i ) )
        stop( sprintf( "No gate named '%s'.", old.name ) )
    if ( new.name %in% nm[ setdiff( seq_along( nm ), i ) ] )
        stop( sprintf( "A gate named '%s' already exists.", new.name ) )
    if ( has.value( gate.def$parent.gate ) ) {
        others <- defs[ setdiff( seq_along( defs ), i ) ]
        parent <- gates.by.name( others )[[ gate.def$parent.gate ]]
        if ( is.null( parent ) )
            stop( sprintf( "The parent gate '%s' is not in the model.", gate.def$parent.gate ) )
        parent.pops <- names( effective.gate.def( parent, gates.by.name( others ) )$populations )
        if ( !( gate.def$parent.popul %in% parent.pops ) )
            stop( sprintf( "Gate '%s' has no population '%s'.", gate.def$parent.gate, gate.def$parent.popul ) )
    }

    if ( is.na( i ) ) {
        numbers <- vapply( defs, function( g ) suppressWarnings( as.integer( g$gate.number %||% NA ) ), integer( 1 ) )
        gate.def$gate.number <- if ( any( !is.na( numbers ) ) ) max( numbers, na.rm = TRUE ) + 1L else 1L
        model$gate.definitions <- c( defs, list( gate.def ) )
        model$application$mode.per.gate[[ new.name ]] <- mode %||% model$application$default.mode %||% "fixed"
        return( model )
    }

    old <- defs[[ i ]]
    gate.def$gate.number <- old$gate.number
    old.regions <- population.regions( old )
    new.regions <- population.regions( gate.def )
    pop.map <- character( 0 )
    for ( p in names( old.regions ) ) {
        q <- names( new.regions )[ new.regions == old.regions[[ p ]] ]
        if ( length( q ) == 1L && q != p ) pop.map[ p ] <- q
    }
    rename.key <- function( key ) {
        if ( !has.value( key ) ) return( key )
        k <- parse.population.key( key )
        if ( k[ 1L ] != old.name ) return( key )
        pop <- if ( k[ 2L ] %in% names( pop.map ) ) pop.map[[ k[ 2L ] ]] else k[ 2L ]
        if ( nzchar( k[ 2L ] ) ) population.key( new.name, pop ) else new.name
    }
    for ( j in seq_along( defs ) ) {
        if ( j == i ) next
        g <- defs[[ j ]]
        if ( identical( g$parent.gate, old.name ) ) {
            g$parent.gate <- new.name
            if ( has.value( g$parent.popul ) && g$parent.popul %in% names( pop.map ) )
                g$parent.popul <- pop.map[[ g$parent.popul ]]
        }
        if ( identical( g$origin.gate, old.name ) ) g$origin.gate <- new.name
        for ( axis in c( "x", "y" ) ) {
            src <- g$gate.param[[ paste0( "replicate.", axis ) ]]
            if ( !is.list( src ) ) next
            if ( identical( src$gate, old.name ) ) src$gate <- new.name
            if ( has.value( src$population ) ) src$population <- rename.key( src$population )
            g$gate.param[[ paste0( "replicate.", axis ) ]] <- src
        }
        if ( length( g$stats.parent ) )
            g$stats.parent <- lapply( g$stats.parent, rename.key )
        defs[[ j ]] <- g
    }

    strip <- function( g ) {
        gp <- g$gate.param %||% list()
        gp <- gp[ lengths( gp ) > 0L ]
        empty <- function( x ) if ( length( x ) == 0L ) NULL else x
        list( type = g$gate.type, x = empty( g$gate.marker.x ), y = empty( g$gate.marker.y ),
              sep.axis = if ( identical( g$gate.type, "1dsep" ) && has.value( g$gate.marker.y ) )
                  g$sep.axis %||% "x" else NULL,
              parent = if ( has.value( g$parent.gate ) ) c( g$parent.gate, g$parent.popul ) else NULL,
              algorithm = empty( g$algorithm ), origin = empty( g$origin.gate ),
              gate.param = empty( gp[ order( names( gp ) ) ] ),
              regions = sort( unname( population.regions( g ) ) ), template = empty( g$template ) )
    }
    boundary <- model$trained.boundaries[[ old.name ]]
    model$trained.boundaries[[ old.name ]] <- NULL
    if ( !is.null( boundary ) && identical( strip( old ), strip( gate.def ) ) ) {
        names( boundary ) <- vapply( names( boundary ), function( p )
            if ( p %in% names( pop.map ) ) pop.map[[ p ]] else p, character( 1 ) )
        for ( p in names( boundary ) )
            if ( !is.null( boundary[[ p ]]$label ) ) boundary[[ p ]]$label <- p
        model$trained.boundaries[[ new.name ]] <- boundary
    }
    current.mode <- model$application$mode.per.gate[[ old.name ]]
    model$application$mode.per.gate[[ old.name ]] <- NULL
    model$application$mode.per.gate[[ new.name ]] <- mode %||% current.mode %||%
        model$application$default.mode %||% "fixed"
    if ( length( gate.def$stats.parent ) && length( pop.map ) )
        names( gate.def$stats.parent ) <- vapply( names( gate.def$stats.parent ), function( p )
            if ( p %in% names( pop.map ) ) pop.map[[ p ]] else p, character( 1 ) )
    defs[[ i ]] <- gate.def
    ordered.gate.defs( defs )
    model$gate.definitions <- defs
    model
}


#' Remove a gate and its descendants
#'
#' @title Drop Gate
#'
#' @param model Gating model.
#' @param name Gate name.
#'
#' @return List with \code{model}, \code{removed} (gate names) and
#'   \code{broken} (remaining gates that replicate a removed gate's
#'   threshold or population).
#'
#' @keywords internal
drop.gate <- function( model, name ) {
    defs <- model$gate.definitions
    removed <- name
    repeat {
        more <- vapply( defs, function( g ) has.value( g$parent.gate ) && g$parent.gate %in% removed &&
                                            !( g$gate.name %in% removed ), logical( 1 ) )
        if ( !any( more ) ) break
        removed <- c( removed, vapply( defs[ more ], function( g ) g$gate.name, character( 1 ) ) )
    }
    keep <- !vapply( defs, function( g ) g$gate.name %in% removed, logical( 1 ) )
    model$gate.definitions <- defs[ keep ]
    model$trained.boundaries <- model$trained.boundaries[ setdiff( names( model$trained.boundaries ), removed ) ]
    model$application$mode.per.gate <-
        model$application$mode.per.gate[ setdiff( names( model$application$mode.per.gate ), removed ) ]
    broken <- character( 0 )
    for ( g in model$gate.definitions )
        if ( any( gate.dependencies( g ) %in% removed ) ) broken <- c( broken, g$gate.name )
    list( model = model, removed = removed, broken = broken )
}


#' Gates of a model as a table
#'
#' @title Gate Table
#'
#' @param model Gating model.
#' @param result Gating result on the events shown, or \code{NULL}.
#'
#' @return \code{data.frame} with number, gate, type, parent, channels,
#'   algorithm ("manual" for a manual threshold), populations (with \% of
#'   parent), mode and status, in hierarchy order.
#'
#' @keywords internal
gate.table <- function( model, result = NULL ) {
    defs <- model$gate.definitions
    if ( length( defs ) == 0L )
        return( data.frame( number = integer( 0 ), gate = character( 0 ), type = character( 0 ),
                            parent = character( 0 ), channels = character( 0 ),
                            algorithm = character( 0 ), populations = character( 0 ),
                            mode = character( 0 ), status = character( 0 ), stringsAsFactors = FALSE ) )
    by.name <- gates.by.name( defs )
    rows <- lapply( ordered.gate.defs( defs ), function( g ) {
        eff <- effective.gate.def( g, by.name )
        b <- result$boundaries[[ g$gate.name ]]
        pops <- vapply( names( eff$populations ), function( p ) {
            f <- b[[ p ]]$fraction
            if ( is.numeric( f ) && length( f ) == 1L ) sprintf( "%s %.1f%%", p, 100 * f ) else p
        }, character( 1 ) )
        shown <- function( axis ) {
            if ( has.value( eff$gate.param[[ paste0( "fixed.threshold.", axis ) ]] ) ) "manual"
            else axis.algorithm( eff, axis )
        }
        algo <- if ( identical( eff$gate.type, "2dsep" ) ) paste( shown( "x" ), shown( "y" ), sep = " / " )
            else if ( eff$gate.type %in% ag.threshold.types ) shown( "x" ) else eff$algorithm %||% ""
        st <- result$status[[ g$gate.name ]]
        status <- if ( is.null( st ) ) "" else if ( has.value( st$flag ) )
            paste0( st$flag, if ( nzchar( st$message %||% "" ) ) paste0( ": ", st$message ) else "" )
            else st$source %||% ""
        data.frame(
            number = suppressWarnings( as.integer( g$gate.number %||% NA ) ), gate = g$gate.name,
            type = g$gate.type %||% "",
            parent = if ( has.value( g$parent.gate ) ) population.key( g$parent.gate, g$parent.popul %||% "" ) else "",
            channels = paste( gate.channels( eff ), collapse = " / " ), algorithm = algo,
            populations = paste( pops, collapse = "; " ),
            mode = model$application$mode.per.gate[[ g$gate.name ]] %||% model$application$default.mode %||% "fixed",
            status = status, stringsAsFactors = FALSE )
    } )
    do.call( rbind, rows )
}


#' Change the transform of one channel
#'
#' @title Retransform Channel
#'
#' @description
#' Puts the channel's events back in raw units with the current inverse
#' transform and transforms them with the new parameters; the forward,
#' inverse and parameter lists in \code{fcs.setup$ap} are updated.
#' Boundaries of a model on the channel are converted with
#' \code{\link{remap.model}}.
#'
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}.
#' @param channel Channel (antigen) name.
#' @param params Transform parameters (see \code{\link{make.transform}}).
#'
#' @return The updated \code{fcs.setup}.
#'
#' @keywords internal
retransform.channel <- function( fcs.setup, channel, params ) {
    if ( !( channel %in% colnames( fcs.setup$fcs.expr.data ) ) )
        stop( sprintf( "No channel '%s' in the data.", channel ) )
    tr <- make.transform( params )
    old.inverse <- fcs.setup$ap$transform.inv[[ channel ]] %||% function( x ) x
    fcs.setup$fcs.expr.data[ , channel ] <- tr$forward( old.inverse( fcs.setup$fcs.expr.data[ , channel ] ) )
    fcs.setup$ap$transform[[ channel ]] <- tr$forward
    fcs.setup$ap$transform.inv[[ channel ]] <- tr$inverse
    fcs.setup$ap$transform.params[[ channel ]] <- params
    fcs.setup
}


#' Write the transforms of an experiment as a CSV file
#'
#' @title Write Transformation CSV
#'
#' @description
#' One row per scatter and fluorescence channel, keyed by dye, with the
#' Gating-ML columns read by \code{\link{read.transformation}}
#' (\code{transformation.filename} of \code{\link{get.autogating.param}}):
#' \code{id} (0 linear, 1 logicle, 2 log), \code{scale.t},
#' \code{linear.a}, \code{logicle.w}, \code{logicle.m}, \code{logicle.a},
#' \code{log.m}.
#'
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}.
#' @param file Output CSV path.
#'
#' @return The table, invisibly.
#'
#' @export
write.transformation.csv <- function( fcs.setup, file ) {
    panel <- fcs.setup$fcs.panel
    params <- fcs.setup$ap$transform.params
    cols <- c( "id", "scale.t", "linear.a", "logicle.w", "logicle.m", "logicle.a", "log.m" )
    keep <- panel$antigen %in% names( params )[ lengths( params ) > 0L ]
    tbl <- data.frame( dye = panel$dye[ keep ], stringsAsFactors = FALSE )
    for ( col in cols )
        tbl[[ col ]] <- vapply( panel$antigen[ keep ], function( ch ) {
            v <- params[[ ch ]][[ col ]]
            if ( is.null( v ) ) NA_real_ else as.numeric( v )
        }, numeric( 1 ) )
    utils::write.csv( tbl, file, row.names = FALSE )
    invisible( tbl )
}
