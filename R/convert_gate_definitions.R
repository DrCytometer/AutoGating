# convert_gate_definitions.R
#
# Converts gate definition lists in the earlier AutoGating format (the
# population and activation gate files of scripts/ST1) to the gate model of
# gate_model.R.
#
# Earlier format:
#   populations     named by position: a 1dsep gate's first population is
#                   below the threshold and its second above; a 2dsep gate's
#                   four are x-y-, x+y-, x-y+, x+y+; "*" fills a position
#                   that is not kept
#   gate.marker.x/y plot axes; a 1dsep gate separates on the axis given by
#                   gate.param$sep.axis (1 = x, 2 = y, default 2)
#   algorithm       one name, or one per axis for 2dsep; "replicate" copies
#                   the boundaries of gate.param$origin.gate by position
#   gate.param      join = "1:2:3": the first population is the union of
#                   those quadrants, the second the remaining one;
#                   region.min/max (2dsep: region.xmin/xmax/ymin/ymax): the
#                   window of the range the threshold is calculated in;
#                   threshold.param: tail quantile (0.99 = a 1% tail)
#   stats.parent    list(gate, popul) of parallel vectors; "parent.N" is the
#                   N-th ancestor population; popul "*" names the gate
#   root gate       algorithm "all" without a parent
#   activation      one entry per parent population with a set of markers;
#                   each marker becomes a 1dsep gate on that population
#
# Conversion:
#   - The root gate is dropped; its children start from all events.
#   - A gate whose parent is several populations of one threshold gate gets
#     a parent population added to that gate: the union of their regions,
#     named by joining their names with ".or.".
#   - Population regions follow the positions; "*" positions are dropped.
#   - A 1dsep gate's separation channel becomes gate.marker.x and its other
#     axis gate.marker.y, with sep.axis "x" or "y" giving the plot axis; a
#     gate whose two axes are the same channel becomes a histogram.
#   - Replicate gates copy their origin's thresholds (replicate.x/y = the
#     origin gate) and take the regions of the origin's kept populations in
#     order.  A replicate without origin.gate copies the gate on its
#     separation channel whose parent population corresponds to its own
#     (the same population of the gate its parent gate replicates), or else
#     the only gate on that channel.
#   - Window and tail settings are kept where they differ from the earlier
#     defaults (window 0.01 to 0.99; tail.fraction = 1 - threshold.param);
#     quantile trims, plot limits, breaks and label sizes are dropped, as
#     the gate engine and the plots use their own.
#   - The statistics parent of each population is the first of the gate's
#     stats.parent entries other than its own parent, whose fraction is
#     always reported.
#
# Gates whose populations already have regions are taken to be in the
# current format and are kept as they are.


earlier.param.keys <- c(
    "quantile.x", "quantile.y", "density.quantile.x", "density.quantile.y",
    "density.trimneg.x", "density.trimneg.y", "limit.xmin", "limit.xmax",
    "limit.ymin", "limit.ymax", "break.x", "break.x.factor", "break.y",
    "break.y.factor", "label.size", "sep.axis", "join", "origin.gate",
    "region.min", "region.max", "region.xmin", "region.xmax", "region.ymin",
    "region.ymax", "threshold.param" )


#' Whether a gate definition is in the current format
#' @title Is Current Gate
#' @param g Gate definition.
#' @return \code{TRUE} for a "replicate" gate or when any population has a
#'   region.
#' @keywords internal
is.current.gate <- function( g ) {
    if ( identical( g$gate.type, "replicate" ) )
        return( TRUE )
    any( vapply( g$populations, function( p ) has.value( p$region ), logical( 1 ) ) )
}


#' Separation channel of an earlier-format 1dsep gate
#' @title Earlier Separation Axes
#' @param g Earlier-format gate definition.
#' @return List with \code{x} (separation channel), \code{y} (other axis or
#'   \code{NULL} for a histogram) and \code{sep.axis}.
#' @keywords internal
earlier.separation <- function( g ) {
    sep <- suppressWarnings( as.integer( g$gate.param[[ "sep.axis" ]] %||% 2L ) )
    x <- g$gate.marker.x
    y <- g$gate.marker.y
    if ( identical( sep, 2L ) && has.value( y ) ) {
        s <- y
        o <- x
        side <- "y"
    } else {
        s <- x
        o <- y
        side <- "x"
    }
    if ( !has.value( o ) || identical( o, s ) )
        return( list( x = s, y = NULL, sep.axis = "x" ) )
    list( x = s, y = o, sep.axis = side )
}


#' Activation gate entries as earlier-format gates
#' @title Activation Entries as Gates
#' @description
#' Each marker of an activation entry becomes a 1dsep gate named
#' \code{<parent.gate>.<parent.popul>.<marker>} with the marker's positive
#' population "pos" above the threshold.  Entries with an empty
#' \code{gate.type} or no markers are display-only and skipped.
#' @param activation.gates List of activation entries.
#' @param start.number Added to \code{popul.number} to number the gates
#'   after the population gates.
#' @return List of earlier-format gate definitions.
#' @keywords internal
activation.as.gates <- function( activation.gates, start.number = 0L ) {
    out <- list()
    for ( ad in activation.gates ) {
        if ( !identical( ad$gate.type, "1dsep" ) || length( ad$markers ) == 0L )
            next
        for ( marker in ad$markers ) {
            alg <- ad$marker.algorithm.override[[ marker ]] %||% ad$algorithm %||% "tail"
            out[[ length( out ) + 1L ]] <- list(
                gate.number   = start.number + as.integer( ad$popul.number %||% 0L ),
                gate.name     = sprintf( "%s.%s.%s", ad$parent.gate, ad$parent.popul, marker ),
                gate.type     = "1dsep",
                gate.marker.x = if ( identical( ad$gate.marker.x, "*" ) ) marker else ad$gate.marker.x,
                gate.marker.y = if ( identical( ad$gate.marker.y, "*" ) ) marker else ad$gate.marker.y,
                parent.gate   = ad$parent.gate,
                parent.popul  = ad$parent.popul,
                populations   = list(
                    `*` = list( label = "*" ),
                    pos = list( label = sprintf( "%s+ %s", marker, ad$popul.label ),
                                label.pos = ad$label.pos ) ),
                algorithm     = alg,
                gate.param    = ad$gate.param %||% list() )
        }
    }
    out
}


#' Populations kept from earlier-format positions
#' @title Kept Populations
#' @param pops Earlier-format populations (named list, "*" for unused
#'   positions).
#' @param regions Region of each position (\code{NA} drops it).
#' @return Named list of populations with label, label.pos and region.
#' @keywords internal
kept.populations <- function( pops, regions ) {
    out <- list()
    for ( i in seq_along( pops ) ) {
        p <- names( pops )[ i ]
        if ( !has.value( p ) || p == "*" || is.na( regions[ i ] ) )
            next
        lab <- pops[[ i ]]$label
        pos <- suppressWarnings( as.integer( pops[[ i ]]$label.pos ) )
        e <- list( label = if ( has.value( lab ) ) as.character( lab ) else p )
        if ( length( pos ) == 1L && !is.na( pos ) )
            e$label.pos <- pos
        e$region <- regions[ i ]
        out[[ p ]] <- e
    }
    json.object( out )
}


#' Convert earlier-format gate definitions
#'
#' @title Upgrade Gate Definitions
#'
#' @description
#' Converts population and activation gate lists written for earlier
#' versions of AutoGating (positional population names, \code{"*"} fillers,
#' \code{join}, \code{origin.gate}, \code{sep.axis} numbers,
#' \code{list(gate, popul)} statistics parents, a root "all" gate) to gate
#' definitions of the current model, which Honeychrome reads as well.  The
#' rules are listed in the file header of \code{convert_gate_definitions.R}.
#'
#' @param population.gates List of population gate definitions.
#' @param activation.gates List of activation entries (optional).
#'
#' @return List with \code{gate.definitions} and \code{messages} (a
#'   character vector describing choices made and settings dropped).
#'
#' @seealso \code{\link{define.gates}}, \code{\link{load.gate.definitions}}
#'
#' @export
convert.gate.definitions <- function( population.gates, activation.gates = list() ) {
    numbers <- vapply( population.gates, function( g )
        suppressWarnings( as.integer( g$gate.number %||% NA ) ), integer( 1 ) )
    start <- if ( any( !is.na( numbers ) ) ) max( numbers, na.rm = TRUE ) else 0L
    old <- c( population.gates, activation.as.gates( activation.gates, start ) )
    if ( length( old ) == 0L )
        return( list( gate.definitions = list(), messages = character( 0 ) ) )
    nm <- vapply( old, function( g ) as.character( g$gate.name ), character( 1 ) )
    if ( anyDuplicated( nm ) )
        stop( "Duplicate gate names: ", paste( unique( nm[ duplicated( nm ) ] ), collapse = ", " ) )
    by.name <- stats::setNames( old, nm )
    current <- stats::setNames( vapply( old, is.current.gate, logical( 1 ) ), nm )
    is.root <- vapply( old, function( g ) identical( as.character( g$algorithm ), "all" ) &&
                           !has.value( g$parent.gate ), logical( 1 ) )
    root <- nm[ !current & is.root ]

    messages <- character( 0 )
    note <- function( fmt, ... ) messages <<- c( messages, sprintf( fmt, ... ) )
    done <- list()
    busy <- character( 0 )
    unions <- list()

    converted <- function( name ) {
        if ( !is.null( done[[ name ]] ) )
            return( done[[ name ]] )
        if ( name %in% busy )
            stop( sprintf( "Replicate gates refer to each other in a cycle at '%s'.", name ) )
        busy <<- c( busy, name )
        out <- if ( current[[ name ]] ) by.name[[ name ]] else convert( by.name[[ name ]] )
        busy <<- setdiff( busy, name )
        done[[ name ]] <<- out
        out
    }

    is.earlier.replicate <- function( n ) {
        !current[[ n ]] && identical( as.character( by.name[[ n ]]$algorithm ), "replicate" )
    }

    parent.of <- function( g ) {
        pg <- unique( as.character( g$parent.gate ) )
        if ( length( pg ) > 1L )
            stop( sprintf( "Gate '%s' has parent populations in %d gates; a gate has one parent.",
                           g$gate.name, length( pg ) ) )
        if ( !has.value( pg ) || pg %in% root )
            return( list( parent.gate = NULL, parent.popul = "root" ) )
        pp <- as.character( g$parent.popul )
        if ( length( pp ) > 1L ) {
            union.name <- paste( pp, collapse = ".or." )
            unions[[ pg ]][[ union.name ]] <<- pp
            return( list( parent.gate = pg, parent.popul = union.name ) )
        }
        list( parent.gate = pg, parent.popul = pp )
    }

    add.union.populations <- function() {
        for ( pg in names( unions ) ) {
            p <- done[[ pg ]]
            if ( !isTRUE( p$gate.type %in% ag.threshold.types ) )
                stop( sprintf( "Gate '%s' is the parent of a union of its populations but is not a threshold gate.", pg ) )
            regions <- population.regions( p )
            for ( union.name in names( unions[[ pg ]] ) ) {
                members <- unions[[ pg ]][[ union.name ]]
                if ( !all( members %in% names( regions ) ) )
                    stop( sprintf( "Gate '%s' has no population %s.", pg,
                                   paste( setdiff( members, names( regions ) ), collapse = ", " ) ) )
                quads <- unique( unlist( lapply( regions[ members ], region.quadrants ) ) )
                region <- if ( p$gate.type == "1dsep" ) NA_character_ else
                    paste( ag.quadrant.regions[ ag.quadrant.regions %in% quads ], collapse = "|" )
                if ( is.na( region ) || !nzchar( region ) )
                    stop( sprintf( "The union of populations %s of gate '%s' is not a region.",
                                   paste( members, collapse = ", " ), pg ) )
                labels <- vapply( members, function( m ) p$populations[[ m ]]$label %||% m, character( 1 ) )
                p$populations[[ union.name ]] <- list( label = paste( labels, collapse = " or " ),
                                                       region = region )
                note( "Gate '%s': population '%s' (%s) added as the parent of gates on %s.",
                      pg, union.name, region, paste( members, collapse = " and " ) )
            }
            done[[ pg ]] <<- p
        }
    }

    axes.of <- function( g, alg ) {
        if ( identical( g$gate.type, "1dsep" ) ) {
            sep <- earlier.separation( g )
            out <- list( gate.type = "1dsep", gate.marker.x = sep$x )
            out[ "gate.marker.y" ] <- list( sep$y )
            if ( !is.null( sep$y ) )
                out$sep.axis <- sep$sep.axis
            return( out )
        }
        type <- if ( identical( g$gate.type, "free" ) && any( alg %in% c( "singlets", "fs.singlets" ) ) )
            "singlets" else g$gate.type
        out <- list( gate.type = type, gate.marker.x = g$gate.marker.x )
        out[ "gate.marker.y" ] <- list( g$gate.marker.y )
        out
    }

    position.regions <- function( g, type ) {
        n <- length( g$populations )
        if ( type == "1dsep" ) {
            if ( n > 2L )
                stop( sprintf( "1dsep gate '%s' has %d populations; at most 2.", g$gate.name, n ) )
            return( c( "below", "above" )[ seq_len( n ) ] )
        }
        if ( type == "2dsep" ) {
            join <- g$gate.param[[ "join" ]]
            if ( has.value( join ) ) {
                idx <- sort( as.integer( strsplit( as.character( join ), ":", fixed = TRUE )[[ 1L ]] ) )
                rest <- setdiff( 1:4, idx )
                regions <- c( paste( ag.quadrant.regions[ idx ], collapse = "|" ),
                              if ( length( rest ) ) paste( ag.quadrant.regions[ rest ], collapse = "|" ) )
                return( regions[ seq_len( n ) ] )
            }
            if ( n > 4L )
                stop( sprintf( "2dsep gate '%s' has %d populations; at most 4.", g$gate.name, n ) )
            return( ag.quadrant.regions[ seq_len( n ) ] )
        }
        rep( "inside", n )
    }

    sep.channels <- function( n, seen = character( 0 ) ) {
        h <- by.name[[ n ]]
        if ( current[[ n ]] )
            return( c( h$gate.marker.x, if ( identical( h$gate.type, "2dsep" ) ) h$gate.marker.y ) )
        o <- h$gate.param[[ "origin.gate" ]]
        if ( is.earlier.replicate( n ) && has.value( o ) && o %in% nm && !( o %in% seen ) )
            return( sep.channels( o, c( seen, n ) ) )
        if ( identical( h$gate.type, "1dsep" ) )
            return( earlier.separation( h )$x )
        if ( identical( h$gate.type, "2dsep" ) )
            return( c( h$gate.marker.x, h$gate.marker.y ) )
        character( 0 )
    }

    inferred.origin <- function( g ) {
        if ( !( g$gate.type %in% ag.threshold.types ) )
            return( NULL )
        own <- if ( g$gate.type == "1dsep" ) earlier.separation( g )$x else
            c( g$gate.marker.x, g$gate.marker.y )
        cands <- nm[ vapply( nm, function( n ) {
            if ( n == g$gate.name || n %in% root )
                return( FALSE )
            if ( is.earlier.replicate( n ) && !has.value( by.name[[ n ]]$gate.param[[ "origin.gate" ]] ) )
                return( FALSE )
            all( own %in% sep.channels( n ) )
        }, logical( 1 ) ) ]
        if ( length( cands ) == 0L )
            return( NULL )
        pg <- g$parent.gate
        parent.origin <- if ( has.value( pg ) && pg %in% nm && is.earlier.replicate( pg ) )
            by.name[[ pg ]]$gate.param[[ "origin.gate" ]] else NULL
        same.place <- cands[ vapply( cands, function( n ) {
            h <- by.name[[ n ]]
            has.value( parent.origin ) && identical( h$parent.gate, parent.origin ) &&
                identical( h$parent.popul, g$parent.popul )
        }, logical( 1 ) ) ]
        if ( length( same.place ) == 1L ) {
            note( "Gate '%s' had no origin gate; it replicates '%s', the gate on %s in population %s of '%s'.",
                  g$gate.name, same.place, paste( own, collapse = "/" ), g$parent.popul, parent.origin )
            return( same.place )
        }
        if ( length( cands ) == 1L ) {
            note( "Gate '%s' had no origin gate; it replicates '%s', the only gate on %s.",
                  g$gate.name, cands, paste( own, collapse = "/" ) )
            return( cands )
        }
        NULL
    }

    convert.param <- function( g, type, alg ) {
        gp <- g$gate.param %||% list()
        out <- list()
        if ( type == "2dsep" && length( alg ) == 2L && alg[ 2 ] != alg[ 1 ] )
            out$algorithm.y <- alg[ 2 ]
        window <- function( key, new, default ) {
            v <- suppressWarnings( as.numeric( gp[[ key ]] ) )
            if ( length( v ) == 1L && !is.na( v ) && abs( v - default ) > 1e-12 )
                out[[ new ]] <<- v
        }
        if ( type == "1dsep" ) {
            window( "region.min", "axis.lo", 0.01 )
            window( "region.max", "axis.hi", 0.99 )
        }
        if ( type == "2dsep" )
            for ( a in c( "x", "y" ) ) {
                window( paste0( "region.", a, "min" ), paste0( "axis.lo.", a ), 0.01 )
                window( paste0( "region.", a, "max" ), paste0( "axis.hi.", a ), 0.99 )
            }
        tp <- suppressWarnings( as.numeric( gp[[ "threshold.param" ]] ) )
        if ( length( tp ) == 1L && !is.na( tp ) ) {
            axes <- if ( type == "2dsep" ) c( "x", "y" ) else "x"
            tail.axes <- axes[ rep_len( alg, length( axes ) ) == "tail" ]
            if ( length( tail.axes ) == length( axes ) )
                out$tail.fraction <- 1 - tp
            else if ( length( tail.axes ) )
                out[[ paste0( "tail.fraction.", tail.axes ) ]] <- 1 - tp
            else
                note( "Gate '%s': threshold.param %g applies to no tail axis and is dropped.",
                      g$gate.name, tp )
        }
        unknown <- setdiff( names( gp ), earlier.param.keys )
        if ( length( unknown ) )
            note( "Gate '%s': parameters dropped: %s.", g$gate.name, paste( unknown, collapse = ", " ) )
        json.object( out )
    }

    ancestor.key <- function( name, k ) {
        cur <- name
        for ( step in seq_len( k ) ) {
            h <- by.name[[ cur ]]
            pg <- unique( as.character( h$parent.gate ) )
            pp <- paste( h$parent.popul, collapse = ".or." )
            if ( !has.value( pg ) || pg %in% root )
                return( "root" )
            cur <- pg
        }
        population.key( pg, pp )
    }

    convert.stats.parent <- function( g, new ) {
        sp <- g$stats.parent
        sg <- as.character( unlist( sp$gate ) )
        so <- as.character( unlist( sp$popul ) )
        pops <- names( if ( identical( new$gate.type, "replicate" ) )
            converted( new$origin.gate )$populations else new$populations )
        if ( length( sg ) == 0L || length( pops ) == 0L )
            return( json.object( NULL ) )
        specs <- character( 0 )
        for ( i in seq_along( sg ) ) {
            a <- sg[ i ]
            b <- if ( length( so ) >= i ) so[ i ] else "*"
            if ( a == "*" )
                next
            m <- regmatches( a, regexec( "^parent\\.([0-9]+)$", a ) )[[ 1L ]]
            if ( length( m ) == 2L && identical( a, b ) ) {
                k <- as.integer( m[ 2 ] )
                if ( k > 1L )
                    specs <- c( specs, ancestor.key( g$gate.name, k ) )
                next
            }
            if ( a %in% root ) {
                specs <- c( specs, "root" )
                next
            }
            if ( !( a %in% nm ) ) {
                note( "Gate '%s': statistics parent gate '%s' does not exist.", g$gate.name, a )
                next
            }
            specs <- c( specs, if ( b == "*" ) a else population.key( a, b ) )
        }
        if ( has.value( new$parent.gate ) )
            specs <- setdiff( specs, population.key( new$parent.gate, new$parent.popul ) )
        specs <- unique( specs )
        if ( length( specs ) == 0L )
            return( json.object( NULL ) )
        if ( length( specs ) > 1L )
            note( "Gate '%s': statistics parents %s; the first is kept.",
                  g$gate.name, paste( specs, collapse = ", " ) )
        stats::setNames( as.list( rep( specs[ 1L ], length( pops ) ) ), pops )
    }

    convert.replicate <- function( g ) {
        name <- g$gate.name
        origin <- g$gate.param[[ "origin.gate" ]]
        if ( has.value( origin ) && !( origin %in% nm ) ) {
            note( "Gate '%s': origin gate '%s' does not exist.", name, origin )
            origin <- NULL
        }
        if ( !has.value( origin ) )
            origin <- inferred.origin( g )
        if ( is.null( origin ) ) {
            note( paste( "Gate '%s' replicates no gate (no origin.gate, and no single gate on its",
                         "channel); set a source in gate.param$replicate.x." ), name )
            axes <- axes.of( g, "replicate" )
            return( c( axes, parent.of( g ), list(
                populations = kept.populations( g$populations, position.regions( g, axes$gate.type ) ),
                algorithm = "replicate", gate.param = json.object( NULL ) ) ) )
        }
        o <- converted( origin )
        if ( isTRUE( o$gate.type %in% ag.threshold.types ) ) {
            axes <- list( gate.type = o$gate.type, gate.marker.x = o$gate.marker.x )
            axes[ "gate.marker.y" ] <- list( o$gate.marker.y )
            if ( !is.null( o$sep.axis ) )
                axes$sep.axis <- o$sep.axis
            regions <- unname( population.regions( o ) )
            pops <- kept.populations( g$populations, regions[ seq_along( g$populations ) ] )
            if ( length( pops ) < sum( names( g$populations ) != "*" ) )
                note( "Gate '%s': populations beyond the %d of origin '%s' are dropped.",
                      name, length( regions ), origin )
            own <- sort( c( g$gate.marker.x, g$gate.marker.y ) )
            theirs <- sort( c( o$gate.marker.x, o$gate.marker.y ) )
            if ( !identical( own, theirs ) )
                note( "Gate '%s' is drawn on %s but replicates '%s' on %s; it now uses the origin's channels.",
                      name, paste( own, collapse = "/" ), origin, paste( theirs, collapse = "/" ) )
            gp <- list( replicate.x = list( gate = origin ) )
            if ( o$gate.type == "2dsep" )
                gp$replicate.y <- list( gate = origin )
            return( c( axes, parent.of( g ), list( populations = pops, algorithm = "replicate",
                                                    gate.param = gp ) ) )
        }
        note( "Gate '%s' replicates the %s gate '%s' and takes its population names.",
              name, o$gate.type %||% "free", origin )
        c( list( gate.type = "replicate", origin.gate = origin ), parent.of( g ),
           list( populations = json.object( NULL ), gate.param = json.object( NULL ) ) )
    }

    convert <- function( g ) {
        name <- g$gate.name
        alg <- as.character( g$algorithm %||% "tail" )
        alg[ alg == "kde.min" ] <- "kde_min"
        out <- list( gate.name = name,
                     gate.number = suppressWarnings( as.integer( g$gate.number %||% NA ) ) )
        if ( identical( alg, "replicate" ) ) {
            out <- c( out, convert.replicate( g ) )
        } else {
            axes <- axes.of( g, alg )
            type <- axes$gate.type
            out <- c( out, axes, parent.of( g ), list(
                populations = kept.populations( g$populations, position.regions( g, type ) ),
                algorithm   = if ( type == "singlets" ) "singlets" else alg[ 1L ],
                gate.param  = convert.param( g, type, alg ) ) )
        }
        out$stats.parent <- convert.stats.parent( g, out )
        out
    }

    kept <- nm[ !( nm %in% root ) ]
    for ( n in kept )
        converted( n )
    add.union.populations()
    list( gate.definitions = unname( done[ kept ] ), messages = messages )
}
