# gate_application.R
#
# Applies gates to events and runs a gate hierarchy on one event matrix
# (one sample, or pooled training events), as Honeychrome Automated Gating
# does.
#
# A gating result is a list:
#   boundaries     gate -> population -> boundary entry, as applied, with
#                  n.events and fraction of the parent
#   events         gate -> population -> integer row indices of the data
#   parent.events  gate -> integer row indices the gate was applied to
#   status         gate -> list(source, drift, flag, message); source is
#                  "trained", "fixed", "recalculated" or "reference"; flag
#                  is NA, "drift", "failed" or "skipped"
#
# Populations are stored as row indices rather than full-length logical
# masks, to keep memory proportional to population sizes.
#
# Application modes: "fixed" applies the stored boundary; "recalculate"
# calculates the gate on the sample, keeps the stored boundary as a
# reference, and flags (or, with drift.action = "reference", replaces)
# results that move further than drift.limit (transformed units: the
# largest threshold shift, or the centroid distance of a polygon).


#' Points inside a polygon
#'
#' @title Points in Polygon
#'
#' @description
#' Even-odd ray casting; points exactly on an edge may fall either side.
#'
#' @param x,y Point coordinates.
#' @param poly Polygon vertices (n x 2 matrix, closed or open).
#'
#' @return Logical vector.
#'
#' @keywords internal
inside.polygon <- function( x, y, poly ) {
    px <- poly[ , 1 ]
    py <- poly[ , 2 ]
    n <- length( px )
    inside <- logical( length( x ) )
    j <- n
    for ( i in seq_len( n ) ) {
        yi <- py[ i ]; yj <- py[ j ]
        if ( yi != yj ) {
            crosses <- ( yi > y ) != ( yj > y )
            if ( any( crosses ) ) {
                xc <- ( px[ j ] - px[ i ] ) * ( y[ crosses ] - yi ) / ( yj - yi ) + px[ i ]
                hit <- x[ crosses ] < xc
                inside[ crosses ][ hit ] <- !inside[ crosses ][ hit ]
            }
        }
        j <- i
    }
    inside
}


#' Values moved just inside display limits
#'
#' @title Edge Clamp
#'
#' @param values Numeric vector.
#' @param limits \code{c(lo, hi)} or \code{NULL}.
#'
#' @return Values clamped to \code{(lo, hi)}, where the plots draw them.
#'
#' @keywords internal
edge.clamped <- function( values, limits ) {
    if ( length( limits ) < 2L ) return( values )
    lo <- limits[ 1 ]
    hi <- limits[ 2 ]
    if ( !( hi > lo ) ) return( values )
    eps <- 1e-6 * ( hi - lo )
    pmin( pmax( values, lo + eps ), hi - eps )
}


#' Population masks of one gate
#'
#' @title Population Masks
#'
#' @description
#' Threshold gates are open-ended (\code{x >= t} is the positive side);
#' polygons and ranges contain the events inside them, with events beyond
#' the display limits tested on the axis edge.
#'
#' @param gate.def Effective gate definition.
#' @param boundary Named list of boundary entries.
#' @param data Event matrix (columns named by channel).
#' @param axis.limits Named list of display limits, or \code{NULL}.
#'
#' @return Named list population -> logical vector over the rows of
#'   \code{data}; empty when a channel is missing or the boundary is empty.
#'
#' @export
population.masks <- function( gate.def, boundary, data, axis.limits = NULL ) {
    if ( length( boundary ) == 0L ) return( list() )
    chans <- gate.channels( gate.def )
    if ( !all( chans %in% colnames( data ) ) ) return( list() )
    if ( is.shape.boundary( boundary ) )
        return( shape.masks( boundary, data, axis.limits ) )
    n <- nrow( data )
    gate.type <- gate.type.of( gate.def )
    regions <- population.regions( gate.def, boundary )
    masks <- list()

    if ( gate.type %in% ag.threshold.types ) {
        tx <- boundary[[ 1L ]]$threshold.x
        if ( !has.value( tx ) ) return( list() )
        above.x <- data[ , chans[ 1 ] ] >= tx
        if ( gate.type == "1dsep" ) {
            for ( p in names( boundary ) ) {
                r <- if ( p %in% names( regions ) ) regions[[ p ]] else "above"
                masks[[ p ]] <- if ( r == "above" ) above.x else !above.x
            }
            return( masks )
        }
        ty <- boundary[[ 1L ]]$threshold.y
        if ( !has.value( ty ) || length( chans ) < 2L ) return( list() )
        above.y <- data[ , chans[ 2 ] ] >= ty
        for ( p in names( boundary ) ) {
            r <- if ( p %in% names( regions ) ) regions[[ p ]] else "x+y+"
            masks[[ p ]] <- quadrant.mask( r, above.x, above.y )
        }
        return( masks )
    }

    x <- edge.clamped( data[ , chans[ 1 ] ], axis.limits[[ chans[ 1 ] ]] )
    y <- if ( length( chans ) > 1L ) edge.clamped( data[ , chans[ 2 ] ], axis.limits[[ chans[ 2 ] ]] ) else NULL
    for ( p in names( boundary ) ) {
        e <- boundary[[ p ]]
        if ( !is.null( e$range ) && length( chans ) == 1L ) {
            masks[[ p ]] <- x >= e$range[ 1 ] & x <= e$range[ 2 ]
            next
        }
        poly <- e$boundary
        if ( is.null( y ) || !is.matrix( poly ) || nrow( poly ) < 3L ) {
            masks[[ p ]] <- logical( n )
            next
        }
        masks[[ p ]] <- inside.polygon( x, y, poly )
    }
    masks
}


#' Masks of shapes as drawn
#'
#' @title Shape Masks
#'
#' @param boundary Named list of shape entries (with \code{shape.channels}).
#' @param data Event matrix.
#' @param axis.limits Named list of display limits, or \code{NULL}.
#'
#' @return Named list of logical vectors; empty when a channel is missing.
#'
#' @keywords internal
shape.masks <- function( boundary, data, axis.limits = NULL ) {
    n <- nrow( data )
    masks <- list()
    for ( p in names( boundary ) ) {
        e <- boundary[[ p ]]
        chans <- unlist( e$shape.channels )
        chans <- chans[ vapply( chans, has.value, logical( 1 ) ) ]
        if ( !all( chans %in% colnames( data ) ) ) return( list() )
        x <- edge.clamped( data[ , chans[ 1 ] ], axis.limits[[ chans[ 1 ] ]] )
        if ( length( chans ) == 1L ) {
            masks[[ p ]] <- if ( !is.null( e$range ) ) x >= e$range[ 1 ] & x <= e$range[ 2 ] else logical( n )
            next
        }
        if ( !is.matrix( e$boundary ) || nrow( e$boundary ) < 3L ) {
            masks[[ p ]] <- logical( n )
            next
        }
        y <- edge.clamped( data[ , chans[ 2 ] ], axis.limits[[ chans[ 2 ] ]] )
        masks[[ p ]] <- inside.polygon( x, y, e$boundary )
    }
    masks
}


#' Threshold a gate has on a channel
#'
#' @title Threshold on Channel
#'
#' @param gate.def Effective gate definition.
#' @param entries Boundary entries of the gate.
#' @param channel Channel name.
#'
#' @return Numeric threshold, or \code{NULL}.
#'
#' @keywords internal
threshold.on.channel <- function( gate.def, entries, channel ) {
    if ( length( entries ) == 0L || is.shape.boundary( entries ) ||
         !( gate.type.of( gate.def ) %in% ag.threshold.types ) )
        return( NULL )
    first <- entries[[ 1L ]]
    if ( identical( gate.def$gate.marker.x, channel ) && has.value( first$threshold.x ) )
        return( first$threshold.x )
    if ( identical( gate.def$gate.type, "2dsep" ) && identical( gate.def$gate.marker.y, channel ) &&
         has.value( first$threshold.y ) )
        return( first$threshold.y )
    NULL
}


#' Fix replicated thresholds of a gate
#'
#' @title Resolve Replicated Thresholds
#'
#' @description
#' Copies the gate definition with the thresholds of its replicated axes in
#' \code{gate.param$fixed.threshold.x/y}: copied from another gate's
#' threshold on the same channel, or calculated on a reference population's
#' events (multiplied by the axis's threshold factor).
#'
#' @param gate.def Effective gate definition.
#' @param data Event matrix the hierarchy is run on.
#' @param result Gating result so far.
#' @param by.name Named list of gate definitions.
#' @param ctx Calculation context (notes are added to it).
#'
#' @return List with \code{gate.def} (\code{NULL} when a source is missing)
#'   and \code{reason}.
#'
#' @keywords internal
resolve.replicates <- function( gate.def, data, result, by.name, ctx ) {
    sources <- replicate.sources( gate.def )
    if ( length( sources ) == 0L )
        return( list( gate.def = gate.def, reason = "" ) )
    out <- gate.def
    gp <- out$gate.param %||% list()
    fail <- function( reason ) list( gate.def = NULL, reason = reason )
    for ( axis in names( sources ) ) {
        src <- sources[[ axis ]]
        channel <- axis.channel( gate.def, axis )
        if ( is.null( src ) )
            return( fail( sprintf( "no replicate source chosen for %s", channel ) ) )
        if ( has.value( src$gate ) ) {
            origin <- by.name[[ src$gate ]]
            if ( is.null( origin ) )
                return( fail( sprintf( "replicate source gate '%s' does not exist", src$gate ) ) )
            eff <- effective.gate.def( origin, by.name )
            t <- threshold.on.channel( eff, result$boundaries[[ src$gate ]], channel )
            if ( is.null( t ) )
                return( fail( sprintf( "gate '%s' has no threshold on %s to replicate",
                                       src$gate, channel ) ) )
        } else {
            key <- parse.population.key( src$population )
            idx <- result$events[[ key[ 1 ] ]][[ key[ 2 ] ]]
            if ( is.null( idx ) || !( channel %in% colnames( data ) ) )
                return( fail( sprintf( "reference population %s is not available", src$population ) ) )
            values <- trim.events( data[ idx, channel ],
                as.numeric( gp[[ "quantile.trim" ]] %||% ag.defaults$quantile.trim ) )
            if ( length( values ) < ag.defaults$min.events )
                return( fail( sprintf( "reference population %s has too few events", src$population ) ) )
            suffix <- if ( identical( gate.def$gate.type, "2dsep" ) ) paste0( ".", axis ) else ""
            algo <- if ( has.value( src$algorithm ) && src$algorithm %in% ag.threshold.algorithms )
                src$algorithm else "tail"
            calc <- crop.to.region( values,
                as.numeric( gp[[ paste0( "axis.lo", suffix ) ]] %||% default.axis.lo( algo ) ),
                as.numeric( gp[[ paste0( "axis.hi", suffix ) ]] %||% ag.defaults$axis.hi ) )
            if ( length( calc ) < ag.defaults$min.events ) calc <- values
            params <- list(
                tail.fraction = gp[[ paste0( "tail.fraction", suffix ) ]] %||%
                    gp[[ "tail.fraction" ]] %||% ag.defaults$tail.fraction,
                mixture.components = gp[[ paste0( "mixture.components", suffix ) ]] %||%
                    gp[[ "mixture.components" ]] %||% 2L,
                mixture.cut = gp[[ paste0( "mixture.cut", suffix ) ]] %||%
                    gp[[ "mixture.cut" ]] %||% "outer" )
            t <- scaled.threshold( ctx, compute.threshold( calc, algo, params ), gp, suffix, channel )
        }
        gp[[ paste0( "fixed.threshold.", axis ) ]] <- as.numeric( t )
    }
    out$gate.param <- gp
    list( gate.def = out, reason = "" )
}


#' Distance a recalculated boundary moved
#'
#' @title Boundary Drift
#'
#' @param gate.def Effective gate definition.
#' @param boundary,reference Boundary entries.
#'
#' @return Largest threshold shift (threshold gates) or centroid distance
#'   (polygon gates), in transformed units; \code{NA} when not comparable.
#'
#' @keywords internal
boundary.drift <- function( gate.def, boundary, reference ) {
    if ( length( boundary ) == 0L || length( reference ) == 0L ) return( NA_real_ )
    if ( gate.type.of( gate.def ) %in% ag.threshold.types ) {
        a <- boundary[[ 1L ]]
        b <- reference[[ 1L ]]
        shifts <- c( if ( has.value( a$threshold.x ) && has.value( b$threshold.x ) )
                         abs( a$threshold.x - b$threshold.x ),
                     if ( has.value( a$threshold.y ) && has.value( b$threshold.y ) )
                         abs( a$threshold.y - b$threshold.y ) )
        return( if ( length( shifts ) ) max( shifts ) else NA_real_ )
    }
    shifts <- numeric( 0 )
    for ( p in names( boundary ) ) {
        pa <- boundary[[ p ]]$boundary
        pb <- reference[[ p ]]$boundary
        if ( is.matrix( pa ) && is.matrix( pb ) && nrow( pa ) && nrow( pb ) )
            shifts <- c( shifts, sqrt( sum( ( colMeans( pa ) - colMeans( pb ) )^2 ) ) )
    }
    if ( length( shifts ) ) max( shifts ) else NA_real_
}


#' Run a gate hierarchy on one event matrix
#'
#' @title Run Gating
#'
#' @description
#' Without \code{reference} every gate is calculated from the data
#' (training).  With \code{reference}, each gate uses its mode from
#' \code{modes} (default \code{default.mode}): "fixed" applies the
#' reference boundary; "recalculate" calculates a new one, falls back to the
#' reference when that fails, and flags a shift larger than
#' \code{drift.limit}.  Each gate is applied to its parent population only,
#' in hierarchy order.
#'
#' @param gate.defs List of gate definitions.
#' @param data Event matrix (transformed), columns named by channel.
#' @param reference Named list gate -> boundary entries, or \code{NULL}.
#' @param modes Named list gate -> "fixed" or "recalculate".
#' @param default.mode Mode of gates not in \code{modes}.
#' @param drift.limit Largest drift not flagged.
#' @param drift.action "flag", or "reference" to use the reference boundary
#'   for flagged gates.
#' @param axis.limits Named list of display limits, or \code{NULL}.
#' @param transforms Named list of transform parameters (for threshold
#'   factors).
#'
#' @return Gating result (see the file header of \code{gate_application.R}).
#'
#' @export
run.gating <- function( gate.defs, data, reference = NULL, modes = list(),
                        default.mode = "fixed", drift.limit = ag.defaults$drift.limit,
                        drift.action = "flag", axis.limits = NULL, transforms = NULL ) {
    n <- nrow( data )
    by.name <- gates.by.name( gate.defs )
    ctx <- gate.context( axis.limits, transforms )
    result <- list( boundaries = list(), events = list(), parent.events = list(),
                    status = list() )
    all.rows <- seq_len( n )

    for ( gate.def in ordered.gate.defs( gate.defs ) ) {
        name <- gate.def$gate.name
        eff  <- effective.gate.def( gate.def, by.name )
        status <- list( source = NA_character_, drift = NA_real_, flag = NA_character_, message = "" )

        if ( has.value( gate.def$parent.gate ) ) {
            parent.idx <- result$events[[ gate.def$parent.gate ]][[ gate.def$parent.popul %||% "" ]]
            if ( is.null( parent.idx ) ) {
                status$flag <- "skipped"
                status$message <- sprintf( "parent population %s/%s is not available",
                                           gate.def$parent.gate, gate.def$parent.popul )
                result$status[[ name ]] <- status
                next
            }
        } else parent.idx <- all.rows
        parent.data <- data[ parent.idx, , drop = FALSE ]

        ref  <- reference[[ name ]]
        mode <- if ( is.null( reference ) ) "recalculate" else ( modes[[ name ]] %||% default.mode )
        boundary <- list()

        if ( mode == "fixed" ) {
            boundary <- ref %||% list()
            status$source <- "fixed"
            if ( length( boundary ) == 0L ) {
                status$flag <- "failed"
                status$message <- "no stored boundary"
            }
        } else {
            boundary <- tryCatch( {
                res <- resolve.replicates( eff, data, result, by.name, ctx )
                if ( is.null( res$gate.def ) ) {
                    status$message <- res$reason
                    list()
                } else {
                    notes <- ctx$notes
                    b <- calculate.gate( parent.data, res$gate.def, ctx, trained = result$boundaries )
                    notes <- c( notes, ctx$notes )
                    if ( length( notes ) ) status$message <- paste( notes, collapse = "; " )
                    b
                }
            }, error = function( e ) {
                status$message <<- conditionMessage( e )
                list()
            } )
            ctx$notes <- character( 0 )
            status$source <- if ( is.null( reference ) ) "trained" else "recalculated"
            if ( length( boundary ) == 0L ) {
                if ( length( ref ) ) {
                    boundary <- ref
                    status$source <- "reference"
                }
                status$flag <- "failed"
                if ( !nzchar( status$message ) ) status$message <- "too few events to calculate"
            } else if ( length( ref ) ) {
                drift <- boundary.drift( eff, boundary, ref )
                status$drift <- drift
                if ( !is.na( drift ) && drift > drift.limit ) {
                    status$flag <- "drift"
                    status$message <- sprintf( "moved %.3g from the reference", drift )
                    if ( drift.action == "reference" ) {
                        boundary <- ref
                        status$source <- "reference"
                    }
                }
            }
        }

        result$status[[ name ]] <- status
        if ( length( boundary ) == 0L ) next
        masks <- population.masks( eff, boundary, parent.data, ctx$axis.limits )
        if ( length( masks ) == 0L ) {
            if ( is.na( status$flag ) ) status$flag <- "failed"
            if ( !nzchar( status$message ) ) status$message <- "gate channels not in the data"
            result$status[[ name ]] <- status
            next
        }
        n.parent <- length( parent.idx )
        events <- list()
        for ( p in names( masks ) ) {
            events[[ p ]] <- parent.idx[ masks[[ p ]] ]
            if ( !is.null( boundary[[ p ]] ) ) {
                k <- sum( masks[[ p ]] )
                boundary[[ p ]]$n.events <- k
                boundary[[ p ]]$fraction <- if ( n.parent ) k / n.parent else 0
            }
        }
        result$events[[ name ]] <- events
        result$parent.events[[ name ]] <- parent.idx
        result$boundaries[[ name ]] <- boundary
    }
    result
}


#' Train gate boundaries on pooled events
#'
#' @title Train Boundaries
#'
#' @description
#' Calculates every gate on pooled training events (see
#' \code{\link{run.gating}}); gates in \code{done} are applied as they are.
#'
#' @inheritParams run.gating
#' @param done Named list gate -> boundary entries already trained.
#'
#' @return Gating result; status source "trained" for every calculated gate.
#'
#' @export
train.boundaries <- function( gate.defs, data, axis.limits = NULL, done = NULL,
                              transforms = NULL ) {
    if ( length( done ) == 0L )
        return( run.gating( gate.defs, data, reference = NULL, axis.limits = axis.limits,
                            transforms = transforms ) )
    nm <- vapply( gate.defs, function( g ) g$gate.name, character( 1 ) )
    modes <- stats::setNames( as.list( ifelse( nm %in% names( done ), "fixed", "recalculate" ) ), nm )
    result <- run.gating( gate.defs, data, reference = done, modes = modes,
                          axis.limits = axis.limits, transforms = transforms )
    for ( g in names( result$status ) )
        if ( result$status[[ g ]]$source %in% c( "recalculated", "fixed" ) )
            result$status[[ g ]]$source <- "trained"
    result
}


#' Population counts and fractions
#'
#' @title Population Table
#'
#' @description
#' One row per population: key ("gate/pop"), gate, population, label,
#' n.events, parent.key, n.parent, fraction.of.parent (NaN for an empty
#' parent) and, when the gate's \code{stats.parent} names a reference for
#' the population ("root", "gate/pop" or a gate name), stats.parent,
#' n.stats.parent and fraction.of.stats.parent (NaN with
#' stats.parent.missing when the reference is not in this sample).
#'
#' @param gate.defs List of gate definitions.
#' @param result Gating result.
#' @param n.total Number of events in the data (for "root").
#'
#' @return \code{data.frame}.
#'
#' @export
population.table <- function( gate.defs, result, n.total ) {
    rows <- list()
    for ( g in ordered.gate.defs( gate.defs ) ) {
        name <- g$gate.name
        events <- result$events[[ name ]]
        if ( length( events ) == 0L ) next
        n.parent <- length( result$parent.events[[ name ]] )
        parent.key <- if ( has.value( g$parent.gate ) )
            population.key( g$parent.gate, g$parent.popul %||% "" ) else "root"
        for ( p in names( events ) ) {
            k <- length( events[[ p ]] )
            lab <- g$populations[[ p ]][[ "label" ]]
            row <- list( key = population.key( name, p ), gate = name, population = p,
                         label = if ( has.value( lab ) ) lab else p, n.events = k,
                         parent.key = parent.key, n.parent = n.parent,
                         fraction.of.parent = if ( n.parent ) k / n.parent else NaN,
                         stats.parent = NA_character_, n.stats.parent = NA_integer_,
                         fraction.of.stats.parent = NA_real_, stats.parent.missing = FALSE )
            spec <- if ( is.list( g$stats.parent ) ) g$stats.parent[[ p ]] else NULL
            if ( has.value( spec ) ) {
                row$stats.parent <- spec
                sp.idx <- if ( spec == "root" ) seq_len( n.total ) else if ( grepl( "/", spec, fixed = TRUE ) ) {
                    key <- parse.population.key( spec )
                    result$events[[ key[ 1 ] ]][[ key[ 2 ] ]]
                } else result$parent.events[[ spec ]]
                if ( is.null( sp.idx ) ) {
                    row$fraction.of.stats.parent <- NaN
                    row$stats.parent.missing <- TRUE
                } else {
                    row$n.stats.parent <- length( sp.idx )
                    row$fraction.of.stats.parent <- if ( length( sp.idx ) ) k / length( sp.idx ) else NaN
                }
            }
            rows[[ length( rows ) + 1L ]] <- row
        }
    }
    if ( length( rows ) == 0L ) return( data.frame() )
    do.call( rbind, lapply( rows, as.data.frame, stringsAsFactors = FALSE ) )
}


#' Median of each channel in each population
#'
#' @title Population Medians
#'
#' @param gate.defs List of gate definitions.
#' @param result Gating result.
#' @param data Event matrix the result was computed on.
#' @param channels Channels to summarise.
#' @param min.events Populations with fewer events are left out.
#'
#' @return \code{data.frame} with key, channel and median.
#'
#' @export
population.medians <- function( gate.defs, result, data, channels, min.events = 1L ) {
    channels <- intersect( channels, colnames( data ) )
    rows <- list()
    for ( g in gate.defs ) {
        for ( p in names( result$events[[ g$gate.name ]] ) ) {
            idx <- result$events[[ g$gate.name ]][[ p ]]
            if ( length( idx ) < min.events || length( channels ) == 0L ) next
            med <- apply( data[ idx, channels, drop = FALSE ], 2L, stats::median )
            rows[[ length( rows ) + 1L ]] <- data.frame(
                key = population.key( g$gate.name, p ), channel = channels,
                median = unname( med ), stringsAsFactors = FALSE )
        }
    }
    if ( length( rows ) == 0L ) return( data.frame() )
    do.call( rbind, rows )
}
