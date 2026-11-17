# gate_model.R
#
# Gate definitions, population regions and the gate hierarchy, following
# the Honeychrome Automated Gating model (.agmodel) so that models move
# between the two.
#
# A gate definition is a named list:
#   gate.name, gate.number
#   gate.type       "singlets", "1dsep", "2dsep", "free" or "replicate"
#   gate.marker.x   separation channel (1dsep) or x channel
#   gate.marker.y   y channel; for 1dsep the channel it is viewed against
#                   (NULL: a histogram)
#   sep.axis        "x" or "y": plot axis of a 1dsep separation channel
#   parent.gate, parent.popul
#   populations     named list; each population has label, label.pos,
#                   region and optionally reference / anchor (TRUE)
#   algorithm       see ag.algorithms.by.type
#   gate.param      named list of algorithm settings
#   stats.parent    named list: population -> "root", "gate/pop" or gate
#   origin.gate     source of a "replicate" gate
#   template        shapes as drawn, for algorithm "imported"
#
# Regions: "below" / "above" the threshold of a 1dsep gate; quadrants
# "x-y-", "x+y-", "x-y+", "x+y+" of a 2dsep gate, joined with "|" for
# several (three make an L); "inside" for polygon gates.
#
# Threshold gates are open-ended: a population is every event on its side
# of the thresholds, however far outside the display range.


ag.gate.types <- c( "singlets", "1dsep", "2dsep", "free", "replicate" )
ag.threshold.types <- c( "1dsep", "2dsep" )
ag.polygon.types <- c( "singlets", "free" )
ag.threshold.algorithms <- c( "mixture", "tail", "otsu", "bimodal", "kde_min" )
ag.algorithms.by.type <- list(
    singlets  = c( "singlets", "imported" ),
    `1dsep`   = c( "mixture", "tail", "replicate", "otsu", "bimodal", "kde_min",
                   "imported" ),
    `2dsep`   = c( "mixture", "tail", "replicate", "otsu", "bimodal", "kde_min",
                   "imported" ),
    free      = c( "cells", "all", "imported" ),
    replicate = "replicate"
)
ag.quadrant.regions <- c( "x-y-", "x+y-", "x-y+", "x+y+" )

ag.defaults <- list(
    tail.fraction         = 1e-4,
    threshold.factor      = 1,
    edge.density          = 0.002,
    debris.algorithm      = "mixture",
    kde.min.fallback.tail = 0.01,
    drift.limit           = 0.05,
    min.events            = 10L,
    kde.exact.max.events  = 20000L,
    quantile.trim         = 0.001,
    axis.hi               = 0.99
)


#' Test whether a value is set
#'
#' @title Value Is Set
#'
#' @param x Any value.
#'
#' @return \code{TRUE} for a length-1 value that is not \code{NULL},
#'   \code{NA} or an empty string.
#'
#' @keywords internal
has.value <- function( x ) {
    !is.null( x ) && length( x ) == 1L && !is.na( x ) && !identical( x, "" )
}


#' Type of a gate
#'
#' @title Gate Type
#'
#' @param gate.def Gate definition.
#'
#' @return \code{gate.type}, or \code{"free"} when unset.
#'
#' @keywords internal
gate.type.of <- function( gate.def ) {
    if ( has.value( gate.def$gate.type ) ) gate.def$gate.type else "free"
}


#' Quadrants of a 2dsep region
#'
#' @title Region Quadrants
#'
#' @param region Character.  One quadrant (\code{"x+y-"}) or several joined
#'   with \code{"|"}.
#'
#' @return Character vector of quadrants in \code{ag.quadrant.regions}
#'   order.
#'
#' @keywords internal
region.quadrants <- function( region ) {
    if ( !has.value( region ) )
        return( character( 0 ) )
    parts <- strsplit( region, "|", fixed = TRUE )[[ 1 ]]
    ag.quadrant.regions[ ag.quadrant.regions %in% parts ]
}


#' Region of a set of quadrants
#'
#' @title Join Quadrants
#'
#' @param quadrants Character vector of quadrants.
#'
#' @return The quadrants in \code{ag.quadrant.regions} order joined with
#'   \code{"|"}.
#'
#' @keywords internal
join.quadrants <- function( quadrants ) {
    paste( ag.quadrant.regions[ ag.quadrant.regions %in% quadrants ], collapse = "|" )
}


#' Side of an axis threshold a 2dsep region lies on
#'
#' @title Region Sign
#'
#' @param region Character region.
#' @param axis \code{"x"} or \code{"y"}.
#'
#' @return \code{"+"} or \code{"-"} when every quadrant of the region lies
#'   on that side, else \code{NULL}.
#'
#' @keywords internal
region.sign <- function( region, axis ) {
    k <- if ( axis == "x" ) 2L else 4L
    signs <- unique( substr( region.quadrants( region ), k, k ) )
    if ( length( signs ) == 1L ) signs else NULL
}


#' Events of a 2dsep region
#'
#' @title Quadrant Mask
#'
#' @param region Character region.
#' @param above.x,above.y Logical vectors: events at or above each threshold.
#'
#' @return Logical vector.
#'
#' @keywords internal
quadrant.mask <- function( region, above.x, above.y ) {
    out <- logical( length( above.x ) )
    for ( q in region.quadrants( region ) ) {
        mx <- if ( substr( q, 2L, 2L ) == "+" ) above.x else !above.x
        my <- if ( substr( q, 4L, 4L ) == "+" ) above.y else !above.y
        out <- out | ( mx & my )
    }
    out
}


#' Closed rectangle
#'
#' @title Rectangle Polygon
#'
#' @param x.lo,x.hi,y.lo,y.hi Numeric limits.
#'
#' @return 5 x 2 matrix, first vertex repeated at the end.
#'
#' @keywords internal
rect.polygon <- function( x.lo, x.hi, y.lo = 0, y.hi = 1 ) {
    cbind( c( x.lo, x.hi, x.hi, x.lo, x.lo ), c( y.lo, y.lo, y.hi, y.hi, y.lo ) )
}


#' Outline of a 2dsep region within display limits
#'
#' @title Region Polygon
#'
#' @description
#' A rectangle for one quadrant or a whole side, an L for three quadrants,
#' two rectangles meeting at the threshold crossing for opposite quadrants.
#'
#' @param region Character region.
#' @param tx,ty Thresholds.
#' @param x.lo,x.hi,y.lo,y.hi Display limits.
#'
#' @return Closed polygon (n x 2 matrix), or an empty matrix.
#'
#' @keywords internal
region.polygon <- function( region, tx, ty, x.lo, x.hi, y.lo, y.hi ) {
    quads <- region.quadrants( region )
    xs <- list( `-` = c( x.lo, tx ), `+` = c( tx, x.hi ) )
    ys <- list( `-` = c( y.lo, ty ), `+` = c( ty, y.hi ) )
    sx <- substr( quads, 2L, 2L )
    sy <- substr( quads, 4L, 4L )

    if ( length( quads ) == 4L )
        return( rect.polygon( x.lo, x.hi, y.lo, y.hi ) )
    if ( length( quads ) == 1L )
        return( rect.polygon( xs[[ sx ]][ 1 ], xs[[ sx ]][ 2 ], ys[[ sy ]][ 1 ], ys[[ sy ]][ 2 ] ) )
    if ( length( quads ) == 2L ) {
        if ( length( unique( sx ) ) == 1L )
            return( rect.polygon( xs[[ sx[ 1 ] ]][ 1 ], xs[[ sx[ 1 ] ]][ 2 ], y.lo, y.hi ) )
        if ( length( unique( sy ) ) == 1L )
            return( rect.polygon( x.lo, x.hi, ys[[ sy[ 1 ] ]][ 1 ], ys[[ sy[ 1 ] ]][ 2 ] ) )
        qa <- sort( quads, method = "radix" )
        ra <- rect.polygon( xs[[ substr( qa[ 1 ], 2, 2 ) ]][ 1 ], xs[[ substr( qa[ 1 ], 2, 2 ) ]][ 2 ],
                            ys[[ substr( qa[ 1 ], 4, 4 ) ]][ 1 ], ys[[ substr( qa[ 1 ], 4, 4 ) ]][ 2 ] )
        rb <- rect.polygon( xs[[ substr( qa[ 2 ], 2, 2 ) ]][ 1 ], xs[[ substr( qa[ 2 ], 2, 2 ) ]][ 2 ],
                            ys[[ substr( qa[ 2 ], 4, 4 ) ]][ 1 ], ys[[ substr( qa[ 2 ], 4, 4 ) ]][ 2 ] )
        from.corner <- function( r ) {
            pts <- r[ -nrow( r ), , drop = FALSE ]
            k <- which.min( abs( pts[ , 1 ] - tx ) + abs( pts[ , 2 ] - ty ) )
            pts[ c( k:nrow( pts ), seq_len( k - 1L ) ), , drop = FALSE ]
        }
        a <- from.corner( ra )
        b <- from.corner( rb )
        return( rbind( a, a[ 1, ], b[ -1, , drop = FALSE ], b[ 1, ] ) )
    }
    if ( length( quads ) != 3L )
        return( matrix( numeric( 0 ), ncol = 2L ) )
    missing <- setdiff( ag.quadrant.regions, quads )
    mx <- substr( missing, 2L, 2L )
    my <- substr( missing, 4L, 4L )
    cx <- if ( mx == "-" ) x.lo else x.hi
    cy <- if ( my == "-" ) y.lo else y.hi
    ox <- if ( mx == "-" ) x.hi else x.lo
    oy <- if ( my == "-" ) y.hi else y.lo
    cbind( c( tx, ox, ox, cx, cx, tx, tx ), c( cy, cy, oy, oy, ty, ty, cy ) )
}


#' Default populations of a new gate
#'
#' @title Default Populations for a Gate Type
#'
#' @param gate.type Character gate type.
#'
#' @return Named list of populations with label, label.pos and region.
#'
#' @keywords internal
default.populations.for.type <- function( gate.type ) {
    pop <- function( label, pos, region ) list( label = label, label.pos = pos, region = region )
    if ( gate.type == "1dsep" )
        return( list( neg = pop( "neg", 1L, "below" ), pos = pop( "pos", 2L, "above" ) ) )
    if ( gate.type == "2dsep" )
        return( list( DN = pop( "DN", 3L, "x-y-" ), `X+` = pop( "X+", 4L, "x+y-" ),
                      `Y+` = pop( "Y+", 1L, "x-y+" ), DP = pop( "DP", 2L, "x+y+" ) ) )
    if ( gate.type %in% ag.polygon.types )
        return( list( `in` = pop( "in", 1L, "inside" ) ) )
    list()
}


#' Default lower crop of a threshold calculation
#'
#' @title Default Axis Crop
#'
#' @description
#' \code{tail} mirrors the events below the negative population's mode, so
#' it starts uncropped; other algorithms drop the lowest 1\% of the range.
#'
#' @param algorithm Character algorithm name.
#'
#' @return Numeric fraction of the range.
#'
#' @keywords internal
default.axis.lo <- function( algorithm ) {
    if ( identical( algorithm, "tail" ) ) 0 else 0.01
}


#' Algorithm of one axis of a threshold gate
#'
#' @title Axis Algorithm
#'
#' @param gate.def Gate definition.
#' @param axis \code{"x"} or \code{"y"}.
#'
#' @return Character algorithm: \code{algorithm.x}/\code{algorithm.y} of a
#'   2dsep gate's parameters, else the gate's algorithm.
#'
#' @keywords internal
axis.algorithm <- function( gate.def, axis ) {
    default <- if ( has.value( gate.def$algorithm ) ) gate.def$algorithm else "tail"
    if ( !identical( gate.def$gate.type, "2dsep" ) )
        return( default )
    a <- gate.def$gate.param[[ paste0( "algorithm.", axis ) ]]
    if ( has.value( a ) ) a else default
}


#' Separation axes of a threshold gate
#'
#' @title Threshold Axes
#'
#' @param gate.def Gate definition.
#'
#' @return \code{"x"}, or \code{c("x", "y")} for a 2dsep gate.
#'
#' @keywords internal
threshold.axes <- function( gate.def ) {
    if ( identical( gate.def$gate.type, "2dsep" ) ) c( "x", "y" ) else "x"
}


#' Channel of a gate axis
#'
#' @title Axis Channel
#'
#' @param gate.def Gate definition.
#' @param axis \code{"x"} or \code{"y"}.
#'
#' @return Channel name or \code{NULL}.
#'
#' @keywords internal
axis.channel <- function( gate.def, axis ) {
    if ( axis == "x" ) gate.def$gate.marker.x else gate.def$gate.marker.y
}


#' Replicated threshold sources of a gate
#'
#' @title Replicate Sources
#'
#' @description
#' For each axis of a threshold gate whose algorithm is \code{"replicate"},
#' its source from \code{gate.param$replicate.x} / \code{replicate.y}:
#' \code{list(gate = name)} copies another gate's threshold on the same
#' channel; \code{list(population = "gate/pop", algorithm = a)} calculates
#' it with algorithm \code{a} on the events of that (reference) population.
#'
#' @param gate.def Gate definition.
#'
#' @return Named list axis -> source list, or \code{NULL} for a missing
#'   source.
#'
#' @keywords internal
replicate.sources <- function( gate.def ) {
    if ( !( gate.type.of( gate.def ) %in% ag.threshold.types ) ||
         identical( gate.def$algorithm, "imported" ) )
        return( list() )
    out <- list()
    for ( axis in threshold.axes( gate.def ) ) {
        if ( axis.algorithm( gate.def, axis ) != "replicate" )
            next
        src <- gate.def$gate.param[[ paste0( "replicate.", axis ) ]]
        ok  <- is.list( src ) && ( has.value( src$gate ) || has.value( src$population ) )
        out[ axis ] <- list( if ( ok ) src else NULL )
    }
    out
}


#' Population key
#'
#' @title Population Key
#'
#' @param gate.name,pop.name Character.
#'
#' @return \code{"gate/pop"}.
#'
#' @keywords internal
population.key <- function( gate.name, pop.name ) paste0( gate.name, "/", pop.name )


#' Gate and population of a population key
#'
#' @title Parse Population Key
#'
#' @param key Character \code{"gate/pop"} (split at the first \code{"/"}).
#'
#' @return Character vector \code{c(gate, pop)}.
#'
#' @keywords internal
parse.population.key <- function( key ) {
    key <- as.character( key )
    k <- regexpr( "/", key, fixed = TRUE )
    if ( k < 0 ) return( c( key, "" ) )
    c( substr( key, 1L, k - 1L ), substr( key, k + 1L, nchar( key ) ) )
}


#' Gates a gate depends on
#'
#' @title Gate Dependencies
#'
#' @param gate.def Gate definition.
#'
#' @return Character vector: parent, replicate origin, and the gates or
#'   reference populations its replicated thresholds come from.
#'
#' @keywords internal
gate.dependencies <- function( gate.def ) {
    deps <- c( gate.def$parent.gate, gate.def$origin.gate )
    for ( src in replicate.sources( gate.def ) ) {
        if ( is.null( src ) ) next
        if ( has.value( src$gate ) ) deps <- c( deps, src$gate )
        else if ( has.value( src$population ) )
            deps <- c( deps, parse.population.key( src$population )[ 1 ] )
    }
    unique( deps[ vapply( deps, has.value, logical( 1 ) ) ] )
}


#' Gate definitions keyed by name
#'
#' @title Gates by Name
#'
#' @param gate.defs List of gate definitions.
#'
#' @return Named list.
#'
#' @keywords internal
gates.by.name <- function( gate.defs ) {
    nm <- vapply( gate.defs, function( g ) if ( has.value( g$gate.name ) ) g$gate.name else "",
                  character( 1 ) )
    stats::setNames( gate.defs[ nm != "" ], nm[ nm != "" ] )
}


#' Gate definitions in dependency order
#'
#' @title Ordered Gate Definitions
#'
#' @description
#' Orders the gates so every gate follows the gates it depends on (see
#' \code{\link{gate.dependencies}}); ties are broken by \code{gate.number}
#' then name.  Stops on a cycle.
#'
#' @param gate.defs List of gate definitions.
#'
#' @return Reordered list.
#'
#' @keywords internal
ordered.gate.defs <- function( gate.defs ) {
    if ( length( gate.defs ) == 0L )
        return( list() )
    by.name <- gates.by.name( gate.defs )
    num <- vapply( gate.defs, function( g ) {
        v <- suppressWarnings( as.integer( g$gate.number ) )
        if ( length( v ) == 1L && !is.na( v ) ) v else 1000000L
    }, integer( 1 ) )
    nm  <- vapply( gate.defs, function( g ) as.character( g$gate.name %||% "" ), character( 1 ) )
    rank.order   <- order( num, nm, method = "radix" )
    ranked       <- gate.defs[ rank.order ]
    names.ranked <- nm[ rank.order ]

    deps <- lapply( ranked, function( g ) intersect( gate.dependencies( g ), names( by.name ) ) )
    n.unmet <- lengths( deps )
    dependents <- if ( sum( n.unmet ) > 0L )
        split( rep( seq_along( deps ), n.unmet ), unlist( deps ) ) else list()

    ready  <- which( n.unmet == 0L )
    done   <- logical( length( ranked ) )
    placed <- character( 0 )
    out    <- vector( "list", 0L )
    while ( length( ready ) > 0L ) {
        r <- min( ready )
        ready <- ready[ ready != r ]
        done[ r ] <- TRUE
        out[[ length( out ) + 1L ]] <- ranked[[ r ]]
        nm.r <- names.ranked[ r ]
        if ( nm.r %in% placed ) next
        placed <- c( placed, nm.r )
        for ( d in dependents[[ nm.r ]] ) {
            n.unmet[ d ] <- n.unmet[ d ] - 1L
            if ( n.unmet[ d ] == 0L ) ready <- c( ready, d )
        }
    }
    if ( !all( done ) )
        stop( sprintf( "Gate hierarchy has a cycle among: %s",
                       paste( names.ranked[ !done ], collapse = ", " ) ) )
    out
}


#' Effective definition of a replicate gate
#'
#' @title Effective Gate Definition
#'
#' @description
#' For a \code{"replicate"} gate, the origin's type, axes and populations
#' under the replicate's own name and parent (following chains of
#' replicates); other gates are returned unchanged.
#'
#' @param gate.def Gate definition.
#' @param by.name Named list of all gate definitions.
#'
#' @return Gate definition.
#'
#' @keywords internal
effective.gate.def <- function( gate.def, by.name ) {
    if ( !identical( gate.def$gate.type, "replicate" ) )
        return( gate.def )
    origin <- by.name[[ gate.def$origin.gate %||% "" ]]
    seen <- gate.def$gate.name
    while ( !is.null( origin ) && identical( origin$gate.type, "replicate" ) ) {
        if ( origin$gate.name %in% seen )
            return( gate.def )
        seen <- c( seen, origin$gate.name )
        origin <- by.name[[ origin$origin.gate %||% "" ]]
    }
    if ( is.null( origin ) )
        return( gate.def )
    merged <- origin
    for ( key in c( "gate.name", "gate.number", "parent.gate", "parent.popul",
                    "stats.parent", "origin.gate" ) )
        merged[ key ] <- list( gate.def[[ key ]] )
    merged$replicate.of <- origin$gate.name
    merged
}


#' Channels a gate reads
#'
#' @title Gate Channels
#'
#' @param gate.def Gate definition.
#'
#' @return Character vector: x, plus y for 2dsep, singlets and free gates.
#'
#' @keywords internal
gate.channels <- function( gate.def ) {
    chans <- gate.def$gate.marker.x
    if ( gate.type.of( gate.def ) %in% c( "2dsep", "singlets", "free" ) &&
         has.value( gate.def$gate.marker.y ) )
        chans <- c( chans, gate.def$gate.marker.y )
    chans[ vapply( chans, has.value, logical( 1 ) ) ]
}


#' Plot axis of a 1dsep separation channel
#'
#' @title Separation Axis
#'
#' @param gate.def Gate definition.
#'
#' @return \code{"x"} or \code{"y"}; always \code{"x"} for other gate
#'   types and histograms.
#'
#' @keywords internal
separation.axis <- function( gate.def ) {
    if ( !identical( gate.def$gate.type, "1dsep" ) || !has.value( gate.def$gate.marker.y ) )
        return( "x" )
    if ( identical( gate.def$sep.axis, "y" ) ) "y" else "x"
}


#' Whether boundary entries are shapes as drawn
#'
#' @title Is Shape Boundary
#'
#' @param boundary Named list of boundary entries.
#'
#' @return Logical.
#'
#' @keywords internal
is.shape.boundary <- function( boundary ) {
    length( boundary ) > 0L &&
        all( vapply( boundary, function( e ) is.list( e ) && length( e$shape.channels ) > 0L,
                     logical( 1 ) ) )
}


#' Channels of an imported gate's shapes
#'
#' @title Shape Channels
#'
#' @param gate.def Gate definition.
#'
#' @return Character vector (empty when the template holds no shapes).
#'
#' @keywords internal
shape.channels <- function( gate.def ) {
    template <- gate.def$template
    if ( !is.shape.boundary( template ) )
        return( character( 0 ) )
    ch <- unlist( template[[ 1L ]]$shape.channels )
    ch[ vapply( ch, has.value, logical( 1 ) ) ]
}


#' Every channel a gate reads or is drawn on
#'
#' @title Referenced Channels
#'
#' @param gate.def Gate definition.
#'
#' @return Character vector.
#'
#' @keywords internal
referenced.channels <- function( gate.def ) {
    out <- gate.channels( gate.def )
    for ( ch in c( gate.def$gate.marker.y, shape.channels( gate.def ) ) )
        if ( has.value( ch ) && !( ch %in% out ) ) out <- c( out, ch )
    out
}


#' Region of each population of a gate
#'
#' @title Population Regions
#'
#' @description
#' Regions come from the boundary entries or population definitions when
#' set; otherwise they are inferred: a 1dsep population takes the side the
#' others leave free, else population order decides (first below, the rest
#' above); 2dsep populations follow DN, X+, Y+, DP order; polygon gates are
#' "inside".
#'
#' @param gate.def Gate definition.
#' @param boundary Optional named list of boundary entries.
#'
#' @return Named character vector population -> region.
#'
#' @keywords internal
population.regions <- function( gate.def, boundary = NULL ) {
    pops <- names( gate.def$populations )
    if ( length( pops ) == 0L ) pops <- names( boundary )
    regions <- character( 0 )
    for ( p in pops ) {
        r <- boundary[[ p ]]$region
        if ( !has.value( r ) ) r <- gate.def$populations[[ p ]]$region
        if ( has.value( r ) ) regions[ p ] <- r
    }
    sides <- unique( regions[ regions %in% c( "below", "above" ) ] )
    free.side <- if ( length( sides ) == 1L ) setdiff( c( "below", "above" ), sides ) else NULL
    for ( i in seq_along( pops ) ) {
        p <- pops[ i ]
        if ( p %in% names( regions ) ) next
        r <- if ( identical( gate.def$gate.type, "1dsep" ) ) {
            if ( !is.null( free.side ) ) free.side else if ( i == 1L ) "below" else "above"
        } else if ( identical( gate.def$gate.type, "2dsep" ) ) {
            if ( i <= 4L ) ag.quadrant.regions[ i ] else NA_character_
        } else "inside"
        if ( !is.na( r ) ) regions[ p ] <- r
    }
    regions[ pops[ pops %in% names( regions ) ] ]
}


#' Viability dye channels
#'
#' @title Viability Channels
#'
#' @description
#' Channels whose fluorophore is a viability dye in AutoSpectral's
#' fluorophore database (\code{is.viability}), matched by name and
#' synonyms.  With a panel, each antigen is looked up by its dye as well.
#'
#' @param channels Character vector of channel (antigen) names.
#' @param panel Optional panel \code{data.frame} with \code{dye} and
#'   \code{antigen}.
#'
#' @return Character vector, a subset of \code{channels}.
#'
#' @importFrom AutoSpectral match.fluorophores
#'
#' @keywords internal
viability.channels <- function( channels, panel = NULL ) {
    channels <- unique( channels[ vapply( channels, has.value, logical( 1 ) ) ] )
    if ( length( channels ) == 0L )
        return( character( 0 ) )
    db.file <- system.file( "extdata", "fluorophore_database.csv", package = "AutoSpectral" )
    if ( db.file == "" )
        return( character( 0 ) )
    db <- utils::read.csv( db.file, stringsAsFactors = FALSE )
    viable <- db$fluorophore[ toupper( trimws( as.character( db$is.viability ) ) ) == "TRUE" ]
    lookup <- channels
    if ( !is.null( panel ) ) {
        dye <- panel$dye[ match( channels, panel$antigen ) ]
        lookup <- ifelse( is.na( dye ), channels, dye )
    }
    hit <- vapply( seq_along( channels ), function( i ) {
        m <- c( match.fluorophores( lookup[ i ], db, verbose = FALSE ),
                match.fluorophores( channels[ i ], db, verbose = FALSE ) )
        any( m %in% viable )
    }, logical( 1 ) )
    channels[ hit ]
}


#' Anchor population of a gate hierarchy
#'
#' @title Anchor Population
#'
#' @description
#' The population marked \code{anchor = TRUE}; else the viable population
#' of the deepest gate on a viability dye (its population on the dye's
#' negative side, or its only population); else the deepest singlets
#' population; else the population labelled "Cells".  The anchor is the
#' default reference population for replicated thresholds and the default
#' abundance denominator.
#'
#' @param gate.defs List of gate definitions.
#' @param panel Optional panel, for dye lookup (see
#'   \code{\link{viability.channels}}).
#'
#' @return Population key \code{"gate/pop"}, or \code{NULL}.
#'
#' @export
anchor.population <- function( gate.defs, panel = NULL ) {
    for ( g in gate.defs )
        for ( p in names( g$populations ) )
            if ( isTRUE( g$populations[[ p ]]$anchor ) )
                return( population.key( g$gate.name, p ) )

    ordered <- ordered.gate.defs( gate.defs )
    by.name <- gates.by.name( gate.defs )
    depth <- integer( 0 )
    for ( g in ordered ) {
        parent <- g$parent.gate
        depth[ g$gate.name ] <- if ( has.value( parent ) && parent %in% names( by.name ) )
            ( if ( parent %in% names( depth ) ) depth[[ parent ]] else -1L ) + 1L else 0L
    }
    deepest <- function( gates ) {
        if ( length( gates ) == 0L ) return( NULL )
        d <- vapply( gates, function( g ) depth[[ g$gate.name ]], integer( 1 ) )
        gates[[ which.max( d ) ]]
    }
    used <- unique( unlist( lapply( gate.defs, referenced.channels ) ) )
    viability <- viability.channels( used, panel )

    via <- deepest( Filter( function( g ) any( referenced.channels( g ) %in% viability ), ordered ) )
    if ( !is.null( via ) ) {
        regions <- population.regions( via )
        if ( length( regions ) == 1L )
            return( population.key( via$gate.name, names( regions ) ) )
        dye <- referenced.channels( via )[ referenced.channels( via ) %in% viability ][ 1L ]
        for ( p in names( regions ) ) {
            if ( identical( via$gate.type, "1dsep" ) && identical( via$gate.marker.x, dye ) &&
                 regions[[ p ]] == "below" )
                return( population.key( via$gate.name, p ) )
            if ( identical( via$gate.type, "2dsep" ) ) {
                axis <- if ( identical( via$gate.marker.x, dye ) ) "x" else "y"
                if ( identical( region.sign( regions[[ p ]], axis ), "-" ) )
                    return( population.key( via$gate.name, p ) )
            }
        }
        return( population.key( via$gate.name, names( regions )[ 1L ] ) )
    }
    sing <- deepest( Filter( function( g ) identical( g$gate.type, "singlets" ), ordered ) )
    if ( !is.null( sing ) )
        return( population.key( sing$gate.name, names( sing$populations )[ 1L ] %||% "" ) )
    for ( g in ordered )
        for ( p in names( g$populations ) ) {
            lab <- g$populations[[ p ]][[ "label" ]]
            if ( tolower( trimws( if ( has.value( lab ) ) lab else p ) ) == "cells" )
                return( population.key( g$gate.name, p ) )
        }
    NULL
}
