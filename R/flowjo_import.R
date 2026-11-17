# flowjo_import.R
#
# Converts the gates of a FlowJo gating source (read by
# read.flowjo.workspace) into gate definitions, following Honeychrome's
# Gating Model Builder (ag_flowjo).
#
# Every gate keeps its shapes as drawn in template (algorithm "imported",
# applied as they are) and template.original, and records the trainable
# gate it converts to in its type, channels, regions and
# proposed.algorithm:
#   - a shape on the forward-scatter singlet pair (area against height or
#     width of one scatter detector) -> "singlets"; the first one is moved
#     to the top of the hierarchy, with the gates above it placed directly
#     below it;
#   - shapes drawn only on scatter parameters -> "free", proposed "cells"
#     (a debris cut on forward scatter);
#   - sibling gates drawn on the same parameters, with at least one
#     fluorescence parameter, become one gate with one population each;
#     shapes on two fluorescence parameters take the quadrant layout that
#     best reproduces them, and a parameter separating no two shapes is
#     dropped (1dsep on the other);
#   - a histogram gate -> 1dsep with no view channel, proposed "tail", or
#     "mixture" when its positive population held more than 25% of the
#     parent in the recorded counts; other threshold gates are proposed
#     "mixture" with proposal.rule "fraction" (decided on the data by
#     convert.imported.gates);
#   - the side of a fluorescence parameter comes from the drawn
#     coordinates: entirely above the channel's raw zero is the positive
#     side, reaching zero or below the negative side.
#
# Coordinates: FlowJo stores raw units and draws polygon and ellipse edges
# straight on its own display scale.  Rectangle bounds convert exactly;
# polygon edges are densified on FlowJo's display scale first (32 points
# per edge); ellipses become 10-vertex polygons of the same area.  FlowJo
# biex scales are AutoSpectral::biexp.transform with FlowJo's own
# parameters, logicle and log the Gating-ML forms of transform.R, arcsinh
# the Gating-ML formula; hyperlog and anything else are treated as linear
# (Honeychrome uses FlowKit's hyperlog).  Open rectangle sides extend one
# display span beyond the 0 to 1 display range.
#
# Gate names are unique: gates sharing a name are told apart by the
# populations they are drawn on, "TIGIT+ (Tem)", then "TIGIT+ (CD4+ > Tem)"
# (with a single right-pointing angle quotation mark); the shared name is
# kept in name.group.  "/" in FlowJo names becomes a division slash.


fj.ellipse.vertices <- 10L
fj.histogram.mixture.frequency <- 0.25
fj.layout <- list( grid = 80L, min.cover = 0.25, tol.1d = 0.05, empty.cover = 0.1,
                   tie = 0.01, max.extra.1d = 0.09, undrawn.weight = 0.25 )
fj.name.sep <- " \u203a "


#' Display transform of a FlowJo axis
#' @title FlowJo Display Transform
#' @param fjt \code{list(kind, params)} from the workspace, or \code{NULL}.
#' @return List with forward and inverse, or \code{NULL} (linear or not
#'   supported).
#' @keywords internal
fj.display.transform <- function( fjt ) {
    if ( is.null( fjt ) )
        return( NULL )
    p <- vapply( fjt$params, function( v ) fj.num( v ), numeric( 1 ) )
    val <- function( key, default ) {
        v <- if ( key %in% names( p ) ) p[[ key ]] else NA_real_
        if ( is.na( v ) || v == 0 ) default else v
    }
    tryCatch( switch( fjt$kind,
        biex = {
            a <- list( channelRange = 1, maxValue = val( "maxRange", 262144 ),
                       pos = val( "pos", 4.41854 ), neg = val( "neg", 0 ),
                       widthBasis = val( "width", -10 ) )
            list( forward = do.call( biexp.transform, c( a, inverse = FALSE ) ),
                  inverse = do.call( biexp.transform, c( a, inverse = TRUE ) ) )
        },
        logicle = make.transform( list( id = 1L, scale.t = val( "T", 262144 ),
            logicle.w = val( "W", 0.5 ), logicle.m = val( "M", 4.5 ), logicle.a = val( "A", 0 ) ) ),
        log = make.transform( list( id = 2L, scale.t = val( "offset", 1 ),
                                    log.m = val( "decades", 4 ) ) ),
        fasinh = {
            t <- val( "T", 262144 ); m <- val( "M", 4.5 ); a <- val( "A", 0 )
            l10 <- log( 10 )
            list( forward = function( x ) ( asinh( x * sinh( m * l10 ) / t ) + a * l10 ) / ( ( m + a ) * l10 ),
                  inverse = function( y ) t * sinh( y * ( m + a ) * l10 - a * l10 ) / sinh( m * l10 ) )
        },
        NULL ), error = function( e ) NULL )
}


#' Raw and display conversions of one gate axis
#' @title FlowJo Axis
#' @param fj.param FlowJo parameter.
#' @param channel Experiment channel it maps to.
#' @param source Gating source.
#' @param exp Experiment channels.
#' @return List with fj (FlowJo display transform or \code{NULL}),
#'   to.display (raw -> display), linear and zero (raw zero in display
#'   units).
#' @keywords internal
fj.axis <- function( fj.param, channel, source, exp ) {
    params <- exp$transforms[[ channel ]]
    hc <- if ( length( params ) ) make.transform( params ) else NULL
    to.display <- if ( is.null( hc ) ) function( v ) as.numeric( v ) else hc$forward
    list( fj = fj.display.transform( source$transforms[[ fj.param ]] ),
          to.display = to.display,
          linear = length( params ) == 0L || identical( as.integer( params$id ), 0L ),
          zero = to.display( 0 ) )
}


#' Raw values onto FlowJo's display scale and back
#' @title FlowJo Scale
#' @param ax FlowJo axis.
#' @param v Values.
#' @return Values.
#' @keywords internal
fj.to.scale <- function( ax, v ) if ( is.null( ax$fj ) ) v else ax$fj$forward( v )

#' @rdname fj.to.scale
#' @keywords internal
fj.from.scale <- function( ax, v ) if ( is.null( ax$fj ) ) v else ax$fj$inverse( v )


#' A FlowJo polygon on the display scale
#' @title FlowJo Polygon to Display
#' @param raw Raw vertices (n x 2).
#' @param ax,ay FlowJo axes.
#' @param points.per.edge Densification on FlowJo's scale.
#' @return Closed polygon (matrix) in display units.
#' @keywords internal
fj.polygon.display <- function( raw, ax, ay, points.per.edge = 32L ) {
    n <- nrow( raw )
    if ( n && any( abs( raw[ 1L, ] - raw[ n, ] ) > 1e-8 * pmax( 1, abs( raw[ 1L, ] ) ) ) )
        raw <- rbind( raw, raw[ 1L, ] )
    if ( !( is.null( ax$fj ) && is.null( ay$fj ) && ax$linear && ay$linear ) ) {
        fj <- densify.polygon( cbind( fj.to.scale( ax, raw[ , 1L ] ), fj.to.scale( ay, raw[ , 2L ] ) ),
                               points.per.edge )
        raw <- cbind( fj.from.scale( ax, fj[ , 1L ] ), fj.from.scale( ay, fj[ , 2L ] ) )
    }
    unname( cbind( ax$to.display( raw[ , 1L ] ), ay$to.display( raw[ , 2L ] ) ) )
}


#' A FlowJo ellipse as a polygon
#' @title FlowJo Ellipse
#' @description
#' FlowJo's four edge points are the ends of the ellipse's two axes on its
#' display scale; the polygon has the ellipse's area there.
#' @param gate FlowJo gate.
#' @param ax,ay FlowJo axes.
#' @param n.points Vertices.
#' @return Raw vertices (n x 2).
#' @keywords internal
fj.ellipse.raw <- function( gate, ax, ay, n.points = fj.ellipse.vertices ) {
    e <- gate$edge[ 1:4, , drop = FALSE ]
    fj <- cbind( fj.to.scale( ax, e[ , 1L ] ), fj.to.scale( ay, e[ , 2L ] ) )
    centre <- colMeans( fj )
    u <- ( fj[ 1L, ] - fj[ 2L, ] ) / 2
    v <- ( fj[ 3L, ] - fj[ 4L, ] ) / 2
    t <- seq( 0, 2 * pi, length.out = n.points + 1L )[ seq_len( n.points ) ]
    k <- sqrt( 2 * pi / ( n.points * sin( 2 * pi / n.points ) ) )
    pts <- matrix( centre, n.points, 2L, byrow = TRUE ) +
        k * ( outer( cos( t ), u ) + outer( sin( t ), v ) )
    cbind( fj.from.scale( ax, pts[ , 1L ] ), fj.from.scale( ay, pts[ , 2L ] ) )
}


#' Whether two channels are a scatter singlet pair
#' @title Is Singlet Pair
#' @param chans Channel names.
#' @param exp Experiment channels.
#' @return \code{TRUE} for scatter area against the height or width of the
#'   same detector.
#' @keywords internal
fj.is.singlet.pair <- function( chans, exp ) {
    if ( length( chans ) != 2L )
        return( FALSE )
    stem.of <- function( x ) sub( "-[^-]*$", "", x )
    pair <- exp$singlet.pair
    if ( length( pair ) == 2L ) {
        area <- pair[ 1L ]
        alternatives <- setdiff( c( pair[ 2L ], paste0( stem.of( area ), c( "-H", "-W" ) ) ), area )
        if ( area %in% chans && any( chans %in% alternatives ) )
            return( TRUE )
    }
    if ( !all( chans %in% exp$scatter ) || !all( grepl( "-", chans ) ) )
        return( FALSE )
    if ( stem.of( chans[ 1L ] ) != stem.of( chans[ 2L ] ) )
        return( FALSE )
    suffixes <- sort( sub( ".*-", "", chans ) )
    identical( suffixes, c( "A", "H" ) ) || identical( suffixes, c( "A", "W" ) )
}


#' One FlowJo population on the display scale
#' @title FlowJo Piece
#' @param pop FlowJo population.
#' @param source Gating source.
#' @param mapping FlowJo parameter -> channel.
#' @param exp Experiment channels.
#' @return List with pop, channels, shape (template entry without
#'   population fields), extent (per channel c(lo, hi), NA open-ended),
#'   signs (per channel "+", "-", or NA for a scatter channel of a 2-D
#'   gate), notes and singlets.
#' @keywords internal
fj.piece <- function( pop, source, mapping, exp ) {
    g <- pop$gate
    chans <- unname( mapping[ g$channels ] )
    axes <- Map( fj.axis, g$channels, chans, MoreArgs = list( source = source, exp = exp ) )
    names( axes ) <- NULL
    notes <- character( 0 )
    fill <- function( lohi ) c( if ( is.na( lohi[ 1L ] ) ) -1 else lohi[ 1L ],
                                if ( is.na( lohi[ 2L ] ) ) 2 else lohi[ 2L ] )
    poly <- NULL
    if ( g$kind == "polygon" ) {
        poly <- fj.polygon.display( g$vertices, axes[[ 1L ]], axes[[ 2L ]] )
    } else if ( g$kind == "ellipse" ) {
        raw <- fj.ellipse.raw( g, axes[[ 1L ]], axes[[ 2L ]] )
        poly <- cbind( axes[[ 1L ]]$to.display( raw[ , 1L ] ), axes[[ 2L ]]$to.display( raw[ , 2L ] ) )
        poly <- unname( rbind( poly, poly[ 1L, ] ) )
        notes <- sprintf( "Ellipse traced as a %d-vertex polygon of the same area.", fj.ellipse.vertices )
    }
    if ( !is.null( poly ) ) {
        extent <- list( range( poly[ , 1L ] ), range( poly[ , 2L ] ) )
        shape <- list( shape.channels = chans, boundary = poly )
    } else {
        extent <- Map( function( b, ax ) c( if ( is.na( b[ 1L ] ) ) NA_real_ else ax$to.display( b[ 1L ] ),
                                            if ( is.na( b[ 2L ] ) ) NA_real_ else ax$to.display( b[ 2L ] ) ),
                       g$bounds, axes )
        if ( length( chans ) == 1L ) {
            lh <- fill( extent[[ 1L ]] )
            shape <- list( shape.channels = chans, range = lh, boundary = rect.polygon( lh[ 1L ], lh[ 2L ] ) )
        } else {
            x <- fill( extent[[ 1L ]] )
            y <- fill( extent[[ 2L ]] )
            shape <- list( shape.channels = chans, boundary = rect.polygon( x[ 1L ], x[ 2L ], y[ 1L ], y[ 2L ] ) )
        }
    }
    signs <- vapply( seq_along( chans ), function( k ) {
        if ( chans[ k ] %in% exp$scatter && length( chans ) == 2L ) return( NA_character_ )
        lo <- extent[[ k ]][ 1L ]
        if ( !is.na( lo ) && lo > axes[[ k ]]$zero ) "+" else "-"
    }, character( 1 ) )
    list( pop = pop, channels = chans, shape = shape, extent = extent, signs = signs,
          notes = notes, singlets = fj.is.singlet.pair( chans, exp ) )
}


#' Region of one side of a threshold
#' @title Region from Sign
#' @param sign "+" or "-".
#' @return "above" or "below".
#' @keywords internal
fj.region.1d <- function( sign ) if ( identical( sign, "+" ) ) "above" else "below"


#' Trainable form of sibling pieces on one plot
#' @title FlowJo Gate Plan
#' @param grp List of pieces.
#' @return List with gate.type, sep (separation channel index of a 1dsep
#'   gate), regions and algorithm.
#' @keywords internal
fj.plan.gate <- function( grp ) {
    first <- grp[[ 1L ]]
    n <- length( grp )
    if ( first$singlets )
        return( list( gate.type = "singlets", sep = NA_integer_, regions = rep( "inside", n ),
                      algorithm = "singlets" ) )
    if ( length( first$channels ) == 1L )
        return( list( gate.type = "1dsep", sep = 1L, algorithm = "tail",
                      regions = vapply( grp, function( p ) fj.region.1d( p$signs[ 1L ] ), "" ) ) )
    fluor <- which( !is.na( first$signs ) )
    if ( length( fluor ) == 0L )
        return( list( gate.type = "free", sep = NA_integer_, regions = rep( "inside", n ),
                      algorithm = "cells" ) )
    if ( length( fluor ) == 1L )
        return( list( gate.type = "1dsep", sep = fluor, algorithm = "mixture",
                      regions = vapply( grp, function( p ) fj.region.1d( p$signs[ fluor ] ), "" ) ) )
    layout <- if ( n > 1L ) fj.quadrant.layout( grp ) else NULL
    if ( !is.null( layout ) )
        return( c( layout, list( algorithm = "mixture" ) ) )
    list( gate.type = "2dsep", sep = NA_integer_, algorithm = "mixture",
          regions = vapply( grp, function( p ) paste0( "x", p$signs[ 1L ], "y", p$signs[ 2L ] ), "" ) )
}


#' Quadrant layout reproducing sibling shapes
#'
#' @title FlowJo Quadrant Layout
#'
#' @description
#' Port of Honeychrome's quadrant layout.  Over a grid on the shapes' joint
#' extent, each pair of thresholds splits the plot into quadrants; a
#' quadrant goes to the shape covering most of it when that shape covers at
#' least a quarter of it and most of what is drawn there.  A shape's region
#' is the union of its quadrants.  The layout with the fewest mismatched
#' cells wins (drawn but given to another shape or none, plus a quarter of
#' those given but not drawn).  The gate separates on one parameter only
#' when one threshold alone does nearly as well and the other threshold
#' separates no two shapes.
#'
#' @param grp List of pieces drawn on two fluorescence channels.
#'
#' @return List with gate.type, sep and regions, or \code{NULL} when no
#'   layout gives every shape a region of its own.
#'
#' @keywords internal
fj.quadrant.layout <- function( grp ) {
    polys <- lapply( grp, function( p ) p$shape$boundary )
    if ( any( vapply( polys, function( pl ) !is.matrix( pl ) || nrow( pl ) < 3L, logical( 1 ) ) ) )
        return( NULL )
    pts <- do.call( rbind, polys )
    lo <- apply( pts, 2L, min )
    hi <- apply( pts, 2L, max )
    if ( !all( hi > lo ) )
        return( NULL )
    n <- fj.layout$grid
    cx <- lo[ 1L ] + ( seq_len( n ) - 0.5 ) * ( hi[ 1L ] - lo[ 1L ] ) / n
    cy <- lo[ 2L ] + ( seq_len( n ) - 0.5 ) * ( hi[ 2L ] - lo[ 2L ] ) / n
    gx <- rep( cx, times = n )
    gy <- rep( cy, each = n )
    cover <- lapply( polys, function( pl ) matrix( inside.polygon( gx, gy, pl ), n, n, byrow = TRUE ) )
    drawn <- Reduce( `|`, cover )
    if ( !any( drawn ) || !all( vapply( cover, any, logical( 1 ) ) ) )
        return( NULL )
    total <- sum( drawn )
    cum <- function( a ) {
        m <- t( apply( apply( a * 1, 2L, cumsum ), 1L, cumsum ) )
        rbind( 0, cbind( 0, m ) )
    }
    sums <- lapply( cover, cum )
    sum.drawn <- cum( drawn )
    box <- function( C, r0, r1, c0, c1 )
        C[ cbind( r1 + 1, c1 + 1 ) ] - C[ cbind( r0 + 1, c1 + 1 ) ] -
            C[ cbind( r1 + 1, c0 + 1 ) ] + C[ cbind( r0 + 1, c0 + 1 ) ]
    w <- fj.layout$undrawn.weight
    mc <- fj.layout$min.cover

    assign.cells <- function( cells ) {
        cost <- 0
        owner <- integer( 0 )
        for ( lab in names( cells ) ) {
            b <- cells[[ lab ]]
            area <- ( b[ 2 ] - b[ 1 ] ) * ( b[ 4 ] - b[ 3 ] )
            if ( area <= 0 ) next
            got <- vapply( sums, box, numeric( 1 ), b[ 1 ], b[ 2 ], b[ 3 ], b[ 4 ] )
            k <- which.max( got )
            in.drawn <- box( sum.drawn, b[ 1 ], b[ 2 ], b[ 3 ], b[ 4 ] )
            if ( got[ k ] >= mc * area && got[ k ] >= 0.5 * in.drawn ) {
                owner[ lab ] <- k
                cost <- cost + ( in.drawn - got[ k ] ) + w * ( area - got[ k ] )
            } else cost <- cost + in.drawn
        }
        list( cost = cost, owner = owner, cells = cells )
    }
    quadrants <- function( i, j ) list( `x-y-` = c( 0, j, 0, i ), `x+y-` = c( 0, j, i, n ),
                                        `x-y+` = c( j, n, 0, i ), `x+y+` = c( j, n, i, n ) )

    ii <- rep( seq_len( n - 1L ), each = n - 1L )
    jj <- rep( seq_len( n - 1L ), times = n - 1L )
    zero <- rep( 0, length( ii ) )
    full <- rep( n, length( ii ) )
    bounds <- list( list( zero, jj, zero, ii ), list( zero, jj, ii, full ),
                    list( jj, full, zero, ii ), list( jj, full, ii, full ) )
    total.cost <- numeric( length( ii ) )
    claimed <- matrix( FALSE, length( grp ), length( ii ) )
    for ( b in bounds ) {
        area <- ( b[[ 2 ]] - b[[ 1 ]] ) * ( b[[ 4 ]] - b[[ 3 ]] )
        got <- do.call( rbind, lapply( sums, box, b[[ 1 ]], b[[ 2 ]], b[[ 3 ]], b[[ 4 ]] ) )
        k <- max.col( t( got ), ties.method = "first" )
        best <- got[ cbind( k, seq_along( k ) ) ]
        in.drawn <- box( sum.drawn, b[[ 1 ]], b[[ 2 ]], b[[ 3 ]], b[[ 4 ]] )
        owned <- area > 0 & best >= mc * area & best >= 0.5 * in.drawn
        total.cost <- total.cost + ifelse( owned, ( in.drawn - best ) + w * ( area - best ), in.drawn )
        for ( p in seq_along( grp ) )
            claimed[ p, ] <- claimed[ p, ] | ( owned & k == p )
    }
    valid <- apply( claimed, 2L, all )
    layouts <- list()
    if ( any( valid ) ) {
        floor <- min( total.cost[ valid ] )
        for ( s in which( valid & total.cost <= floor + fj.layout$tie * total ) )
            layouts[[ length( layouts ) + 1L ]] <- assign.cells( quadrants( ii[ s ], jj[ s ] ) )
    }
    best2 <- if ( length( layouts ) )
        layouts[[ which.min( vapply( layouts, `[[`, numeric( 1 ), "cost" ) ) ]] else NULL
    best1 <- list()
    for ( k.axis in 0:1 ) {
        for ( i in seq_len( n - 1L ) ) {
            cells <- if ( k.axis == 0L ) list( `-` = c( 0, n, 0, i ), `+` = c( 0, n, i, n ) )
                     else list( `-` = c( 0, i, 0, n ), `+` = c( i, n, 0, n ) )
            a <- assign.cells( cells )
            key <- as.character( k.axis )
            if ( length( unique( a$owner ) ) == length( grp ) &&
                 ( is.null( best1[[ key ]] ) || a$cost < best1[[ key ]]$cost ) )
                best1[[ key ]] <- c( a, list( axis = k.axis ) )
        }
    }
    flip <- function( q, pos ) {
        substr( q, pos, pos ) <- if ( substr( q, pos, pos ) == "-" ) "+" else "-"
        q
    }
    splits <- function( layout, pos ) {
        sides <- tapply( substr( names( layout$owner ), pos, pos ), layout$owner,
                         function( s ) unique( s ), simplify = FALSE )
        single <- unlist( sides[ lengths( sides ) == 1L ] )
        setequal( unique( single ), c( "+", "-" ) )
    }
    redundant <- function( layout, pos ) {
        if ( splits( layout, pos ) ) return( FALSE )
        for ( q in names( layout$owner ) ) {
            k <- layout$owner[[ q ]]
            m <- flip( q, pos )
            if ( m %in% names( layout$owner ) && layout$owner[[ m ]] != k ) return( FALSE )
            if ( !( m %in% names( layout$owner ) ) ) {
                b <- layout$cells[[ m ]]
                area <- ( b[ 2 ] - b[ 1 ] ) * ( b[ 4 ] - b[ 3 ] )
                if ( area && box( sum.drawn, b[ 1 ], b[ 2 ], b[ 3 ], b[ 4 ] ) >=
                     fj.layout$empty.cover * area ) return( FALSE )
            }
        }
        TRUE
    }
    one.d <- unname( best1 )
    if ( !is.null( best2 ) ) {
        allowed <- list()
        close <- list()
        for ( kp in list( c( 0L, 4L ), c( 1L, 2L ) ) ) {
            b1 <- best1[[ as.character( kp[ 1 ] ) ]]
            if ( is.null( b1 ) ) next
            if ( b1$cost <= best2$cost + fj.layout$max.extra.1d * total &&
                 any( vapply( layouts, redundant, logical( 1 ), kp[ 2 ] ) ) )
                allowed[[ length( allowed ) + 1L ]] <- b1
        }
        for ( kp in list( c( 0L, 4L ), c( 1L, 2L ) ) ) {
            b1 <- best1[[ as.character( kp[ 1 ] ) ]]
            if ( is.null( b1 ) ) next
            if ( !splits( best2, kp[ 2 ] ) && b1$cost <= best2$cost + fj.layout$tol.1d * total )
                close[[ length( close ) + 1L ]] <- b1
        }
        one.d <- c( allowed, close )
    }
    if ( length( one.d ) ) {
        b1 <- one.d[[ which.min( vapply( one.d, `[[`, numeric( 1 ), "cost" ) ) ]]
        regions <- rep( NA_character_, length( grp ) )
        for ( side in names( b1$owner ) )
            regions[ b1$owner[[ side ]] ] <- fj.region.1d( side )
        return( list( gate.type = "1dsep", sep = b1$axis + 1L, regions = regions ) )
    }
    if ( is.null( best2 ) )
        return( NULL )
    regions <- vapply( seq_along( grp ), function( p )
        join.quadrants( names( best2$owner )[ best2$owner == p ] ), character( 1 ) )
    list( gate.type = "2dsep", sep = NA_integer_, regions = regions )
}


#' Group sibling pieces into gates
#' @title FlowJo Sibling Groups
#' @param pieces List of pieces with one parent.
#' @param merge Merge siblings drawn on the same channels.
#' @return List of groups (lists of pieces).
#' @keywords internal
fj.group.siblings <- function( pieces, merge ) {
    groups <- list()
    for ( p in pieces ) {
        placed <- FALSE
        if ( merge && !p$singlets && any( !is.na( p$signs ) ) ) {
            for ( gi in seq_along( groups ) ) {
                q <- groups[[ gi ]][[ 1L ]]
                if ( !identical( q$channels, p$channels ) || q$singlets ) next
                regions <- fj.plan.gate( c( groups[[ gi ]], list( p ) ) )$regions
                if ( !anyDuplicated( regions ) ) {
                    groups[[ gi ]] <- c( groups[[ gi ]], list( p ) )
                    placed <- TRUE
                    break
                }
            }
        }
        if ( !placed )
            groups[[ length( groups ) + 1L ]] <- list( p )
    }
    groups
}


#' Algorithm proposed for a histogram gate
#' @title Histogram Algorithm
#' @param positive.fraction Fraction of the parent on the positive side, or
#'   \code{NULL}.
#' @return "mixture" above 25\%, else "tail".
#' @keywords internal
fj.histogram.algorithm <- function( positive.fraction ) {
    if ( length( positive.fraction ) == 1L && !is.na( positive.fraction ) &&
         positive.fraction > fj.histogram.mixture.frequency ) "mixture" else "tail"
}


#' Name usable for a gate or population
#' @title Safe Name
#' @param name Character.
#' @return Name with "/" replaced (it separates gate and population).
#' @keywords internal
fj.safe.name <- function( name ) {
    out <- trimws( gsub( "/", "\u2215", as.character( name ), fixed = TRUE ) )
    if ( nzchar( out ) ) out else "gate"
}


#' Gate definition of a group of pieces
#' @title FlowJo Gate Definition
#' @param grp Pieces.
#' @param plan Result of \code{fj.plan.gate}.
#' @param name Gate name.
#' @param parent.gate,parent.popul Parent.
#' @return Gate definition.
#' @keywords internal
fj.gate.def <- function( grp, plan, name, parent.gate, parent.popul ) {
    chans <- grp[[ 1L ]]$channels
    type <- plan$gate.type
    sep.axis <- "x"
    if ( type == "1dsep" ) {
        k <- plan$sep
        mx <- chans[ k ]
        my <- if ( length( chans ) == 2L ) chans[ 3L - k ] else NULL
        if ( k == 2L ) sep.axis <- "y"
    } else {
        mx <- chans[ 1L ]
        my <- if ( length( chans ) > 1L ) chans[ 2L ] else NULL
    }
    defaults <- default.populations.for.type( type )
    pos.by.region <- stats::setNames( vapply( defaults, `[[`, numeric( 1 ), "label.pos" ),
                                      vapply( defaults, `[[`, character( 1 ), "region" ) )
    populations <- list()
    template <- list()
    for ( i in seq_along( grp ) ) {
        p <- grp[[ i ]]
        region <- plan$regions[ i ]
        pname <- fj.safe.name( p$pop$name )
        pos <- if ( region %in% names( pos.by.region ) ) as.integer( pos.by.region[[ region ]] ) else 1L
        populations[[ pname ]] <- list( label = pname, label.pos = pos, region = region )
        template[[ pname ]] <- c( p$shape, list( threshold.x = NA_real_, threshold.y = NA_real_,
                                                 region = region, label = pname ) )
    }
    g <- list( gate.name = name, gate.type = type, gate.marker.x = mx )
    g[ "gate.marker.y" ] <- list( my )
    g$sep.axis <- sep.axis
    g[ "parent.gate" ] <- list( parent.gate )
    g$parent.popul <- parent.popul
    g$algorithm <- "imported"
    g$proposed.algorithm <- plan$algorithm
    g$gate.param <- json.object( NULL )
    g$stats.parent <- json.object( NULL )
    g$populations <- populations
    g$template <- template
    g
}


#' Unique gate names
#' @title FlowJo Final Names
#' @param bases Named character vector temporary name -> name.
#' @param parent.paths Named list temporary name -> population path.
#' @param taken Names already used.
#' @return Named character vector temporary name -> unique name.
#' @keywords internal
fj.final.names <- function( bases, parent.paths, taken = character( 0 ) ) {
    out <- character( 0 )
    used <- taken
    groups <- split( names( bases ), factor( bases, levels = unique( bases ) ) )
    for ( base in names( groups ) ) {
        m <- groups[[ base ]]
        if ( length( m ) == 1L && !( base %in% used ) ) {
            out[ m ] <- base
            used <- c( used, base )
        }
    }
    for ( base in names( groups ) ) {
        m <- groups[[ base ]]
        if ( m[ 1L ] %in% names( out ) ) next
        paths <- lapply( m, function( k ) { p <- parent.paths[[ k ]]; if ( length( p ) ) p else "All events" } )
        depth <- max( lengths( paths ) )
        chosen <- NULL
        for ( k in seq_len( depth ) ) {
            cands <- vapply( paths, function( p ) sprintf( "%s (%s)", base,
                paste( vapply( utils::tail( p, k ), fj.safe.name, "" ), collapse = fj.name.sep ) ), "" )
            if ( !anyDuplicated( cands ) && !any( cands %in% used ) ) {
                chosen <- cands
                break
            }
        }
        if ( is.null( chosen ) ) {
            chosen <- character( 0 )
            for ( p in paths ) {
                cand <- sprintf( "%s (%s)", base, paste( vapply( p, fj.safe.name, "" ), collapse = fj.name.sep ) )
                stem <- cand
                k <- 2L
                while ( cand %in% c( used, chosen ) ) {
                    cand <- sprintf( "%s %d", stem, k )
                    k <- k + 1L
                }
                chosen <- c( chosen, cand )
            }
        }
        out[ m ] <- chosen
        used <- c( used, chosen )
    }
    out
}


#' Move the first singlet gate to the top of the hierarchy
#' @title FlowJo Singlets First
#' @description
#' The gates on the path from the top to the singlet gate are placed
#' directly below it, in their original order, and the singlet gate's own
#' daughters move to the last of them; the recorded parent counts of the
#' moved gates are cleared.
#' @param gate.defs Gate definitions.
#' @param forward.area Forward-scatter area channel, or \code{NULL}.
#' @return List with gate.defs, moved (names, singlet gate first) and below
#'   (gate and population the singlet gate's daughters now sit in).
#' @keywords internal
fj.singlets.first <- function( gate.defs, forward.area = NULL ) {
    names.all <- vapply( gate.defs, `[[`, character( 1 ), "gate.name" )
    is.singlets <- vapply( gate.defs, function( g ) identical( g$gate.type, "singlets" ), logical( 1 ) )
    pick <- if ( !is.null( forward.area ) )
        which( is.singlets & vapply( gate.defs, function( g )
            forward.area %in% c( g$gate.marker.x, g$gate.marker.y ), logical( 1 ) ) )[ 1L ]
    else which( is.singlets )[ 1L ]
    none <- list( gate.defs = gate.defs, moved = character( 0 ), below = NULL )
    if ( is.na( pick ) || !has.value( gate.defs[[ pick ]]$parent.gate ) )
        return( none )
    single <- gate.defs[[ pick ]]
    chain <- character( 0 )
    name <- single$parent.gate
    while ( has.value( name ) && name %in% names.all && !( name %in% chain ) ) {
        chain <- c( chain, name )
        name <- gate.defs[[ match( name, names.all ) ]]$parent.gate
    }
    chain <- rev( chain )
    last <- c( single$parent.gate, single$parent.popul )
    single.pop <- names( single$populations )[ 1L ] %||% "root"
    for ( i in seq_along( gate.defs ) )
        if ( identical( gate.defs[[ i ]]$parent.gate, single$gate.name ) ) {
            gate.defs[[ i ]]$parent.gate <- last[ 1L ]
            gate.defs[[ i ]]$parent.popul <- last[ 2L ]
        }
    gate.defs[[ pick ]][ "parent.gate" ] <- list( NULL )
    gate.defs[[ pick ]]$parent.popul <- "root"
    top <- match( chain[ 1L ], names.all )
    gate.defs[[ top ]]$parent.gate <- single$gate.name
    gate.defs[[ top ]]$parent.popul <- single.pop
    for ( nm in c( single$gate.name, chain ) ) {
        i <- match( nm, names.all )
        if ( !is.null( gate.defs[[ i ]]$source ) )
            gate.defs[[ i ]]$source[ "parent_count" ] <- list( NULL )
    }
    gate.defs <- ordered.gate.defs( gate.defs )
    for ( k in seq_along( gate.defs ) )
        gate.defs[[ k ]]$gate.number <- k
    list( gate.defs = gate.defs, moved = c( single$gate.name, chain ), below = last )
}


#' Convert a FlowJo gating source to gate definitions
#'
#' @title FlowJo Gate Definitions
#'
#' @description
#' Converts the importable populations of a gating source into gate
#' definitions (rules in the file header of \code{flowjo_import.R}).
#'
#' @param source Gating source of a parsed workspace.
#' @param mapping Named character vector FlowJo parameter -> channel
#'   (\code{NA}: none); populations using an unmapped parameter are dropped.
#' @param exp Experiment channels (\code{fj.experiment}).
#' @param merge.siblings Merge sibling gates drawn on the same channels.
#' @param existing.names Gate names already in use.
#' @param workspace.name Recorded in each gate's \code{source}.
#'
#' @return List with gate.definitions, plan (\code{fj.plan.import}), notes
#'   (gate -> character) and populations (FlowJo population key -> "gate/pop").
#'
#' @keywords internal
fj.build.gate.definitions <- function( source, mapping, exp, merge.siblings = TRUE,
                                       existing.names = character( 0 ), workspace.name = "" ) {
    plan <- fj.plan.import( source, mapping )
    keep <- plan$keep
    gate.defs <- list()
    notes <- list()
    pop.map <- list()
    bases <- character( 0 )
    parent.paths <- list()
    temp.names <- character( 0 )
    title <- sprintf( "%s: %s (%d populations)", if ( source$kind == "group" ) "Group" else "Sample",
                      source$name, length( fj.populations( source ) ) )

    unique.name <- function( name, parent.path ) {
        base <- fj.safe.name( name )
        cand <- base
        k <- 2L
        while ( cand %in% temp.names ) {
            cand <- paste0( base, "\u001f", k )
            k <- k + 1L
        }
        temp.names <<- c( temp.names, cand )
        bases[ cand ] <<- base
        parent.paths[[ cand ]] <<- parent.path
        cand
    }
    label.of <- function( ch ) ch

    emit <- function( parent, children ) {
        kids <- Filter( function( k ) k$key %in% keep, children )
        if ( length( kids ) == 0L )
            return( invisible( NULL ) )
        pg <- if ( is.null( parent ) ) list( NULL, "root" ) else as.list( pop.map[[ parent$key ]] )
        pieces <- lapply( kids, fj.piece, source = source, mapping = mapping, exp = exp )
        for ( grp in fj.group.siblings( pieces, merge.siblings ) ) {
            nm <- vapply( grp, function( p ) p$pop$name, "" )
            gname <- unique.name( if ( length( grp ) == 1L ) nm else paste( nm, collapse = " & " ),
                                  if ( is.null( parent ) ) character( 0 ) else parent$path )
            plan.g <- fj.plan.gate( grp )
            g <- fj.gate.def( grp, plan.g, gname, pg[[ 1L ]], pg[[ 2L ]] )
            g$gate.number <- length( gate.defs ) + 1L
            g$template.original <- g$template
            parent.count <- if ( is.null( parent ) ) source$count else parent$count
            if ( g$gate.type == "1dsep" && length( grp[[ 1L ]]$channels ) == 1L ) {
                g$proposal.rule <- "histogram"
                g$proposed.algorithm <- fj.histogram.algorithm(
                    fj.histogram.positive.fraction( grp, plan.g, parent.count ) )
            } else if ( g$gate.type %in% ag.threshold.types ) {
                g$proposal.rule <- "fraction"
            } else if ( g$gate.type == "free" ) {
                g$gate.param <- list( debris.axis = if ( length( exp$singlet.pair ) &&
                    identical( grp[[ 1L ]]$channels[ 2L ], exp$singlet.pair[ 1L ] ) ) "y" else "x" )
            }
            pnames <- names( g$populations )
            g$source <- list(
                origin = "FlowJo", workspace = workspace.name, gating_source = title,
                sample = if ( source$kind == "sample" ) source$name else source$sample.name,
                parent_count = parent.count,
                populations = stats::setNames( lapply( grp, function( p ) p$pop$key ), pnames ),
                counts = stats::setNames( lapply( grp, function( p ) p$pop$count ), pnames ),
                parameters = as.list( stats::setNames(
                    unlist( lapply( grp, function( p ) p$pop$gate$channels ) ),
                    unlist( lapply( grp, `[[`, "channels" ) ) ) ) )
            gate.defs[[ length( gate.defs ) + 1L ]] <<- g
            notes[[ gname ]] <<- c( fj.describe( g ), unlist( lapply( grp, `[[`, "notes" ) ) )
            for ( i in seq_along( grp ) )
                pop.map[[ grp[[ i ]]$pop$key ]] <<- c( gname, pnames[ i ] )
        }
        for ( k in kids ) emit( k, k$children )
    }
    emit( NULL, source$roots )

    renamed <- fj.final.names( bases, parent.paths, existing.names )
    for ( i in seq_along( gate.defs ) ) {
        g <- gate.defs[[ i ]]
        g$name.group <- bases[[ g$gate.name ]]
        g$gate.name <- renamed[[ g$gate.name ]]
        if ( has.value( g$parent.gate ) ) g$parent.gate <- renamed[[ g$parent.gate ]]
        gate.defs[[ i ]] <- g
    }
    names( notes ) <- unname( renamed[ names( notes ) ] )
    pop.map <- lapply( pop.map, function( v ) c( renamed[[ v[ 1L ] ]], v[ 2L ] ) )

    moved <- fj.singlets.first( gate.defs, exp$singlet.pair[ 1L ] )
    gate.defs <- moved$gate.defs
    if ( length( moved$moved ) ) {
        singlets <- moved$moved[ 1L ]
        single.key <- names( pop.map )[ vapply( pop.map, function( v ) v[ 1L ] == singlets, logical( 1 ) ) ][ 1L ]
        for ( nm in moved$moved[ -1L ] )
            for ( key in names( pop.map ) )
                if ( identical( pop.map[[ key ]][ 1L ], nm ) ) pop.map[[ key ]] <- NULL
        pop.map[[ single.key ]] <- moved$below
        notes[[ singlets ]] <- c( notes[[ singlets ]], sprintf(
            "Moved to the top of the hierarchy; %s now sit%s directly below it.",
            paste0( "'", moved$moved[ -1L ], "'", collapse = ", " ),
            if ( length( moved$moved ) == 2L ) "s" else "" ) )
    }
    list( gate.definitions = gate.defs, plan = plan, notes = notes,
          populations = vapply( pop.map, function( v ) population.key( v[ 1L ], v[ 2L ] ), "" ) )
}


#' Positive fraction of a histogram gate from recorded counts
#' @title Histogram Positive Fraction
#' @param grp Pieces.
#' @param plan Gate plan.
#' @param parent.count Parent count recorded in the workspace.
#' @return Fraction, or \code{NULL}.
#' @keywords internal
fj.histogram.positive.fraction <- function( grp, plan, parent.count ) {
    if ( length( parent.count ) != 1L || is.na( parent.count ) || parent.count == 0 )
        return( NULL )
    counts <- vapply( grp, function( p ) as.numeric( p$pop$count ), numeric( 1 ) )
    above <- plan$regions == "above"
    if ( any( above ) && !anyNA( counts[ above ] ) )
        return( sum( counts[ above ] ) / parent.count )
    if ( any( !above ) && !anyNA( counts[ !above ] ) )
        return( 1 - sum( counts[ !above ] ) / parent.count )
    NULL
}


#' One line describing a converted gate
#' @title Describe FlowJo Gate
#' @param g Gate definition.
#' @return Character.
#' @keywords internal
fj.describe <- function( g ) {
    pops <- paste( names( g$populations ), vapply( g$populations, `[[`, "", "region" ), collapse = ", " )
    tail <- sprintf( "; proposed algorithm '%s'", g$proposed.algorithm )
    switch( g$gate.type,
        `1dsep` = sprintf( "1dsep separating on %s (%s axis), view %s: %s%s.", g$gate.marker.x,
                           g$sep.axis, g$gate.marker.y %||% "histogram", pops, tail ),
        `2dsep` = sprintf( "2dsep on %s x %s: %s%s.", g$gate.marker.x, g$gate.marker.y, pops, tail ),
        singlets = sprintf( "Singlet gate%s.", tail ),
        sprintf( "Scatter gate on %s x %s: cell gate cutting debris on %s%s.", g$gate.marker.x,
                 g$gate.marker.y, if ( identical( g$gate.param$debris.axis, "y" ) ) g$gate.marker.y
                 else g$gate.marker.x, tail ) )
}
