# gate_calculation.R
#
# Calculates one gate's boundary from its parent population's events, as
# Honeychrome Automated Gating does.
#
# Boundary entries (one per population) are named lists:
#   boundary     closed polygon (n x 2 matrix) for drawing and for polygon
#                gates
#   range        c(lo, hi) of a one-dimensional shape
#   threshold.x, threshold.y   thresholds of threshold gates (NA otherwise)
#   region       side of the thresholds, or "inside"
#   label        population label
#   n.events, fraction   added when the gate is applied
#
# Gate types:
#   singlets  line fit of height on area with three rounds of outlier
#             removal (residual SD estimated from the median absolute
#             deviation, so a shoulder of doublets does not widen the band);
#             the gate is a band of +/- width.sd residual SDs along the whole
#             axis; saturated piles are left out of the fit
#   1dsep     one threshold on gate.marker.x
#   2dsep     one threshold per axis; populations are quadrants or unions
#   free      "cells": a rectangle from the debris cut on forward scatter to
#             where the cells' density falls to edge.density of its peak;
#             "all": convex hull of the events
#   imported  shapes as drawn (e.g. from FlowJo), applied as they are


#' New calculation context
#'
#' @title Gate Calculation Context
#'
#' @param axis.limits Named list channel -> \code{c(lo, hi)} display limits,
#'   used to draw threshold rectangles (missing channels use the data
#'   range).
#' @param transforms Named list channel -> transform parameters (see
#'   \code{\link{make.transform}}), needed for threshold factors.
#'
#' @return Environment with \code{axis.limits}, \code{transforms} and
#'   \code{notes}.
#'
#' @keywords internal
gate.context <- function( axis.limits = NULL, transforms = NULL ) {
    ctx <- new.env( parent = emptyenv() )
    ctx$axis.limits <- axis.limits %||% list()
    ctx$transforms  <- transforms %||% list()
    ctx$notes       <- character( 0 )
    ctx
}


#' Threshold multiplied by the gate's threshold factor
#'
#' @title Scaled Threshold
#'
#' @param ctx Calculation context.
#' @param t Threshold.
#' @param gate.param Gate parameters.
#' @param suffix \code{""}, \code{".x"} or \code{".y"}.
#' @param channel Channel name.
#'
#' @return Numeric threshold; a note is added to \code{ctx} when the factor
#'   cannot be applied.
#'
#' @keywords internal
scaled.threshold <- function( ctx, t, gate.param, suffix, channel ) {
    res <- nudge.threshold( t, threshold.factor( gate.param, suffix ),
                            ctx$transforms[[ channel ]] )
    if ( nzchar( res$note ) )
        ctx$notes <- c( ctx$notes, paste0( channel, ": ", res$note ) )
    res$threshold
}


#' Display limits of a channel
#'
#' @title Axis Limits
#'
#' @param ctx Calculation context.
#' @param channel Channel name.
#' @param values Events on the channel.
#'
#' @return \code{c(lo, hi)}.
#'
#' @keywords internal
axis.limits.of <- function( ctx, channel, values ) {
    lim <- ctx$axis.limits[[ channel ]]
    if ( length( lim ) >= 2L ) return( as.numeric( lim[ 1:2 ] ) )
    if ( length( values ) ) return( range( values ) )
    c( 0, 1 )
}


#' Calculate a gate's boundary
#'
#' @title Calculate Gate
#'
#' @param data Matrix of the parent population's events (transformed),
#'   columns named by channel.
#' @param gate.def Gate definition (effective definition for replicate
#'   gates, see \code{\link{effective.gate.def}}).
#' @param ctx Calculation context (\code{\link{gate.context}}).
#' @param trained Boundaries calculated earlier in the same run; replicate
#'   gates copy their origin's boundary from here.
#'
#' @return Named list of boundary entries, empty when the gate cannot be
#'   calculated (missing channels, too few events, no origin boundary).
#'
#' @export
calculate.gate <- function( data, gate.def, ctx = gate.context(), trained = NULL ) {
    ctx$notes <- character( 0 )
    gate.type <- gate.type.of( gate.def )
    origin <- gate.def$replicate.of
    if ( !has.value( origin ) && gate.type == "replicate" ) origin <- gate.def$origin.gate
    if ( has.value( origin ) )
        return( trained[[ origin ]] %||% list() )
    if ( !all( gate.channels( gate.def ) %in% colnames( data ) ) )
        return( list() )

    result <- if ( identical( gate.def$algorithm, "imported" ) )
        calc.imported( data, gate.def, ctx )
    else switch( gate.type,
        singlets = calc.singlets( data, gate.def, ctx ),
        `1dsep`  = calc.1dsep( data, gate.def, ctx ),
        `2dsep`  = calc.2dsep( data, gate.def, ctx ),
        calc.free( data, gate.def, ctx ) )

    for ( p in names( result ) ) {
        lab <- gate.def$populations[[ p ]][[ "label" ]]
        result[[ p ]]$label <- if ( has.value( lab ) ) lab else p
    }
    result
}


#' Boundary of an imported gate
#'
#' @title Imported Gate Boundary
#'
#' @description
#' Shapes as drawn are returned as they are; imported threshold gates keep
#' their thresholds and have their rectangles redrawn to the axis limits.
#'
#' @inheritParams calculate.gate
#'
#' @return Named list of boundary entries.
#'
#' @keywords internal
calc.imported <- function( data, gate.def, ctx ) {
    template <- gate.def$template
    if ( length( template ) == 0L )
        return( list() )
    if ( is.shape.boundary( template ) ) {
        if ( !all( shape.channels( gate.def ) %in% colnames( data ) ) )
            return( list() )
        return( template )
    }
    if ( gate.def$gate.type %in% ag.threshold.types ) {
        tx <- template[[ 1L ]]$threshold.x
        ty <- template[[ 1L ]]$threshold.y
        if ( !has.value( tx ) || ( gate.def$gate.type == "2dsep" && !has.value( ty ) ) )
            return( list() )
        y.values <- if ( gate.def$gate.type == "2dsep" ) data[ , gate.def$gate.marker.y ] else NULL
        return( threshold.boundaries( gate.def, tx, if ( has.value( ty ) ) ty else NA_real_,
                                      data[ , gate.def$gate.marker.x ], y.values, ctx ) )
    }
    keys <- c( "boundary", "range", "threshold.x", "threshold.y", "region" )
    lapply( template, function( e ) e[ intersect( keys, names( e ) ) ] )
}


#' Singlet gate
#'
#' @title Singlet Gate Boundary
#'
#' @inheritParams calculate.gate
#'
#' @return Named list with one boundary entry.
#'
#' @keywords internal
calc.singlets <- function( data, gate.def, ctx ) {
    gp <- gate.def$gate.param
    q <- as.numeric( gp[[ "quantile.trim" ]] %||% ag.defaults$quantile.trim )
    width.sd <- as.numeric( gp[[ "width.sd" ]] %||% 3 )
    sat <- as.numeric( gp[[ "saturation.fraction" ]] %||% 5e-4 )
    mx <- gate.def$gate.marker.x
    xy <- data[ , c( mx, gate.def$gate.marker.y ), drop = FALSE ]
    unsat <- rep( TRUE, nrow( xy ) )
    for ( k in 1:2 ) {
        pile <- saturation.cut( xy[ , k ], sat )
        if ( !is.null( pile ) ) unsat <- unsat & xy[ , k ] < pile
    }
    xy <- trim.events( xy[ unsat, , drop = FALSE ], q )
    if ( nrow( xy ) < 50L )
        return( list() )
    x <- xy[ , 1L ]
    y <- xy[ , 2L ]
    line.fit <- function( keep ) {
        cf <- stats::.lm.fit( cbind( 1, x[ keep ] ), y[ keep ] )$coefficients
        c( slope = cf[ 2 ], intercept = cf[ 1 ] )
    }
    robust.sd <- function( r ) stats::mad( r, constant = 1 / stats::qnorm( 0.75 ) )
    keep <- rep( TRUE, length( x ) )
    for ( it in 1:3 ) {
        cf <- line.fit( keep )
        resid <- y[ keep ] - ( cf[ 1 ] * x[ keep ] + cf[ 2 ] )
        sigma <- robust.sd( resid )
        if ( sigma == 0 ) break
        keep[ keep ] <- abs( resid ) < width.sd * sigma
    }
    cf <- line.fit( keep )
    sigma <- robust.sd( y[ keep ] - ( cf[ 1 ] * x[ keep ] + cf[ 2 ] ) )
    if ( sigma == 0 ) sigma <- 1
    offset <- width.sd * sigma
    full.x <- data[ , mx ]
    lim <- axis.limits.of( ctx, mx, full.x )
    x.lo <- min( min( full.x ), lim[ 1 ] )
    x.hi <- max( max( full.x ), lim[ 2 ] )
    line <- function( v ) cf[ 1 ] * v + cf[ 2 ]
    poly <- unname( cbind(
        c( x.lo, x.hi, x.hi, x.lo, x.lo ),
        c( line( x.lo ) - offset, line( x.hi ) - offset, line( x.hi ) + offset,
           line( x.lo ) + offset, line( x.lo ) - offset ) ) )
    pop <- names( gate.def$populations )[ 1L ] %||% "in"
    stats::setNames( list( list( boundary = poly, threshold.x = NA_real_,
                                 threshold.y = NA_real_, region = "inside" ) ), pop )
}


#' One-threshold gate
#'
#' @title 1dsep Gate Boundary
#'
#' @inheritParams calculate.gate
#'
#' @return Named list of boundary entries.
#'
#' @keywords internal
calc.1dsep <- function( data, gate.def, ctx ) {
    gp <- gate.def$gate.param
    q  <- as.numeric( gp[[ "quantile.trim" ]] %||% ag.defaults$quantile.trim )
    ch <- gate.def$gate.marker.x
    full <- data[ , ch ]
    fixed <- gp[[ "fixed.threshold.x" ]]
    if ( has.value( fixed ) )
        return( threshold.boundaries( gate.def, as.numeric( fixed ), NA_real_, full, NULL, ctx ) )
    if ( identical( gate.def$algorithm, "replicate" ) )
        return( list() )
    col <- trim.events( full, q )
    if ( length( col ) < ag.defaults$min.events )
        return( list() )
    algo <- if ( has.value( gate.def$algorithm ) ) gate.def$algorithm else "tail"
    calc <- crop.to.region( col,
        as.numeric( gp[[ "axis.lo" ]] %||% default.axis.lo( algo ) ),
        as.numeric( gp[[ "axis.hi" ]] %||% ag.defaults$axis.hi ) )
    if ( length( calc ) < ag.defaults$min.events ) calc <- col
    t <- scaled.threshold( ctx, compute.threshold( calc, algo, gp ), gp, "", ch )
    threshold.boundaries( gate.def, t, NA_real_, full, NULL, ctx )
}


#' Two-threshold gate
#'
#' @title 2dsep Gate Boundary
#'
#' @inheritParams calculate.gate
#'
#' @return Named list of boundary entries.
#'
#' @keywords internal
calc.2dsep <- function( data, gate.def, ctx ) {
    gp <- gate.def$gate.param
    q  <- as.numeric( gp[[ "quantile.trim" ]] %||% ag.defaults$quantile.trim )
    ch <- c( gate.def$gate.marker.x, gate.def$gate.marker.y )
    xy <- trim.events( data[ , ch, drop = FALSE ], q )
    both.fixed <- has.value( gp[[ "fixed.threshold.x" ]] ) && has.value( gp[[ "fixed.threshold.y" ]] )
    if ( nrow( xy ) < 20L && !both.fixed )
        return( list() )
    default.algo <- if ( has.value( gate.def$algorithm ) ) gate.def$algorithm else "tail"

    axis.params <- function( axis ) {
        p <- gp
        p$tail.fraction <- gp[[ paste0( "tail.fraction.", axis ) ]] %||%
            gp[[ "tail.fraction" ]] %||% ag.defaults$tail.fraction
        p$mixture.components <- gp[[ paste0( "mixture.components.", axis ) ]] %||%
            gp[[ "mixture.components" ]] %||% 2L
        p$mixture.cut <- gp[[ paste0( "mixture.cut.", axis ) ]] %||%
            gp[[ "mixture.cut" ]] %||% "outer"
        p
    }

    thresholds <- c( NA_real_, NA_real_ )
    for ( k in 1:2 ) {
        axis <- c( "x", "y" )[ k ]
        fixed <- gp[[ paste0( "fixed.threshold.", axis ) ]]
        if ( has.value( fixed ) ) {
            thresholds[ k ] <- as.numeric( fixed )
            next
        }
        algo <- gp[[ paste0( "algorithm.", axis ) ]]
        if ( !has.value( algo ) ) algo <- default.algo
        if ( algo == "replicate" )
            return( list() )
        col <- xy[ , k ]
        calc <- crop.to.region( col,
            as.numeric( gp[[ paste0( "axis.lo.", axis ) ]] %||% default.axis.lo( algo ) ),
            as.numeric( gp[[ paste0( "axis.hi.", axis ) ]] %||% ag.defaults$axis.hi ) )
        if ( length( calc ) < ag.defaults$min.events ) calc <- col
        t <- compute.threshold( calc, algo, axis.params( axis ) )
        thresholds[ k ] <- scaled.threshold( ctx, t, gp, paste0( ".", axis ),
                                             axis.channel( gate.def, axis ) )
    }
    threshold.boundaries( gate.def, thresholds[ 1 ], thresholds[ 2 ],
                          data[ , ch[ 1 ] ], data[ , ch[ 2 ] ], ctx )
}


#' Boundary entries of a threshold gate
#'
#' @title Threshold Boundaries
#'
#' @param gate.def Gate definition.
#' @param tx,ty Thresholds (\code{ty} \code{NA} for 1dsep).
#' @param x.values,y.values Events on each axis (for the drawing limits).
#' @param ctx Calculation context.
#'
#' @return Named list of boundary entries; rectangles are drawn to the
#'   display limits, thresholds are open-ended.
#'
#' @keywords internal
threshold.boundaries <- function( gate.def, tx, ty, x.values = NULL, y.values = NULL,
                                  ctx = gate.context() ) {
    xl <- axis.limits.of( ctx, gate.def$gate.marker.x, x.values )
    x.lo <- min( xl[ 1 ], tx )
    x.hi <- max( xl[ 2 ], tx )
    regions <- population.regions( gate.def )
    out <- list()
    if ( identical( gate.def$gate.type, "1dsep" ) ) {
        for ( p in names( regions ) ) {
            box <- if ( regions[[ p ]] == "below" ) rect.polygon( x.lo, tx ) else rect.polygon( tx, x.hi )
            out[[ p ]] <- list( boundary = box, threshold.x = tx, threshold.y = NA_real_,
                                region = regions[[ p ]] )
        }
        return( out )
    }
    yl <- axis.limits.of( ctx, gate.def$gate.marker.y, y.values )
    y.lo <- min( yl[ 1 ], ty )
    y.hi <- max( yl[ 2 ], ty )
    for ( p in names( regions ) )
        out[[ p ]] <- list(
            boundary = region.polygon( regions[[ p ]], tx, ty, x.lo, x.hi, y.lo, y.hi ),
            threshold.x = tx, threshold.y = ty, region = regions[[ p ]] )
    out
}


#' Polygon gate (cells, all, ellipse)
#'
#' @title Free Gate Boundary
#'
#' @inheritParams calculate.gate
#'
#' @return Named list with one boundary entry.
#'
#' @keywords internal
calc.free <- function( data, gate.def, ctx ) {
    gp  <- gate.def$gate.param
    q   <- as.numeric( gp[[ "quantile.trim" ]] %||% ag.defaults$quantile.trim )
    pop <- names( gate.def$populations )[ 1L ] %||% "in"
    x.ch <- gate.def$gate.marker.x
    y.ch <- gate.def$gate.marker.y
    entry <- function( poly, rng = NULL ) {
        e <- list( boundary = poly, threshold.x = NA_real_, threshold.y = NA_real_,
                   region = "inside" )
        if ( !is.null( rng ) ) e$range <- rng
        stats::setNames( list( e ), pop )
    }
    if ( identical( gate.def$algorithm, "cells" ) && has.value( y.ch ) )
        return( calc.cells( data, gate.def, ctx, pop ) )
    if ( !has.value( y.ch ) ) {
        col <- trim.events( data[ , x.ch ], q )
        if ( length( col ) < ag.defaults$min.events )
            return( list() )
        return( entry( rect.polygon( min( col ), max( col ) ), c( min( col ), max( col ) ) ) )
    }
    cols <- trim.events( data[ , c( x.ch, y.ch ), drop = FALSE ], q )
    if ( nrow( cols ) < 4L )
        return( list() )
    poly <- if ( identical( gate.def$algorithm, "ellipse" ) ) ellipse.polygon( cols )
            else convex.hull.polygon( cols )
    if ( is.null( poly ) ) list() else entry( poly )
}


#' Cell gate on two scatter channels
#'
#' @title Cells Gate Boundary
#'
#' @description
#' A rectangle whose lower side on the forward-scatter axis
#' (\code{debris.axis}, "x" or "y") is the debris cut (\code{\link{debris.cut}},
#' multiplied by the threshold factor); on the events above the cut, the
#' other sides are where each axis's density falls to \code{edge.density}
#' (per axis \code{edge.density.x} / \code{edge.density.y}, default 0.002)
#' of its peak.  \code{cell.side = "debris"} keeps the events below the cut
#' instead.  A saturated pile at the top of either axis is cut off.
#'
#' @inheritParams calculate.gate
#' @param pop Population name.
#'
#' @return Named list with one boundary entry.
#'
#' @keywords internal
calc.cells <- function( data, gate.def, ctx, pop ) {
    gp <- gate.def$gate.param
    q  <- as.numeric( gp[[ "quantile.trim" ]] %||% ag.defaults$quantile.trim )
    chans <- c( gate.def$gate.marker.x, gate.def$gate.marker.y )
    xy <- data[ , chans, drop = FALSE ]
    if ( nrow( xy ) < 50L )
        return( list() )
    d <- if ( identical( gp[[ "debris.axis" ]], "y" ) ) 2L else 1L
    cut <- scaled.threshold( ctx, debris.cut( xy[ , d ], gp ), gp, "", chans[ d ] )
    debris <- identical( gp[[ "cell.side" ]], "debris" )
    kept <- if ( debris ) xy[ xy[ , d ] < cut, , drop = FALSE ] else xy[ xy[ , d ] >= cut, , drop = FALSE ]
    if ( nrow( kept ) < 50L ) kept <- xy
    sat <- as.numeric( gp[[ "saturation.fraction" ]] %||% 5e-4 )
    bounds <- matrix( NA_real_, 2L, 2L )
    for ( k in 1:2 ) {
        edge <- as.numeric( gp[[ paste0( "edge.density.", c( "x", "y" )[ k ] ) ]] %||%
                            gp[[ "edge.density" ]] %||% ag.defaults$edge.density )
        lh <- kde.extent( trim.events( kept[ , k ], q ), edge )
        if ( k == d ) {
            if ( debris ) lh[ 2 ] <- cut else lh[ 1 ] <- cut
        }
        pile <- saturation.cut( xy[ , k ], sat )
        if ( !is.null( pile ) ) lh[ 2 ] <- min( lh[ 2 ], pile )
        bounds[ k, ] <- lh
    }
    stats::setNames( list( list(
        boundary = rect.polygon( bounds[ 1, 1 ], bounds[ 1, 2 ], bounds[ 2, 1 ], bounds[ 2, 2 ] ),
        threshold.x = NA_real_, threshold.y = NA_real_, region = "inside" ) ), pop )
}


#' Debris cut of a cells gate
#'
#' @title Debris Cut
#'
#' @description
#' Forward-scatter value separating debris and dead cells from the cells,
#' before any threshold factor: a saturated pile is left out, the rest is
#' quantile-trimmed and cropped to the \code{axis.lo}/\code{axis.hi} window,
#' and the gate's \code{debris.algorithm} (default "mixture") places the
#' threshold with the gate's own parameters.
#'
#' @param x Forward-scatter values.
#' @param gate.param Gate parameters.
#'
#' @return Numeric cut.
#'
#' @keywords internal
debris.cut <- function( x, gate.param = list() ) {
    gp <- gate.param
    x <- as.numeric( x )
    pile <- saturation.cut( x, as.numeric( gp[[ "saturation.fraction" ]] %||% 5e-4 ) )
    if ( !is.null( pile ) ) x <- x[ x < pile ]
    col <- trim.events( x, as.numeric( gp[[ "quantile.trim" ]] %||% ag.defaults$quantile.trim ) )
    if ( length( col ) < ag.defaults$min.events )
        return( if ( length( col ) ) min( col ) else 0 )
    algo <- gp[[ "debris.algorithm" ]]
    if ( !has.value( algo ) || !( algo %in% ag.threshold.algorithms ) )
        algo <- ag.defaults$debris.algorithm
    calc <- crop.to.region( col,
        as.numeric( gp[[ "axis.lo" ]] %||% default.axis.lo( algo ) ),
        as.numeric( gp[[ "axis.hi" ]] %||% ag.defaults$axis.hi ) )
    if ( length( calc ) < ag.defaults$min.events ) calc <- col
    compute.threshold( calc, algo, gp )
}


#' Extent of the dense part of a distribution
#'
#' @title Density Extent
#'
#' @param x Numeric vector.
#' @param rel.density Fraction of the peak density.
#'
#' @return \code{c(lo, hi)}: the lowest and highest points where the
#'   density is at least \code{rel.density} of its peak.
#'
#' @keywords internal
kde.extent <- function( x, rel.density = ag.defaults$edge.density ) {
    x <- as.numeric( x )
    if ( length( x ) < 20L || max( x ) == min( x ) )
        return( if ( length( x ) ) range( x ) else c( 0, 0 ) )
    k <- kde.grid( x )
    keep <- which( k$density >= rel.density * max( k$density ) )
    k$grid[ c( keep[ 1L ], keep[ length( keep ) ] ) ]
}


#' Cut below a saturated pile
#'
#' @title Saturation Cut
#'
#' @param values Numeric vector.
#' @param min.fraction Fraction of events at the maximum that counts as a
#'   pile.
#'
#' @return Value just below the pile, or \code{NULL} when there is none.
#'
#' @keywords internal
saturation.cut <- function( values, min.fraction = 5e-4 ) {
    values <- as.numeric( values )
    if ( length( values ) == 0L ) return( NULL )
    v.max <- max( values )
    span <- v.max - min( values )
    if ( span == 0 ) span <- 1
    edge <- v.max - 1e-3 * span
    if ( mean( values >= edge ) <= min.fraction ) return( NULL )
    edge - 1e-3 * span
}


#' Convex hull of events
#'
#' @title Convex Hull Polygon
#'
#' @param xy Two-column matrix.
#'
#' @return Closed polygon, or \code{NULL}.
#'
#' @keywords internal
convex.hull.polygon <- function( xy ) {
    h <- grDevices::chull( xy[ , 1 ], xy[ , 2 ] )
    if ( length( h ) < 3L ) return( NULL )
    unname( xy[ c( h, h[ 1 ] ), , drop = FALSE ] )
}


#' Ellipse of events as a polygon
#'
#' @title Ellipse Polygon
#'
#' @description
#' The \code{n.sd} standard-deviation ellipse of the events as an
#' \code{n.points}-vertex polygon of the same area (models saved with the
#' retired "ellipse" algorithm).
#'
#' @param xy Two-column matrix.
#' @param n.points Number of vertices.
#' @param n.sd Size in standard deviations.
#'
#' @return Closed polygon, or \code{NULL}.
#'
#' @keywords internal
ellipse.polygon <- function( xy, n.points = 10L, n.sd = 3 ) {
    e <- tryCatch( eigen( stats::cov( xy ), symmetric = TRUE ), error = function( e ) NULL )
    if ( is.null( e ) ) return( NULL )
    vals <- rev( e$values )
    vecs <- e$vectors[ , 2:1 ]
    k <- sqrt( 2 * pi / ( n.points * sin( 2 * pi / n.points ) ) )
    a <- k * n.sd * sqrt( abs( vals[ 1 ] ) )
    b <- k * n.sd * sqrt( abs( vals[ 2 ] ) )
    angle <- atan2( vecs[ 2, 1 ], vecs[ 1, 1 ] )
    t <- seq( 0, 2 * pi, length.out = n.points + 1L )[ seq_len( n.points ) ]
    xr <- a * cos( t ) * cos( angle ) - b * sin( t ) * sin( angle ) + mean( xy[ , 1 ] )
    yr <- a * cos( t ) * sin( angle ) + b * sin( t ) * cos( angle ) + mean( xy[ , 2 ] )
    cbind( c( xr, xr[ 1 ] ), c( yr, yr[ 1 ] ) )
}
