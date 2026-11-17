# gate_figures.R
#
# Figures of the gates as applied.  Each gate is drawn on the events of its
# parent population, coloured by how crowded each part of the plot is (blue
# sparse, red dense), with the thresholds or polygons that were used and,
# for each population, its label and its share of the parent.  A gate on a
# single channel with no second channel to view it against is drawn as a
# histogram.
#
# Axes are drawn on the display scale the gates are calculated on (0 to 1
# for the top of the instrument range), with tick labels in raw units:
#   - linear channels (scatter, time, linear fluorescence) get evenly spaced
#     ticks labelled with a K / M / G suffix, whatever the instrument range;
#   - logicle and log channels get decade ticks (0, 10^2, 10^3, and the
#     negative decades below zero); decades too close to a neighbour to
#     read are left as unlabelled ticks.
# The axis range of a channel comes from its panel type: "scatter" channels
# use ap$figure.limits$scatter.min/max, fluorescence channels
# ap$figure.limits$expr.min/max (raw units, converted through the channel's
# transform), the time channel its data range.  Events beyond the axis range
# are drawn on its edge, as FlowJo does.
#
# Colour shows the rank of each event's density (a binned two-dimensional
# kernel density estimate over the plot), so the scale adapts to every plot
# rather than to the single densest point.


#' Axis ticks of a channel
#'
#' @title Axis Ticks
#'
#' @description
#' Decade ticks for logicle and log transforms, evenly spaced ticks in raw
#' units for linear transforms and the identity; see the file header of
#' \code{gate_figures.R}.
#'
#' @param params Transform parameters of the channel (\code{NULL} or empty
#'   for the identity).
#' @param limits Display limits \code{c(lo, hi)}.
#' @param forward,inverse Raw to display and display to raw functions.
#' @param min.gap Smallest distance between labelled decades, as a fraction
#'   of the axis.
#'
#' @return List with \code{breaks} (display units), \code{labels}
#'   (character or expression vector, one per break) and \code{minor}
#'   (display units).
#'
#' @keywords internal
axis.ticks <- function( params, limits, forward, inverse, min.gap = 0.05 ) {
    id <- if ( length( params ) ) params$id else NULL
    inside <- function( v ) is.finite( v ) & v >= limits[ 1 ] - 1e-9 & v <= limits[ 2 ] + 1e-9

    if ( !is.null( id ) && id %in% c( 1L, 2L ) ) {
        t.top <- as.numeric( params$scale.t )
        k.max <- floor( log10( t.top ) )
        if ( id == 1L ) {
            k <- seq_len( k.max )
            raw <- c( 0, 10^rev( k ), -10^rev( k ) )
        } else {
            k.min <- ceiling( log10( t.top ) - as.numeric( params$log.m ) )
            raw <- 10^rev( seq( k.min, k.max ) )
        }
        pos <- suppressWarnings( forward( raw ) )
        keep <- inside( pos )
        raw <- raw[ keep ]
        pos <- pos[ keep ]
        gap <- min.gap * diff( limits )
        labelled <- logical( length( pos ) )
        for ( i in seq_along( pos ) )
            labelled[ i ] <- all( abs( pos[ i ] - pos[ labelled ] ) >= gap )
        decade.label <- function( v ) {
            if ( v == 0 ) return( "0" )
            sprintf( "%s10^%d", if ( v < 0 ) "-" else "", as.integer( round( log10( abs( v ) ) ) ) )
        }
        ord <- order( pos[ labelled ] )
        breaks <- pos[ labelled ][ ord ]
        labels <- parse( text = vapply( raw[ labelled ][ ord ], decade.label, character( 1 ) ) )
        k.all <- seq( if ( id == 1L ) 0L else k.min, k.max )
        minor.raw <- as.vector( outer( 2:9, 10^k.all ) )
        if ( id == 1L ) minor.raw <- c( minor.raw, -minor.raw )
        minor <- suppressWarnings( forward( minor.raw ) )
        minor <- sort( c( minor[ inside( minor ) ], pos[ !labelled ] ) )
        return( list( breaks = breaks, labels = labels, minor = minor ) )
    }

    raw.lim <- sort( inverse( limits ) )
    raw <- pretty( raw.lim, n = 5L )
    raw <- raw[ raw >= raw.lim[ 1 ] - 1e-9 * diff( raw.lim ) & raw <= raw.lim[ 2 ] + 1e-9 * diff( raw.lim ) ]
    big <- if ( length( raw ) ) max( abs( raw ) ) else 0
    e <- if ( big >= 1000 ) min( 12, 3 * floor( log10( big ) / 3 ) ) else 0
    suffix <- c( "", "K", "M", "G", "T" )[ e / 3 + 1 ]
    labels <- vapply( raw / 10^e, function( v )
        format( signif( v, 4 ), scientific = FALSE, drop0trailing = TRUE, trim = TRUE ), character( 1 ) )
    labels <- ifelse( raw == 0, "0", paste0( labels, suffix ) )
    step <- if ( length( raw ) > 1L ) diff( raw[ 1:2 ] ) else diff( raw.lim )
    minor.raw <- seq( min( raw ) - step / 2, max( raw ) + step / 2, by = step / 2 )
    minor <- forward( minor.raw )
    minor <- minor[ inside( minor ) ]
    list( breaks = forward( raw ), labels = labels, minor = minor )
}


#' Axes of the gate figures
#'
#' @title Figure Axes
#'
#' @description
#' Title, display limits and ticks of every channel, from the panel type
#' and transform of each channel (see the file header of
#' \code{gate_figures.R}).
#'
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}.
#' @param ap Autogating parameter object (\code{figure.limits}).
#' @param data Event matrix, for the range of the time channel.
#'
#' @return Named list channel -> list(title, kind, limits, breaks, labels,
#'   minor); \code{kind} is "scatter", "time" or "fluorescence".
#'
#' @keywords internal
figure.axes <- function( fcs.setup, ap, data = fcs.setup$fcs.expr.data ) {
    panel <- fcs.setup$fcs.panel
    tr <- fcs.setup$ap
    lims <- ap$figure.limits
    titles <- fcs.setup$fcs.marker.label.figure
    identity.fn <- function( x ) x
    out <- list()
    for ( i in seq_len( nrow( panel ) ) ) {
        ch <- panel$antigen[ i ]
        if ( !( ch %in% colnames( data ) ) ) next
        type <- as.character( panel$type[ i ] )
        kind <- if ( type %in% c( "scatter", "time" ) ) type else "fluorescence"
        params <- tr$transform.params[[ ch ]]
        forward <- tr$transform[[ ch ]] %||% identity.fn
        inverse <- tr$transform.inv[[ ch ]] %||% identity.fn
        raw <- switch( kind,
                       scatter      = c( lims$scatter.min, lims$scatter.max ),
                       fluorescence = c( lims$expr.min, lims$expr.max ),
                       NULL )
        limits <- if ( length( raw ) == 2L ) sort( forward( raw ) ) else
            suppressWarnings( range( data[ , ch ], finite = TRUE ) )
        if ( !all( is.finite( limits ) ) || diff( limits ) <= 0 )
            limits <- c( 0, 1 )
        ticks <- axis.ticks( params, limits, forward, inverse )
        title <- if ( !is.null( titles ) && ch %in% names( titles ) ) titles[[ ch ]] else ch
        out[[ ch ]] <- c( list( title = title, kind = kind, limits = limits ), ticks )
    }
    out
}


#' Density of each event in a scatter plot
#'
#' @title Event Density
#'
#' @description
#' Binned two-dimensional kernel density estimate: events are counted on an
#' \code{n} by \code{n} grid over the plot and the counts smoothed with a
#' Gaussian kernel on each axis (standard deviation the normal-reference
#' bandwidth of \code{\link[MASS]{bandwidth.nrd}} / 4, as
#' \code{\link[MASS]{kde2d}}, times \code{bandwidth.factor}); each event
#' takes the density interpolated between the four nearest cell centres.
#'
#' @param x,y Event coordinates (inside the limits).
#' @param lim.x,lim.y Plot limits.
#' @param n Grid points per axis.
#' @param bandwidth.factor Multiplies the bandwidth.
#'
#' @return Numeric vector of densities (arbitrary units).
#'
#' @keywords internal
event.density <- function( x, y, lim.x, lim.y, n = 128L, bandwidth.factor = 1 ) {
    k <- length( x )
    if ( k < 3L )
        return( rep( 1, k ) )
    cell <- function( v, lim ) pmin( pmax( floor( ( v - lim[ 1 ] ) / diff( lim ) * n ) + 1L, 1L ), n )
    ix <- cell( x, lim.x )
    iy <- cell( y, lim.y )
    counts <- matrix( tabulate( ix + n * ( iy - 1L ), n * n ), n, n )
    smoother <- function( v, lim ) {
        h <- MASS::bandwidth.nrd( v ) / 4
        if ( !is.finite( h ) || h <= 0 ) h <- diff( lim ) / 100
        sd.cells <- max( h * bandwidth.factor / diff( lim ) * n, 0.5 )
        g <- seq_len( n )
        stats::dnorm( outer( g, g, "-" ), sd = sd.cells )
    }
    z <- smoother( x, lim.x ) %*% counts %*% t( smoother( y, lim.y ) )
    corner <- function( v, lim ) {
        u <- ( v - lim[ 1 ] ) / diff( lim ) * n + 0.5
        i <- pmin( pmax( floor( u ), 1L ), n - 1L )
        list( i = i, w = pmin( pmax( u - i, 0 ), 1 ) )
    }
    cx <- corner( x, lim.x )
    cy <- corner( y, lim.y )
    ( 1 - cx$w ) * ( 1 - cy$w ) * z[ cbind( cx$i, cy$i ) ] +
        cx$w * ( 1 - cy$w ) * z[ cbind( cx$i + 1L, cy$i ) ] +
        ( 1 - cx$w ) * cy$w * z[ cbind( cx$i, cy$i + 1L ) ] +
        cx$w * cy$w * z[ cbind( cx$i + 1L, cy$i + 1L ) ]
}


#' Channels a gate is drawn on
#'
#' @title Figure Channels
#'
#' @param gate.def Effective gate definition.
#' @param boundary Boundary entries of the gate.
#'
#' @return Character vector: x channel, then y channel (absent for a
#'   histogram).
#'
#' @keywords internal
figure.channels <- function( gate.def, boundary ) {
    if ( is.shape.boundary( boundary ) ) {
        ch <- unlist( boundary[[ 1L ]]$shape.channels )
        return( ch[ vapply( ch, has.value, logical( 1 ) ) ] )
    }
    if ( identical( gate.type.of( gate.def ), "1dsep" ) ) {
        sep <- gate.def$gate.marker.x
        view <- gate.def$gate.marker.y
        if ( !has.value( view ) ) return( sep )
        return( if ( separation.axis( gate.def ) == "y" ) c( view, sep ) else c( sep, view ) )
    }
    gate.channels( gate.def )
}


#' Lines, outlines and labels of a gate figure
#'
#' @title Gate Annotations
#'
#' @description
#' Thresholds become lines across the plot; polygons and drawn shapes
#' become outlines.  Each population label shows its share of the parent:
#' threshold populations are labelled in the outer corner of their region
#' (for three quadrants, the corner opposite the missing one), polygon
#' populations above their outline (below it when the outline reaches the
#' top of the plot).
#'
#' @param gate.def Effective gate definition.
#' @param boundary Boundary entries as applied (with \code{fraction}).
#' @param chans Channels drawn (\code{\link{figure.channels}}).
#' @param lim.x,lim.y Plot limits (\code{lim.y} \code{NULL} for a
#'   histogram).
#'
#' @return List with \code{lines} (\code{data.frame} axis, at),
#'   \code{paths} (\code{data.frame} x, y, group) and \code{labels}
#'   (\code{data.frame} x, y, text, hjust, vjust).
#'
#' @keywords internal
gate.annotations <- function( gate.def, boundary, chans, lim.x, lim.y = NULL ) {
    histogram <- length( chans ) < 2L
    inset.x <- 0.02 * diff( lim.x )
    inset.y <- if ( histogram ) 0 else 0.02 * diff( lim.y )
    clamp <- function( v, lim ) pmin( pmax( v, lim[ 1 ] ), lim[ 2 ] )
    lines <- data.frame( axis = character( 0 ), at = numeric( 0 ), stringsAsFactors = FALSE )
    paths <- data.frame( x = numeric( 0 ), y = numeric( 0 ), group = character( 0 ),
                         stringsAsFactors = FALSE )
    labels <- data.frame( x = numeric( 0 ), y = numeric( 0 ), text = character( 0 ),
                          hjust = numeric( 0 ), vjust = numeric( 0 ), stringsAsFactors = FALSE )
    add.line <- function( axis, at, lim )
        lines[ nrow( lines ) + 1L, ] <<- list( axis, clamp( at, lim ) )
    add.label <- function( p, x, y, hjust, vjust ) {
        lab <- gate.def$populations[[ p ]]$label
        if ( !has.value( lab ) ) lab <- p
        f <- boundary[[ p ]]$fraction
        text <- if ( is.numeric( f ) && length( f ) == 1L && is.finite( f ) )
            sprintf( "%s\n%.1f%%", lab, 100 * f ) else lab
        labels[ nrow( labels ) + 1L, ] <<- list( x, y, text, hjust, vjust )
    }
    top.y <- if ( histogram ) Inf else lim.y[ 2 ] - inset.y
    top.vjust <- if ( histogram ) 1.2 else 1
    outline.label <- function( p, poly ) {
        if ( histogram ) return( invisible() )
        top <- max( poly[ , 2 ] )
        x <- clamp( mean( range( poly[ , 1 ] ) ), lim.x )
        if ( top + 0.12 * diff( lim.y ) < lim.y[ 2 ] )
            add.label( p, x, top + inset.y, 0.5, 0 )
        else
            add.label( p, x, max( min( poly[ , 2 ] ) - inset.y, lim.y[ 1 ] + inset.y ), 0.5,
                       if ( min( poly[ , 2 ] ) - inset.y > lim.y[ 1 ] + inset.y ) 1 else 0 )
    }

    if ( is.shape.boundary( boundary ) || !( gate.type.of( gate.def ) %in% ag.threshold.types ) ) {
        for ( p in names( boundary ) ) {
            e <- boundary[[ p ]]
            if ( histogram && length( e$range ) == 2L ) {
                for ( v in e$range ) add.line( "x", v, lim.x )
                add.label( p, clamp( mean( e$range ), lim.x ), top.y, 0.5, top.vjust )
                next
            }
            poly <- e$boundary
            if ( histogram || !is.matrix( poly ) || nrow( poly ) < 2L ) next
            if ( any( poly[ 1, ] != poly[ nrow( poly ), ] ) ) poly <- rbind( poly, poly[ 1, ] )
            paths <- rbind( paths, data.frame( x = clamp( poly[ , 1 ], lim.x ), y = clamp( poly[ , 2 ], lim.y ),
                                               group = p, stringsAsFactors = FALSE ) )
            outline.label( p, poly )
        }
        return( list( lines = lines, paths = paths, labels = labels ) )
    }

    first <- boundary[[ 1L ]]
    regions <- population.regions( gate.def, boundary )
    if ( identical( gate.type.of( gate.def ), "1dsep" ) ) {
        t <- first$threshold.x
        if ( !has.value( t ) ) return( list( lines = lines, paths = paths, labels = labels ) )
        on.y <- !histogram && separation.axis( gate.def ) == "y"
        if ( on.y ) add.line( "y", t, lim.y ) else add.line( "x", t, lim.x )
        for ( p in names( boundary ) ) {
            below <- identical( regions[ p ][[ 1 ]], "below" )
            if ( on.y )
                add.label( p, lim.x[ 2 ] - inset.x,
                           if ( below ) lim.y[ 1 ] + inset.y else lim.y[ 2 ] - inset.y,
                           1, if ( below ) 0 else 1 )
            else
                add.label( p, if ( below ) lim.x[ 1 ] + inset.x else lim.x[ 2 ] - inset.x,
                           top.y, if ( below ) 0 else 1, top.vjust )
        }
        return( list( lines = lines, paths = paths, labels = labels ) )
    }

    tx <- first$threshold.x
    ty <- first$threshold.y
    if ( has.value( tx ) ) add.line( "x", tx, lim.x )
    if ( has.value( ty ) ) add.line( "y", ty, lim.y )
    flip <- function( s ) if ( s == "-" ) "+" else "-"
    for ( p in names( boundary ) ) {
        q <- region.quadrants( regions[ p ][[ 1 ]] )
        if ( length( q ) == 0L ) next
        if ( length( q ) == 3L ) {
            m <- setdiff( ag.quadrant.regions, q )
            q <- sprintf( "x%sy%s", flip( substr( m, 2L, 2L ) ), flip( substr( m, 4L, 4L ) ) )
        }
        sx <- substr( q[ 1 ], 2L, 2L )
        sy <- substr( q[ 1 ], 4L, 4L )
        add.label( p, if ( sx == "-" ) lim.x[ 1 ] + inset.x else lim.x[ 2 ] - inset.x,
                   if ( sy == "-" ) lim.y[ 1 ] + inset.y else lim.y[ 2 ] - inset.y,
                   if ( sx == "-" ) 0 else 1, if ( sy == "-" ) 0 else 1 )
    }
    list( lines = lines, paths = paths, labels = labels )
}


#' Figure of one gate
#'
#' @title Gate Figure
#'
#' @description
#' Draws a gate on its parent population's events as applied in a gating
#' result; see the file header of \code{gate_figures.R}.  At most
#' \code{ap$figure.max.events} events are drawn (a seeded subsample, the
#' caller's random number stream is left as it was).
#'
#' @param gate.def Gate definition (a replicate gate is drawn with its
#'   origin's channels; see \code{\link{effective.gate.def}}).
#' @param result Gating result (\code{\link{run.gating}}).
#' @param data Event matrix the result refers to.
#' @param axes Figure axes (\code{\link{figure.axes}}).
#' @param ap Autogating parameter object (\code{figure.*},
#'   \code{density.grid.n}, \code{bandwidth.factor}).
#' @param title Optional text put before the gate name (e.g. the sample).
#' @param seed Seed of the subsample.
#'
#' @return A \pkg{ggplot2} object, or \code{NULL} when the gate was not
#'   applied or a channel has no axis.
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline geom_hline geom_path
#' @importFrom ggplot2 geom_label scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 scale_colour_gradientn labs theme_bw theme element_blank
#' @importFrom ggplot2 element_text guides guide_axis expansion
#' @importFrom scattermore geom_scattermore
#'
#' @export
gate.figure <- function( gate.def, result, data, axes, ap, title = NULL, seed = 0L ) {
    name <- gate.def$gate.name
    boundary <- result$boundaries[[ name ]]
    rows <- result$parent.events[[ name ]]
    if ( length( boundary ) == 0L || is.null( rows ) )
        return( NULL )
    chans <- figure.channels( gate.def, boundary )
    if ( length( chans ) == 0L || !all( chans %in% names( axes ) ) || !all( chans %in% colnames( data ) ) )
        return( NULL )
    histogram <- length( chans ) < 2L
    ax <- axes[[ chans[ 1 ] ]]
    ay <- if ( histogram ) NULL else axes[[ chans[ 2 ] ]]

    n.max <- as.integer( ap$figure.max.events %||% 100000L )
    if ( length( rows ) > n.max )
        rows <- sort( rows[ seeded.call( seed, function() sample.int( length( rows ), n.max ) ) ] )
    squish <- function( v, lim ) {
        v[ !is.finite( v ) ] <- lim[ 1 ]
        pmin( pmax( v, lim[ 1 ] ), lim[ 2 ] )
    }
    x <- squish( data[ rows, chans[ 1 ] ], ax$limits )
    ann <- gate.annotations( gate.def, boundary, chans, ax$limits, ay$limits )

    scale.for <- function( a, axis ) {
        f <- if ( axis == "x" ) scale_x_continuous else scale_y_continuous
        f( name = a$title, breaks = a$breaks, labels = a$labels, minor_breaks = a$minor,
           limits = a$limits, expand = expansion( mult = 0.01 ) )
    }

    if ( histogram ) {
        fig <- ggplot( data.frame( x = x ), aes( x = x ) ) +
            geom_histogram( bins = 128L, boundary = ax$limits[ 1 ], fill = "grey55", colour = NA,
                            na.rm = TRUE ) +
            scale.for( ax, "x" ) +
            scale_y_continuous( name = "Events", expand = expansion( mult = c( 0, 0.12 ) ) )
    } else {
        y <- squish( data[ rows, chans[ 2 ] ], ay$limits )
        z <- event.density( x, y, ax$limits, ay$limits, n = as.integer( ap$density.grid.n %||% 128L ),
                            bandwidth.factor = ap$bandwidth.factor %||% 1 )
        z <- if ( length( z ) ) rank( z ) / length( z ) else z
        raster <- as.integer( round( 0.7 * min( ap$figure.width %||% 6, ap$figure.height %||% 6 ) *
                                     ( ap$figure.resolution %||% 150 ) ) )
        o <- order( z )
        fig <- ggplot( data.frame( x = x[ o ], y = y[ o ], z = z[ o ] ), aes( x = x, y = y, colour = z ) ) +
            geom_scattermore( pointsize = ap$figure.point.size %||% 1, pixels = rep( raster, 2L ) ) +
            scale_colour_gradientn( colours = c( "#0000FF", "#00FFFF", "#00FF00", "#FFFF00", "#FF0000" ),
                                    limits = c( 0, 1 ), guide = "none" ) +
            scale.for( ax, "x" ) + scale.for( ay, "y" )
    }

    v <- ann$lines[ ann$lines$axis == "x", , drop = FALSE ]
    h <- ann$lines[ ann$lines$axis == "y", , drop = FALSE ]
    if ( nrow( v ) ) fig <- fig + geom_vline( data = v, aes( xintercept = at ), colour = "black", linewidth = 0.6 )
    if ( nrow( h ) ) fig <- fig + geom_hline( data = h, aes( yintercept = at ), colour = "black", linewidth = 0.6 )
    if ( nrow( ann$paths ) )
        fig <- fig + geom_path( data = ann$paths, aes( x = x, y = y, group = group ),
                                inherit.aes = FALSE, colour = "black", linewidth = 0.6 )
    if ( nrow( ann$labels ) )
        fig <- fig + geom_label( data = ann$labels,
                                 aes( x = x, y = y, label = text, hjust = hjust, vjust = vjust ),
                                 inherit.aes = FALSE, size = ap$figure.label.size %||% 3.5,
                                 lineheight = 0.9,
                                 fill = sprintf( "#FFFFFF%02X", as.integer( round( 255 * ap$figure.label.alpha %||% 0.7 ) ) ) )

    parent <- if ( has.value( gate.def$parent.gate ) )
        sprintf( "%s/%s", gate.def$parent.gate, gate.def$parent.popul %||% "" ) else "all events"
    status <- result$status[[ name ]]
    sub <- sprintf( "Parent %s, %s events", parent,
                    format( length( result$parent.events[[ name ]] ), big.mark = "," ) )
    if ( has.value( status$source ) ) sub <- paste0( sub, "; boundary ", status$source )
    if ( has.value( status$flag ) )
        sub <- paste0( sub, "; flagged ", status$flag,
                       if ( nzchar( status$message %||% "" ) ) paste0( ": ", status$message ) else "" )
    text.size <- ap$figure.text.size %||% 11
    wrap <- function( s, size ) paste( strwrap( s, width = max( 20L, floor( 0.8 * 72 * ( ap$figure.width %||% 6 ) /
                                                                          ( 0.6 * size ) ) ) ),
                                       collapse = "\n" )
    fig <- fig + labs( title = wrap( if ( is.null( title ) ) name else paste0( title, ": ", name ),
                                     1.2 * text.size ),
                       subtitle = wrap( sub, 0.75 * text.size ) ) +
        theme_bw( base_size = text.size ) +
        theme( aspect.ratio = 1, panel.grid = element_blank(),
               plot.subtitle = element_text( size = 0.75 * text.size ) )
    if ( utils::packageVersion( "ggplot2" ) >= "3.5.0" ) {
        minor.guide <- function() guide_axis( minor.ticks = TRUE )
        fig <- fig + guides( x = minor.guide(), y = if ( histogram ) "axis" else minor.guide() )
    }
    fig
}


#' Write a figure of every gate
#'
#' @title Write Gate Figures
#'
#' @description
#' One PNG per applied gate (\code{\link{gate.figure}}), named by gate
#' number and gate name, written with \pkg{ragg}.
#'
#' @param gate.defs List of gate definitions.
#' @param result Gating result (\code{\link{run.gating}}).
#' @param data Event matrix the result refers to.
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}
#'   (panel types, transforms and axis titles).
#' @param ap Autogating parameter object (\code{figure.*}).
#' @param dir Output directory (created when missing).
#' @param title Optional text put before each gate name.
#' @param seed Seed of the event subsamples.
#'
#' @return Character vector of the files written, invisibly.
#'
#' @export
write.gate.figures <- function( gate.defs, result, data, fcs.setup, ap, dir, title = NULL,
                                seed = ap$seed %||% 0L ) {
    axes <- figure.axes( fcs.setup, ap, data )
    by.name <- gates.by.name( gate.defs )
    dir.create( dir, recursive = TRUE, showWarnings = FALSE )
    width <- as.integer( ap$gate.number.width %||% 3L )
    files <- character( 0 )
    ordered <- ordered.gate.defs( gate.defs )
    for ( i in seq_along( ordered ) ) {
        g <- ordered[[ i ]]
        fig <- gate.figure( effective.gate.def( g, by.name ), result, data, axes, ap,
                            title = title, seed = seed )
        if ( is.null( fig ) ) next
        number <- if ( has.value( g$gate.number ) ) as.integer( g$gate.number ) else i
        file <- file.path( dir, sprintf( "%0*d_%s.png", width, number,
                                         gsub( "[^A-Za-z0-9._+-]+", "_", g$gate.name ) ) )
        ggplot2::ggsave( file, fig, device = ragg::agg_png,
                         width = ap$figure.width %||% 6, height = ap$figure.height %||% 6,
                         units = "in", dpi = ap$figure.resolution %||% 150 )
        files <- c( files, file )
    }
    invisible( files )
}
