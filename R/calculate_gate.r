# calculate_gate.r
#
# Gate boundary calculation functions.
#
# Each function receives the data for the parent population and returns a named
# list of boundary matrices (one per non-ignored population).
#
# Shared helpers trim.quantiles() and make.rect.boundary() are defined in
# analyze_aux.r and expected to be available in the calling environment.

# ==============================================================================
# calculate.gate.free
# ==============================================================================

#' Calculate gate boundaries for a free-draw (all-events or singlets) gate
#'
#' @title Calculate Free-Draw Gate Boundary
#'
#' @description
#' Computes the boundary polygon for a single-population free gate.  Two
#' algorithms are supported:
#' \describe{
#'   \item{\code{"replicate"}}{Copies the boundary from a previously computed
#'     gate identified by \code{p$origin.gate}.}
#'   \item{any other string}{Dispatches to
#'     \code{calculate.boundary.<algorithm>()} (e.g.
#'     \code{\link{calculate.boundary.all}} or
#'     \code{\link{calculate.boundary.singlets}}).}
#' }
#'
#' @param gate.data Numeric matrix (n x 2).  Events from the parent
#'   population.
#' @param popul.names Character vector of length 1.  Name of the single
#'   population to gate.
#' @param algorithm Character scalar.  Boundary algorithm name (e.g.
#'   \code{"all"}, \code{"singlets"}, \code{"replicate"}).
#' @param p Named list.  Resolved gate parameters from
#'   \code{\link{resolve.gate.param}}.
#' @param gate.number Integer.  Gate index used for figure filenames.
#' @param gate.name Character.  Gate identifier used for figure filenames.
#' @param fcs.gates Named list.  Current gate results; needed only when
#'   \code{algorithm == "replicate"}.
#' @param figure.dir Character.  Output directory (currently unused for free
#'   gates but kept for interface consistency).
#' @param ap Autogating parameter object.
#'
#' @return Named list with one element (named by \code{popul.names[1]})
#'   containing the boundary matrix (n x 2).
#'
#' @seealso \code{\link{calculate.boundary.all}},
#'   \code{\link{calculate.boundary.singlets}},
#'   \code{\link{do.gate.generic}}
#'
#' @keywords internal
calculate.gate.free <- function( gate.data, popul.names, algorithm, p,
    gate.number, gate.name, fcs.gates, figure.dir, ap ) {

    stopifnot( length( algorithm ) == 1L && length( popul.names ) == 1L )

    boundary.list <- list()
    pop <- popul.names[ 1L ]

    if ( algorithm == "replicate" ) {

        origin <- fcs.gates[[ p$origin.gate ]]
        stopifnot( length( origin ) == 1L )
        boundary.list[[ pop ]] <- origin[[ 1L ]]$boundary

    } else {

        data.trim <- trim.quantiles( gate.data, p$quantile.x, p$quantile.y )

        calc.fn <- get( paste0( "calculate.boundary.", algorithm ) )
        message( sprintf( "  -> %s()", paste0( "calculate.boundary.", algorithm ) ) )
        boundary.list[[ pop ]] <- calc.fn( data.trim, p, ap )
    }

    boundary.list
}


# ==============================================================================
# calculate.gate.1dsep
# ==============================================================================

#' Calculate gate boundaries for a one-dimensional threshold gate
#'
#' @title Calculate 1-D Separating Gate Boundary
#'
#' @description
#' Computes axis-aligned rectangular boundary polygons that separate events
#' into a negative (low) and an optional positive (high) population along one
#' axis.  The threshold is determined by the selected algorithm (e.g.
#' \code{"tail"}, \code{"bimodal"}, \code{"mixture"}, \code{"kde.min"}).  A
#' diagnostic density plot is saved to \code{figure.dir}.
#'
#' Either or both population names may be the ignore sentinel (\code{"*"}),
#' in which case the corresponding boundary is omitted from the output.
#'
#' @param gate.data Numeric matrix (n x 2).  Events from the parent
#'   population.
#' @param popul.names Character vector of length 1 or 2.  Names of the
#'   negative (\code{[1]}) and positive (\code{[2]}) populations.
#' @param algorithm Character scalar.  Threshold algorithm name (one of
#'   \code{"tail"}, \code{"bimodal"}, \code{"mixture"}, \code{"kde.min"},
#'   or \code{"replicate"}).
#' @param p Named list.  Resolved gate parameters, including \code{sep.axis}
#'   (1 = separate on x, 2 = separate on y), \code{region.min},
#'   \code{region.max}, \code{limit.xmin}, \code{limit.xmax},
#'   \code{limit.ymin}, \code{limit.ymax}, and
#'   \code{default.tail.probability}.
#' @param gate.number Integer.  Gate index for figure filenames.
#' @param gate.name Character.  Gate identifier for figure filenames.
#' @param fcs.gates Named list.  Current gate results (needed only when
#'   \code{algorithm == "replicate"}).
#' @param figure.dir Character.  Output directory for the diagnostic PNG.
#' @param ap Autogating parameter object.
#'
#' @return Named list of 4 x 2 boundary matrices, one per non-ignored
#'   population.
#'
#' @seealso \code{\link{calculate.threshold.tail}},
#'   \code{\link{calculate.threshold.bimodal}},
#'   \code{\link{calculate.threshold.mixture}},
#'   \code{\link{calculate.threshold.kde.min}},
#'   \code{\link{do.gate.generic}}
#'
#' @keywords internal
calculate.gate.1dsep <- function( gate.data, popul.names, algorithm, p,
    gate.number, gate.name, fcs.gates, figure.dir, ap ) {

    stopifnot( length( algorithm ) == 1L )

    boundary.list <- list()

    if ( algorithm == "replicate" ) {

        origin <- fcs.gates[[ p$origin.gate ]]
        stopifnot( length( popul.names ) == length( origin ) )

        for ( i in seq_along( popul.names ) ) {
            pop <- popul.names[ i ]
            if ( pop != "*" )
                boundary.list[[ pop ]] <- origin[[ i ]]$boundary
        }

        return( boundary.list )
    }

    stopifnot( length( popul.names ) >= 1L && length( popul.names ) <= 2L )

    sep.axis   <- as.integer( p$sep.axis )
    region.min <- p$region.min
    region.max <- p$region.max
    lim.xmin   <- p$limit.xmin
    lim.xmax   <- p$limit.xmax
    lim.ymin   <- p$limit.ymin
    lim.ymax   <- p$limit.ymax

    data.trim <- trim.quantiles( gate.data, p$quantile.x, p$quantile.y )

    # Work on x-axis; transpose temporarily if separating on y
    if ( sep.axis == 2L )
        data.trim <- data.trim[ , c( 2L, 1L ) ]

    x.rng  <- range( data.trim[ , 1L ] )
    y.rng  <- range( data.trim[ , 2L ] )

    data.calc <- data.trim[
        data.trim[ , 1L ] >= ( 1 - region.min ) * x.rng[ 1 ] +
                              region.min        * x.rng[ 2 ] &
        data.trim[ , 1L ] <= ( 1 - region.max ) * x.rng[ 1 ] +
                              region.max        * x.rng[ 2 ], , drop = FALSE ]

    calc.fn  <- get( paste0( "calculate.threshold.", algorithm ) )
    message( sprintf( "  -> %s()", paste0( "calculate.threshold.", algorithm ) ) )
    x.thr    <- calc.fn( data.calc[ , 1L ], p$default.tail.probability )

    # Diagnostic density plot
    grDevices::png( filename = file.path( figure.dir,
            sprintf( "%0*d - %s - calculation.png",
                ap$gate.number.width, gate.number, gate.name ) ),
        width = 1024, height = 768 )
    graphics::par( mar = c( 5, 5, 4, 1 ) )
    plot( stats::density( data.trim[ , 1L ] ),
        xlab  = colnames( data.trim )[ 1L ],
        main  = sprintf( "%s  -  1dsep %s", gate.name, algorithm ),
        lwd   = 3, cex.lab = 2, cex.axis = 2, cex.main = 2.5, font.main = 1 )
    graphics::abline( v = range( data.calc[ , 1L ] ), lwd = 2, lty = 2 )
    graphics::abline( v = x.thr, lwd = 2, col = "blue3" )
    grDevices::dev.off()

    # Negative population (popul.names[1])
    if ( length( popul.names ) >= 1L && popul.names[ 1L ] != "*" ) {
        xmin <- ( 1 - lim.xmin ) * x.rng[ 1 ] + lim.xmin * x.rng[ 2 ]
        boundary.list[[ popul.names[ 1L ] ]] <- make.rect.boundary(
            xmin, x.thr,
            ( 1 - lim.ymin ) * y.rng[ 1 ] + lim.ymin * y.rng[ 2 ],
            ( 1 - lim.ymax ) * y.rng[ 1 ] + lim.ymax * y.rng[ 2 ],
            swap = ( sep.axis == 2L )
        )
    }

    # Positive population (popul.names[2])
    if ( length( popul.names ) >= 2L && popul.names[ 2L ] != "*" ) {
        xmax <- ( 1 - lim.xmax ) * x.rng[ 1 ] + lim.xmax * x.rng[ 2 ]
        boundary.list[[ popul.names[ 2L ] ]] <- make.rect.boundary(
            x.thr, xmax,
            ( 1 - lim.ymin ) * y.rng[ 1 ] + lim.ymin * y.rng[ 2 ],
            ( 1 - lim.ymax ) * y.rng[ 1 ] + lim.ymax * y.rng[ 2 ],
            swap = ( sep.axis == 2L )
        )
    }

    boundary.list
}


# ==============================================================================
# calculate.gate.2dsep
# ==============================================================================

#' Calculate gate boundaries for a two-dimensional threshold gate
#'
#' @title Calculate 2-D Separating Gate Boundary
#'
#' @description
#' Computes axis-aligned rectangular boundary polygons separating events into
#' up to four quadrants (double-negative, x-positive, y-positive,
#' double-positive) using independent thresholds on each axis.  Optionally,
#' a \code{join} parameter can merge adjacent quadrants into a single L-shaped
#' boundary (currently \code{"1:2:3"} is supported, combining DN, x-positive,
#' and y-positive).  A two-panel diagnostic density plot is saved to
#' \code{figure.dir}.
#'
#' @param gate.data Numeric matrix (n x 2).  Events from the parent
#' @param popul.names Character vector of length 4 (when \code{join} is
#'   \code{NULL}), length 1 (joined polygon only), or length 2 (joined polygon
#'   plus named complement population) when \code{join} is set.
#' @param algorithm Character vector of length 2.  Threshold algorithm for the
#'   x-axis (\code{[1]}) and y-axis (\code{[2]}), or a length-1 vector
#'   \code{"replicate"}.
#' @param p Named list.  Resolved gate parameters, including
#'   \code{region.xmin}, \code{region.xmax}, \code{region.ymin},
#'   \code{region.ymax}, limit fields, \code{default.tail.probability}, and
#'   optionally \code{join} (character or \code{NULL}).
#' @param gate.number Integer.  Gate index for figure filenames.
#' @param gate.name Character.  Gate identifier for figure filenames.
#' @param fcs.gates Named list.  Current gate results (needed only when
#'   \code{algorithm == "replicate"}).
#' @param figure.dir Character.  Output directory for the diagnostic PNG.
#' @param ap Autogating parameter object.
#'
#' @return Named list of boundary matrices (4 x 2 each), one per non-ignored
#'   population, or a single L-shaped polygon if \code{join} is set.
#'
#' @seealso \code{\link{calculate.gate.1dsep}}, \code{\link{do.gate.generic}}
#'
#' @keywords internal
calculate.gate.2dsep <- function( gate.data, popul.names, algorithm, p,
    gate.number, gate.name, fcs.gates, figure.dir, ap ) {

    join <- p$join    # NULL or e.g. "1:2:3"

    stopifnot( length( algorithm ) == 2L ||
        ( length( algorithm ) == 1L && algorithm[ 1L ] == "replicate" ) )

    boundary.list <- list()

    if ( algorithm[ 1L ] == "replicate" ) {

        origin <- fcs.gates[[ p$origin.gate ]]
        stopifnot( is.null( join ) &&
            length( popul.names ) == length( origin ) )

        for ( i in seq_along( popul.names ) ) {
            pop <- popul.names[ i ]
            if ( pop != "*" )
                boundary.list[[ pop ]] <- origin[[ i ]]$boundary
        }

        return( boundary.list )
    }

    # ---- standard 2-D separation ---------------------------------------------

    stopifnot(
        ( is.null( join ) && length( popul.names ) == 4L ) ||
        ( !is.null( join ) && length( popul.names ) %in% c(1L, 2L) )
    )

    lim.xmin    <- p$limit.xmin
    lim.xmax    <- p$limit.xmax
    lim.ymin    <- p$limit.ymin
    lim.ymax    <- p$limit.ymax
    region.xmin <- p$region.xmin
    region.xmax <- p$region.xmax
    region.ymin <- p$region.ymin
    region.ymax <- p$region.ymax

    data.trim <- trim.quantiles( gate.data, p$quantile.x, p$quantile.y )

    x.rng <- range( data.trim[ , 1L ] )
    y.rng <- range( data.trim[ , 2L ] )

    data.calc <- data.trim[
        data.trim[ , 1L ] >= ( 1 - region.xmin ) * x.rng[ 1 ] +
                              region.xmin         * x.rng[ 2 ] &
        data.trim[ , 1L ] <= ( 1 - region.xmax ) * x.rng[ 1 ] +
                              region.xmax         * x.rng[ 2 ] &
        data.trim[ , 2L ] >= ( 1 - region.ymin ) * y.rng[ 1 ] +
                              region.ymin         * y.rng[ 2 ] &
        data.trim[ , 2L ] <= ( 1 - region.ymax ) * y.rng[ 1 ] +
                              region.ymax         * y.rng[ 2 ], , drop = FALSE ]

    calc.fn.x <- get( paste0( "calculate.threshold.", algorithm[ 1L ] ) )
    calc.fn.y <- get( paste0( "calculate.threshold.", algorithm[ 2L ] ) )

    message( sprintf( "  -> %s() / %s()",
                      paste0( "calculate.threshold.", algorithm[ 1L ] ),
                      paste0( "calculate.threshold.", algorithm[ 2L ] ) ) )

    x.thr <- calc.fn.x( data.calc[ , 1L ], p$default.tail.probability )
    y.thr <- calc.fn.y( data.calc[ , 2L ], p$default.tail.probability )

    # Diagnostic density plot
    grDevices::png( filename = file.path( figure.dir,
            sprintf( "%0*d - %s - calculation.png",
                ap$gate.number.width, gate.number, gate.name ) ),
        width = 1024, height = 1536 )
    graphics::par( mfrow = c( 2L, 1L ), mar = c( 5, 5, 1, 1 ), oma = c( 0, 0, 3, 0 ) )
    plot( stats::density( data.trim[ , 1L ] ),
        xlab = colnames( data.trim )[ 1L ], main = "",
        lwd = 3, cex.lab = 2, cex.axis = 2 )
    graphics::abline( v = range( data.calc[ , 1L ] ), lwd = 2, lty = 2 )
    graphics::abline( v = x.thr, lwd = 2, col = "blue3" )
    plot( stats::density( data.trim[ , 2L ] ),
        xlab = colnames( data.trim )[ 2L ], main = "",
        lwd = 3, cex.lab = 2, cex.axis = 2 )
    graphics::abline( v = range( data.calc[ , 2L ] ), lwd = 2, lty = 2 )
    graphics::abline( v = y.thr, lwd = 2, col = "blue3" )
    graphics::mtext( sprintf( "%s  -  2dsep %s", gate.name,
            paste0( algorithm, collapse = ":" ) ),
        outer = TRUE, cex = 2.5 )
    grDevices::dev.off()

    xmin.full <- ( 1 - lim.xmin ) * x.rng[ 1 ] + lim.xmin * x.rng[ 2 ]
    xmax.full <- ( 1 - lim.xmax ) * x.rng[ 1 ] + lim.xmax * x.rng[ 2 ]
    ymin.full <- ( 1 - lim.ymin ) * y.rng[ 1 ] + lim.ymin * y.rng[ 2 ]
    ymax.full <- ( 1 - lim.ymax ) * y.rng[ 1 ] + lim.ymax * y.rng[ 2 ]

    if ( is.null( join ) ) {

        # Four-quadrant layout: DN, Xpos, Ypos, DP
        if ( popul.names[ 1L ] != "*" )
            boundary.list[[ popul.names[ 1L ] ]] <-
                make.rect.boundary( xmin.full, x.thr, ymin.full, y.thr )

        if ( popul.names[ 2L ] != "*" )
            boundary.list[[ popul.names[ 2L ] ]] <-
                make.rect.boundary( x.thr, xmax.full, ymin.full, y.thr )

        if ( popul.names[ 3L ] != "*" )
            boundary.list[[ popul.names[ 3L ] ]] <-
                make.rect.boundary( xmin.full, x.thr, y.thr, ymax.full )

        if ( popul.names[ 4L ] != "*" )
            boundary.list[[ popul.names[ 4L ] ]] <-
                make.rect.boundary( x.thr, xmax.full, y.thr, ymax.full )

    } else {

      # ---- joined L-shaped boundary ----------------------------------------
      #
      # Quadrant index convention:
      #   1 = DN   [xmin, xthr] x [ymin, ythr]
      #   2 = Xpos [xthr, xmax] x [ymin, ythr]
      #   3 = Ypos [xmin, xthr] x [ythr, ymax]
      #   4 = DP   [xthr, xmax] x [ythr, ymax]
      #
      # join is a sorted integer vector of the three quadrant indices to merge.
      # The complement quadrant (if popul.names[2] is given and not "*") gets
      # a standard rectangular boundary.

      complement <- setdiff( 1L:4L, join )   # single integer

      # Rectangle lookup: list indexed 1..4 -> make.rect.boundary() args
      quad.rect <- list(
        make.rect.boundary( xmin.full, x.thr,     ymin.full, y.thr     ),  # 1 DN
        make.rect.boundary( x.thr,     xmax.full, ymin.full, y.thr     ),  # 2 Xpos
        make.rect.boundary( xmin.full, x.thr,     y.thr,     ymax.full ),  # 3 Ypos
        make.rect.boundary( x.thr,     xmax.full, y.thr,     ymax.full )   # 4 DP
      )

      # L-shaped polygon vertices (clockwise) for every possible complement.
      # Each shape is the full bounding box minus one corner rectangle.
      l.shape <- switch( complement,
                         # complement = 1 (DN missing): Xpos + Ypos + DP
                         rbind( c( x.thr,     ymin.full ),
                                c( xmax.full, ymin.full ),
                                c( xmax.full, ymax.full ),
                                c( xmin.full, ymax.full ),
                                c( xmin.full, y.thr     ),
                                c( x.thr,     y.thr     ) ),
                         # complement = 2 (Xpos missing): DN + Ypos + DP
                         rbind( c( xmin.full, ymin.full ),
                                c( x.thr,     ymin.full ),
                                c( x.thr,     y.thr     ),
                                c( xmax.full, y.thr     ),
                                c( xmax.full, ymax.full ),
                                c( xmin.full, ymax.full ) ),
                         # complement = 3 (Ypos missing): DN + Xpos + DP
                         rbind( c( xmin.full, ymin.full ),
                                c( xmax.full, ymin.full ),
                                c( xmax.full, ymax.full ),
                                c( x.thr,     ymax.full ),
                                c( x.thr,     y.thr     ),
                                c( xmin.full, y.thr     ) ),
                         # complement = 4 (DP missing): DN + Xpos + Ypos
                         rbind( c( xmin.full, ymin.full ),
                                c( xmax.full, ymin.full ),
                                c( xmax.full, y.thr     ),
                                c( x.thr,     y.thr     ),
                                c( x.thr,     ymax.full ),
                                c( xmin.full, ymax.full ) )
      )

      if ( popul.names[ 1L ] != "*" )
        boundary.list[[ popul.names[ 1L ] ]] <- l.shape

      if ( length( popul.names ) == 2L && popul.names[ 2L ] != "*" )
        boundary.list[[ popul.names[ 2L ] ]] <- quad.rect[[ complement ]]
    }

    return(boundary.list)
}
