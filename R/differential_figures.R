# differential_figures.R
#
# Figures of the differential statistics (differential.statistics()):
#
#   volcano plots  one point per population (frequencies) or per population
#                  and marker (medians) and comparison: the difference
#                  between the groups across, the strength of the evidence
#                  (-log10 p-value) up.  Significant results are coloured by
#                  direction; dashed lines mark the effect-size threshold.
#   effect maps    the differences of every result significant in at least
#                  one comparison, one column per comparison (red higher,
#                  blue lower in the second group), starred where
#                  significant.
#   group plots    each sample's value for the most significant results,
#                  by group, with the group median as a bar: the
#                  frequencies as a percentage of the denominator, the
#                  marker medians on the display scale.
#
# Points within a group are spread sideways by rank rather than at random,
# so a figure drawn twice is the same.


#' Labels of statistics results for figures
#'
#' @title Result Labels
#'
#' @param tests Result table (feature, label).
#'
#' @return Character vector: the label, else the feature; labels shared by
#'   several features are followed by the feature.
#'
#' @keywords internal
result.labels <- function( tests ) {
    lab <- if ( "label" %in% names( tests ) ) as.character( tests$label ) else tests$feature
    lab[ is.na( lab ) | !nzchar( lab ) ] <- tests$feature[ is.na( lab ) | !nzchar( lab ) ]
    if ( "channel" %in% names( tests ) )
        lab <- paste( lab, tests$channel, sep = " | " )
    by.label <- tapply( tests$feature, lab, function( f ) length( unique( f ) ) )
    shared <- lab %in% names( by.label )[ by.label > 1L ]
    lab[ shared ] <- sprintf( "%s (%s)", lab[ shared ], tests$feature[ shared ] )
    lab
}


#' Volcano plot of statistics results
#'
#' @title Volcano Figure
#'
#' @param tests Result table (feature, comparison, logFC, P.Value,
#'   significant, threshold).
#' @param x.title Title of the difference axis.
#' @param n.labels Significant results labelled per comparison (smallest
#'   p-values first).
#'
#' @return A \pkg{ggplot2} object, or \code{NULL} without results.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_text
#' @importFrom ggplot2 scale_colour_manual scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 expansion facet_wrap vars labs theme_bw theme
#' @importFrom ggplot2 element_blank
#'
#' @keywords internal
volcano.figure <- function( tests, x.title, n.labels = 8L ) {
    tests <- tests[ is.finite( tests$logFC ) & is.finite( tests$P.Value ), , drop = FALSE ]
    if ( nrow( tests ) == 0L )
        return( NULL )
    tests$display <- result.labels( tests )
    tests$neglog10p <- -log10( pmax( tests$P.Value, 1e-300 ) )
    tests$significance <- factor(
        ifelse( !tests$significant, "not significant", ifelse( tests$logFC > 0, "higher", "lower" ) ),
        levels = c( "higher", "lower", "not significant" ) )
    sig <- tests[ tests$significant, , drop = FALSE ]
    sig <- sig[ order( sig$comparison, sig$P.Value ), , drop = FALSE ]
    sig <- do.call( rbind, lapply( split( sig, sig$comparison ), utils::head, n.labels ) )
    threshold <- unique( tests$threshold[ is.finite( tests$threshold ) ] )
    fig <- ggplot( tests, aes( x = logFC, y = neglog10p, colour = significance ) ) +
        geom_point( size = 1.6, alpha = 0.8 ) +
        scale_colour_manual( values = c( higher = "#C0392B", lower = "#2E6FB7",
                                         `not significant` = "grey65" ),
                             name = NULL, drop = FALSE ) +
        facet_wrap( vars( comparison ) ) +
        scale_x_continuous( expand = expansion( mult = 0.25 ) ) +
        scale_y_continuous( expand = expansion( mult = c( 0.03, 0.12 ) ) ) +
        labs( x = x.title, y = "-log10 p-value" ) +
        theme_bw() + theme( panel.grid.minor = element_blank() )
    if ( length( threshold ) == 1L && threshold > 0 )
        fig <- fig + geom_vline( xintercept = c( -threshold, threshold ), linetype = "dashed",
                                 colour = "grey40" )
    if ( !is.null( sig ) && nrow( sig ) )
        fig <- fig + geom_text( data = sig, aes( label = display ), size = 2.8, vjust = -0.7,
                                check_overlap = TRUE, show.legend = FALSE )
    fig
}


#' Effect map of statistics results
#'
#' @title Effect Map
#'
#' @param tests Result table (feature, comparison, logFC, P.Value,
#'   significant).
#' @param fill.title Title of the colour scale.
#' @param max.rows Most results shown (smallest p-values first).
#'
#' @return A \pkg{ggplot2} object, or \code{NULL} when nothing is
#'   significant.
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient2 labs
#' @importFrom ggplot2 theme_bw theme element_blank element_text
#'
#' @keywords internal
effect.map <- function( tests, fill.title, max.rows = 60L ) {
    tests$display <- result.labels( tests )
    hit <- unique( tests$feature[ tests$significant ] )
    if ( length( hit ) == 0L )
        return( NULL )
    tests <- tests[ tests$feature %in% hit, , drop = FALSE ]
    best <- tapply( tests$P.Value, tests$feature, min, na.rm = TRUE )
    keep <- names( sort( best ) )[ seq_len( min( max.rows, length( best ) ) ) ]
    tests <- tests[ tests$feature %in% keep, , drop = FALSE ]
    order.display <- unique( tests$display[ match( keep, tests$feature ) ] )
    tests$display <- factor( tests$display, levels = rev( order.display ) )
    tests$comparison <- factor( tests$comparison, levels = unique( tests$comparison ) )
    tests$star <- ifelse( tests$significant, "*", "" )
    lim <- max( abs( tests$logFC ), na.rm = TRUE )
    ggplot( tests, aes( x = comparison, y = display, fill = logFC ) ) +
        geom_tile( colour = "white" ) +
        geom_text( aes( label = star ), size = 4, vjust = 0.75 ) +
        scale_fill_gradient2( low = "#2E6FB7", mid = "white", high = "#C0392B", midpoint = 0,
                              limits = c( -lim, lim ), name = fill.title ) +
        labs( x = NULL, y = NULL ) +
        theme_bw() +
        theme( panel.grid = element_blank(), axis.text.x = element_text( angle = 30, hjust = 1 ) )
}


#' Values of the most significant results by group
#'
#' @title Group Figure
#'
#' @param values Matrix samples x features.
#' @param groups Named character vector sample -> group.
#' @param group.levels Group order.
#' @param tests Result table (feature, P.Value, significant).
#' @param y.title Title of the value axis.
#' @param max.panels Most results shown (smallest p-values first).
#'
#' @return A \pkg{ggplot2} object, or \code{NULL} when nothing is
#'   significant.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar facet_wrap vars
#' @importFrom ggplot2 scale_x_continuous labs theme_bw theme element_blank
#' @importFrom ggplot2 element_text
#'
#' @keywords internal
group.figure <- function( values, groups, group.levels, tests, y.title, max.panels = 24L ) {
    tests <- tests[ tests$feature %in% colnames( values ), , drop = FALSE ]
    tests$display <- result.labels( tests )
    sig <- tests[ tests$significant, , drop = FALSE ]
    if ( nrow( sig ) == 0L )
        return( NULL )
    best <- tapply( sig$P.Value, sig$feature, min, na.rm = TRUE )
    keep <- names( sort( best ) )[ seq_len( min( max.panels, length( best ) ) ) ]
    display <- tests$display[ match( keep, tests$feature ) ]
    samples <- intersect( rownames( values ), names( groups ) )
    long <- do.call( rbind, lapply( seq_along( keep ), function( i ) data.frame(
        display = display[ i ], sample = samples, group = groups[ samples ],
        value = values[ samples, keep[ i ] ], stringsAsFactors = FALSE ) ) )
    long <- long[ is.finite( long$value ), , drop = FALSE ]
    long$display <- factor( long$display, levels = display )
    long$group <- factor( long$group, levels = group.levels )
    cell <- interaction( long$display, long$group, drop = TRUE )
    rank.in <- stats::ave( long$value, cell, FUN = function( v ) rank( v, ties.method = "first" ) )
    size.in <- stats::ave( long$value, cell, FUN = length )
    long$group.position <- as.integer( long$group ) +
        ifelse( size.in > 1, 0.5 * ( ( rank.in - 1 ) / ( size.in - 1 ) - 0.5 ), 0 )
    med <- stats::aggregate( value ~ display + group, data = long, FUN = stats::median )
    med$group.position <- as.integer( med$group )
    ggplot( long, aes( x = group.position, y = value ) ) +
        geom_errorbar( data = med, aes( ymin = value, ymax = value ), width = 0.6,
                       colour = "grey30", linewidth = 0.7 ) +
        geom_point( aes( colour = group ), size = 1.8, show.legend = FALSE ) +
        facet_wrap( vars( display ), scales = "free_y" ) +
        scale_x_continuous( breaks = seq_along( group.levels ), labels = group.levels ) +
        labs( x = NULL, y = y.title ) +
        theme_bw() +
        theme( panel.grid = element_blank(), strip.text = element_text( size = 8 ),
               axis.text.x = element_text( angle = 30, hjust = 1 ) )
}


#' Write the figures of differential statistics
#'
#' @title Write Differential Figures
#'
#' @description
#' Volcano plots, effect maps and group plots of the frequency and marker
#' tests (see the file header of \code{differential_figures.R}), as PNG
#' files written with \pkg{ragg}: \code{frequency_volcano.png},
#' \code{frequency_effects.png}, \code{frequency_groups.png} and the same
#' for \code{marker_}.  Effect maps and group plots are written only when
#' something is significant.
#'
#' @param stats Result of \code{\link{differential.statistics}}.
#' @param dir Output directory (created when missing).
#' @param resolution Dots per inch.
#'
#' @return Character vector of the files written, invisibly.
#'
#' @export
write.differential.figures <- function( stats, dir, resolution = 150 ) {
    dir.create( dir, recursive = TRUE, showWarnings = FALSE )
    files <- character( 0 )
    save <- function( fig, name, width, height ) {
        if ( is.null( fig ) ) return( invisible() )
        file <- file.path( dir, name )
        ggplot2::ggsave( file, fig, device = ragg::agg_png, width = width, height = height,
                         units = "in", dpi = resolution, limitsize = FALSE )
        files <<- c( files, file )
    }
    n.comp <- length( stats$comparisons )
    levels <- stats$group.levels %||% unique( stats$groups )
    panels <- function( n ) ceiling( sqrt( n ) )
    sets <- list(
        list( prefix = "frequency", tests = stats$frequencies, values = stats$frequency.values,
              x.title = "Difference (log2 odds ratio)", y.title = "% of parent",
              fill.title = "log2 odds\nratio" ),
        list( prefix = "marker", tests = stats$markers, values = stats$marker.values,
              x.title = "Difference in median (display scale)", y.title = "Median (display scale)",
              fill.title = "Difference" ) )
    for ( s in sets ) {
        if ( is.null( s$tests ) || nrow( s$tests ) == 0L ) next
        k <- panels( n.comp )
        save( volcano.figure( s$tests, s$x.title ), paste0( s$prefix, "_volcano.png" ),
              3.5 * k + 1.5, 3.2 * ceiling( n.comp / k ) + 0.5 )
        n.rows <- min( 60L, length( unique( s$tests$feature[ s$tests$significant ] ) ) )
        save( effect.map( s$tests, s$fill.title ), paste0( s$prefix, "_effects.png" ),
              4 + 0.6 * n.comp, 1.5 + 0.2 * n.rows )
        n.panels <- min( 24L, length( unique( s$tests$feature[ s$tests$significant ] ) ) )
        k <- panels( n.panels )
        save( group.figure( s$values, stats$groups, levels, s$tests, s$y.title ),
              paste0( s$prefix, "_groups.png" ), max( 3, 2.2 * k + 0.6 ),
              max( 3, 2 * ceiling( n.panels / max( k, 1 ) ) + 0.6 ) )
    }
    invisible( files )
}
