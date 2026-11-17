# mixture_impact_test.R
#
# Where and how much the threshold calculation differs between AutoGating and
# Honeychrome on the same inputs (UCL IMMPROVE IC controls, IMMPROVE FlowJo
# template, hc_harness scenario ag_r_match).
#
# Two models with the same gate definitions:
#   hc  Honeychrome's trained model
#   r   the same gates trained here on Honeychrome's training events
# Both are applied here to every event of both samples, so only the trained
# thresholds differ (the gating engines agree on Honeychrome's model to 2 of
# 3,244 population counts).
#
# Each gate is also recalculated here on exactly the parent events
# Honeychrome's model gives it (recalculate mode, Honeychrome's boundary kept
# for the gates below, drift = the threshold change). Each population is
# classed by why it can differ:
#   own threshold   its gate's threshold differs by more than 0.001 (display
#                   units) on identical parent events
#   inherited       its own does not, but a gate above it does
#   unchanged       neither (it can still move a little when its parent does)
#
# Outputs in file.path( hc.dir, "r", "impact" ):
#   populations.csv  every population: counts, % of parent, differences,
#                    class, parent size, depth, gate algorithm and marker
#   gates.csv        every threshold gate: threshold in each model, change,
#                    parent size on the training events
#   by_size.csv      summary by parent size and class
#   f1_parent_pct.png, f2_live_pct.png, f3_thresholds.png,
#   f4_largest.png, f5_hladr.png, f6_inherited.png
#   summary.txt
#
# Run after ag_r_match_test.R (uses its panel, experiment and transformation
# files in file.path( hc.dir, "r" )).

library( AutoGating )
library( ggplot2 )

hc.dir  <- "~/dev/claude-sandboxes/hc_harness/results/ag_r_match/20261117T102742Z"
fcs.dir <- "~/Experiments/UCL IMMPROVE/Unmixed/IC controls"

r.dir   <- file.path( hc.dir, "r" )
out.dir <- file.path( r.dir, "impact" )
dir.create( out.dir, showWarnings = FALSE, recursive = TRUE )
summary.lines <- character( 0 )
say <- function( ... ) {
    line <- sprintf( ... )
    summary.lines <<- c( summary.lines, line )
    cat( line, "\n" )
}
class.colours <- c( "own threshold" = "#E41A1C", "inherited" = "#377EB8", "unchanged" = "grey60" )
size.breaks <- c( 0, 100, 1000, 10000, Inf )
size.labels <- c( "< 100", "100-1,000", "1,000-10,000", ">= 10,000" )


# ---- Data and models ----------------------------------------------------------

ap <- get.autogating.param(
    cytometer               = "aurora",
    data.dir                = fcs.dir,
    panel.filename          = file.path( r.dir, "panel.csv" ),
    experiment.filename     = file.path( r.dir, "experiment.csv" ),
    transformation.filename = file.path( r.dir, "transformation.csv" ),
    figure.training         = FALSE )
fcs.setup <- setup.autogating( ap )
transforms <- fcs.setup$ap$transform.params
panel <- fcs.setup$fcs.panel

events <- jsonlite::read_json( file.path( hc.dir, "events.json" ), simplifyVector = TRUE )
training <- matrix( readBin( file.path( hc.dir, "training_events.f32" ), "numeric",
                             n = prod( events$training_shape ), size = 4L ),
                    nrow = events$training_shape[ 1 ], byrow = TRUE )
colnames( training ) <- events$channels
training <- training[ , panel$dye, drop = FALSE ]
colnames( training ) <- panel$antigen
training <- training[ , colnames( fcs.setup$fcs.expr.data ), drop = FALSE ]

model.hc <- read.agmodel( file.path( hc.dir, "model.agmodel" ), panel )
defs <- model.hc$gate.definitions
names( defs ) <- vapply( defs, `[[`, character( 1 ), "gate.name" )
boundaries <- list(
    hc = model.hc$trained.boundaries,
    r  = train.boundaries( defs, training, transforms = transforms )$boundaries )
isolated <- run.gating( defs, training, reference = boundaries$hc,
                        modes = stats::setNames( as.list( rep( "recalculate", length( defs ) ) ), names( defs ) ),
                        drift.limit = 0, drift.action = "reference", transforms = transforms )
training.parent <- isolated$parent.events
isolated.change <- vapply( names( defs ), function( g ) {
    v <- isolated$status[[ g ]]$drift
    if ( is.null( v ) ) NA_real_ else as.numeric( v )
}, numeric( 1 ) )


# ---- Threshold of every gate ----------------------------------------------------

threshold.of <- function( b, axis ) {
    if ( length( b ) == 0L ) return( NA_real_ )
    v <- b[[ 1 ]][[ paste0( "threshold.", axis ) ]]
    if ( is.null( v ) ) NA_real_ else as.numeric( v )
}
axis.algorithm <- function( g, axis ) {
    a <- g$gate.param[[ paste0( "algorithm.", axis ) ]]
    if ( is.null( a ) ) g$algorithm else a
}
gates <- do.call( rbind, lapply( defs, function( g ) {
    do.call( rbind, lapply( c( "x", "y" ), function( axis ) {
        hc <- threshold.of( boundaries$hc[[ g$gate.name ]], axis )
        if ( is.na( hc ) ) return( NULL )
        data.frame( gate = g$gate.name, axis = axis, gate.type = g$gate.type,
                    algorithm = axis.algorithm( g, axis ),
                    channel = if ( axis == "x" ) g$gate.marker.x else g$gate.marker.y,
                    n.parent.training = length( training.parent[[ g$gate.name ]] ),
                    hc = hc, r = threshold.of( boundaries$r[[ g$gate.name ]], axis ),
                    stringsAsFactors = FALSE )
    } ) )
} ) )
gates$change <- gates$r - gates$hc
gates$isolated.change <- isolated.change[ gates$gate ]
utils::write.csv( gates, file.path( out.dir, "gates.csv" ), row.names = FALSE )
moved <- names( isolated.change )[ !is.na( isolated.change ) & isolated.change > 1e-3 ]
say( "Threshold gates: %d axes on %d gates", nrow( gates ), length( unique( gates$gate ) ) )
for ( a in sort( unique( gates$algorithm ) ) ) {
    s <- gates$algorithm == a
    g <- unique( gates$gate[ s ] )
    say( "   %-9s %4d axes; trained model: moved > 0.001 in %d, > 0.01 in %d; on identical parents: %d of %d gates differ by > 0.001",
         a, sum( s ), sum( abs( gates$change[ s ] ) > 1e-3, na.rm = TRUE ),
         sum( abs( gates$change[ s ] ) > 1e-2, na.rm = TRUE ), sum( g %in% moved ), length( g ) )
}


# ---- Every population in every sample -----------------------------------------------

ancestors <- function( name ) {
    out <- character( 0 )
    p <- defs[[ name ]]$parent.gate
    while ( !is.null( p ) && length( p ) && nzchar( p ) && p %in% names( defs ) ) {
        out <- c( out, p )
        p <- defs[[ p ]]$parent.gate
    }
    out
}
anc <- lapply( names( defs ), ancestors )
names( anc ) <- names( defs )

gated <- list()
pops <- do.call( rbind, lapply( fcs.setup$fcs.sample, function( s ) {
    rows <- which( fcs.setup$fcs.event.sample == s )
    d <- fcs.setup$fcs.expr.data[ rows, , drop = FALSE ]
    tabs <- lapply( c( "hc", "r" ), function( m ) {
        res <- run.gating( defs, d, reference = boundaries[[ m ]], transforms = transforms )
        gated[[ paste( s, m ) ]] <<- res
        population.table( defs, res, nrow( d ) )
    } )
    t <- merge( tabs[[ 1 ]][ , c( "key", "gate", "population", "n.events", "n.parent", "fraction.of.parent" ) ],
                tabs[[ 2 ]][ , c( "key", "n.events", "n.parent", "fraction.of.parent" ) ],
                by = "key", suffixes = c( ".hc", ".r" ) )
    live <- tabs[[ 1 ]]$n.events[ tabs[[ 1 ]]$key == "FSC-A, viability subset/FSC-A, viability subset" ]
    cbind( sample = s, t, live.hc = live, stringsAsFactors = FALSE )
} ) )
pops$own <- pops$gate %in% moved
pops$inherited <- vapply( pops$gate, function( g ) any( anc[[ g ]] %in% moved ), logical( 1 ) )
pops$class <- ifelse( pops$own, "own threshold", ifelse( pops$inherited, "inherited", "unchanged" ) )
pops$depth <- lengths( anc[ pops$gate ] )
pops$algorithm <- vapply( pops$gate, function( g ) axis.algorithm( defs[[ g ]], "x" ), character( 1 ) )
pops$marker <- vapply( pops$gate, function( g ) defs[[ g ]]$gate.marker.x, character( 1 ) )
pops$d.pct.parent <- 100 * ( pops$fraction.of.parent.r - pops$fraction.of.parent.hc )
pops$d.count <- pops$n.events.r - pops$n.events.hc
pops$d.pct.live <- 100 * pops$d.count / pops$live.hc
pops$d.count.relative <- pops$d.count / pmax( pops$n.events.hc, 1 )
pops$parent.size <- cut( pops$n.parent.hc, size.breaks, size.labels, right = FALSE )
pops$class <- factor( pops$class, names( class.colours ) )
utils::write.csv( pops, file.path( out.dir, "populations.csv" ), row.names = FALSE )


# ---- Summary by parent size -------------------------------------------------------

by.size <- do.call( rbind, lapply( split( pops, list( pops$parent.size, pops$class ), drop = TRUE ), function( p )
    data.frame( parent.size = p$parent.size[ 1 ], class = p$class[ 1 ], populations = nrow( p ),
                median.abs.d.pct.parent = stats::median( abs( p$d.pct.parent ), na.rm = TRUE ),
                over.1pp.parent = sum( abs( p$d.pct.parent ) > 1, na.rm = TRUE ),
                over.5pp.parent = sum( abs( p$d.pct.parent ) > 5, na.rm = TRUE ),
                max.abs.d.pct.parent = max( abs( p$d.pct.parent ), na.rm = TRUE ),
                over.0.1pct.live = sum( abs( p$d.pct.live ) > 0.1, na.rm = TRUE ),
                max.abs.d.pct.live = max( abs( p$d.pct.live ), na.rm = TRUE ) ) ) )
by.size <- by.size[ order( by.size$parent.size, by.size$class ), ]
utils::write.csv( by.size, file.path( out.dir, "by_size.csv" ), row.names = FALSE )

say( "Populations: %d (both samples); own threshold %d, inherited %d, unchanged %d",
     nrow( pops ), sum( pops$class == "own threshold" ), sum( pops$class == "inherited" ),
     sum( pops$class == "unchanged" ) )
big.pp <- abs( pops$d.pct.parent ) > 5 & !is.na( pops$d.pct.parent )
say( "Populations whose %% of parent moves by > 5 points: %d; parent < 1,000 events in %d; own threshold %d, inherited only %d",
     sum( big.pp ), sum( big.pp & pops$n.parent.hc < 1000 ),
     sum( big.pp & pops$class == "own threshold" ), sum( big.pp & pops$class == "inherited" ) )
big.live <- abs( pops$d.pct.live ) > 0.1 & !is.na( pops$d.pct.live )
say( "Populations whose count moves by > 0.1%% of live cells: %d; parent >= 10,000 events in %d",
     sum( big.live ), sum( big.live & pops$n.parent.hc >= 10000 ) )
say( "Largest count change: %.2f%% of live cells; median over all populations %.4f%%",
     max( abs( pops$d.pct.live ), na.rm = TRUE ), stats::median( abs( pops$d.pct.live ), na.rm = TRUE ) )
say( "Spearman correlation of |change in %% of parent| with parent size: %.2f; of |change as %% of live cells| with parent size: %.2f",
     stats::cor( abs( pops$d.pct.parent ), pops$n.parent.hc, method = "spearman", use = "complete.obs" ),
     stats::cor( abs( pops$d.pct.live ), pops$n.parent.hc, method = "spearman", use = "complete.obs" ) )


# ---- Figures ---------------------------------------------------------------------------

theme.set <- theme_bw( base_size = 10 ) + theme( panel.grid.minor = element_blank(), legend.position = "bottom" )
pops.plot <- pops[ pops$n.parent.hc > 0 & !is.na( pops$d.pct.parent ), ]

f1 <- ggplot( pops.plot, aes( n.parent.hc, abs( d.pct.parent ), colour = class ) ) +
    geom_point( size = 0.9, alpha = 0.7 ) +
    scale_x_log10( labels = scales::comma ) +
    scale_colour_manual( values = class.colours, drop = FALSE ) +
    geom_hline( yintercept = c( 1, 5 ), linetype = "dashed", colour = "grey40" ) +
    facet_wrap( ~ sample ) +
    labs( x = "Parent population (events in the sample)", y = "|R - Honeychrome|, % of parent",
          colour = NULL, title = "Change in % of parent against parent size" ) + theme.set
ggsave( file.path( out.dir, "f1_parent_pct.png" ), f1, width = 9, height = 4.8, dpi = 150 )

f2 <- ggplot( pops.plot, aes( n.parent.hc, pmax( abs( d.pct.live ), 1e-4 ), colour = class ) ) +
    geom_point( size = 0.9, alpha = 0.7 ) +
    scale_x_log10( labels = scales::comma ) + scale_y_log10() +
    scale_colour_manual( values = class.colours, drop = FALSE ) +
    geom_hline( yintercept = 0.1, linetype = "dashed", colour = "grey40" ) +
    facet_wrap( ~ sample ) +
    labs( x = "Parent population (events in the sample)",
          y = "|R - Honeychrome| count, % of live cells (floor 1e-4)",
          colour = NULL, title = "Change in cell count against parent size" ) + theme.set
ggsave( file.path( out.dir, "f2_live_pct.png" ), f2, width = 9, height = 4.8, dpi = 150 )

g.plot <- gates[ !duplicated( gates$gate ) & !is.na( gates$isolated.change ), ]
g.plot$label <- ifelse( g.plot$isolated.change > 0.02 & g.plot$n.parent.training >= 2000, g.plot$gate, NA )
f3 <- ggplot( g.plot, aes( n.parent.training, isolated.change + 1e-7, colour = algorithm ) ) +
    geom_point( size = 0.9, alpha = 0.7 ) +
    geom_text( aes( label = label ), size = 2.3, hjust = -0.05, vjust = 0.5, na.rm = TRUE, show.legend = FALSE ) +
    scale_x_log10( labels = scales::comma ) + scale_y_log10() +
    scale_colour_brewer( palette = "Dark2" ) +
    geom_hline( yintercept = 1e-3, linetype = "dashed", colour = "grey40" ) +
    labs( x = "Parent population (training events)",
          y = "|R - Honeychrome threshold| on identical parents (display units, + 1e-7)",
          colour = NULL, title = "Threshold difference of each gate, calculated on identical parent events" ) + theme.set
ggsave( file.path( out.dir, "f3_thresholds.png" ), f3, width = 8, height = 5.5, dpi = 150 )

# The gates whose own threshold moves the most cells, drawn on their parent
# in the sample where they move most: Honeychrome's threshold solid, R's dashed.
# Axes are named by fluorophore: the exported panel names only the first
# marker of a stacked fluorophore.
dye.of <- stats::setNames( panel$dye, panel$antigen )
own <- pops[ pops$class == "own threshold" & !is.na( pops$d.pct.live ), ]
own <- own[ order( -abs( own$d.pct.live ) ), ]
own <- own[ !duplicated( own$gate ), ][ seq_len( min( 8L, sum( !duplicated( own$gate ) ) ) ), ]
panels <- lapply( seq_len( nrow( own ) ), function( i ) {
    g <- defs[[ own$gate[ i ] ]]
    s <- own$sample[ i ]
    rows <- which( fcs.setup$fcs.event.sample == s )
    parent <- gated[[ paste( s, "hc" ) ]]$parent.events[[ g$gate.name ]]
    d <- fcs.setup$fcs.expr.data[ rows, , drop = FALSE ][ parent, , drop = FALSE ]
    tx <- c( hc = threshold.of( boundaries$hc[[ g$gate.name ]], "x" ), r = threshold.of( boundaries$r[[ g$gate.name ]], "x" ) )
    ty <- c( hc = threshold.of( boundaries$hc[[ g$gate.name ]], "y" ), r = threshold.of( boundaries$r[[ g$gate.name ]], "y" ) )
    title <- sprintf( "%s (%s)\n%s: parent %s; %s %.1f%% vs %.1f%%; %+.2f%% of live",
                      g$gate.name, axis.algorithm( g, "x" ), s, format( length( parent ), big.mark = "," ),
                      own$population[ i ], 100 * own$fraction.of.parent.r[ i ], 100 * own$fraction.of.parent.hc[ i ],
                      own$d.pct.live[ i ] )
    title <- paste( strwrap( title, 60 ), collapse = "\n" )
    if ( identical( g$gate.type, "2dsep" ) ) {
        p <- ggplot( data.frame( x = d[ , g$gate.marker.x ], y = d[ , g$gate.marker.y ] ), aes( x, y ) ) +
            geom_point( size = 0.2, alpha = 0.3, colour = "grey30" ) +
            geom_vline( xintercept = tx, linetype = c( "solid", "dashed" ), colour = c( "black", "#E41A1C" ) ) +
            geom_hline( yintercept = ty[ !is.na( ty ) ], linetype = c( "solid", "dashed" )[ !is.na( ty ) ],
                        colour = c( "black", "#E41A1C" )[ !is.na( ty ) ] ) +
            labs( x = dye.of[[ g$gate.marker.x ]], y = dye.of[[ g$gate.marker.y ]] )
    } else {
        p <- ggplot( data.frame( x = d[ , g$gate.marker.x ] ), aes( x ) ) +
            geom_histogram( bins = 128, fill = "grey55" ) +
            geom_vline( xintercept = tx, linetype = c( "solid", "dashed" ), colour = c( "black", "#E41A1C" ) ) +
            labs( x = dye.of[[ g$gate.marker.x ]], y = "Events" )
    }
    p + labs( title = title ) + theme_bw( base_size = 8 ) + theme( plot.title = element_text( size = 7 ) )
} )
ggsave( file.path( out.dir, "f4_largest.png" ), patchwork::wrap_plots( panels, ncol = 4 ),
        width = 14, height = 7.5, dpi = 150 )

# The populations that change most without their own threshold differing:
# the parent as gated by each model (Honeychrome filled grey, R red outline)
# and each model's threshold (Honeychrome solid, R dashed).
inh <- pops[ pops$class == "inherited" & !is.na( pops$d.pct.parent ) & pops$n.parent.hc >= 100, ]
inh <- inh[ order( -abs( inh$d.pct.parent ) ), ]
inh <- inh[ !duplicated( inh$gate ), ][ seq_len( min( 4L, sum( !duplicated( inh$gate ) ) ) ), ]
panels <- lapply( seq_len( nrow( inh ) ), function( i ) {
    g <- defs[[ inh$gate[ i ] ]]
    s <- inh$sample[ i ]
    rows <- which( fcs.setup$fcs.event.sample == s )
    d <- fcs.setup$fcs.expr.data[ rows, g$gate.marker.x ]
    both <- rbind( data.frame( model = "Honeychrome", x = d[ gated[[ paste( s, "hc" ) ]]$parent.events[[ g$gate.name ]] ] ),
                   data.frame( model = "R", x = d[ gated[[ paste( s, "r" ) ]]$parent.events[[ g$gate.name ]] ] ) )
    tx <- c( threshold.of( boundaries$hc[[ g$gate.name ]], "x" ), threshold.of( boundaries$r[[ g$gate.name ]], "x" ) )
    title <- sprintf( "%s (%s, %s)\n%s: parent %s vs %s; %s %.1f%% vs %.1f%%",
                      g$gate.name, axis.algorithm( g, "x" ), g$gate.type, s,
                      format( inh$n.parent.r[ i ], big.mark = "," ), format( inh$n.parent.hc[ i ], big.mark = "," ),
                      inh$population[ i ], 100 * inh$fraction.of.parent.r[ i ], 100 * inh$fraction.of.parent.hc[ i ] )
    ggplot( both, aes( x ) ) +
        geom_histogram( data = both[ both$model == "Honeychrome", ], bins = 64, fill = "grey70" ) +
        geom_histogram( data = both[ both$model == "R", ], bins = 64, fill = NA, colour = "#E41A1C" ) +
        geom_vline( xintercept = tx, linetype = c( "solid", "dashed" ), colour = c( "black", "#E41A1C" ) ) +
        labs( x = dye.of[[ g$gate.marker.x ]], y = "Events", title = paste( strwrap( title, 60 ), collapse = "\n" ) ) +
        theme_bw( base_size = 8 ) + theme( plot.title = element_text( size = 7 ) )
} )
ggsave( file.path( out.dir, "f6_inherited.png" ), patchwork::wrap_plots( panels, ncol = 4 ),
        width = 14, height = 4, dpi = 150 )

hla <- pops.plot[ grepl( "^HLADR", pops.plot$gate ), ]
f5 <- ggplot( hla, aes( n.parent.hc, d.pct.parent, colour = class ) ) +
    geom_point( size = 1.4 ) +
    scale_x_log10( labels = scales::comma ) +
    scale_colour_manual( values = class.colours, drop = FALSE ) +
    geom_hline( yintercept = 0, colour = "grey40" ) +
    facet_wrap( ~ sample ) +
    labs( x = "Parent population (events in the sample)", y = "R - Honeychrome, % of parent",
          colour = NULL, title = "HLA-DR+ gates" ) + theme.set
ggsave( file.path( out.dir, "f5_hladr.png" ), f5, width = 9, height = 4.5, dpi = 150 )
say( "HLA-DR+ gates: %d populations; |change in %% of parent| > 1 point in %d (parent >= 10,000 events in %d)",
     nrow( hla ), sum( abs( hla$d.pct.parent ) > 1 ), sum( abs( hla$d.pct.parent ) > 1 & hla$n.parent.hc >= 10000 ) )

writeLines( summary.lines, file.path( out.dir, "summary.txt" ) )
