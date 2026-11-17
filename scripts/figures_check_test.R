# figures_check_test.R
#
# Draws the gate figures and statistics figures on simulated data, for a
# visual check after installing the package:
#   - an Aurora-range experiment (scatter to 4.2M, logicle fluorescence)
#     with singlets, cells, a CD3 threshold viewed against SSC, CD4 x CD8
#     quadrants, a CD69 histogram, a threshold on the y axis, a three-
#     quadrant population and a time gate;
#   - training figures from define.gates() and per-sample figures from
#     run.autogating();
#   - volcano plots, effect maps and group plots from
#     differential.statistics() on three simulated groups.
# Figures are written to out.dir and listed at the end.

library( AutoGating )

out.dir <- file.path( tempdir(), "autogating_figures" )
unlink( out.dir, recursive = TRUE )

set.seed( 1 )
n <- 60000
fsc <- c( rnorm( n * 0.85, 1.2e6, 2e5 ), runif( n * 0.15, 5e4, 3e6 ) )
fsch <- fsc * 0.8 + rnorm( n, 0, 4e4 )
doublets <- sample( n, n * 0.05 )
fsch[ doublets ] <- fsc[ doublets ] * 0.5
ssc <- abs( c( rnorm( n * 0.85, 6e5, 1.5e5 ), runif( n * 0.15, 5e4, 3e6 ) ) )
marker <- function( frac, m ) ifelse( runif( n ) < frac, 10^rnorm( n, m, 0.25 ), rnorm( n, 0, 300 ) )
raw <- cbind( `FSC-A` = fsc, `FSC-H` = fsch, `SSC-A` = ssc, Time = sort( runif( n, 0, 500 ) ),
              CD3 = marker( 0.6, 4.3 ), CD4 = marker( 0.55, 4.6 ), CD8 = marker( 0.3, 4.0 ),
              CD69 = marker( 0.1, 3.5 ) )
panel <- data.frame(
    dye = c( "FSC-A", "FSC-H", "SSC-A", "Time", "BV421-A", "FITC-A", "PE-A", "APC-A" ),
    antigen = colnames( raw ),
    type = factor( c( "scatter", "scatter", "scatter", "time", "lineage", "lineage", "lineage", "activation" ),
                   levels = c( "scatter", "time", "lineage", "activation" ) ),
    stringsAsFactors = FALSE )

ap <- get.autogating.param( "aurora", data.dir = ".", panel.filename = "panel.csv",
                            experiment.filename = "experiment.csv",
                            statistics.dir = file.path( out.dir, "statistics" ),
                            figure.dir = file.path( out.dir, "figures" ), figure.samples = TRUE )
tr <- AutoGating:::read.transformation( NULL, panel, which( panel$type == "scatter" ),
                                        which( panel$type == "lineage" ), which( panel$type == "activation" ),
                                        ap, which( panel$type == "time" ) )
ap$transform <- tr[[ 1 ]]
ap$transform.inv <- tr[[ 2 ]]
ap$transform.params <- tr[[ 3 ]]
data <- raw
for ( ch in colnames( data ) ) data[ , ch ] <- ap$transform[[ ch ]]( data[ , ch ] )
fcs.setup <- list(
    fcs.expr.data = data, fcs.sample = c( "s1", "s2" ), fcs.sample.n = 2L, fcs.sample.gate.calc = "s1",
    fcs.event.n = n, fcs.event.sample = factor( rep( c( "s1", "s2" ), each = n / 2 ), levels = c( "s1", "s2" ) ),
    fcs.panel = panel, fcs.marker = colnames( data ), fcs.marker.scatter = c( "FSC-A", "FSC-H", "SSC-A" ),
    fcs.marker.lineage = c( "CD3", "CD4", "CD8" ), fcs.marker.activation = "CD69",
    fcs.marker.label.figure = stats::setNames(
        c( "FSC-A", "FSC-H", "SSC-A", "Time", "BV421 - CD3", "FITC - CD4", "PE - CD8", "APC - CD69" ),
        colnames( data ) ),
    ap = ap )

populations <- function( type ) AutoGating:::default.populations.for.type( type )
gates <- list(
    list( gate.name = "singlets", gate.number = 1L, gate.type = "singlets", gate.marker.x = "FSC-A",
          gate.marker.y = "FSC-H", algorithm = "singlets", populations = populations( "singlets" ) ),
    list( gate.name = "cells", gate.number = 2L, gate.type = "free", gate.marker.x = "FSC-A",
          gate.marker.y = "SSC-A", parent.gate = "singlets", parent.popul = "in", algorithm = "cells",
          populations = populations( "free" ) ),
    list( gate.name = "T cells", gate.number = 3L, gate.type = "1dsep", gate.marker.x = "CD3",
          gate.marker.y = "SSC-A", sep.axis = "x", parent.gate = "cells", parent.popul = "in",
          algorithm = "mixture",
          populations = list( `CD3-` = list( label = "CD3-", region = "below" ),
                              T = list( label = "T", region = "above" ) ) ),
    list( gate.name = "CD4 CD8", gate.number = 4L, gate.type = "2dsep", gate.marker.x = "CD4",
          gate.marker.y = "CD8", parent.gate = "T cells", parent.popul = "T", algorithm = "mixture",
          populations = populations( "2dsep" ) ),
    list( gate.name = "CD69", gate.number = 5L, gate.type = "1dsep", gate.marker.x = "CD69",
          parent.gate = "T cells", parent.popul = "T", algorithm = "tail", populations = populations( "1dsep" ) ),
    list( gate.name = "CD8 on y", gate.number = 6L, gate.type = "1dsep", gate.marker.x = "CD8",
          gate.marker.y = "CD4", sep.axis = "y", parent.gate = "T cells", parent.popul = "T",
          algorithm = "mixture", populations = populations( "1dsep" ) ),
    list( gate.name = "L gate", gate.number = 7L, gate.type = "2dsep", gate.marker.x = "CD4",
          gate.marker.y = "CD8", parent.gate = "T cells", parent.popul = "T", algorithm = "mixture",
          populations = list( `not DP` = list( label = "not DP", region = "x-y-|x+y-|x-y+" ),
                              DP = list( label = "DP", region = "x+y+" ) ) ),
    list( gate.name = "time", gate.number = 8L, gate.type = "free", gate.marker.x = "Time",
          gate.marker.y = "SSC-A", parent.gate = "cells", parent.popul = "in", algorithm = "all",
          populations = populations( "free" ) ) )

t0 <- Sys.time()
model <- define.gates( fcs.setup, ap, new.agmodel( gates ) )
before <- .Random.seed
gating <- run.autogating( fcs.setup, model, ap )
cat( sprintf( "Gating and figures: %.1f s; random number stream unchanged: %s\n",
              as.numeric( difftime( Sys.time(), t0, units = "secs" ) ), identical( before, .Random.seed ) ) )

# Statistics figures: 18 simulated samples in three groups, 20 populations.
samples <- sprintf( "s%02d", 1:18 )
group <- rep( c( "ctrl", "drug A", "drug B" ), each = 6 )
keys <- sprintf( "g%d/p%d", rep( 1:4, each = 5 ), 1:5 )
labels <- sprintf( "Pop %d.%d", rep( 1:4, each = 5 ), 1:5 )
effect.a <- c( 2, 1, rep( 0, 18 ) )
effect.b <- c( 0, 1.6, 0, 0, -1.5, rep( 0, 15 ) )
population.rows <- do.call( rbind, lapply( seq_along( samples ), function( i ) {
    parent <- round( runif( 1, 3000, 6000 ) )
    shift <- if ( group[ i ] == "drug A" ) 2^effect.a else if ( group[ i ] == "drug B" ) 2^effect.b else 1
    p <- pmin( 0.9, seq( 0.05, 0.3, length.out = 20 ) * shift * exp( rnorm( 20, 0, 0.1 ) ) )
    data.frame( sample = samples[ i ], key = keys, gate = sub( "/.*", "", keys ),
                population = sub( ".*/", "", keys ), label = labels, n.events = rbinom( 20, parent, p ),
                parent.key = "root", n.parent = parent, n.stats.parent = NA_real_,
                stats.parent.missing = FALSE, stringsAsFactors = FALSE )
} ) )
median.rows <- do.call( rbind, lapply( seq_along( samples ), function( i ) do.call( rbind, lapply( keys, function( k )
    data.frame( sample = samples[ i ], key = k, channel = c( "CD25", "CD69", "PD1" ),
                median = 0.4 + 0.02 * rnorm( 3 ) + ( group[ i ] == "drug A" & k == "g1/p1" ) * c( 0.1, 0, 0 ),
                n.events = 500, stringsAsFactors = FALSE ) ) ) ) )
stats <- differential.statistics(
    list( population.table = population.rows, medians = median.rows, model = list( gate.definitions = list() ) ),
    data.frame( sample = samples, group = group ), output.dir = file.path( out.dir, "differential" ) )

print( list.files( out.dir, pattern = "png$", recursive = TRUE ) )
cat( "Figures in", out.dir, "\n" )
if ( Sys.info()[[ "sysname" ]] == "Darwin" ) system2( "open", shQuote( out.dir ) )
