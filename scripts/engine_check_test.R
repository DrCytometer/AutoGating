# engine_check_test.R
#
# Compares the R gating engine with the Honeychrome engine (ag_core).
#
# Before running, in the Honeychrome virtual environment:
#   python ~/dev/claude-sandboxes/AutoGating-sandbox/scripts/engine_reference.py
# Then:
#   source( "~/dev/claude-sandboxes/AutoGating-sandbox/scripts/engine_check_test.R" )
#
# Differences are on the 0-1 display scale (thresholds) or in events.

pkg.dir <- "~/dev/claude-sandboxes/AutoGating-sandbox"
suppressPackageStartupMessages( devtools::load_all( pkg.dir, quiet = TRUE ) )
ref.dir <- file.path( pkg.dir, "scripts", "engine_reference" )


# ---- 1. Threshold algorithms ------------------------------------------------

vec <- utils::read.csv( file.path( ref.dir, "vectors.csv" ), stringsAsFactors = FALSE )
thr <- utils::read.csv( file.path( ref.dir, "thresholds.csv" ), stringsAsFactors = FALSE )

r.value <- function( v, quantity ) {
    switch( quantity,
        `tail_1e-4`    = compute.threshold( v, "tail", list( tail.fraction = 1e-4 ) ),
        `tail_1e-2`    = compute.threshold( v, "tail", list( tail.fraction = 1e-2 ) ),
        kde_min        = compute.threshold( v, "kde_min" ),
        bimodal        = compute.threshold( v, "bimodal" ),
        otsu           = compute.threshold( v, "otsu" ),
        mixture2       = compute.threshold( v, "mixture" ),
        mixture3_outer = compute.threshold( v, "mixture", list( mixture.components = 3 ) ),
        mixture3_lower = compute.threshold( v, "mixture", list( mixture.components = 3, mixture.cut = "lower" ) ),
        mixture3_upper = compute.threshold( v, "mixture", list( mixture.components = 3, mixture.cut = "upper" ) ),
        kde_mode       = kde.mode( v ),
        delta_bic      = bimodality.delta.bic( v ),
        debris_cut     = debris.cut( v, list() ),
        extent_lo      = kde.extent( v )[ 1 ],
        extent_hi      = kde.extent( v )[ 2 ] )
}

cat( "1. Threshold algorithms (|R - Python|; delta_bic compared only up to 5,000 events)\n" )
thr$r <- NA_real_
for ( i in seq_len( nrow( thr ) ) ) {
    v <- vec$value[ vec$vector == thr$vector[ i ] ]
    thr$r[ i ] <- r.value( v, thr$quantity[ i ] )
}
thr$diff <- abs( thr$r - thr$value )
thr$n <- as.integer( table( vec$vector )[ thr$vector ] )
thr$diff[ thr$quantity == "delta_bic" & thr$n > 5000 ] <- NA
wide <- stats::reshape( thr[ , c( "vector", "quantity", "diff" ) ], idvar = "quantity",
                        timevar = "vector", direction = "wide" )
names( wide ) <- sub( "^diff\\.", "", names( wide ) )
print( format( wide, digits = 2 ), row.names = FALSE )


# ---- 2. Gate hierarchy -----------------------------------------------------

to.r <- function( x, parent = "" ) {
    if ( !is.list( x ) ) return( x )
    nm <- names( x )
    out <- lapply( seq_along( x ), function( i ) to.r( x[[ i ]], if ( is.null( nm ) ) "" else nm[ i ] ) )
    if ( !is.null( nm ) )
        names( out ) <- if ( parent %in% c( "populations", "stats_parent" ) ) nm else gsub( "_", ".", nm )
    out
}
gates <- to.r( jsonlite::read_json( file.path( ref.dir, "gates.json" ), simplifyVector = FALSE ) )
read.events <- function( f ) as.matrix( utils::read.csv( file.path( ref.dir, f ), check.names = FALSE ) )
data.a <- read.events( "data_a.csv" )
data.b <- read.events( "data_b.csv" )
tf <- list( CD3 = list( id = 1, scale.t = 262144, logicle.w = 0.5, logicle.m = 4.5, logicle.a = 0 ) )

el <- system.time( {
    train <- train.boundaries( gates, data.a, transforms = tf )
    applied <- run.gating( gates, data.b, reference = train$boundaries,
                           modes = list( cd3 = "recalculate", live = "recalculate" ),
                           transforms = tf )
} )[ "elapsed" ]

py <- utils::read.csv( file.path( ref.dir, "gating.csv" ), stringsAsFactors = FALSE )
cat( sprintf( "\n2. Gate hierarchy (train on sample A, apply to B; %.1f s)\n", el ) )
rows <- lapply( seq_len( nrow( py ) ), function( i ) {
    res <- if ( py$run[ i ] == "train" ) train else applied
    e <- res$boundaries[[ py$gate[ i ] ]][[ py$population[ i ] ]]
    st <- res$status[[ py$gate[ i ] ]]
    b <- if ( is.null( e ) ) matrix( NA_real_, 1, 2 ) else e$boundary
    data.frame(
        run = py$run[ i ], gate = py$gate[ i ], pop = py$population[ i ],
        n.py = py$n_events[ i ], n.r = if ( is.null( e ) ) -1L else e$n.events,
        dtx = abs( ( if ( is.null( e ) ) NA else e$threshold.x ) - py$threshold_x[ i ] ),
        dty = abs( ( if ( is.null( e ) ) NA else e$threshold.y ) - py$threshold_y[ i ] ),
        dbox = suppressWarnings( max( abs( c( min( b[ , 1 ] ), max( b[ , 1 ] ), min( b[ , 2 ] ), max( b[ , 2 ] ) ) -
                         c( py$x_lo[ i ], py$x_hi[ i ], py$y_lo[ i ], py$y_hi[ i ] ) ), na.rm = TRUE ) ),
        src.py = py$source[ i ], src.r = st$source %||% "",
        flag.py = py$flag[ i ], flag.r = if ( is.na( st$flag ) ) "" else st$flag,
        stringsAsFactors = FALSE )
} )
cmp <- do.call( rbind, rows )
cmp$dbox[ !is.finite( cmp$dbox ) ] <- NA
print( format( cmp, digits = 3 ), row.names = FALSE )


# ---- 3. Population table and anchor -----------------------------------------

pt.py <- utils::read.csv( file.path( ref.dir, "population_table.csv" ), stringsAsFactors = FALSE )
pt.r  <- population.table( gates, applied, nrow( data.b ) )
m <- merge( pt.py, pt.r, by = "key", all = TRUE, suffixes = c( ".py", ".r" ) )
cat( "\n3. Population table (sample B)\n" )
print( format( data.frame( key = m$key, n.py = m$n_events, n.r = m$n.events,
    parent.py = m$n_parent, parent.r = m$n.parent,
    stats.parent.py = m$n_stats_parent, stats.parent.r = m$n.stats.parent ), digits = 3 ),
    row.names = FALSE )
cat( sprintf( "\n   anchor: Python %s, R %s\n",
    trimws( readLines( file.path( ref.dir, "anchor.txt" ) ) ),
    anchor.population( gates ) %||% "None" ) )
