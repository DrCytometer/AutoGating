# agmodel_check_test.R
#
# Round trip of gating models between Honeychrome (ag_core snapshot) and
# AutoGating, on the synthetic events of engine_reference.py.
#
# In the Honeychrome virtual environment, first:
#   python ~/dev/claude-sandboxes/AutoGating-sandbox/scripts/agmodel_reference.py make
# Then:
#   source( "~/dev/claude-sandboxes/AutoGating-sandbox/scripts/agmodel_check_test.R" )
# Then, in the Honeychrome virtual environment:
#   python ~/dev/claude-sandboxes/AutoGating-sandbox/scripts/agmodel_reference.py apply \
#       ~/dev/claude-sandboxes/AutoGating-sandbox/scripts/engine_reference/model_r.agmodel
# and source this file again for part 3.

pkg.dir <- "~/dev/claude-sandboxes/AutoGating-sandbox"
if ( !exists( "read.agmodel" ) )
    suppressPackageStartupMessages( devtools::load_all( pkg.dir, quiet = TRUE ) )
ref.dir <- file.path( pkg.dir, "scripts", "engine_reference" )
read.events <- function( f ) as.matrix( utils::read.csv( file.path( ref.dir, f ), check.names = FALSE ) )
data.b <- read.events( "data_b.csv" )

apply.model <- function( model ) {
    app <- model$application
    run.gating( model$gate.definitions, data.b, reference = model$trained.boundaries,
                modes = app$mode.per.gate, default.mode = app$default.mode,
                drift.limit = app$drift.limit, drift.action = app$drift.action,
                transforms = model$transforms )
}
counts <- function( res ) {
    do.call( rbind, lapply( names( res$boundaries ), function( g )
        do.call( rbind, lapply( names( res$boundaries[[ g ]] ), function( p ) {
            e <- res$boundaries[[ g ]][[ p ]]
            data.frame( gate = g, population = p, n.events = e$n.events,
                        threshold.x = e$threshold.x %||% NA_real_,
                        threshold.y = e$threshold.y %||% NA_real_,
                        source = res$status[[ g ]]$source, stringsAsFactors = FALSE )
        } ) ) ) )
}
compare <- function( r, py, title ) {
    m <- merge( py, r, by = c( "gate", "population" ), all = TRUE, suffixes = c( ".py", ".r" ) )
    cat( sprintf( "\n%s\n", title ) )
    print( format( data.frame( gate = m$gate, pop = m$population,
        n.py = m$n_events, n.r = m$n.events,
        dtx = abs( m$threshold_x - m$threshold.x ), dty = abs( m$threshold_y - m$threshold.y ),
        src.py = m$source.py, src.r = m$source.r ), digits = 3 ), row.names = FALSE )
}


# ---- 1. Model saved by Honeychrome, applied in R ----------------------------

model.py <- read.agmodel( file.path( ref.dir, "model_py.agmodel" ) )
res.r <- apply.model( model.py )
compare( counts( res.r ), utils::read.csv( file.path( ref.dir, "counts_py.csv" ), stringsAsFactors = FALSE ),
         "1. Honeychrome model applied to sample B: Python vs R" )


# ---- 2. Written by R and read back --------------------------------------------

path.r <- file.path( ref.dir, "model_r.agmodel" )
f <- file.path( ref.dir, "counts_r_model.csv" )
applied.by.python <- file.exists( f ) && file.exists( path.r ) && file.mtime( f ) >= file.mtime( path.r )
if ( !applied.by.python )
    write.agmodel( model.py, path.r )
model.back <- read.agmodel( path.r )
same.defs <- identical( model.back$gate.definitions, model.py$gate.definitions )
d <- unlist( Map( function( a, b ) Map( function( x, y ) {
    if ( is.matrix( x$boundary ) ) max( abs( x$boundary - y$boundary ) ) else 0
}, a, b ), model.py$trained.boundaries, model.back$trained.boundaries ) )
cat( sprintf( "\n2. R write and read back: gate definitions identical %s; largest boundary difference %.2g\n",
              same.defs, max( d ) ) )


# ---- 3. Written by R, applied by Honeychrome ------------------------------------

if ( applied.by.python ) {
    compare( counts( apply.model( model.back ) ), utils::read.csv( f, stringsAsFactors = FALSE ),
             "3. Model written by R, applied to sample B: Honeychrome vs R" )
} else {
    cat( "\n3. Run agmodel_reference.py apply on model_r.agmodel, then source this file again.\n" )
}
