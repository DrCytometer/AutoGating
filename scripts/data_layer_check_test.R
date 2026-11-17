# data_layer_check_test.R
#
# Checks the data layer (AutoSpectral::readFCS, compensation, Gating-ML
# transforms) of AutoGating.
#
# Before running, in the Honeychrome virtual environment:
#   python ~/dev/claude-sandboxes/AutoGating-sandbox/scripts/transform_reference.py
# Then:
#   source( "~/dev/claude-sandboxes/AutoGating-sandbox/scripts/data_layer_check_test.R" )
#
# Section 2 reads the first three Aurora spleen files in
# claude-sandboxes/FlowCode context/FlowCodeAnalysis/SpectroFlo unmixed/Spleen,
# with a panel made by create.panel.file() and a synthetic spillover matrix
# (the files are already unmixed; the matrix only exercises compensation).
# flowCore is used here only as a reference.

data.dir  <- "~/dev/claude-sandboxes/FlowCode context/FlowCodeAnalysis/SpectroFlo unmixed/Spleen"
n.samples <- 3L

pkg.dir <- "~/dev/claude-sandboxes/AutoGating-sandbox"
suppressPackageStartupMessages( devtools::load_all( pkg.dir, quiet = TRUE ) )


# ---- 1. Transforms against FlowKit ------------------------------------------

ref <- utils::read.csv( file.path( pkg.dir, "scripts", "transform_reference.csv" ) )
key <- with( ref, paste( id, scale_t, linear_a, logicle_w, logicle_m, logicle_a, log_m ) )
cat( "1. Transforms against FlowKit (max |difference| on the 0-1 scale)\n" )
for ( k in unique( key ) ) {
    r  <- ref[ key == k, ]
    p  <- list( id = r$id[ 1 ], scale.t = r$scale_t[ 1 ], linear.a = r$linear_a[ 1 ],
                logicle.w = r$logicle_w[ 1 ], logicle.m = r$logicle_m[ 1 ],
                logicle.a = r$logicle_a[ 1 ], log.m = r$log_m[ 1 ] )
    tf <- tryCatch( make.transform( p ), error = function( e ) e )
    if ( inherits( tf, "error" ) ) {
        cat( sprintf( "   id %d  %-48s not built: %s\n", p$id, sub( "^[0-9]+ ", "", k ),
            conditionMessage( tf ) ) )
        next
    }
    y  <- tf$forward( r$x )
    x2 <- tf$inverse( y )
    cat( sprintf( "   id %d  %-48s forward %.2e   round trip (rel) %.2e\n",
        p$id, sub( "^[0-9]+ ", "", k ), max( abs( y - r$y ) ),
        max( abs( x2 - r$x ) / pmax( 1, abs( r$x ) ) ) ) )
}


suppressPackageStartupMessages( library( flowCore ) )

cat( sprintf( "\n2. Data layer on %d Aurora files\n", n.samples ) )
work.dir <- file.path( tempdir(), "ag_data_layer" )
dir.create( work.dir, showWarnings = FALSE )
fcs.files <- sort( list.files( path.expand( data.dir ), pattern = "\\.fcs$" ) )[ seq_len( n.samples ) ]

panel.file <- file.path( work.dir, "panel" )
panel <- create.panel.file( file.path( path.expand( data.dir ), fcs.files[ 1 ] ), panel.file )
panel$antigen <- make.unique( panel$antigen )
utils::write.csv( panel, paste0( panel.file, ".csv" ), row.names = FALSE )
cat( sprintf( "   panel: %d channels (%s)\n", nrow( panel ),
    paste( names( table( panel$type ) ), table( panel$type ), sep = " ", collapse = ", " ) ) )

expt.file <- file.path( work.dir, "experiment.csv" )
utils::write.csv( data.frame( filename = fcs.files, sample = paste0( "S", seq_len( n.samples ) ),
    gate.calculation = TRUE ), expt.file, row.names = FALSE )

fluor.dye <- panel$dye[ panel$type == "lineage" & panel$antigen != "AutoF" ]
set.seed( 1L )
spill <- diag( length( fluor.dye ) )
spill[ spill == 0 ] <- stats::runif( sum( spill == 0 ), 0, 0.05 )
dimnames( spill ) <- list( fluor.dye, fluor.dye )
comp.file <- file.path( work.dir, "spillover.csv" )
utils::write.csv( spill, comp.file )

ap <- get.autogating.param(
    cytometer             = "aurora",
    data.dir              = path.expand( data.dir ),
    panel.filename        = paste0( panel.file, ".csv" ),
    experiment.filename   = expt.file,
    compensation.filename = comp.file )

el <- system.time( fcs.setup <- setup.autogating( ap ) )[ "elapsed" ]
cat( sprintf( "   setup.autogating: %.1f s, %d events, %.0f MB matrix\n", el,
    fcs.setup$fcs.event.n, as.numeric( utils::object.size( fcs.setup$fcs.expr.data ) ) / 2^20 ) )

panel <- fcs.setup$fcs.panel
comp  <- read.compensation.csv( ap$compensation.filename, panel,
    fcs.setup$fcs.panel.lineage.idx, fcs.setup$fcs.panel.activation.idx )
tf    <- fcs.setup$ap$transform

for ( i in seq_len( n.samples ) ) {
    path <- file.path( ap$data.dir, fcs.files[ i ] )
    ff <- suppressWarnings( flowCore::read.FCS( path, transformation = FALSE,
        truncate_max_range = FALSE ) )
    ff <- ff[ , panel$dye ]
    colnames( ff ) <- panel$antigen
    ff <- flowCore::compensate( ff, flowCore::compensation( comp ) )
    raw <- flowCore::exprs( ff )

    ref.t <- raw
    for ( ant in panel$antigen ) ref.t[ , ant ] <- tf[[ ant ]]( raw[ , ant ] )
    new <- fcs.setup$fcs.expr.data[ fcs.setup$fcs.event.sample == paste0( "S", i ), , drop = FALSE ]
    ref.t <- ref.t[ , colnames( new ), drop = FALSE ]

    trunc <- flowCore::exprs( suppressWarnings( flowCore::read.FCS( path,
        transformation = FALSE, truncate_max_range = TRUE ) ) )[ , panel$dye ]
    untrunc <- flowCore::exprs( suppressWarnings( flowCore::read.FCS( path,
        transformation = FALSE, truncate_max_range = FALSE ) ) )[ , panel$dye ]

    cat( sprintf( "   S%d: events %d vs %d, max |difference| %.2e, values the old truncation changed %d\n",
        i, nrow( new ), nrow( ref.t ), max( abs( new - ref.t ) ), sum( trunc != untrunc ) ) )
}

cat( "\n   Display range per channel type (S1):\n" )
s1 <- fcs.setup$fcs.expr.data[ fcs.setup$fcs.event.sample == "S1", ]
for ( ty in c( "scatter", "time", "lineage", "activation" ) ) {
    ch <- panel$antigen[ panel$type == ty ]
    if ( length( ch ) )
        cat( sprintf( "   %-10s [%9.4g, %9.4g]\n", ty, min( s1[ , ch ] ), max( s1[ , ch ] ) ) )
}
