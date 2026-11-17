# stats_check_test.R
#
# Compares AutoGating's differential.statistics with Honeychrome's
# ag_results.run_statistics on one synthetic experiment and five designs.
#
# Before running, in the Honeychrome virtual environment:
#   python ~/dev/claude-sandboxes/AutoGating-sandbox/scripts/stats_reference.py
# Then:
#   source( "~/dev/claude-sandboxes/AutoGating-sandbox/scripts/stats_check_test.R" )

pkg.dir <- "~/dev/claude-sandboxes/AutoGating-sandbox"
if ( !exists( "differential.statistics" ) )
    suppressPackageStartupMessages( devtools::load_all( pkg.dir, quiet = TRUE ) )
ref.dir <- file.path( pkg.dir, "scripts", "stats_reference" )
rd <- function( f ) utils::read.csv( file.path( ref.dir, f ), stringsAsFactors = FALSE, check.names = FALSE )

gates <- lapply( jsonlite::read_json( file.path( ref.dir, "gates.json" ), simplifyVector = FALSE ), gate.from.json )
pt <- rd( "population_table.csv" )
pt$stats.parent[ pt$stats.parent == "" ] <- NA
gating <- list( population.table = pt, medians = rd( "medians.csv" ), model = list( gate.definitions = gates ) )
st <- rd( "samples.csv" )
designs <- list(
    reference = list( contrast.mode = "reference", reference = "A" ),
    pairwise  = list( contrast.mode = "pairwise" ),
    fixed     = list( contrast.mode = "reference", reference = "A", pairing.column = "donor" ),
    random    = list( contrast.mode = "reference", reference = "A", pairing.column = "mouse" ),
    covariate = list( contrast.mode = "reference", reference = "A", covariates = "age" ) )

compare <- function( py, r, by ) {
    m <- merge( py, r, by = by, suffixes = c( ".py", ".r" ) )
    d <- function( col ) suppressWarnings( max( abs( m[[ paste0( col, ".py" ) ]] - m[[ paste0( col, ".r" ) ]] ), na.rm = TRUE ) )
    sprintf( "%d/%d rows; max |diff| logFC %.1e, SE %.1e, P %.1e, adj %.1e, adj.global %.1e; significance differs in %d",
             nrow( m ), nrow( py ), d( "logFC" ), d( "SE" ), d( "P.Value" ), d( "adj.P.Val" ), d( "adj.P.Val.global" ),
             sum( as.logical( m$significant.py ) != m$significant.r ) )
}

for ( name in names( designs ) ) {
    res <- suppressMessages( do.call( differential.statistics,
        c( list( gating = gating, sample.table = st, group.order = c( "A", "B", "C" ) ), designs[[ name ]] ) ) )
    f.py <- rd( sprintf( "freq_%s.csv", name ) )
    m.py <- rd( sprintf( "markers_%s.csv", name ) )
    m.py$population <- m.py$cluster_id
    cat( sprintf( "\n%s (pairing: %s)\n  frequencies: %s\n  markers:     %s\n", name,
                  if ( nzchar( res$pairing.mode ) ) res$pairing.mode else "none",
                  compare( f.py, res$frequencies, c( "feature", "comparison" ) ),
                  compare( m.py, res$markers, c( "feature", "comparison" ) ) ) )
    mm <- merge( m.py, res$markers, by = c( "feature", "comparison" ) )
    cat( sprintf( "  stage-wise:  max |diff| %.1e\n",
                  max( abs( mm$stagewise.adj.P.Val.x - mm$stagewise.adj.P.Val.y ), na.rm = TRUE ) ) )
}
