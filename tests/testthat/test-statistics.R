test_that( "the share prior gives an empty population the same log odds in every sample", {
    counts <- cbind( a = c( 0, 0, 5 ), b = c( 10, 40, 20 ) )
    denom  <- cbind( a = c( 100, 1000, 500 ), b = c( 100, 1000, 500 ) )
    y <- population.log2.odds( counts, denom )
    expect_equal( y[ 1, "a" ], y[ 2, "a" ] )
    expect_equal( unname( y[ 1, "a" ] ), -log2( 1 + 2 * mean( denom[ , "a" ] ) ) )
} )

test_that( "frequencies are compared on the samples where the parent is large enough", {
    set.seed( 2 )
    samples <- sprintf( "s%02d", 1:12 )
    group <- rep( c( "ctrl", "treated arm" ), each = 6 )
    rows <- do.call( rbind, lapply( seq_along( samples ), function( i ) {
        parent <- c( 2000, 15, 3000, 2500, 1800, 2200, 2100, 2400, 1900, 2600, 2300, 2000 )[ i ]
        p <- c( 0.2, 0.4, 0.1 ) * if ( group[ i ] == "ctrl" ) 1 else 1.8
        n <- stats::rbinom( 3, parent, p )
        data.frame( sample = samples[ i ], key = c( "g/a", "g/b", "g/c" ), gate = "g",
                    population = c( "a", "b", "c" ), label = c( "a", "b", "c" ), n.events = n,
                    parent.key = "root", n.parent = parent, n.stats.parent = NA_real_,
                    stats.parent.missing = FALSE, stringsAsFactors = FALSE )
    } ) )
    gating <- list( population.table = rows, medians = data.frame(), model = list( gate.definitions = list() ) )
    res <- suppressMessages( differential.statistics( gating, data.frame( sample = samples, group = group ),
                                                      run.markers = FALSE ) )
    expect_true( all( is.na( res$frequency.values[ "s02", ] ) ) )
    expect_equal( unique( res$frequencies$comparison ), "treated arm vs ctrl" )
    expect_true( all( res$frequencies$logFC > 0 ) )
} )
