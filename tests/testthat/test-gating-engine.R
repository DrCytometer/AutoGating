make.events <- function( n, pos.mean, seed ) {
    set.seed( seed )
    pos <- stats::runif( n ) < 0.3
    cbind( A = ifelse( pos, stats::rnorm( n, pos.mean, 0.05 ), stats::rnorm( n, 0.25, 0.04 ) ),
           B = ifelse( stats::runif( n ) < 0.5, stats::rnorm( n, 0.6, 0.05 ), stats::rnorm( n, 0.2, 0.04 ) ) )
}

gate.1d <- list( gate.name = "a", gate.number = 1L, gate.type = "1dsep",
    gate.marker.x = "A", populations = list(
        neg = list( label = "neg", region = "below" ),
        pos = list( label = "pos", region = "above" ) ),
    algorithm = "mixture", gate.param = list() )

test_that( "threshold populations cover the parent on brighter data", {
    train <- train.boundaries( list( gate.1d ), make.events( 5000, 0.6, 1 ) )
    bright <- make.events( 5000, 0.95, 2 )
    res <- run.gating( list( gate.1d ), bright, reference = train$boundaries )
    expect_equal( length( res$events$a$neg ) + length( res$events$a$pos ), nrow( bright ) )
    expect_length( intersect( res$events$a$neg, res$events$a$pos ), 0L )
} )

test_that( "joined quadrant regions are the union of their quadrants", {
    g2 <- list( gate.name = "q", gate.number = 1L, gate.type = "2dsep",
        gate.marker.x = "A", gate.marker.y = "B",
        populations = list( dn = list( region = "x-y-" ), rest = list( region = "x+y-|x-y+|x+y+" ) ),
        algorithm = "tail", gate.param = list() )
    res <- train.boundaries( list( g2 ), make.events( 4000, 0.6, 3 ) )
    expect_equal( sort( c( res$events$q$dn, res$events$q$rest ) ), seq_len( 4000 ) )
    poly <- res$boundaries$q$rest$boundary
    expect_equal( nrow( poly ), 7L )
} )

test_that( "a threshold is replicated by gate name and by reference population", {
    rep.gate <- list( gate.name = "r", gate.number = 2L, gate.type = "1dsep",
        gate.marker.x = "A", parent.gate = "a", parent.popul = "neg",
        populations = list( lo = list( region = "below" ), hi = list( region = "above" ) ),
        algorithm = "replicate", gate.param = list( replicate.x = list( gate = "a" ) ) )
    ref.gate <- rep.gate
    ref.gate$gate.name <- "s"
    ref.gate$gate.param <- list( replicate.x = list( population = "a/pos", algorithm = "tail" ) )
    res <- train.boundaries( list( gate.1d, rep.gate, ref.gate ), make.events( 5000, 0.6, 4 ) )
    expect_equal( res$boundaries$r$lo$threshold.x, res$boundaries$a$neg$threshold.x )
    expect_true( is.finite( res$boundaries$s$lo$threshold.x ) )
    expect_true( is.na( res$status$s$flag ) )
} )

test_that( "gates are ordered after the gates they depend on", {
    g <- list( list( gate.name = "child", gate.number = 1L, parent.gate = "parent" ),
               list( gate.name = "parent", gate.number = 2L ) )
    expect_equal( vapply( ordered.gate.defs( g ), `[[`, "", "gate.name" ), c( "parent", "child" ) )
} )

test_that( "the mixture crossover uses unequal variances", {
    set.seed( 5 )
    x <- c( stats::rnorm( 8000, 0, 10 ), stats::rnorm( 2000, 60, 20 ) )
    fit <- fit.mixture( x, 2L )
    expect_gt( max( fit$sd ) / min( fit$sd ), 1.5 )
    t <- threshold.mixture( x )
    expect_gt( t, 20 )
    expect_lt( t, 45 )
} )

test_that( "Otsu splits two separated populations", {
    set.seed( 6 )
    x <- c( stats::rnorm( 3000, 0.2, 0.03 ), stats::rnorm( 3000, 0.7, 0.03 ) )
    t <- threshold.otsu( x )
    expect_gt( t, 0.3 )
    expect_lt( t, 0.6 )
} )
