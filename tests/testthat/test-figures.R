linear.ticks <- function( top ) {
    tr <- make.transform( list( id = 0L, scale.t = top, linear.a = 0 ) )
    axis.ticks( list( id = 0L, scale.t = top ), tr$forward( c( 0, top ) ), tr$forward, tr$inverse )
}

test_that( "linear axes are labelled in raw units for any instrument range", {
    aurora <- linear.ticks( 4194304 )
    expect_equal( aurora$labels, c( "0", "1M", "2M", "3M", "4M" ) )
    expect_equal( aurora$breaks, c( 0, 1e6, 2e6, 3e6, 4e6 ) / 4194304 )
    symphony <- linear.ticks( 262144 )
    expect_equal( symphony$labels, c( "0", "50K", "100K", "150K", "200K", "250K" ) )
    expect_equal( linear.ticks( 900 )$labels, c( "0", "200", "400", "600", "800" ) )
} )

test_that( "logicle axes get decade labels and leave crowded decades unlabelled", {
    params <- list( id = 1L, scale.t = 262144, logicle.w = 0.5, logicle.m = 4.418539922, logicle.a = 0 )
    tr <- make.transform( params )
    limits <- tr$forward( c( -2000, 262144 ) )
    ticks <- axis.ticks( params, limits, tr$forward, tr$inverse )
    text <- vapply( ticks$labels, function( l ) paste( deparse( l ), collapse = "" ), character( 1 ) )
    expect_true( all( c( "0", "10^3", "10^4", "10^5" ) %in% text ) )
    expect_true( all( diff( ticks$breaks ) >= 0.05 * diff( limits ) - 1e-12 ) )
    expect_equal( ticks$breaks[ text == "10^4" ], tr$forward( 1e4 ) )
    expect_true( any( abs( ticks$minor - tr$forward( 200 ) ) < 1e-12 ) )
} )

test_that( "axis ranges follow the panel type of each channel", {
    panel <- data.frame( dye = c( "SSC-A", "B2-A", "Time" ), antigen = c( "SSC-A", "CD4", "Time" ),
                         type = c( "scatter", "lineage", "time" ), stringsAsFactors = FALSE )
    lin <- make.transform( list( id = 0L, scale.t = 4194304 ) )
    lgc <- make.transform( list( id = 1L, scale.t = 4194304, logicle.w = 1.5, logicle.m = 5.62, logicle.a = 0 ) )
    id <- function( x ) x
    setup <- list( fcs.panel = panel,
                   ap = list( transform = list( `SSC-A` = lin$forward, CD4 = lgc$forward, Time = id ),
                              transform.inv = list( `SSC-A` = lin$inverse, CD4 = lgc$inverse, Time = id ),
                              transform.params = list( `SSC-A` = list( id = 0L, scale.t = 4194304 ),
                                                       CD4 = list( id = 1L, scale.t = 4194304, logicle.w = 1.5,
                                                                   logicle.m = 5.62, logicle.a = 0 ),
                                                       Time = list() ) ) )
    data <- cbind( `SSC-A` = c( 0.1, 0.5 ), CD4 = c( 0.3, 0.6 ), Time = c( 2, 90 ) )
    ap <- list( figure.limits = list( scatter.min = 0, scatter.max = 4200000, expr.min = -5000, expr.max = 4200000 ) )
    axes <- figure.axes( setup, ap, data )
    expect_equal( axes$`SSC-A`$kind, "scatter" )
    expect_equal( axes$`SSC-A`$limits, lin$forward( c( 0, 4200000 ) ) )
    expect_equal( axes$`SSC-A`$labels[ 2 ], "1M" )
    expect_equal( axes$CD4$kind, "fluorescence" )
    expect_equal( axes$CD4$limits, lgc$forward( c( -5000, 4200000 ) ) )
    expect_equal( axes$Time$limits, c( 2, 90 ) )
} )

test_that( "population labels sit in the outer corner of their region", {
    gate <- list( gate.name = "q", gate.type = "2dsep", gate.marker.x = "a", gate.marker.y = "b",
                  populations = list( rest = list( region = "x-y-|x+y-|x-y+" ), DP = list( region = "x+y+" ) ) )
    boundary <- list( rest = list( threshold.x = 0.5, threshold.y = 0.4, region = "x-y-|x+y-|x-y+", fraction = 0.75 ),
                      DP = list( threshold.x = 0.5, threshold.y = 0.4, region = "x+y+", fraction = 0.25 ) )
    ann <- gate.annotations( gate, boundary, c( "a", "b" ), c( 0, 1 ), c( 0, 1 ) )
    expect_equal( sort( ann$lines$at ), c( 0.4, 0.5 ) )
    expect_equal( ann$labels$hjust, c( 0, 1 ) )
    expect_equal( ann$labels$vjust, c( 0, 1 ) )
    expect_equal( ann$labels$text, c( "rest\n75.0%", "DP\n25.0%" ) )
} )

test_that( "subsamples leave the caller's random numbers as they were", {
    set.seed( 7 )
    before <- .Random.seed
    a <- seeded.call( 1L, function() sample.int( 100, 5 ) )
    expect_identical( .Random.seed, before )
    expect_identical( a, seeded.call( 1L, function() sample.int( 100, 5 ) ) )
} )

test_that( "a gate figure is drawn from the gating result", {
    set.seed( 4 )
    data <- cbind( a = c( stats::rnorm( 300, 0.3, 0.03 ), stats::rnorm( 300, 0.7, 0.03 ) ),
                   b = stats::runif( 600 ) )
    gate <- list( gate.name = "g", gate.number = 1L, gate.type = "1dsep", gate.marker.x = "a",
                  gate.marker.y = "b", algorithm = "mixture",
                  populations = default.populations.for.type( "1dsep" ) )
    res <- run.gating( list( gate ), data )
    axes <- lapply( c( a = "a", b = "b" ), function( ch ) list( title = ch, kind = "fluorescence",
        limits = c( 0, 1 ), breaks = c( 0, 0.5, 1 ), labels = c( "0", "0.5", "1" ), minor = numeric( 0 ) ) )
    fig <- gate.figure( gate, res, data, axes, list( figure.max.events = 200L ) )
    expect_s3_class( fig, "ggplot" )
    expect_null( gate.figure( modifyList( gate, list( gate.name = "missing" ) ), res, data, axes, list() ) )
} )
