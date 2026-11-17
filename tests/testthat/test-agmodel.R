make.events <- function( n, pos.mean, seed ) {
    set.seed( seed )
    pos <- stats::runif( n ) < 0.3
    cbind( A = ifelse( pos, stats::rnorm( n, pos.mean, 0.05 ), stats::rnorm( n, 0.25, 0.04 ) ),
           B = ifelse( stats::runif( n ) < 0.5, stats::rnorm( n, 0.6, 0.05 ), stats::rnorm( n, 0.2, 0.04 ) ) )
}

gates <- list(
    list( gate.name = "a", gate.number = 1L, gate.type = "1dsep", gate.marker.x = "A",
          populations = list( neg = list( label = "neg", region = "below" ),
                              pos = list( label = "pos", region = "above" ) ),
          algorithm = "mixture", gate.param = list( mixture.components = 2L ) ),
    list( gate.name = "q", gate.number = 2L, gate.type = "2dsep", gate.marker.x = "A",
          gate.marker.y = "B", parent.gate = "a", parent.popul = "pos",
          populations = list( dn = list( region = "x-y-" ), rest = list( region = "x+y-|x-y+|x+y+" ) ),
          algorithm = "tail", gate.param = list( tail.fraction.y = 0.01 ),
          stats.parent = list( dn = "root", rest = "a" ) ) )

test_that( "a model written and read back applies identically", {
    train <- train.boundaries( gates, make.events( 4000, 0.6, 1 ) )
    model <- new.agmodel( gates, train$boundaries,
                          transforms = list( A = list( id = 0L, scale.t = 262144, linear.a = 0 ) ),
                          application = list( mode.per.gate = list( a = "fixed", q = "recalculate" ) ) )
    path <- tempfile( fileext = ".agmodel" )
    write.agmodel( model, path )
    back <- read.agmodel( path )
    for ( i in seq_along( gates ) )
        expect_equal( back$gate.definitions[[ i ]][ names( gates[[ i ]] ) ], gates[[ i ]] )
    expect_equal( back$trained.boundaries, model$trained.boundaries, tolerance = 1e-14 )
    expect_equal( back$application$mode.per.gate, model$application$mode.per.gate )
    data <- make.events( 4000, 0.8, 2 )
    run <- function( m ) run.gating( m$gate.definitions, data, reference = m$trained.boundaries,
                                     modes = m$application$mode.per.gate )
    expect_identical( run( back )$events, run( model )$events )
} )

test_that( "earlier gate lists are converted", {
    old <- list(
        list( gate.number = 0L, gate.name = "all", gate.type = "free", gate.marker.x = "A",
              gate.marker.y = "B", populations = list( all = list( label = "All" ) ),
              algorithm = "all", stats.parent = list( gate = "*", popul = "*" ) ),
        list( gate.number = 1L, gate.name = "sep", gate.type = "1dsep", gate.marker.x = "B",
              gate.marker.y = "A", parent.gate = "all", parent.popul = "all",
              populations = list( `*` = list( label = "*" ), hi = list( label = "Hi" ) ),
              algorithm = "kde.min", gate.param = list( sep.axis = 2L, region.min = 0.2 ) ),
        list( gate.number = 2L, gate.name = "quad", gate.type = "2dsep", gate.marker.x = "A",
              gate.marker.y = "B", parent.gate = "sep", parent.popul = "hi",
              populations = list( l = list( label = "L" ) ), algorithm = c( "mixture", "tail" ),
              gate.param = list( join = "1:2:3", threshold.param = 0.99 ),
              stats.parent = list( gate = c( "parent.1", "parent.2" ), popul = c( "parent.1", "parent.2" ) ) ),
        list( gate.number = 3L, gate.name = "copy", gate.type = "free", gate.marker.x = "A",
              gate.marker.y = "B", parent.gate = "sep", parent.popul = "hi",
              populations = list( l2 = list( label = "L2" ) ), algorithm = "replicate",
              gate.param = list( origin.gate = "quad" ) ) )
    up <- convert.gate.definitions( old )
    g <- gates.by.name( up$gate.definitions )
    expect_named( g, c( "sep", "quad", "copy" ) )
    expect_equal( g$sep$gate.marker.x, "A" )
    expect_equal( g$sep$gate.marker.y, "B" )
    expect_equal( g$sep$sep.axis, "y" )
    expect_equal( g$sep$parent.popul, "root" )
    expect_equal( g$sep$algorithm, "kde_min" )
    expect_equal( g$sep$populations$hi$region, "above" )
    expect_equal( g$sep$gate.param, list( axis.lo = 0.2 ) )
    expect_equal( g$quad$populations$l$region, "x-y-|x+y-|x-y+" )
    expect_equal( g$quad$gate.param, list( algorithm.y = "tail", tail.fraction.y = 0.01 ) )
    expect_equal( g$quad$stats.parent, list( l = "root" ) )
    expect_equal( g$copy$gate.type, "2dsep" )
    expect_equal( g$copy$populations$l2$region, "x-y-|x+y-|x-y+" )
    expect_equal( g$copy$gate.param$replicate.x, list( gate = "quad" ) )
    res <- train.boundaries( up$gate.definitions, make.events( 4000, 0.6, 3 ) )
    expect_identical( res$events$copy$l2, res$events$quad$l )
} )

test_that( "boundaries move with a change of transform", {
    model <- new.agmodel( gates[ 1 ],
        list( a = list( pos = list( boundary = rect.polygon( 0.5, 1 ), threshold.x = 0.5,
                                    threshold.y = NA_real_, region = "above" ) ) ),
        transforms = list( A = list( id = 0L, scale.t = 1000, linear.a = 0 ) ) )
    out <- remap.model( model, list( A = list( id = 0L, scale.t = 2000, linear.a = 0 ) ) )
    expect_equal( out$model$trained.boundaries$a$pos$threshold.x, 0.25 )
    expect_equal( out$model$transforms$A$scale.t, 2000 )
    expect_length( out$messages, 1L )
} )
