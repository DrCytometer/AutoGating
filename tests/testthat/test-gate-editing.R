editing.model <- function() {
    gates <- list(
        list( gate.name = "T", gate.number = 1L, gate.type = "1dsep", gate.marker.x = "a", algorithm = "mixture",
              populations = list( neg = list( label = "neg", region = "below" ),
                                  pos = list( label = "pos", region = "above" ) ) ),
        list( gate.name = "Q", gate.number = 2L, gate.type = "2dsep", gate.marker.x = "b", gate.marker.y = "c",
              parent.gate = "T", parent.popul = "pos", algorithm = "tail",
              gate.param = list( algorithm.x = "tail", algorithm.y = "mixture", axis.lo.x = 0.1 ),
              populations = list( rest = list( label = "rest", region = "x-y-|x+y-|x-y+" ),
                                  DP = list( label = "DP", region = "x+y+" ) ) ),
        list( gate.name = "R", gate.number = 3L, gate.type = "1dsep", gate.marker.x = "b", parent.gate = "Q",
              parent.popul = "DP", algorithm = "replicate",
              gate.param = list( replicate.x = list( gate = "Q" ) ),
              populations = default.populations.for.type( "1dsep" ) ) )
    m <- new.agmodel( gates )
    m$trained.boundaries <- list( T = list( neg = list( threshold.x = 0.4, region = "below" ),
                                            pos = list( threshold.x = 0.4, region = "above" ) ) )
    m
}

test_that( "form fields round-trip a gate definition", {
    m <- editing.model()
    q <- m$gate.definitions[[ 2 ]]
    v <- gate.form.values( q )
    expect_equal( unname( v$populations ), c( "rest", "rest", "rest", "DP" ) )
    expect_equal( v$axis.lo.x, 0.1 )
    expect_equal( v$algorithm.y, "mixture" )
    g <- gate.from.form( v, q )
    expect_equal( population.regions( g ), population.regions( q ) )
    expect_equal( g$gate.param$axis.lo.x, 0.1 )
    r <- gate.from.form( gate.form.values( m$gate.definitions[[ 3 ]] ), m$gate.definitions[[ 3 ]] )
    expect_equal( r$gate.param$replicate.x, list( gate = "Q" ) )
    v$populations[ "x+y-" ] <- ""
    expect_equal( unname( population.regions( gate.from.form( v, q ) ) ), c( "x-y-|x-y+", "x+y+" ) )
    v$x <- ""
    expect_null( gate.from.form( v, q )$gate.marker.x )
} )

test_that( "renaming a gate or population follows every reference", {
    m <- editing.model()
    v <- gate.form.values( m$gate.definitions[[ 1 ]] )
    v$name <- "T cells"
    v$populations[ "above" ] <- "T+"
    m2 <- put.gate( m, gate.from.form( v, m$gate.definitions[[ 1 ]] ), "T" )
    q <- gates.by.name( m2$gate.definitions )[[ "Q" ]]
    expect_equal( c( q$parent.gate, q$parent.popul ), c( "T cells", "T+" ) )
    expect_equal( names( m2$trained.boundaries[[ "T cells" ]] ), c( "neg", "T+" ) )
    expect_null( m2$trained.boundaries[[ "T" ]] )
    v2 <- gate.form.values( q )
    v2$name <- "Quad"
    m3 <- put.gate( m2, gate.from.form( v2, q ), "Q" )
    r <- gates.by.name( m3$gate.definitions )[[ "R" ]]
    expect_equal( r$parent.gate, "Quad" )
    expect_equal( r$gate.param$replicate.x$gate, "Quad" )
    v$algorithm <- "tail"
    m4 <- put.gate( m, gate.from.form( v, m$gate.definitions[[ 1 ]] ), "T" )
    expect_null( m4$trained.boundaries[[ "T cells" ]] )
    expect_error( put.gate( m, m$gate.definitions[[ 2 ]], NULL ), "already exists" )
} )

test_that( "a new gate needs an existing parent population and gets the next number", {
    m <- editing.model()
    v <- gate.form.values( NULL )
    v$name <- "child"
    v$x <- "a"
    v$parent <- "Q/DP"
    m2 <- put.gate( m, gate.from.form( v ), NULL, mode = "recalculate" )
    expect_equal( m2$gate.definitions[[ 4 ]]$gate.number, 4L )
    expect_equal( m2$application$mode.per.gate$child, "recalculate" )
    v$parent <- "Q/missing"
    expect_error( put.gate( m, gate.from.form( v ), NULL ), "no population" )
} )

test_that( "dropping a gate removes its descendants", {
    d <- drop.gate( editing.model(), "Q" )
    expect_equal( sort( d$removed ), c( "Q", "R" ) )
    expect_equal( length( d$model$gate.definitions ), 1L )
} )

test_that( "a channel can be given a new transform and the transforms written as CSV", {
    lin <- list( id = 0L, scale.t = 1000, linear.a = 0 )
    f <- make.transform( lin )
    raw <- c( -50, 0, 250, 900 )
    panel <- data.frame( dye = c( "D1", "D2" ), antigen = c( "a", "b" ), type = c( "lineage", "scatter" ),
                         stringsAsFactors = FALSE )
    setup <- list( fcs.panel = panel, fcs.expr.data = cbind( a = f$forward( raw ), b = f$forward( raw ) ),
                   ap = list( transform = list( a = f$forward, b = f$forward ),
                              transform.inv = list( a = f$inverse, b = f$inverse ),
                              transform.params = list( a = lin, b = lin ) ) )
    lgc <- list( id = 1L, scale.t = 1000, logicle.w = 0.5, logicle.m = 3, logicle.a = 0 )
    s2 <- retransform.channel( setup, "a", lgc )
    expect_equal( s2$ap$transform.inv$a( s2$fcs.expr.data[ , "a" ] ), raw, tolerance = 1e-6 )
    expect_equal( s2$fcs.expr.data[ , "b" ], setup$fcs.expr.data[ , "b" ] )
    file <- tempfile( fileext = ".csv" )
    write.transformation.csv( s2, file )
    back <- read.transformation( file, panel, 2L, 1L, integer( 0 ),
                                 list( transform.defaults = list( fluorescence = lin, scatter = lin ) ) )
    expect_equal( back[[ 3 ]]$a$logicle.w, 0.5 )
    expect_equal( back[[ 3 ]]$a$id, 1 )
    expect_equal( back[[ 3 ]]$b$id, 0 )
} )
