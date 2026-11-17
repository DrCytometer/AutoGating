workspace.xml <- '<?xml version="1.0" encoding="UTF-8"?>
<Workspace version="20.0" flowJoVersion="10.10.0" xmlns:gating="http://www.isac-net.org/std/Gating-ML/v2.0/gating" xmlns:data-type="http://www.isac-net.org/std/Gating-ML/v2.0/datatypes">
<SampleList><Sample>
<Keywords><Keyword name="$P1N" value="FSC-A"/><Keyword name="$P2N" value="SSC-A"/><Keyword name="$P3N" value="Comp-BV421-A"/><Keyword name="$P3S" value="CD3"/></Keywords>
<SampleNode name="s1.fcs" sampleID="1" count="1000"><Subpopulations>
<Population name="Cells" count="800"><Gate gating:id="g1"><gating:PolygonGate eventsInside="1">
<gating:dimension><data-type:fcs-dimension data-type:name="FSC-A"/></gating:dimension>
<gating:dimension><data-type:fcs-dimension data-type:name="SSC-A"/></gating:dimension>
<gating:vertex><gating:coordinate data-type:value="50000"/><gating:coordinate data-type:value="10000"/></gating:vertex>
<gating:vertex><gating:coordinate data-type:value="200000"/><gating:coordinate data-type:value="10000"/></gating:vertex>
<gating:vertex><gating:coordinate data-type:value="200000"/><gating:coordinate data-type:value="150000"/></gating:vertex>
</gating:PolygonGate></Gate><Subpopulations>
<Population name="T/cells" count="300"><Gate gating:id="g2"><gating:RectangleGate eventsInside="1">
<gating:dimension gating:min="2000"><data-type:fcs-dimension data-type:name="Comp-BV421-A"/></gating:dimension>
</gating:RectangleGate></Gate></Population>
</Subpopulations></Population>
</Subpopulations></SampleNode></Sample></SampleList></Workspace>'

test_that( "a FlowJo workspace imports as gates with their shapes as drawn", {
    path <- tempfile( fileext = ".wsp" )
    writeLines( workspace.xml, path )
    ws <- read.flowjo.workspace( path )
    src <- fj.default.source( ws )
    expect_equal( fj.gated.parameters( src ), c( "FSC-A", "SSC-A", "Comp-BV421-A" ) )
    panel <- data.frame( dye = c( "FSC-A", "SSC-A", "BV421-A" ), antigen = c( "FSC-A", "SSC-A", "CD3" ),
                         type = c( "scatter", "scatter", "lineage" ) )
    lin <- list( id = 0L, scale.t = 262144, linear.a = 0 )
    exp <- fj.experiment( panel, list( `FSC-A` = lin, `SSC-A` = lin, CD3 = lin ) )
    al <- align.flowjo.parameters( fj.gated.parameters( src ), src$labels, exp, ws$affixes )
    expect_equal( al$channel, c( "FSC-A", "SSC-A", "CD3" ) )
    res <- fj.build.gate.definitions( src, stats::setNames( al$channel, al$parameter ), exp )
    g <- gates.by.name( res$gate.definitions )
    expect_named( g, c( "Cells", "T∕cells" ) )
    expect_equal( g$Cells$gate.type, "free" )
    expect_equal( g$Cells$proposed.algorithm, "cells" )
    expect_equal( g[[ 2 ]]$gate.type, "1dsep" )
    expect_equal( g[[ 2 ]]$populations[[ 1 ]]$region, "above" )
    expect_equal( g[[ 2 ]]$proposed.algorithm, "mixture" )
    expect_equal( g[[ 2 ]]$template[[ 1 ]]$range[ 1 ], 2000 / 262144 )
    set.seed( 1 )
    data <- cbind( `FSC-A` = stats::runif( 2000, 0.3, 0.7 ), `SSC-A` = stats::runif( 2000, 0.1, 0.4 ),
                   CD3 = c( stats::rnorm( 1400, 0, 0.002 ), stats::rnorm( 600, 0.05, 0.005 ) ) )
    conv <- convert.imported.gates( res$gate.definitions, data, exp$transforms )
    expect_true( all( vapply( conv$gate.definitions, function( x ) x$algorithm != "imported", logical( 1 ) ) ) )
    expect_equal( nrow( conv$proposals ), 2L )
} )
