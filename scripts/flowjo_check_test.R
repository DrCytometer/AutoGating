# flowjo_check_test.R
#
# Compares the R FlowJo import and conversion with Honeychrome's (ag_flowjo,
# ag_review snapshots).
#
# Before running, in the Honeychrome virtual environment:
#   python ~/dev/claude-sandboxes/AutoGating-sandbox/scripts/flowjo_reference.py \
#       ~/dev/claude-sandboxes/FlowJo_workspace/Workspace_test.wsp \
#       "~/dev/claude-sandboxes/FlowJo_workspace/20251212_IMMPROVE panel template_new_panel_v2.wspt"
# Then:
#   source( "~/dev/claude-sandboxes/AutoGating-sandbox/scripts/flowjo_check_test.R" )

pkg.dir <- "~/dev/claude-sandboxes/AutoGating-sandbox"
ws.dir  <- "~/dev/claude-sandboxes/FlowJo_workspace"
if ( !exists( "import.flowjo.workspace" ) )
    suppressPackageStartupMessages( devtools::load_all( pkg.dir, quiet = TRUE ) )
ref.dir <- file.path( pkg.dir, "scripts", "flowjo_reference" )
workspaces <- c( Workspace_test = "Workspace_test.wsp",
                 `20251212_IMMPROVE_panel_template_new_panel_v2` = "20251212_IMMPROVE panel template_new_panel_v2.wspt" )

num <- function( v ) suppressWarnings( as.numeric( v ) )

for ( tag in names( workspaces ) ) {
    cat( sprintf( "\n==== %s ====\n", tag ) )
    ws <- read.flowjo.workspace( file.path( ws.dir, workspaces[[ tag ]] ) )
    src <- fj.default.source( ws )
    params <- fj.gated.parameters( src )
    ex <- utils::read.csv( file.path( ref.dir, paste0( tag, "_experiment.csv" ) ), stringsAsFactors = FALSE,
                           check.names = FALSE )
    panel <- data.frame( dye = ex$channel, antigen = ex$antigen, type = ex$type, stringsAsFactors = FALSE )
    keys <- c( "id", "scale_t", "linear_a", "logicle_w", "logicle_m", "logicle_a" )
    transforms <- stats::setNames( lapply( seq_len( nrow( ex ) ), function( i ) {
        p <- lapply( keys, function( k ) ex[[ k ]][ i ] )
        names( p ) <- gsub( "_", ".", keys )
        p[ !is.na( unlist( p ) ) ]
    } ), panel$antigen )
    pair <- if ( all( c( "FSC-A", "FSC-H" ) %in% panel$antigen ) ) c( "FSC-A", "FSC-H" ) else NULL
    exp <- fj.experiment( panel, transforms, pair )

    # ---- 1. Parameter alignment ---------------------------------------------
    al.py <- utils::read.csv( file.path( ref.dir, paste0( tag, "_alignment.csv" ) ), stringsAsFactors = FALSE,
                              check.names = FALSE, na.strings = character( 0 ) )
    stems <- vapply( params, fj.parameter.stem, "", affixes = ws$affixes )
    named <- ifelse( grepl( "-A$", panel$dye ) & panel$type != "scatter", sub( "-A$", "", panel$dye ), panel$dye )
    lab <- vapply( panel$dye, function( s ) {
        p <- params[ match( s, stems ) ]
        l <- if ( !is.na( p ) && p %in% names( src$labels ) ) src$labels[[ p ]] else ""
        t <- strsplit( l, "_", fixed = TRUE )[[ 1L ]][ 1L ]
        if ( is.na( t ) || !nzchar( t ) ) "" else t
    }, "" )
    panel.al <- data.frame( dye = named, antigen = ifelse( nzchar( lab ), lab, named ), type = panel$type,
                            stringsAsFactors = FALSE )
    al.r <- align.flowjo.parameters( params, src$labels, fj.experiment( panel.al ), ws$affixes )
    m <- merge( al.py, al.r, by = "parameter", suffixes = c( ".py", ".r" ) )
    m$channel.r[ is.na( m$channel.r ) ] <- ""
    m$method.r[ is.na( m$method.r ) ] <- ""
    diff.al <- m[ m$antigen != m$channel.r | m$method.py != m$method.r, c( "parameter", "antigen", "method.py",
                                                                           "channel.r", "method.r" ) ]
    cat( sprintf( "1. Alignment: %d parameters, %d differ\n", nrow( m ), nrow( diff.al ) ) )
    if ( nrow( diff.al ) ) print( diff.al, row.names = FALSE )

    # ---- 2. Gate definitions --------------------------------------------------
    py <- jsonlite::read_json( file.path( ref.dir, paste0( tag, "_gates.json" ) ), simplifyVector = FALSE )
    gp <- lapply( py$gate_definitions, gate.from.json )
    names( gp ) <- vapply( gp, `[[`, "", "gate.name" )
    mapping <- stats::setNames( stems, params )
    el <- system.time( res <- fj.build.gate.definitions( src, mapping, exp, workspace.name = basename( workspaces[[ tag ]] ) ) )[ "elapsed" ]
    gr <- res$gate.definitions
    names( gr ) <- vapply( gr, `[[`, "", "gate.name" )
    cat( sprintf( "2. Gates: Python %d, R %d (%.1f s); names only in Python %d, only in R %d\n",
                  length( gp ), length( gr ), el, length( setdiff( names( gp ), names( gr ) ) ),
                  length( setdiff( names( gr ), names( gp ) ) ) ) )
    common <- intersect( names( gp ), names( gr ) )
    fields <- c( "gate.type", "gate.marker.x", "gate.marker.y", "sep.axis", "parent.gate", "parent.popul",
                 "proposed.algorithm", "proposal.rule", "name.group", "gate.number" )
    field.diff <- sapply( fields, function( f ) sum( !vapply( common, function( n )
        identical( as.character( gp[[ n ]][[ f ]] %||% "" ), as.character( gr[[ n ]][[ f ]] %||% "" ) ), logical( 1 ) ) ) )
    reg.diff <- sum( !vapply( common, function( n ) identical(
        vapply( gp[[ n ]]$populations, `[[`, "", "region" ), vapply( gr[[ n ]]$populations, `[[`, "", "region" ) ),
        logical( 1 ) ) )
    shape.d <- unlist( lapply( common, function( n ) lapply( names( gp[[ n ]]$template ), function( p ) {
        a <- gp[[ n ]]$template[[ p ]]$boundary
        b <- gr[[ n ]]$template[[ p ]]$boundary
        if ( is.matrix( a ) && is.matrix( b ) && all( dim( a ) == dim( b ) ) ) max( abs( a - b ) ) else NA
    } ) ) )
    cat( "   gates differing per field:\n" )
    print( field.diff )
    cat( sprintf( "   population regions differ in %d gates; shapes: %d compared, %d with different vertex counts, largest vertex difference %.2g\n",
                  reg.diff, sum( !is.na( shape.d ) ), sum( is.na( shape.d ) ), max( shape.d, na.rm = TRUE ) ) )
    cat( sprintf( "   dropped: Python %d lines, R %d lines; same %s\n", length( py$dropped ),
                  length( fj.plan.messages( res$plan ) ),
                  setequal( as.character( unlist( py$dropped ) ), fj.plan.messages( res$plan ) ) ) )

    # ---- 3. Conversion of imported gates on synthetic events -----------------
    ev <- as.matrix( utils::read.csv( file.path( ref.dir, paste0( tag, "_events.csv" ) ), check.names = FALSE ) )
    el <- system.time( conv <- convert.imported.gates( unname( gr ), ev, transforms ) )[ "elapsed" ]
    rv.py <- utils::read.csv( file.path( ref.dir, paste0( tag, "_review.csv" ) ), stringsAsFactors = FALSE,
                              na.strings = "" )
    pr <- conv$proposals
    mm <- merge( rv.py, pr, by = "gate" )
    same.alg <- mm$type.x == mm$type.y & mm$algorithm_x == mm$algorithm.x &
        ( is.na( mm$algorithm_y ) & is.na( mm$algorithm.y ) | mm$algorithm_y %in% mm$algorithm.y )
    cat( sprintf( "3. Conversion (%.1f s): %d gates, type and algorithms differ in %d\n", el, nrow( mm ),
                  sum( !same.alg ) ) )
    cat( sprintf( "   largest difference: below zero %.2g, drawn positive %.2g, agreement %.2g\n",
                  max( abs( num( mm$below_zero_x ) - mm$below.zero.x ), na.rm = TRUE ),
                  max( abs( num( mm$positive_x ) - mm$positive.x ), na.rm = TRUE ),
                  suppressWarnings( max( abs( num( mm$agreement.x ) - mm$agreement.y ), na.rm = TRUE ) ) ) )
    if ( any( !same.alg ) )
        print( mm[ !same.alg, c( "gate", "type.x", "algorithm_x", "algorithm_y", "type.y", "algorithm.x",
                                 "algorithm.y", "agreement.x", "agreement.y" ) ], row.names = FALSE )
    cat( sprintf( "   overlapping populations (R): %d\n", nrow( conv$overlaps ) ) )
}
