# flowjo_workspace.R
#
# Reads FlowJo 10 workspaces (.wsp) and workspace templates (.wspt) and
# matches the parameters their gates use to the experiment's channels, as
# Honeychrome's Gating Model Builder does (ag_flowjo).  The gates are
# converted to gate definitions in flowjo_import.R.
#
# A parsed workspace is a list:
#   path, flowjo.version, cytometer
#   affixes   list of c(prefix, suffix) of the compensation or unmixing
#             matrices (e.g. "Comp-"), plus "Comp-", "Unmixed-", "Spectral-"
#   sources   list of gating sources: each group with gates, then each
#             sample; a source has kind ("group" or "sample"), name, roots
#             (top populations), labels ($PnN -> $PnS), transforms ($PnN ->
#             list(kind, params)), count, sample.name (a group's only
#             sample, whose counts the group's counts are)
#
# A population has name, path (names from the top), key (path joined with
# "/"), node.type ("Population", "AndNode", "OrNode", "NotNode"), gate,
# count and children.  A gate has kind ("rectangle", "polygon", "ellipse"
# or the element name of anything else), channels (FlowJo parameters),
# bounds (rectangles: c(min, max) per dimension, NA open-ended), vertices
# (polygons), foci and edge (ellipses), inside (FALSE when the gate selects
# the events outside its shape) and quad.id.  Coordinates are raw
# (untransformed) parameter units, as FlowJo stores them.
#
# Parameter matching (align.flowjo.parameters) differs from Honeychrome in
# two details: fluorophores and markers are recognised with AutoSpectral's
# match.fluorophores and match.markers (same databases, same
# longest-match rule), and the last-resort closest name uses a normalised
# edit distance (1 - distance / longer length >= 0.8) in place of Python's
# difflib ratio.


flowjo.population.nodes <- c( "Population", "AndNode", "OrNode", "NotNode" )
flowjo.default.prefixes <- c( "Comp-", "Unmixed-", "Spectral-" )


#' Local name of an XML tag or attribute
#' @title FlowJo Local Name
#' @param x Character vector.
#' @return Names without namespace prefix.
#' @keywords internal
fj.local <- function( x ) sub( "^.*[:}]", "", x )


#' Attribute of an XML node by local name
#' @title FlowJo Attribute
#' @param node XML node.
#' @param name Local attribute name.
#' @param default Returned when the attribute is absent.
#' @return Character value.
#' @keywords internal
fj.attr <- function( node, name, default = NA_character_ ) {
    a <- xml2::xml_attrs( node )
    hit <- which( fj.local( names( a ) ) == name )
    if ( length( hit ) ) unname( a[ hit[ 1L ] ] ) else default
}


#' Child elements by local name
#' @title FlowJo Children
#' @param node XML node.
#' @param name Local element name.
#' @return Node set.
#' @keywords internal
fj.children <- function( node, name ) {
    kids <- xml2::xml_children( node )
    kids[ xml2::xml_name( kids ) == name ]
}


#' First child element by local name
#' @title FlowJo Child
#' @inheritParams fj.children
#' @return Node, or \code{NULL}.
#' @keywords internal
fj.child <- function( node, name ) {
    if ( is.null( node ) ) return( NULL )
    kids <- fj.children( node, name )
    if ( length( kids ) ) kids[[ 1L ]] else NULL
}


#' Elements by local name at any depth
#' @title FlowJo Descendants
#' @inheritParams fj.children
#' @return Node set, including \code{node} itself when it matches.
#' @keywords internal
fj.descendants <- function( node, name ) {
    xml2::xml_find_all( node, sprintf( "descendant-or-self::*[local-name()='%s']", name ) )
}


#' Number from an attribute value
#' @title FlowJo Number
#' @param v Character value.
#' @return Finite number, or \code{NA}.
#' @keywords internal
fj.num <- function( v ) {
    x <- suppressWarnings( as.numeric( v ) )
    if ( length( x ) == 1L && is.finite( x ) ) x else NA_real_
}


#' Integer from an attribute value
#' @title FlowJo Integer
#' @param v Character value.
#' @return Integer, or \code{NA}.
#' @keywords internal
fj.int <- function( v ) {
    x <- suppressWarnings( as.integer( v ) )
    if ( length( x ) == 1L ) x else NA_integer_
}


#' Read a FlowJo workspace
#'
#' @title Read FlowJo Workspace
#'
#' @description
#' Parses a FlowJo 10 workspace (\code{.wsp}) or workspace template
#' (\code{.wspt}): every group with gates and every sample becomes a gating
#' source, with its parameter labels and FlowJo display transforms.  The
#' structure is described in the file header of \code{flowjo_workspace.R}.
#'
#' @param path Path to the workspace.
#'
#' @return Parsed workspace (named list).
#'
#' @seealso \code{\link{import.flowjo.workspace}}
#'
#' @export
read.flowjo.workspace <- function( path ) {
    root <- xml2::read_xml( path )
    if ( xml2::xml_name( root ) != "Workspace" )
        stop( sprintf( "Not a FlowJo workspace: the root element is <%s>, not <Workspace>.",
                       xml2::xml_name( root ) ) )
    ws <- list( path = path, flowjo.version = fj.attr( root, "flowJoVersion", "" ),
                cytometer = "", affixes = fj.matrix.affixes( root ), sources = list() )

    cyt.transforms <- list()
    for ( cyt in fj.children( fj.child( root, "Cytometers" ) %||% root, "Cytometer" ) ) {
        if ( !nzchar( ws$cytometer ) )
            ws$cytometer <- fj.attr( cyt, "cyt", "" )
        store <- fj.child( cyt, "TransformStore" )
        if ( !is.null( store ) )
            for ( el in fj.descendants( store, "Transforms" ) ) {
                xf <- fj.parse.transforms( el )
                cyt.transforms[ names( xf ) ] <- xf
            }
    }

    samples <- list()
    sample.list <- fj.child( root, "SampleList" )
    for ( s in if ( is.null( sample.list ) ) list() else fj.children( sample.list, "Sample" ) ) {
        node <- fj.child( s, "SampleNode" )
        if ( is.null( node ) )
            next
        sid <- fj.attr( node, "sampleID" )
        if ( is.na( sid ) )
            sid <- fj.attr( fj.child( s, "DataSet" ) %||% node, "sampleID" )
        transforms <- cyt.transforms
        xf.el <- fj.child( s, "Transformations" )
        if ( !is.null( xf.el ) ) {
            xf <- fj.parse.transforms( xf.el )
            transforms[ names( xf ) ] <- xf
        }
        name <- fj.attr( node, "name", "" )
        samples[[ sid ]] <- list(
            kind = "sample", name = if ( has.value( name ) ) name else sprintf( "Sample %s", sid ),
            roots = fj.subpopulations( node, character( 0 ) ), sample.id = sid,
            labels = fj.parameter.labels( fj.child( s, "Keywords" ) ),
            transforms = transforms, count = fj.int( fj.attr( node, "count" ) ),
            sample.name = NA_character_ )
    }

    groups <- list()
    groups.el <- fj.child( root, "Groups" )
    for ( g in if ( is.null( groups.el ) ) list() else fj.children( groups.el, "GroupNode" ) ) {
        roots <- fj.subpopulations( g, character( 0 ) )
        if ( length( roots ) == 0L )
            next
        refs <- character( 0 )
        sr <- fj.child( fj.child( g, "Group" ), "SampleRefs" )
        if ( !is.null( sr ) )
            refs <- vapply( fj.children( sr, "SampleRef" ), fj.attr, character( 1 ), name = "sampleID" )
        members <- samples[ refs[ refs %in% names( samples ) ] ]
        member <- if ( length( members ) ) members[[ 1L ]] else NULL
        single <- if ( length( members ) == 1L ) member else NULL
        name <- fj.attr( g, "name", "" )
        groups[[ length( groups ) + 1L ]] <- list(
            kind = "group", name = if ( has.value( name ) ) name else "Group", roots = roots,
            sample.id = NA_character_, labels = member$labels %||% character( 0 ),
            transforms = if ( is.null( member ) ) cyt.transforms else member$transforms,
            count = if ( is.null( single ) ) NA_integer_ else single$count,
            sample.name = if ( is.null( single ) ) NA_character_ else single$name )
    }
    ws$sources <- c( groups, unname( samples ) )
    ws
}


#' Prefixes and suffixes of the workspace's matrices
#' @title FlowJo Matrix Affixes
#' @param root Workspace root node.
#' @return List of \code{c(prefix, suffix)}.
#' @keywords internal
fj.matrix.affixes <- function( root ) {
    out <- list()
    add <- function( pair ) {
        if ( !any( vapply( out, identical, logical( 1 ), pair ) ) )
            out[[ length( out ) + 1L ]] <<- pair
    }
    for ( el in fj.descendants( root, "spilloverMatrix" ) ) {
        pair <- c( fj.attr( el, "prefix", "" ), fj.attr( el, "suffix", "" ) )
        pair[ is.na( pair ) ] <- ""
        if ( any( nzchar( pair ) ) )
            add( pair )
    }
    for ( p in flowjo.default.prefixes )
        add( c( p, "" ) )
    out
}


#' FlowJo transforms of a Transformations element
#' @title FlowJo Transforms
#' @param el Transformations or Transforms node.
#' @return Named list parameter -> \code{list(kind, params)}.
#' @keywords internal
fj.parse.transforms <- function( el ) {
    out <- list()
    for ( xf in xml2::xml_children( el ) ) {
        a <- xml2::xml_attrs( xf )
        names( a ) <- fj.local( names( a ) )
        for ( p in fj.descendants( xf, "parameter" ) ) {
            name <- fj.attr( p, "name" )
            if ( has.value( name ) )
                out[[ name ]] <- list( kind = xml2::xml_name( xf ), params = as.list( a ) )
        }
    }
    out
}


#' Parameter labels of a sample
#' @title FlowJo Parameter Labels
#' @param keywords Keywords node, or \code{NULL}.
#' @return Named character vector $PnN -> $PnS.
#' @keywords internal
fj.parameter.labels <- function( keywords ) {
    if ( is.null( keywords ) )
        return( character( 0 ) )
    kws <- fj.children( keywords, "Keyword" )
    keys <- vapply( kws, fj.attr, character( 1 ), name = "name", default = "" )
    vals <- vapply( kws, fj.attr, character( 1 ), name = "value", default = "" )
    vals[ is.na( vals ) ] <- ""
    labels <- character( 0 )
    for ( i in grep( "^\\$P[0-9]+N$", keys ) ) {
        s <- vals[ match( sub( "N$", "S", keys[ i ] ), keys ) ]
        labels[ vals[ i ] ] <- trimws( if ( is.na( s ) ) "" else s )
    }
    labels
}


#' Populations below a node
#' @title FlowJo Subpopulations
#' @param node SampleNode, GroupNode or population node.
#' @param path Names of the populations above.
#' @return List of populations.
#' @keywords internal
fj.subpopulations <- function( node, path ) {
    sub <- fj.child( node, "Subpopulations" )
    if ( is.null( sub ) )
        return( list() )
    out <- list()
    for ( child in xml2::xml_children( sub ) ) {
        tag <- xml2::xml_name( child )
        if ( !( tag %in% flowjo.population.nodes ) )
            next
        name <- fj.attr( child, "name", "" )
        if ( is.na( name ) ) name <- ""
        p <- c( path, name )
        pop <- list( name = name, path = p, key = paste( p, collapse = "/" ), node.type = tag,
                     gate = NULL, count = fj.int( fj.attr( child, "count" ) ) )
        if ( tag == "Population" ) {
            g <- fj.child( child, "Gate" )
            if ( !is.null( g ) )
                pop$gate <- fj.parse.gate( g )
        }
        pop$children <- fj.subpopulations( child, p )
        out[[ length( out ) + 1L ]] <- pop
    }
    out
}


#' Vertices of a polygon, foci or edge element
#' @title FlowJo Vertices
#' @param el XML node, or \code{NULL}.
#' @return n x 2 matrix.
#' @keywords internal
fj.vertices <- function( el ) {
    if ( is.null( el ) )
        return( matrix( numeric( 0 ), 0L, 2L ) )
    rows <- lapply( fj.children( el, "vertex" ), function( v ) {
        coords <- vapply( fj.children( v, "coordinate" ),
                          function( k ) fj.num( fj.attr( k, "value" ) ), numeric( 1 ) )
        if ( length( coords ) >= 2L && !anyNA( coords[ 1:2 ] ) ) coords[ 1:2 ] else NULL
    } )
    rows <- rows[ lengths( rows ) > 0L ]
    if ( length( rows ) == 0L )
        return( matrix( numeric( 0 ), 0L, 2L ) )
    do.call( rbind, rows )
}


#' A population's gate
#' @title FlowJo Gate
#' @param gate.el Gate node.
#' @return Gate (see the file header of \code{flowjo_workspace.R}).
#' @keywords internal
fj.parse.gate <- function( gate.el ) {
    kids <- xml2::xml_children( gate.el )
    shapes <- kids[ endsWith( xml2::xml_name( kids ), "Gate" ) ]
    gid <- fj.attr( gate.el, "id", "" )
    if ( length( shapes ) == 0L )
        return( list( kind = "none", channels = character( 0 ), inside = TRUE, gate.id = gid ) )
    shape <- shapes[[ 1L ]]
    tag <- xml2::xml_name( shape )
    dims <- fj.children( shape, "dimension" )
    dim.name <- function( d ) {
        el <- xml2::xml_find_first( d, ".//*[local-name()='fcs-dimension' or local-name()='new-dimension']" )
        if ( inherits( el, "xml_missing" ) ) "" else fj.attr( el, "name", "" )
    }
    qid <- fj.int( fj.attr( shape, "quadId" ) )
    g <- list(
        kind = switch( tag, RectangleGate = "rectangle", PolygonGate = "polygon",
                       EllipsoidGate = "ellipse", tag ),
        channels = vapply( dims, dim.name, character( 1 ) ),
        inside = !identical( fj.attr( shape, "eventsInside", "1" ), "0" ),
        quad.id = if ( !is.na( qid ) && qid >= 0L ) qid else NA_integer_,
        gate.id = gid )
    if ( g$kind == "rectangle" )
        g$bounds <- lapply( dims, function( d ) c( fj.num( fj.attr( d, "min" ) ),
                                                   fj.num( fj.attr( d, "max" ) ) ) )
    else if ( g$kind == "polygon" )
        g$vertices <- fj.vertices( shape )
    else if ( g$kind == "ellipse" ) {
        g$foci <- fj.vertices( fj.child( shape, "foci" ) )
        g$edge <- fj.vertices( fj.child( shape, "edge" ) )
    }
    g
}


#' Populations of a gating source, parents first
#' @title FlowJo Source Populations
#' @param source Gating source.
#' @return Flat list of populations.
#' @keywords internal
fj.populations <- function( source ) {
    out <- list()
    visit <- function( pops ) for ( p in pops ) {
        out[[ length( out ) + 1L ]] <<- p
        visit( p$children )
    }
    visit( source$roots )
    out
}


#' Default gating source of a workspace
#'
#' @title Default FlowJo Source
#'
#' @description
#' The group with the most populations, else the sample with the most.
#'
#' @param ws Parsed workspace.
#'
#' @return Gating source, or \code{NULL}.
#'
#' @keywords internal
fj.default.source <- function( ws ) {
    n <- vapply( ws$sources, function( s ) length( fj.populations( s ) ), integer( 1 ) )
    groups <- which( vapply( ws$sources, `[[`, character( 1 ), "kind" ) == "group" & n > 0L )
    if ( length( groups ) )
        return( ws$sources[[ groups[ which.max( n[ groups ] ) ] ]] )
    if ( any( n > 0L ) ) ws$sources[[ which.max( n ) ]] else NULL
}


#' FlowJo parameters a source's gates use
#' @title FlowJo Gated Parameters
#' @param source Gating source.
#' @return Character vector, in order of first use.
#' @keywords internal
fj.gated.parameters <- function( source ) {
    unique( unlist( lapply( fj.populations( source ), function( p ) p$gate$channels ) ) )
}


#' Parameter name without its matrix prefix or suffix
#' @title FlowJo Parameter Stem
#' @param name Parameter name.
#' @param affixes List of \code{c(prefix, suffix)}.
#' @return Character.
#' @keywords internal
fj.parameter.stem <- function( name, affixes = NULL ) {
    stem <- as.character( name )
    if ( is.null( affixes ) )
        affixes <- lapply( flowjo.default.prefixes, function( p ) c( p, "" ) )
    for ( a in affixes ) {
        if ( nzchar( a[ 1 ] ) && startsWith( stem, a[ 1 ] ) )
            stem <- substring( stem, nchar( a[ 1 ] ) + 1L )
        if ( nzchar( a[ 2 ] ) && endsWith( stem, a[ 2 ] ) )
            stem <- substring( stem, 1L, nchar( stem ) - nchar( a[ 2 ] ) )
    }
    stem
}


#' The experiment side of a FlowJo import
#'
#' @title FlowJo Experiment Channels
#'
#' @description
#' Channels are the panel antigens (the data's column names); each has its
#' dye ($PnN), its type and its display transform.  The singlet pair is the
#' forward-scatter area and height (or width) channels.
#'
#' @param panel Panel \code{data.frame} (\code{dye}, \code{antigen},
#'   \code{type}).
#' @param transforms Named list antigen -> transform parameters.
#' @param singlet.pair Character \code{c(area, height)} of dyes or
#'   antigens, e.g. \code{ap$singlet.scatter.parameter}.
#'
#' @return List with channels, dyes, fluorescence, scatter, transforms and
#'   singlet.pair.
#'
#' @keywords internal
fj.experiment <- function( panel, transforms = list(), singlet.pair = NULL ) {
    type <- as.character( panel$type )
    to.antigen <- function( x ) ifelse( x %in% panel$dye, panel$antigen[ match( x, panel$dye ) ], x )
    pair <- if ( length( singlet.pair ) == 2L ) to.antigen( singlet.pair ) else NULL
    list(
        channels     = panel$antigen,
        dyes         = stats::setNames( panel$dye, panel$antigen ),
        fluorescence = panel$antigen[ type %in% c( "lineage", "activation" ) ],
        scatter      = panel$antigen[ type == "scatter" ],
        transforms   = transforms,
        singlet.pair = if ( !is.null( pair ) && all( pair %in% panel$antigen ) ) pair else NULL )
}


#' Match FlowJo parameters to the experiment's channels
#'
#' @title Align FlowJo Parameters
#'
#' @description
#' Proposes a panel channel (antigen) for each FlowJo parameter, in this
#' order:
#' \enumerate{
#'   \item the same name after removing the matrix prefix (e.g.
#'     \code{"Comp-"}), against the panel dyes and antigens ("name");
#'   \item the same fluorophore on an area parameter ("fluorophore"); the
#'     FlowJo $PnS label is checked against the channel's antigen and a
#'     disagreement is reported ("fluorophore_label_differs");
#'   \item the $PnS label against the antigens, marker by marker for
#'     combined labels such as \code{"CD4_IgM"} ("antigen");
#'   \item the closest dye name ("fuzzy").
#' }
#' Fluorescence height and width parameters are not matched.  Matches other
#' than "name" and "fluorophore" should be checked.
#'
#' @param params FlowJo parameter names.
#' @param labels Named character vector $PnN -> $PnS.
#' @param exp Experiment channels (see \code{fj.experiment}).
#' @param affixes Matrix prefixes and suffixes of the workspace.
#'
#' @return \code{data.frame} with parameter, channel (antigen or
#'   \code{NA}), method, note and label.
#'
#' @importFrom AutoSpectral match.fluorophores match.markers
#'
#' @keywords internal
align.flowjo.parameters <- function( params, labels, exp, affixes = NULL ) {
    fl.db <- autospectral.database( "fluorophore_database.csv" )
    mk.db <- autospectral.database( "marker_database.csv" )
    norm <- function( x ) tolower( gsub( "[[:space:]_.+-]+", "", as.character( x ) ) )
    canon.fluor <- function( x ) {
        if ( is.null( fl.db ) ) return( NA_character_ )
        m <- unname( unlist( match.fluorophores( x, fl.db, verbose = FALSE ) ) )[ 1L ]
        if ( is.null( m ) || is.na( m ) || m %in% c( "No match", "Negative", "AF" ) ) NA_character_ else m
    }
    marker.key <- function( x ) {
        m <- if ( is.null( mk.db ) ) "" else unname( unlist( match.markers( x, mk.db, verbose = FALSE ) ) )[ 1L ]
        norm( if ( is.null( m ) || is.na( m ) || !nzchar( m ) ) x else m )
    }
    tokens.of <- function( label ) {
        t <- trimws( strsplit( as.character( label ), "[_/,;|]+" )[[ 1L ]] )
        t[ nzchar( t ) ]
    }

    exp.fluor <- list()
    for ( ch in exp$fluorescence ) {
        canon <- canon.fluor( exp$dyes[[ ch ]] )
        k <- norm( if ( is.na( canon ) ) exp$dyes[[ ch ]] else canon )
        exp.fluor[[ k ]] <- c( exp.fluor[[ k ]], ch )
    }
    exp.marker <- list()
    for ( ch in exp$fluorescence ) {
        if ( identical( ch, exp$dyes[[ ch ]] ) ) next
        k <- marker.key( ch )
        exp.marker[[ k ]] <- c( exp.marker[[ k ]], ch )
    }
    scatter.stems <- sub( "-[^-]*$", "", exp$dyes[ exp$scatter ] )
    pool <- exp$dyes[ c( exp$fluorescence, exp$scatter ) ]

    rows <- lapply( params, function( p ) {
        label <- if ( p %in% names( labels ) ) labels[[ p ]] else ""
        if ( is.na( label ) ) label <- ""
        stem <- fj.parameter.stem( p, affixes )
        channel <- NA_character_
        method <- NA_character_
        note <- ""
        tokens <- vapply( tokens.of( label ), marker.key, character( 1 ) )
        by.name <- c( stem, p )[ c( stem, p ) %in% c( exp$dyes, exp$channels ) ]
        if ( length( by.name ) ) {
            n <- by.name[ 1L ]
            channel <- if ( n %in% exp$dyes ) names( exp$dyes )[ match( n, exp$dyes ) ] else n
            method <- "name"
        } else if ( grepl( "-[HW]$", stem ) && !any( startsWith( stem, scatter.stems ) ) ) {
            note <- "Height and width parameters are not unmixed."
        } else {
            base <- sub( "-A$", "", stem )
            canon <- canon.fluor( base )
            hits <- if ( is.na( canon ) ) character( 0 ) else exp.fluor[[ norm( canon ) ]]
            if ( length( hits ) == 1L ) {
                channel <- hits
                method <- "fluorophore"
                if ( length( tokens ) && !identical( channel, exp$dyes[[ channel ]] ) &&
                     !( marker.key( channel ) %in% tokens ) ) {
                    method <- "fluorophore_label_differs"
                    note <- sprintf( "FlowJo label '%s', experiment antigen '%s'.", label, channel )
                }
            } else if ( length( hits ) > 1L ) {
                channel <- hits[ 1L ]
                method <- "fuzzy"
                note <- paste( "Fluorophore is on more than one channel:", paste( hits, collapse = ", " ) )
            }
            if ( is.na( channel ) && length( tokens ) ) {
                found <- unique( unlist( exp.marker[ tokens[ tokens %in% names( exp.marker ) ] ] ) )
                if ( length( found ) == 1L ) {
                    channel <- found
                    method <- "antigen"
                    note <- sprintf( "Matched on the label '%s'; fluorophore%s not found.", label,
                                     if ( is.na( canon ) ) "" else sprintf( " (%s)", canon ) )
                } else if ( length( found ) > 1L ) {
                    channel <- found[ 1L ]
                    method <- "fuzzy"
                    note <- sprintf( "Label '%s' matches %s", label, paste( found, collapse = ", " ) )
                }
            }
            if ( is.na( channel ) && length( pool ) ) {
                sim <- 1 - c( utils::adist( base, pool ) ) / pmax( nchar( base ), nchar( pool ) )
                if ( max( sim ) >= 0.8 ) {
                    channel <- names( pool )[ which.max( sim ) ]
                    method <- "fuzzy"
                    note <- "Closest channel name only."
                }
            }
        }
        data.frame( parameter = p, channel = channel, method = method, note = note,
                    label = label, stringsAsFactors = FALSE )
    } )
    if ( length( rows ) == 0L )
        return( data.frame( parameter = character( 0 ), channel = character( 0 ),
                            method = character( 0 ), note = character( 0 ),
                            label = character( 0 ), stringsAsFactors = FALSE ) )
    do.call( rbind, rows )
}


#' An AutoSpectral database
#' @title AutoSpectral Database
#' @param file File name in AutoSpectral's extdata.
#' @return \code{data.frame}, or \code{NULL} when not found.
#' @keywords internal
autospectral.database <- function( file ) {
    path <- system.file( "extdata", file, package = "AutoSpectral" )
    if ( !nzchar( path ) ) return( NULL )
    utils::read.csv( path, stringsAsFactors = FALSE, check.names = FALSE )
}


#' Populations a FlowJo source can import
#'
#' @title Plan FlowJo Import
#'
#' @description
#' A population is dropped when its gate type is not supported (Boolean
#' gates, gates selecting the outside of a shape, incomplete shapes) or a
#' parameter it uses maps to no channel; every population below it is
#' dropped with it.
#'
#' @param source Gating source.
#' @param mapping Named character vector FlowJo parameter -> channel
#'   (\code{NA}: none).
#'
#' @return List with \code{keep} (population keys) and \code{dropped}
#'   (named character vector key -> reason).
#'
#' @keywords internal
fj.plan.import <- function( source, mapping ) {
    keep <- character( 0 )
    dropped <- character( 0 )
    unsupported <- function( pop ) {
        if ( pop$node.type != "Population" )
            return( sprintf( "Boolean gates (%s) are not imported.", sub( "Node$", "", pop$node.type ) ) )
        g <- pop$gate
        if ( is.null( g ) || g$kind == "none" ) return( "The population has no gate." )
        if ( !( g$kind %in% c( "rectangle", "polygon", "ellipse" ) ) )
            return( sprintf( "%s gates are not supported.", g$kind ) )
        if ( !g$inside ) return( "Gates that select events outside the shape are not supported." )
        if ( g$kind == "polygon" && ( nrow( g$vertices ) < 3L || length( g$channels ) != 2L ) )
            return( "The polygon is incomplete." )
        if ( g$kind == "ellipse" && ( nrow( g$edge ) < 4L || length( g$channels ) != 2L ) )
            return( "The ellipse is incomplete." )
        if ( g$kind == "rectangle" && all( is.na( unlist( g$bounds ) ) ) )
            return( "The rectangle has no bounds." )
        NULL
    }
    visit <- function( pop, parent.dropped ) {
        reason <- if ( !is.null( parent.dropped ) ) sprintf( "its parent '%s' was dropped.", parent.dropped )
                  else unsupported( pop )
        if ( is.null( reason ) ) {
            ch <- pop$gate$channels
            missing <- ch[ !( ch %in% names( mapping ) ) | is.na( mapping[ ch ] ) ]
            if ( length( missing ) )
                reason <- paste0( "no channel for ", paste( missing, collapse = ", " ), "." )
        }
        if ( is.null( reason ) ) {
            keep <<- c( keep, pop$key )
            for ( k in pop$children ) visit( k, NULL )
        } else {
            dropped[ pop$key ] <<- reason
            top <- parent.dropped %||% pop$key
            for ( k in pop$children ) visit( k, top )
        }
    }
    for ( r in source$roots ) visit( r, NULL )
    list( keep = keep, dropped = dropped )
}


#' Messages about dropped populations
#' @title FlowJo Import Messages
#' @param plan Result of \code{fj.plan.import}.
#' @return One line per dropped branch.
#' @keywords internal
fj.plan.messages <- function( plan ) {
    keys <- names( plan$dropped )
    tops <- keys[ !vapply( keys, function( k ) any( startsWith( k, paste0( setdiff( keys, k ), "/" ) ) ),
                           logical( 1 ) ) ]
    vapply( tops, function( k ) {
        n <- sum( startsWith( keys, paste0( k, "/" ) ) )
        sprintf( "Dropped '%s'%s: %s", k,
                 if ( n ) sprintf( " and %d population(s) below it", n ) else "", plan$dropped[[ k ]] )
    }, character( 1 ), USE.NAMES = FALSE )
}
