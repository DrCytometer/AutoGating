# transform.R
#
# Gating-ML 2.0 scale transformations used for display and gating:
# linear (flin), logicle (flogicle) and logarithmic (flog). Every transform
# maps raw intensities onto a common display scale on which the top of the
# instrument range (T) is 1, the same scale as Honeychrome and FlowKit, so
# gate boundaries are portable between the two.
#
# Transform parameters are kept as a named list per channel:
#   id         0 = linear, 1 = logicle, 2 = log, NULL = identity (e.g. time)
#   scale.t    T, top of the data range in raw units
#   linear.a   A of flin (offset)
#   logicle.w  W of flogicle (linearisation width, decades)
#   logicle.m  M of flogicle (decades of the full display)
#   logicle.a  A of flogicle (additional negative decades)
#   log.m      M of flog (decades)
#
# The logicle transform is AutoSpectral::biexp.transform() (Parks, Roederer
# & Moore, Cytometry A 2006) with channelRange = 1, which puts it on the
# Gating-ML 0 to 1 scale: T = maxValue, M = pos, W = log10(-widthBasis) / 2.
# biexp.transform() supports A = 0 only; a logicle with A != 0 stops with
# its error message.


#' Build forward and inverse functions for one channel transform
#'
#' @title Make a Channel Transform
#'
#' @description
#' Returns the forward (raw to display) and inverse (display to raw)
#' functions for a Gating-ML 2.0 linear, logicle or log transform, or the
#' identity when \code{params} is \code{NULL} or has no recognised
#' \code{id}.
#'
#' @param params Named list of transform parameters (see the file header of
#'   \code{transform.R}): \code{id}, \code{scale.t}, and \code{linear.a},
#'   \code{logicle.w}/\code{logicle.m}/\code{logicle.a} or \code{log.m}.
#'
#' @return A list with elements \code{forward} and \code{inverse}, each a
#'   vectorised function of a numeric vector.
#'
#' @seealso \code{\link{read.transformation}},
#'   \code{\link[AutoSpectral]{biexp.transform}}
#'
#' @importFrom AutoSpectral biexp.transform
#'
#' @keywords internal
make.transform <- function( params ) {

    identity.fn <- function( x ) x
    id <- if ( is.null( params ) ) NULL else params$id

    if ( is.null( id ) || !( id %in% 0:2 ) )
        return( list( forward = identity.fn, inverse = identity.fn ) )

    t.top <- as.numeric( params$scale.t )
    stopifnot( length( t.top ) == 1L, t.top > 0 )

    if ( id == 0 ) {
        a <- as.numeric( params$linear.a %||% 0 )
        return( list(
            forward = function( x ) ( x + a ) / ( t.top + a ),
            inverse = function( y ) y * ( t.top + a ) - a
        ) )
    }

    if ( id == 2 ) {
        m <- as.numeric( params$log.m )
        return( list(
            forward = function( x ) suppressWarnings( log10( x / t.top ) / m + 1 ),
            inverse = function( y ) t.top * 10^( ( y - 1 ) * m )
        ) )
    }

    w <- as.numeric( params$logicle.w )
    logicle.args <- list( channelRange = 1, maxValue = t.top,
        pos = as.numeric( params$logicle.m ),
        neg = as.numeric( params$logicle.a %||% 0 ),
        widthBasis = -10^( 2 * w ) )
    list(
        forward = do.call( biexp.transform, c( logicle.args, inverse = FALSE ) ),
        inverse = do.call( biexp.transform, c( logicle.args, inverse = TRUE ) )
    )
}


#' Convert FlowJo biexponential settings to Gating-ML logicle parameters
#'
#' @title FlowJo Biex to Logicle Parameters
#'
#' @description
#' Converts FlowJo-style biex settings (the arguments of
#' \code{\link[AutoSpectral]{biexp.transform}}) to Gating-ML logicle
#' parameters: T = \code{max.range}, M = \code{pos}, A = \code{neg} and
#' W = \eqn{\log_{10}(-\mathrm{width}) / 2}, on the 0 to 1 display scale
#' instead of \code{length} channels.
#'
#' @param max.range,pos,neg,width FlowJo biex settings (as in the legacy
#'   transformation CSV).
#'
#' @return Named list of logicle transform parameters.
#'
#' @keywords internal
flowjo.biex.to.logicle <- function( max.range, pos, neg, width ) {
    w <- if ( width < 0 ) log10( -width ) / 2 else 0
    list( id = 1L, scale.t = max.range, logicle.w = w,
          logicle.m = pos, logicle.a = neg )
}


#' Test whether two channel transforms are the same
#'
#' @title Same Transform
#'
#' @param a,b Transform parameter lists (see \code{\link{make.transform}}).
#' @param tol Relative tolerance for numeric parameters.
#'
#' @return \code{TRUE} when every parameter matches (display limits are
#'   ignored); \code{FALSE} when either list is empty.
#'
#' @keywords internal
same.transform <- function( a, b, tol = 1e-9 ) {
    if ( length( a ) == 0L || length( b ) == 0L )
        return( FALSE )
    for ( key in c( "id", "scale.t", "linear.a", "logicle.w", "logicle.m",
                    "logicle.a", "log.m" ) ) {
        va <- a[[ key ]]
        vb <- b[[ key ]]
        if ( is.numeric( va ) && is.numeric( vb ) && length( va ) == 1L && length( vb ) == 1L ) {
            if ( abs( va - vb ) > tol * max( 1, abs( va ) ) )
                return( FALSE )
        } else if ( !identical( va, vb ) ) {
            return( FALSE )
        }
    }
    TRUE
}


#' Add points along polygon edges
#'
#' @title Densify Polygon
#'
#' @description
#' Inserts evenly spaced points along each edge, so that a polygon keeps its
#' shape when its vertices are moved to a different (non-linear) axis
#' transform.
#'
#' @param coords Polygon vertices (n x 2 matrix).
#' @param points.per.edge Points per edge, including its first vertex.
#'
#' @return Matrix with \code{(n - 1) * points.per.edge + 1} rows.
#'
#' @keywords internal
densify.polygon <- function( coords, points.per.edge = 64L ) {
    n <- nrow( coords )
    if ( n < 2L )
        return( coords )
    t <- seq( 0, 1, length.out = points.per.edge + 1L )[ seq_len( points.per.edge ) ]
    idx <- rep( seq_len( n - 1L ), each = points.per.edge )
    a <- coords[ idx, , drop = FALSE ]
    b <- coords[ idx + 1L, , drop = FALSE ]
    unname( rbind( a + rep( t, times = n - 1L ) * ( b - a ), coords[ n, , drop = FALSE ] ) )
}


#' Convert gate boundaries between channel transforms
#'
#' @title Remap Boundaries
#'
#' @description
#' Moves thresholds, ranges and polygon vertices from the transforms they
#' were calculated in to other transforms (raw units in between).
#' Thresholds convert exactly; polygons of polygon gates and shapes as drawn
#' are first densified along each edge.  Channels whose transforms match, or
#' that lack parameters on either side, are left unchanged.
#'
#' @param gate.defs List of gate definitions.
#' @param boundaries Named list gate -> population -> boundary entry.
#' @param src.params,dst.params Named lists channel -> transform parameters.
#' @param points.per.edge Densification of polygon edges.
#'
#' @return List with \code{boundaries} and \code{report} (named character
#'   vector channel -> "same", "remapped" or "unknown").
#'
#' @keywords internal
remap.boundaries <- function( gate.defs, boundaries, src.params, dst.params,
                              points.per.edge = 64L ) {
    report <- character( 0 )
    cache <- list()
    converter <- function( ch ) {
        if ( ch %in% names( report ) )
            return( cache[[ ch ]] )
        sp <- src.params[[ ch ]]
        dp <- dst.params[[ ch ]]
        if ( length( sp ) == 0L || length( dp ) == 0L ) {
            report[ ch ] <<- "unknown"
            cache[ ch ] <<- list( NULL )
        } else if ( same.transform( sp, dp ) ) {
            report[ ch ] <<- "same"
            cache[ ch ] <<- list( NULL )
        } else {
            report[ ch ] <<- "remapped"
            src <- make.transform( sp )
            dst <- make.transform( dp )
            cache[[ ch ]] <<- function( v ) dst$forward( src$inverse( as.numeric( v ) ) )
        }
        cache[[ ch ]]
    }

    by.name <- gates.by.name( gate.defs )
    for ( name in names( boundaries ) ) {
        gdef <- by.name[[ name ]]
        if ( is.null( gdef ) )
            next
        eff <- effective.gate.def( gdef, by.name )
        chans <- gate.channels( eff )
        entries <- boundaries[[ name ]]
        for ( p in names( entries ) ) {
            e <- entries[[ p ]]
            shape <- unlist( e[[ "shape.channels" ]] )
            shape <- shape[ vapply( shape, has.value, logical( 1 ) ) ]
            convs <- lapply( if ( length( shape ) ) shape else chans, converter )
            if ( all( vapply( convs, is.null, logical( 1 ) ) ) )
                next
            for ( k in 1:2 ) {
                key <- c( "threshold.x", "threshold.y" )[ k ]
                if ( k <= length( convs ) && !is.null( convs[[ k ]] ) && has.value( e[[ key ]] ) )
                    e[[ key ]] <- convs[[ k ]]( e[[ key ]] )
            }
            if ( !is.null( e[[ "range" ]] ) && !is.null( convs[[ 1L ]] ) )
                e[[ "range" ]] <- convs[[ 1L ]]( e[[ "range" ]] )
            b <- e[[ "boundary" ]]
            if ( is.matrix( b ) && ncol( b ) == 2L && nrow( b ) > 0L ) {
                if ( isTRUE( eff$gate.type %in% ag.polygon.types ) || length( shape ) )
                    b <- densify.polygon( b, points.per.edge )
                for ( k in seq_len( min( 2L, length( convs ) ) ) )
                    if ( !is.null( convs[[ k ]] ) )
                        b[ , k ] <- convs[[ k ]]( b[ , k ] )
                e[[ "boundary" ]] <- b
            }
            entries[[ p ]] <- e
        }
        boundaries[[ name ]] <- entries
    }
    list( boundaries = boundaries, report = report )
}


#' Convert a model's boundaries to other channel transforms
#'
#' @title Remap Model
#'
#' @description
#' Converts the trained boundaries and the shapes as drawn of a gating model
#' (see \code{\link{read.agmodel}}) from the model's own transforms to
#' \code{dst.params}, the transforms of the data it will be applied to, as
#' Honeychrome does when a model is opened in another experiment.  The
#' returned model's \code{transforms} describe the converted boundaries.
#'
#' @param model Gating model.
#' @param dst.params Named list channel -> transform parameters.
#'
#' @return List with \code{model} and \code{messages} (character).
#'
#' @export
remap.model <- function( model, dst.params ) {
    src <- model$transforms
    if ( length( src ) == 0L )
        return( list( model = model, messages = character( 0 ) ) )
    gate.defs <- model$gate.definitions
    res <- remap.boundaries( gate.defs, model$trained.boundaries, src, dst.params )
    report <- res$report
    model$trained.boundaries <- res$boundaries
    for ( key in c( "template", "template.original" ) ) {
        has <- vapply( gate.defs, function( g ) length( g[[ key ]] ) > 0L, logical( 1 ) )
        if ( !any( has ) )
            next
        templates <- stats::setNames( lapply( gate.defs[ has ], `[[`, key ),
            vapply( gate.defs[ has ], function( g ) g$gate.name, character( 1 ) ) )
        tr <- remap.boundaries( gate.defs, templates, src, dst.params )
        for ( ch in setdiff( names( tr$report ), names( report ) ) )
            report[ ch ] <- tr$report[[ ch ]]
        for ( i in which( has ) )
            model$gate.definitions[[ i ]][[ key ]] <- tr$boundaries[[ gate.defs[[ i ]]$gate.name ]]
    }
    remapped <- sort( names( report )[ report == "remapped" ] )
    for ( ch in remapped )
        model$transforms[[ ch ]] <- dst.params[[ ch ]]
    unknown <- names( report )[ report == "unknown" ]
    no.source <- sort( unknown[ !( unknown %in% names( src ) ) ] )
    no.target <- sort( unknown[ unknown %in% names( src ) & !( unknown %in% names( dst.params ) ) ] )
    messages <- c(
        if ( length( remapped ) )
            paste( "Boundaries were converted to the data's axis transforms on:",
                   paste( remapped, collapse = ", " ) ),
        if ( length( no.source ) )
            paste0( "The model has no stored transform for: ", paste( no.source, collapse = ", " ),
                    ". Boundaries on these channels are used as stored." ),
        if ( length( no.target ) )
            paste0( "The data have no transform for: ", paste( no.target, collapse = ", " ),
                    ". Boundaries on these channels are used as stored." ) )
    list( model = model, messages = messages )
}
