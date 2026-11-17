# flowjo_review.R
#
# Data-driven conversion of imported gates (shapes as drawn in FlowJo) into
# trainable gates, following Honeychrome's review of imported gates
# (ag_review).  The events are gated with the hierarchy as imported; each
# imported gate is then judged on its parent population's events:
#
#   1dsep gates, by the share of the parent below the channel's raw zero:
#     >= 40%  "tail" (the negative population is centred on zero)
#     >= 5%   the candidate whose thresholds best reproduce the shapes as
#             drawn: tail, or a 2- or 3-component mixture
#     < 5% and more than 40% of the parent drawn positive: "replicate", the
#             threshold calculated on the anchor population (viable cells)
#             above the gate, with the 25% rule on the anchor's share above
#             the drawn threshold
#     otherwise the 25% rule: "mixture" when more than 25% of the parent is
#             drawn positive, else "tail"
#   2dsep gates: the best-fitting combination of candidates per axis, else
#     the 25% rule per axis.  A single shape whose edges on one axis leave
#     out at most 1% of the parent's events within its band on the other
#     axis separates nothing there, and the gate becomes 1dsep on the other
#     axis.
#   scatter shapes: "cells", keeping the debris side when most drawn events
#     lie below the debris cut.
#
# Agreement of a candidate is the mean over the drawn populations of the
# events in both the drawn and the calculated population divided by the
# events in the larger of the two.  The agreement reported in the proposals
# is NA when, on any threshold axis, less than 2% of the parent or fewer
# than 50 events are drawn positive (nothing to agree on); it does not
# affect the choice of candidate.  At most 10,000 parent events are used
# per gate (seeded subsample).  Honeychrome's interactive split of 2dsep
# gates into 1-D gates is not part of this conversion.


fj.review <- list( centred.below.zero = 0.40, spans.zero.below = 0.05, replicate.positive = 0.40,
                   separates.min = 0.01, fit.max.events = 10000L, evidence.min.events = 50L,
                   score.min.positive = 0.02,
                   overlap.limit = 0.5, overlap.min.events = 20L, overlap.same.as.ancestor = 0.9 )

fj.fit.candidates <- list(
    tail           = list( algorithm = "tail", params = list() ),
    mixture        = list( algorithm = "mixture",
                           params = list( mixture.components = 2L, mixture.cut = "outer" ) ),
    mixture3_lower = list( algorithm = "mixture",
                           params = list( mixture.components = 3L, mixture.cut = "lower" ) ),
    mixture3_upper = list( algorithm = "mixture",
                           params = list( mixture.components = 3L, mixture.cut = "upper" ) ) )

fj.fit.labels <- c( tail = "tail", mixture = "mixture",
                    mixture3_lower = "mixture (3 components, lower cut)",
                    mixture3_upper = "mixture (3 components, upper cut)" )


#' Trainable algorithms of a gate type
#' @title Trainable Algorithms
#' @param gate.def Gate definition.
#' @return Character vector (without "imported").
#' @keywords internal
fj.trainable.algorithms <- function( gate.def ) {
    setdiff( ag.algorithms.by.type[[ gate.type.of( gate.def ) ]] %||% character( 0 ), "imported" )
}


#' Positive fraction of each axis of a threshold gate
#' @title Positive Fractions
#' @param gate.def Gate definition.
#' @param fractions Named list population -> fraction of the parent.
#' @return Named list axis -> fraction on the positive side.
#' @keywords internal
fj.positive.fractions <- function( gate.def, fractions ) {
    if ( !( gate.type.of( gate.def ) %in% ag.threshold.types ) )
        return( list() )
    regions <- population.regions( gate.def )
    out <- list()
    for ( axis in threshold.axes( gate.def ) ) {
        side <- vapply( regions, function( r ) {
            if ( identical( gate.def$gate.type, "1dsep" ) ) return( if ( r == "above" ) "+" else "-" )
            region.sign( r, axis ) %||% NA_character_
        }, character( 1 ) )
        f <- vapply( names( regions ), function( p ) as.numeric( fractions[[ p ]] %||% NA ), numeric( 1 ) )
        plus <- f[ side %in% "+" ]
        minus <- f[ side %in% "-" ]
        if ( length( plus ) && !anyNA( plus ) ) out[[ axis ]] <- sum( plus )
        else if ( length( minus ) && !anyNA( minus ) ) out[[ axis ]] <- 1 - sum( minus )
    }
    out
}


#' Extent of a drawn shape on a channel
#' @title Shape Extent
#' @param gate.def Gate definition with a template.
#' @param pop Population.
#' @param channel Channel.
#' @return \code{c(lo, hi)}, or \code{NULL}.
#' @keywords internal
fj.shape.extent <- function( gate.def, pop, channel ) {
    e <- gate.def$template[[ pop ]]
    chans <- e$shape.channels
    k <- match( channel, chans )
    if ( is.null( e ) || is.na( k ) ) return( NULL )
    if ( !is.null( e$range ) && k == 1L ) return( as.numeric( e$range[ 1:2 ] ) )
    b <- e$boundary
    if ( !is.matrix( b ) || ncol( b ) < k || nrow( b ) == 0L ) return( NULL )
    range( b[ , k ] )
}


#' Axis of a single drawn shape that separates nothing
#' @title Separating Axis Dropped
#' @param gate.def 2dsep gate definition with one population.
#' @param parent Parent events.
#' @return "x", "y" or \code{NULL}.
#' @keywords internal
fj.separating.axis.dropped <- function( gate.def, parent ) {
    if ( !identical( gate.def$gate.type, "2dsep" ) || length( gate.def$populations ) != 1L ||
         !is.shape.boundary( gate.def$template ) )
        return( NULL )
    pop <- names( gate.def$populations )
    chans <- c( x = gate.def$gate.marker.x, y = gate.def$gate.marker.y )
    if ( !all( chans %in% colnames( parent ) ) ) return( NULL )
    ext <- lapply( chans, function( ch ) fj.shape.extent( gate.def, pop, ch ) )
    if ( any( vapply( ext, is.null, logical( 1 ) ) ) ) return( NULL )
    inside <- lapply( c( x = "x", y = "y" ), function( a ) {
        v <- parent[ , chans[[ a ]] ]
        v >= ext[[ a ]][ 1L ] & v <= ext[[ a ]][ 2L ]
    } )
    cuts <- c()
    for ( a in c( "x", "y" ) ) {
        band <- inside[[ if ( a == "x" ) "y" else "x" ]]
        if ( sum( band ) >= 20L ) cuts[ a ] <- 1 - mean( inside[[ a ]][ band ] )
    }
    if ( length( cuts ) != 2L ) return( NULL )
    axis <- names( cuts )[ which.min( cuts ) ]
    if ( cuts[[ axis ]] <= fj.review$separates.min ) axis else NULL
}


#' A 2dsep gate as a 1dsep gate on its other axis
#' @title Without Axis
#' @param gate.def 2dsep gate definition.
#' @param axis Axis dropped.
#' @return 1dsep gate definition drawn against the dropped channel; stops
#'   when two populations share a side, or one spans both.
#' @keywords internal
fj.without.axis <- function( gate.def, axis ) {
    keep <- if ( axis == "x" ) "y" else "x"
    sep <- axis.channel( gate.def, keep )
    view <- axis.channel( gate.def, axis )
    pops <- gate.def$populations
    regions <- population.regions( gate.def )
    sides <- character( 0 )
    for ( p in names( pops ) ) {
        r <- regions[[ p ]]
        sign <- if ( r == "above" ) "+" else if ( r == "below" ) "-" else region.sign( r, keep )
        if ( is.null( sign ) )
            stop( sprintf( "'%s' lies on both sides of %s.", p, sep ) )
        if ( sign %in% sides )
            stop( sprintf( "'%s' and '%s' are on the same side of %s.", names( sides )[ sides == sign ], p, sep ) )
        sides[ p ] <- sign
        pops[[ p ]]$region <- if ( sign == "+" ) "above" else "below"
        pops[[ p ]]$label.pos <- if ( sign == "+" ) 2L else 1L
    }
    out <- gate.def
    out$gate.type <- "1dsep"
    out$gate.marker.x <- sep
    out$gate.marker.y <- view
    out$sep.axis <- keep
    out$populations <- pops
    gp <- gate.def$gate.param %||% list()
    new <- gp[ !( grepl( "\\.[xy]$", names( gp ) ) | names( gp ) == "algorithm" ) ]
    for ( key in c( "tail.fraction", "mixture.components", "axis.lo", "axis.hi", "threshold.factor" ) )
        if ( !is.null( gp[[ paste0( key, ".", keep ) ]] ) ) new[[ key ]] <- gp[[ paste0( key, ".", keep ) ]]
    if ( is.list( gp[[ paste0( "replicate.", keep ) ]] ) ) new$replicate.x <- gp[[ paste0( "replicate.", keep ) ]]
    out$gate.param <- json.object( new )
    if ( has.value( gp[[ paste0( "algorithm.", keep ) ]] ) && !identical( out$algorithm, "imported" ) )
        out$algorithm <- gp[[ paste0( "algorithm.", keep ) ]]
    out
}


#' Default reference population of a gate
#' @title Default Reference
#' @param gate.defs Gate definitions.
#' @param gate.def Gate definition.
#' @param anchor Anchor population key, or \code{NULL}.
#' @return The anchor when it lies above the population the gate is drawn
#'   on, else \code{NULL}.
#' @keywords internal
fj.default.reference <- function( gate.defs, gate.def, anchor ) {
    if ( is.null( anchor ) ) return( NULL )
    by.name <- gates.by.name( gate.defs )
    g <- by.name[[ gate.def$parent.gate %||% "" ]]
    seen <- character( 0 )
    while ( !is.null( g ) && !( g$gate.name %in% seen ) && has.value( g$parent.gate ) ) {
        seen <- c( seen, g$gate.name )
        if ( identical( population.key( g$parent.gate, g$parent.popul %||% "root" ), anchor ) )
            return( anchor )
        g <- by.name[[ g$parent.gate ]]
    }
    NULL
}


#' Drawn threshold of a 1dsep gate
#' @title Drawn Threshold
#' @param gate.def 1dsep gate definition with a template.
#' @return Lower edge of the shapes on the positive side, else upper edge of
#'   those on the negative side, or \code{NULL}.
#' @keywords internal
fj.drawn.threshold <- function( gate.def ) {
    ch <- gate.def$gate.marker.x
    regions <- population.regions( gate.def )
    above <- lapply( names( regions )[ regions == "above" ], fj.shape.extent, gate.def = gate.def, channel = ch )
    below <- lapply( names( regions )[ regions == "below" ], fj.shape.extent, gate.def = gate.def, channel = ch )
    if ( length( above ) && !any( vapply( above, is.null, logical( 1 ) ) ) )
        return( min( vapply( above, `[`, numeric( 1 ), 1L ) ) )
    if ( length( below ) && !any( vapply( below, is.null, logical( 1 ) ) ) )
        return( max( vapply( below, `[`, numeric( 1 ), 2L ) ) )
    NULL
}


#' Give one axis a fit candidate
#' @title Apply Candidate
#' @param d Gate definition.
#' @param axis "x" or "y".
#' @param label Name in \code{fj.fit.candidates}.
#' @return Gate definition.
#' @keywords internal
fj.apply.candidate <- function( d, axis, label ) {
    cand <- fj.fit.candidates[[ label ]]
    gp <- d$gate.param %||% list()
    suffix <- if ( identical( d$gate.type, "2dsep" ) ) paste0( ".", axis ) else ""
    if ( nzchar( suffix ) ) gp[[ paste0( "algorithm", suffix ) ]] <- cand$algorithm
    else d$algorithm <- cand$algorithm
    for ( key in c( "mixture.components", "mixture.cut" ) ) gp[[ paste0( key, suffix ) ]] <- NULL
    for ( key in names( cand$params ) ) gp[[ paste0( key, suffix ) ]] <- cand$params[[ key ]]
    d$gate.param <- gp
    d
}


#' Agreement of calculated with drawn populations
#' @title Agreement
#' @param drawn,calculated Named lists population -> logical mask.
#' @return Mean over drawn populations of |both| / max(|drawn|, |calculated|).
#' @keywords internal
fj.agreement <- function( drawn, calculated ) {
    if ( length( drawn ) == 0L ) return( NULL )
    mean( vapply( names( drawn ), function( p ) {
        d <- drawn[[ p ]]
        c <- calculated[[ p ]]
        if ( is.null( c ) ) return( 0 )
        m <- max( sum( d ), sum( c ) )
        if ( m ) sum( c & d ) / m else 1
    }, numeric( 1 ) ) )
}


#' Agreement of each candidate combination with the shapes as drawn
#' @title Fit Scores
#' @param definition Gate definition.
#' @param parent Parent events.
#' @param drawn Drawn masks.
#' @param ctx Calculation context.
#' @param candidates Named list axis -> candidate labels (\code{NA}: the
#'   axis keeps its algorithm).
#' @return Named numeric vector, names the per-axis labels joined with "|".
#' @keywords internal
fj.fit.scores <- function( definition, parent, drawn, ctx, candidates ) {
    axes <- names( candidates )
    labels <- sort( unique( stats::na.omit( unlist( candidates ) ) ), method = "radix" )
    thresholds <- list()
    for ( label in if ( length( labels ) ) labels else NA_character_ ) {
        d <- definition
        if ( identical( d$gate.type, "2dsep" ) ) d$algorithm <- "tail"
        for ( a in axes )
            if ( !is.na( label ) && label %in% candidates[[ a ]] ) d <- fj.apply.candidate( d, a, label )
        b <- tryCatch( calculate.gate( parent, d, ctx ), error = function( e ) list() )
        if ( length( b ) )
            thresholds[[ if ( is.na( label ) ) ".keep" else label ]] <-
                vapply( axes, function( a ) as.numeric( b[[ 1L ]][[ paste0( "threshold.", a ) ]] %||% NA ), numeric( 1 ) )
    }
    if ( length( thresholds ) == 0L ) return( numeric( 0 ) )
    any.t <- thresholds[[ 1L ]]
    combos <- expand.grid( rev( lapply( candidates, function( v ) v ) ), stringsAsFactors = FALSE )
    combos <- combos[ , rev( seq_along( axes ) ), drop = FALSE ]
    out <- numeric( 0 )
    for ( r in seq_len( nrow( combos ) ) ) {
        combo <- unlist( combos[ r, ] )
        if ( any( !is.na( combo ) & !( combo %in% names( thresholds ) ) ) ) next
        d <- definition
        d$algorithm <- "tail"
        gp <- d$gate.param %||% list()
        ok <- TRUE
        for ( k in seq_along( axes ) ) {
            t <- ( if ( is.na( combo[ k ] ) ) any.t else thresholds[[ combo[ k ] ]] )[[ axes[ k ] ]]
            if ( is.na( t ) ) { ok <- FALSE; break }
            gp[[ paste0( "fixed.threshold.", axes[ k ] ) ]] <- t
        }
        if ( !ok ) next
        d$gate.param <- gp
        masks <- tryCatch( {
            b <- calculate.gate( parent, d, ctx )
            if ( length( b ) ) population.masks( d, b, parent ) else list()
        }, error = function( e ) list() )
        if ( length( masks ) )
            out[ paste( ifelse( is.na( combo ), "", combo ), collapse = "|" ) ] <- fj.agreement( drawn, masks )
    }
    out
}


#' Evidence on converting one imported gate
#'
#' @title Conversion Evidence
#'
#' @param gate.def Imported gate definition.
#' @param gate.defs All gate definitions.
#' @param data Event matrix the hierarchy was gated on.
#' @param result Gating result of the hierarchy as imported.
#' @param zeros Named numeric channel -> raw zero in display units.
#' @param anchor Anchor population key, or \code{NULL}.
#' @param ctx Calculation context.
#' @param seed Seed of the parent subsample.
#'
#' @return Empty list when the gate has no drawn shape or too few parent
#'   events; else a list with definition (the gate, or as 1dsep when an axis
#'   separates nothing), drop.axis, n, positive, below.zero, fit,
#'   reference, reference.positive and, for a scatter shape, cell.side.
#'
#' @keywords internal
fj.conversion.evidence <- function( gate.def, gate.defs, data, result, zeros, anchor, ctx,
                                    seed = 42L ) {
    name <- gate.def$gate.name
    idx <- result$parent.events[[ name ]]
    if ( is.null( idx ) || !is.shape.boundary( gate.def$template ) ||
         length( idx ) < fj.review$evidence.min.events )
        return( list() )
    n <- length( idx )
    if ( n > fj.review$fit.max.events )
        idx <- sort( idx[ seeded.call( seed, function() sample.int( n, fj.review$fit.max.events ) ) ] )
    if ( !all( shape.channels( gate.def ) %in% colnames( data ) ) )
        return( list() )
    parent <- data[ idx, , drop = FALSE ]
    drawn <- population.masks( gate.def, gate.def$template, parent )
    ev <- list( definition = gate.def, drop.axis = NULL, n = n )
    gp <- gate.def$gate.param %||% list()
    if ( identical( gate.def$gate.type, "free" ) ) {
        chans <- shape.channels( gate.def )
        if ( length( chans ) == 2L && length( drawn ) ) {
            col <- parent[ , chans[ if ( identical( gp$debris.axis, "y" ) ) 2L else 1L ] ]
            cut <- debris.cut( col, gp )
            inside <- Reduce( `|`, drawn )
            if ( any( inside ) )
                ev$cell.side <- if ( mean( col[ inside ] < cut ) > 0.5 ) "debris" else "cells"
        }
        return( ev )
    }
    if ( !( gate.type.of( gate.def ) %in% ag.threshold.types ) )
        return( ev )
    definition <- gate.def
    drop <- fj.separating.axis.dropped( gate.def, parent )
    if ( !is.null( drop ) ) {
        d <- tryCatch( fj.without.axis( gate.def, drop ), error = function( e ) NULL )
        if ( !is.null( d ) ) {
            definition <- d
            ev$drop.axis <- drop
        }
    }
    ev$definition <- definition
    ev$positive <- fj.positive.fractions( definition, lapply( drawn, mean ) )
    below <- list()
    for ( axis in threshold.axes( definition ) ) {
        ch <- axis.channel( definition, axis )
        if ( ch %in% colnames( parent ) && ch %in% names( zeros ) )
            below[[ axis ]] <- mean( parent[ , ch ] < zeros[[ ch ]] )
    }
    ev$below.zero <- below
    dgp <- definition$gate.param %||% list()
    if ( identical( definition$gate.type, "2dsep" ) ) {
        cands <- lapply( c( x = "x", y = "y" ), function( a )
            if ( has.value( dgp[[ paste0( "algorithm.", a ) ]] ) ) NA_character_ else names( fj.fit.candidates ) )
        if ( !( "replicate" %in% c( dgp$algorithm.x, dgp$algorithm.y ) ) )
            ev$fit <- fj.fit.scores( definition, parent, drawn, ctx, cands )
    } else {
        b0 <- below$x
        if ( !is.null( b0 ) && b0 >= fj.review$spans.zero.below && b0 < fj.review$centred.below.zero )
            ev$fit <- fj.fit.scores( definition, parent, drawn, ctx, list( x = names( fj.fit.candidates ) ) )
        ref <- fj.default.reference( gate.defs, definition, anchor )
        thr <- fj.drawn.threshold( definition )
        ch <- definition$gate.marker.x
        if ( !is.null( ref ) && !is.null( thr ) && ch %in% colnames( data ) ) {
            key <- parse.population.key( ref )
            ridx <- result$events[[ key[ 1L ] ]][[ key[ 2L ] ]]
            if ( length( ridx ) >= fj.review$evidence.min.events ) {
                ev$reference <- ref
                ev$reference.positive <- mean( data[ ridx, ch ] >= thr )
            }
        }
    }
    ev
}


#' Algorithms the evidence proposes per axis
#' @title Evidence Axes
#' @param ev Conversion evidence.
#' @return Named list axis -> list(algorithm, reason, params).
#' @keywords internal
fj.evidence.axes <- function( ev ) {
    d <- ev$definition
    if ( is.null( d ) || !( gate.type.of( d ) %in% ag.threshold.types ) )
        return( list() )
    pos <- ev$positive %||% list()
    fit <- ev$fit %||% numeric( 0 )
    by.rule <- function( axis ) {
        f <- pos[[ axis ]]
        if ( is.null( f ) ) return( NULL )
        list( algorithm = fj.histogram.algorithm( f ),
              reason = sprintf( "%.1f%% of the parent drawn positive", 100 * f ), params = list() )
    }
    fitted <- function( label, why ) {
        cand <- fj.fit.candidates[[ label ]]
        list( algorithm = cand$algorithm,
              reason = sprintf( "%sbest fit to the shapes as drawn: %s", why, fj.fit.labels[[ label ]] ),
              params = cand$params )
    }
    best <- function() names( fit )[ which.max( fit ) ]
    if ( identical( d$gate.type, "2dsep" ) ) {
        if ( length( fit ) ) {
            combo <- strsplit( best(), "|", fixed = TRUE )[[ 1L ]]
            combo <- c( combo, rep( "", 2L - length( combo ) ) )
            why <- sprintf( "agreement %.2f; ", max( fit ) )
            out <- list()
            for ( k in 1:2 ) if ( nzchar( combo[ k ] ) ) out[[ c( "x", "y" )[ k ] ]] <- fitted( combo[ k ], why )
            return( out )
        }
        out <- list()
        for ( a in c( "x", "y" ) ) { r <- by.rule( a ); if ( !is.null( r ) ) out[[ a ]] <- r }
        return( out )
    }
    b0 <- ev$below.zero$x
    rule.x <- function() { r <- by.rule( "x" ); if ( is.null( r ) ) list() else list( x = r ) }
    if ( is.null( b0 ) ) return( rule.x() )
    if ( b0 >= fj.review$centred.below.zero )
        return( list( x = list( algorithm = "tail", params = list(),
            reason = sprintf( "%.0f%% of the parent below zero (centred on zero)", 100 * b0 ) ) ) )
    if ( b0 >= fj.review$spans.zero.below ) {
        if ( length( fit ) )
            return( list( x = fitted( best(), sprintf( "%.0f%% of the parent below zero; agreement %.2f; ",
                                                         100 * b0, max( fit ) ) ) ) )
        return( rule.x() )
    }
    f <- pos$x
    if ( !is.null( f ) && f > fj.review$replicate.positive && !is.null( ev$reference ) )
        return( list( x = list( algorithm = "replicate", params = list(), reason = sprintf(
            "%.0f%% of the parent drawn positive and %.1f%% below zero: its own events cannot place the threshold",
            100 * f, 100 * b0 ) ) ) )
    rule.x()
}


#' Algorithm an imported gate converts to
#' @title Convert Target
#' @param gate.def Imported gate definition.
#' @param ev Conversion evidence.
#' @return List with algorithm and reason.
#' @keywords internal
fj.convert.target <- function( gate.def, ev ) {
    allowed <- fj.trainable.algorithms( gate.def )
    if ( length( ev ) && isTRUE( gate.def$proposal.rule %in% c( "histogram", "fraction" ) ) ) {
        r <- fj.evidence.axes( ev )$x
        if ( !is.null( r ) && r$algorithm %in% allowed ) return( r[ c( "algorithm", "reason" ) ] )
    }
    if ( isTRUE( gate.def$proposed.algorithm %in% allowed ) )
        return( list( algorithm = gate.def$proposed.algorithm, reason = "proposed at import" ) )
    type <- gate.type.of( gate.def )
    algo <- if ( type %in% ag.threshold.types ) "tail" else if ( type == "singlets" ) "singlets"
            else if ( "cells" %in% allowed ) "cells" else allowed[ 1L ]
    list( algorithm = algo, reason = "default" )
}


#' Trainable gate an imported gate converts to
#' @title Proposed Definition
#' @param gate.def Imported gate definition.
#' @param ev Conversion evidence.
#' @return List with definition and reasons (named character, per axis).
#' @keywords internal
fj.proposed.definition <- function( gate.def, ev ) {
    base <- if ( length( ev ) ) ev$definition %||% gate.def else gate.def
    ruled <- if ( length( ev ) && isTRUE( gate.def$proposal.rule %in% c( "histogram", "fraction" ) ) )
        fj.evidence.axes( ev ) else list()
    target <- fj.convert.target( gate.def, ev )
    reasons <- c( x = target$reason )
    out <- base
    out$algorithm <- target$algorithm
    gp <- out$gate.param %||% list()
    if ( gate.type.of( base ) %in% ag.threshold.types ) {
        axes <- list( x = target$algorithm )
        if ( identical( base$gate.type, "2dsep" ) ) {
            y <- gate.def$gate.param$algorithm.y %||% ruled$y$algorithm %||% target$algorithm
            axes$y <- if ( y %in% fj.trainable.algorithms( base ) ) y else target$algorithm
            reasons[ "y" ] <- ruled$y$reason %||% reasons[[ "x" ]]
            for ( a in names( axes ) )
                if ( is.null( gp[[ paste0( "algorithm.", a ) ]] ) ) gp[[ paste0( "algorithm.", a ) ]] <- axes[[ a ]]
        }
        suffix <- function( a ) if ( identical( base$gate.type, "2dsep" ) ) paste0( ".", a ) else ""
        for ( a in names( axes ) ) {
            r <- ruled[[ a ]]
            if ( !is.null( r ) && identical( r$algorithm, axes[[ a ]] ) )
                for ( k in names( r$params ) )
                    if ( is.null( gp[[ paste0( k, suffix( a ) ) ]] ) ) gp[[ paste0( k, suffix( a ) ) ]] <- r$params[[ k ]]
            key <- paste0( "replicate.", a )
            if ( identical( axes[[ a ]], "replicate" ) && !is.list( gp[[ key ]] ) && !is.null( ev$reference ) )
                gp[[ key ]] <- list( population = ev$reference,
                                     algorithm = fj.histogram.algorithm( ev$reference.positive ),
                                     default = TRUE )
        }
    }
    if ( !is.null( ev$cell.side ) && identical( base$gate.type, "free" ) && is.null( gp$cell.side ) )
        gp$cell.side <- ev$cell.side
    out$gate.param <- json.object( gp )
    list( definition = out, reasons = reasons )
}


#' Populations defined twice in a hierarchy
#'
#' @title Overlapping Populations
#'
#' @description
#' Pairs of populations that gates are drawn on and that hold largely the
#' same events (intersection over union at least \code{limit}) although
#' neither lies below the other, as Honeychrome reports them.  A population
#' sharing at least 90\% of its events with a population above it is
#' represented by that one.  Populations of one gate, and populations with
#' fewer than \code{min.events} events, are not compared.
#'
#' @param gate.defs Gate definitions.
#' @param result Gating result.
#' @param n Number of events gated.
#' @param limit Smallest share reported.
#' @param min.events Smallest population compared.
#'
#' @return \code{data.frame} with population.a, population.b and share.
#'
#' @keywords internal
fj.overlapping.populations <- function( gate.defs, result, n, limit = fj.review$overlap.limit,
                                        min.events = fj.review$overlap.min.events ) {
    by.name <- gates.by.name( gate.defs )
    ancestors <- function( g ) {
        out <- character( 0 )
        seen <- character( 0 )
        while ( !is.null( g ) && has.value( g$parent.gate ) && !( g$parent.gate %in% seen ) &&
                g$parent.gate %in% names( by.name ) ) {
            seen <- c( seen, g$gate.name )
            out <- c( out, population.key( g$parent.gate, g$parent.popul ) )
            g <- by.name[[ g$parent.gate ]]
        }
        out
    }
    parents <- unique( unlist( lapply( gate.defs, function( g )
        if ( has.value( g$parent.gate ) ) population.key( g$parent.gate, g$parent.popul ) ) ) )
    keys <- character( 0 )
    idx <- list()
    for ( k in parents ) {
        kk <- parse.population.key( k )
        e <- result$events[[ kk[ 1L ] ]][[ kk[ 2L ] ]]
        if ( length( e ) >= min.events ) {
            keys <- c( keys, k )
            idx[[ k ]] <- e
        }
    }
    empty <- data.frame( population.a = character( 0 ), population.b = character( 0 ),
                         share = numeric( 0 ), stringsAsFactors = FALSE )
    if ( length( keys ) < 2L ) return( empty )
    m <- matrix( 0, n, length( keys ) )
    for ( j in seq_along( keys ) ) m[ idx[[ j ]], j ] <- 1
    inter <- crossprod( m )
    sizes <- diag( inter )
    union <- outer( sizes, sizes, `+` ) - inter
    share <- ifelse( union > 0, inter / union, 0 )
    dimnames( share ) <- list( keys, keys )
    anc <- lapply( keys, function( k ) ancestors( by.name[[ parse.population.key( k )[ 1L ] ]] ) )
    names( anc ) <- keys
    distinct <- vapply( keys, function( k ) {
        a <- intersect( anc[[ k ]], keys )
        !any( share[ k, a ] >= fj.review$overlap.same.as.ancestor )
    }, logical( 1 ) )
    out <- list()
    for ( i in seq_along( keys ) ) for ( j in seq_along( keys ) ) {
        if ( j <= i || !distinct[ i ] || !distinct[ j ] || share[ i, j ] < limit ) next
        a <- keys[ i ]; b <- keys[ j ]
        if ( parse.population.key( a )[ 1L ] == parse.population.key( b )[ 1L ] ||
             a %in% anc[[ b ]] || b %in% anc[[ a ]] ) next
        out[[ length( out ) + 1L ]] <- data.frame( population.a = a, population.b = b,
                                                   share = share[ i, j ], stringsAsFactors = FALSE )
    }
    if ( length( out ) == 0L ) return( empty )
    out <- do.call( rbind, out )
    out[ order( -out$share ), , drop = FALSE ]
}


#' Convert imported gates into trainable gates
#'
#' @title Convert Imported Gates
#'
#' @description
#' Gates the events with the hierarchy as imported (shapes as drawn), then
#' turns each imported gate into the trainable gate the events support
#' (rules in the file header of \code{flowjo_review.R}).  The shapes as
#' drawn stay in \code{template}.  Run it on the events the model will be
#' trained on (e.g. the gate-calculation samples).
#'
#' @param gate.defs Gate definitions (e.g. from
#'   \code{\link{import.flowjo.workspace}}).
#' @param data Event matrix (display units), columns named by channel.
#' @param transforms Named list channel -> transform parameters of
#'   \code{data}.
#' @param panel Optional panel, for the viability lookup of the anchor
#'   population.
#' @param gates Names of the gates to convert; default every imported gate.
#'
#' @return List with gate.definitions (converted), proposals (one row per
#'   converted gate: type before and after, algorithm per axis, reason,
#'   share below zero, share drawn positive, best agreement, NA when less
#'   than 2\% of the parent or fewer than 50 events are drawn positive on an
#'   axis) and overlaps
#'   (see \code{fj.overlapping.populations}, on the hierarchy as imported).
#'
#' @export
convert.imported.gates <- function( gate.defs, data, transforms = list(), panel = NULL,
                                    gates = NULL ) {
    imported <- vapply( gate.defs, function( g ) identical( g$algorithm, "imported" ), logical( 1 ) )
    names.all <- vapply( gate.defs, `[[`, character( 1 ), "gate.name" )
    todo <- names.all[ imported & ( if ( is.null( gates ) ) TRUE else names.all %in% gates ) ]
    result <- train.boundaries( gate.defs, data, transforms = transforms )
    zeros <- vapply( transforms[ lengths( transforms ) > 0L ],
                     function( p ) make.transform( p )$forward( 0 ), numeric( 1 ) )
    anchor <- anchor.population( gate.defs, panel )
    ctx <- gate.context( transforms = transforms )
    rows <- list()
    out <- gate.defs
    for ( g in ordered.gate.defs( gate.defs ) ) {
        if ( !( g$gate.name %in% todo ) ) next
        ev <- tryCatch( fj.conversion.evidence( g, gate.defs, data, result, zeros, anchor, ctx ),
                        error = function( e ) list() )
        prop <- fj.proposed.definition( g, ev )
        d <- prop$definition
        i <- match( g$gate.name, names.all )
        for ( key in c( "gate.type", "gate.marker.x", "gate.marker.y", "sep.axis", "populations",
                        "gate.param", "algorithm" ) )
            out[[ i ]][ key ] <- list( d[[ key ]] )
        axis.alg <- if ( identical( d$gate.type, "2dsep" ) )
            c( d$gate.param$algorithm.x %||% d$algorithm, d$gate.param$algorithm.y %||% d$algorithm )
            else c( d$algorithm, NA_character_ )
        num <- function( v ) if ( is.null( v ) ) NA_real_ else as.numeric( v )
        pos <- unlist( ev$positive %||% list() )
        scored <- length( ev$fit ) > 0L && length( pos ) > 0L &&
            all( pos >= fj.review$score.min.positive & pos * ev$n >= fj.review$evidence.min.events )
        rows[[ length( rows ) + 1L ]] <- data.frame(
            gate = g$gate.name, imported.type = g$gate.type, type = d$gate.type,
            algorithm.x = axis.alg[ 1L ], algorithm.y = axis.alg[ 2L ],
            reason.x = prop$reasons[[ "x" ]],
            reason.y = if ( "y" %in% names( prop$reasons ) ) prop$reasons[[ "y" ]] else NA_character_,
            dropped.axis = ev$drop.axis %||% NA_character_, n.parent = num( ev$n ),
            below.zero.x = num( ev$below.zero$x ), below.zero.y = num( ev$below.zero$y ),
            positive.x = num( ev$positive$x ), positive.y = num( ev$positive$y ),
            agreement = if ( scored ) max( ev$fit ) else NA_real_,
            stringsAsFactors = FALSE )
    }
    list( gate.definitions = out,
          proposals = if ( length( rows ) ) do.call( rbind, rows ) else data.frame(),
          overlaps = fj.overlapping.populations( gate.defs, result, nrow( data ) ) )
}


#' Import the gates of a FlowJo workspace
#'
#' @title Import FlowJo Workspace
#'
#' @description
#' Reads a FlowJo workspace or template, matches the parameters its gates
#' use to the panel (\code{align.flowjo.parameters}) and converts the gates
#' into gate definitions whose shapes as drawn are applied as they are
#' (algorithm "imported"), in the transforms of the data.  Convert them to
#' trainable gates with \code{\link{convert.imported.gates}}, then train
#' with \code{\link{define.gates}}.  Choices and dropped populations are
#' reported with \code{message()}.
#'
#' @param path Workspace path.
#' @param fcs.setup Named list returned by \code{\link{setup.autogating}}
#'   (panel and transforms).
#' @param ap Autogating parameter object (singlet scatter pair).
#' @param source Name of the group or sample to import; default the group
#'   with the most populations, else the sample with the most.
#' @param mapping Optional named character vector FlowJo parameter ->
#'   antigen overriding the proposed match (\code{NA} drops a parameter).
#' @param merge.siblings Merge sibling gates drawn on the same channels.
#'
#' @return List with gate.definitions, mapping (\code{data.frame} of the
#'   proposed and final matches), dropped (messages), notes (gate ->
#'   character) and populations (FlowJo population -> "gate/pop").
#'
#' @seealso \code{\link{read.flowjo.workspace}},
#'   \code{\link{convert.imported.gates}}
#'
#' @export
import.flowjo.workspace <- function( path, fcs.setup, ap, source = NULL, mapping = NULL,
                                     merge.siblings = TRUE ) {
    ws <- read.flowjo.workspace( path )
    src <- if ( is.null( source ) ) fj.default.source( ws ) else
        Filter( function( s ) identical( s$name, source ), ws$sources )[ 1L ][[ 1L ]]
    if ( is.null( src ) )
        stop( sprintf( "No gating source%s in '%s'.",
                       if ( is.null( source ) ) "" else sprintf( " named '%s'", source ), path ) )
    transforms <- fcs.setup$ap$transform.params
    exp <- fj.experiment( fcs.setup$fcs.panel, transforms[ lengths( transforms ) > 0L ],
                          ap$singlet.scatter.parameter )
    params <- fj.gated.parameters( src )
    proposal <- align.flowjo.parameters( params, src$labels, exp, ws$affixes )
    final <- stats::setNames( proposal$channel, proposal$parameter )
    if ( length( mapping ) )
        final[ intersect( names( mapping ), names( final ) ) ] <- mapping[ intersect( names( mapping ), names( final ) ) ]
    proposal$final <- unname( final[ proposal$parameter ] )
    uncertain <- proposal[ !is.na( proposal$final ) & !( proposal$method %in% c( "name", "fluorophore" ) ) &
                           !( proposal$parameter %in% names( mapping ) ), , drop = FALSE ]
    for ( i in seq_len( nrow( uncertain ) ) )
        message( sprintf( "Check: %s -> %s (%s). %s", uncertain$parameter[ i ], uncertain$final[ i ],
                          uncertain$method[ i ], uncertain$note[ i ] ) )
    res <- fj.build.gate.definitions( src, final, exp, merge.siblings = merge.siblings,
                                      workspace.name = basename( path ) )
    dropped <- fj.plan.messages( res$plan )
    for ( msg in dropped ) message( msg )
    message( sprintf( "%d gate(s) imported from %s.", length( res$gate.definitions ), basename( path ) ) )
    list( gate.definitions = res$gate.definitions, mapping = proposal, dropped = dropped,
          notes = res$notes, populations = res$populations )
}
