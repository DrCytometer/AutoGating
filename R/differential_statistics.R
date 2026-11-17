# differential_statistics.R
#
# Group comparisons of gated populations, as Honeychrome Automated Gating
# makes them (differential_stats and ag_results.run_statistics).
#
# Frequencies.  Each population is compared as the log2 odds of its events
# within its parent population (or its statistics parent),
#     y = log2( (n + p) / (d - n + p) ),
# n the population's events and d the parent's.  The prior count p keeps
# 0% and 100% finite and is the same share of every sample's parent:
# p = 0.5 * d / mean(d), half an event in a sample of average size (the
# mean over the samples tested; edgeR's scaling of the prior count).  An
# empty population then has the same value in every sample whatever the
# size of its parent, so a difference in parent size between groups is not
# read as a difference in frequency.  Differences are log2 odds ratios.
# A frequency is measured only where the parent has at least
# min.parent.events events (default 20); elsewhere it is missing and the
# population is tested on the samples that have it.  Populations measured
# in fewer than 3 samples, or with events in fewer than 3, are left out.
#
# Marker medians.  The median of each lineage and activation channel in
# each population (display scale), from samples where the population has
# at least min.marker.events events (default 20).  The channels that define
# a population (those of its own gate and the gates above it) are left
# out.  With population.first (default), each population is a family: it
# is screened on all its markers together (Simes) and its markers are
# tested only when it passes (stage-wise FDR, Benjamini and Bogomolov
# 2014); otherwise every population and marker is one list.
#
# Model.  limma's moderated linear model on the values (lmFit, eBayes):
# an intercept, one indicator per group other than the baseline, optional
# pairing and adjustment covariates.  Several samples from one individual
# enter as a fixed pairing term when individuals have samples in more than
# one group, or as a random effect (duplicateCorrelation, then lmFit with
# the consensus correlation) when each individual belongs to one group.
# P-values are TREAT p-values for |logFC| above the threshold (McCarthy and
# Smyth 2009; default 0.5 for frequencies, 0.03 for marker medians), with
# Benjamini-Hochberg FDR per comparison and over all comparisons
# (significance from the latter by default).  Groups need at least 3
# samples.  Comparisons are named by group, "other vs baseline", so group
# names need not be syntactic.
#
# Not here (v1): the counts test, ANOVA variance tables and per-comparison
# refits for missing covariate values (samples with a missing covariate
# are left out, with a message).


ds.defaults <- list( min.samples.per.group = 3L, min.measured.samples = 3L,
                     min.nonzero.samples = 3L, min.parent.events = 20L,
                     min.marker.events = 20L, pval.threshold = 0.05,
                     fc.threshold = 0.5, marker.threshold = 0.03 )


#' Prior counts that are the same share of every denominator
#' @title Share Prior
#' @param d Matrix of denominators (samples x features).
#' @param mean.d Mean denominator per feature, or \code{NULL} for the
#'   column means of the positive, finite denominators.
#' @return Matrix \code{0.5 * d / mean.d}; \code{NA} where \code{d} is not
#'   positive.
#' @keywords internal
share.prior <- function( d, mean.d = NULL ) {
    usable <- is.finite( d ) & d > 0
    if ( is.null( mean.d ) )
        mean.d <- colSums( ifelse( usable, d, 0 ) ) / colSums( usable )
    p <- 0.5 * sweep( d, 2L, mean.d, "/" )
    p[ !usable ] <- NA_real_
    p
}


#' Log2 odds of populations within their denominators
#'
#' @title Log2 Odds
#'
#' @description
#' \eqn{\log_2((n + p) / (d - n + p))} with the prior count \eqn{p} the same
#' share of every sample's denominator (see the file header of
#' \code{differential_statistics.R}), or half an event with
#' \code{share.prior = FALSE}.
#'
#' @param counts,denominators Matrices samples x populations.
#' @param share.prior Logical.
#'
#' @return Matrix of log2 odds; \code{NA} where the denominator is missing
#'   or not positive.
#'
#' @keywords internal
population.log2.odds <- function( counts, denominators, share.prior = TRUE ) {
    if ( any( counts > denominators, na.rm = TRUE ) )
        stop( "A population count exceeds its denominator." )
    p <- if ( share.prior ) share.prior( denominators ) else 0.5
    out <- log2( ( counts + p ) / ( denominators - counts + p ) )
    out[ !is.finite( denominators ) | denominators <= 0 ] <- NA_real_
    out
}


#' Population counts and denominators per sample
#'
#' @title Frequency Tables
#'
#' @param population.table Population table of every sample (as written by
#'   \code{\link{run.autogating}}).
#' @param samples Samples to include.
#' @param use.stats.parent Use a population's statistics parent as its
#'   denominator when it has one.
#' @param min.parent.events Smallest denominator measured.
#'
#' @return List with counts and denominators (samples x populations;
#'   \code{NA} where not measured), dropped and partial (population keys).
#'
#' @keywords internal
population.frequency.tables <- function( population.table, samples, use.stats.parent = TRUE,
                              min.parent.events = ds.defaults$min.parent.events ) {
    pt <- population.table[ population.table$sample %in% samples, , drop = FALSE ]
    keys <- unique( pt$key )
    counts <- matrix( NA_real_, length( samples ), length( keys ), dimnames = list( samples, keys ) )
    denom <- counts
    idx <- cbind( match( pt$sample, samples ), match( pt$key, keys ) )
    use.sp <- use.stats.parent & !is.na( pt$n.stats.parent )
    missing <- use.stats.parent & as.logical( pt$stats.parent.missing ) %in% TRUE
    counts[ idx ] <- as.numeric( pt$n.events )
    denom[ idx ] <- ifelse( missing, NA_real_, ifelse( use.sp, pt$n.stats.parent, pt$n.parent ) )
    unmeasurable <- is.na( counts ) | is.na( denom ) | denom <= 0 | denom < min.parent.events
    counts[ unmeasurable ] <- NA_real_
    denom[ unmeasurable ] <- NA_real_
    measured <- colSums( !unmeasurable )
    nonzero <- colSums( !is.na( counts ) & counts > 0 )
    dropped <- keys[ measured < ds.defaults$min.measured.samples | nonzero < ds.defaults$min.nonzero.samples ]
    keep <- setdiff( keys, dropped )
    list( counts = counts[ , keep, drop = FALSE ], denominators = denom[ , keep, drop = FALSE ],
          dropped = dropped, partial = keep[ colSums( is.na( counts[ , keep, drop = FALSE ] ) ) > 0 ] )
}


#' Channels that define each gate's populations
#' @title Ancestor Channels
#' @param gate.defs Gate definitions.
#' @param gate.name Gate.
#' @return Channels of the gate and every gate above it.
#' @keywords internal
ancestor.channels <- function( gate.defs, gate.name ) {
    by.name <- gates.by.name( gate.defs )
    chans <- character( 0 )
    seen <- character( 0 )
    cur <- by.name[[ gate.name ]]
    while ( !is.null( cur ) && !( cur$gate.name %in% seen ) ) {
        seen <- c( seen, cur$gate.name )
        chans <- union( chans, gate.channels( effective.gate.def( cur, by.name ) ) )
        cur <- if ( has.value( cur$parent.gate ) ) by.name[[ cur$parent.gate ]] else NULL
    }
    chans
}


#' Marker medians per sample
#'
#' @title Marker Median Table
#'
#' @param medians Medians of every sample (as written by
#'   \code{\link{run.autogating}}).
#' @param population.table Population table of every sample.
#' @param samples Samples to include.
#' @param gate.defs Gate definitions.
#' @param exclude.gating.markers Leave out the channels defining each
#'   population.
#' @param min.events Smallest population whose medians are used.
#'
#' @return List with values (samples x features), meta (feature,
#'   population, label, channel) and dropped (population keys).
#'
#' @keywords internal
marker.median.table <- function( medians, population.table, samples, gate.defs,
                                 exclude.gating.markers = TRUE,
                                 min.events = ds.defaults$min.marker.events ) {
    pt <- population.table[ population.table$sample %in% samples, , drop = FALSE ]
    md <- medians[ medians$sample %in% samples, , drop = FALSE ]
    n <- stats::setNames( pt$n.events, paste( pt$sample, pt$key, sep = "\r" ) )
    labels <- stats::setNames( pt$label, pt$key )
    keys <- unique( pt$key[ pt$key %in% md$key ] )
    values <- list()
    meta <- list()
    dropped <- character( 0 )
    for ( key in keys ) {
        k <- n[ paste( samples, key, sep = "\r" ) ]
        measured <- samples[ !is.na( k ) & k >= min.events & samples %in% md$sample[ md$key == key ] ]
        if ( length( measured ) < ds.defaults$min.measured.samples ) {
            dropped <- c( dropped, key )
            next
        }
        gate <- parse.population.key( key )[ 1L ]
        skip <- if ( exclude.gating.markers ) ancestor.channels( gate.defs, gate ) else character( 0 )
        rows <- md[ md$key == key & md$sample %in% measured, , drop = FALSE ]
        chans <- Reduce( intersect, split( rows$channel, rows$sample ) )
        for ( ch in sort( setdiff( chans, skip ), method = "radix" ) ) {
            feature <- paste( key, ch, sep = " | " )
            v <- stats::setNames( rep( NA_real_, length( samples ) ), samples )
            r <- rows[ rows$channel == ch, , drop = FALSE ]
            v[ r$sample ] <- r$median
            values[[ feature ]] <- v
            meta[[ length( meta ) + 1L ]] <- data.frame( feature = feature, population = key,
                label = labels[[ key ]], channel = ch, stringsAsFactors = FALSE )
        }
    }
    list( values = if ( length( values ) ) do.call( cbind, values ) else
              matrix( numeric( 0 ), length( samples ), 0L, dimnames = list( samples, NULL ) ),
          meta = if ( length( meta ) ) do.call( rbind, meta ) else
              data.frame( feature = character( 0 ), population = character( 0 ), label = character( 0 ),
                          channel = character( 0 ), stringsAsFactors = FALSE ),
          dropped = dropped )
}


#' How the pairing variable enters the model
#' @title Resolve Pairing
#' @param group Group of each sample.
#' @param pairing Individual of each sample, or \code{NULL}.
#' @param mode "auto", "fixed" or "random".
#' @return List with pairing (fixed term), block (random effect) and mode.
#' @keywords internal
resolve.pairing <- function( group, pairing, mode = "auto" ) {
    if ( is.null( pairing ) )
        return( list( pairing = NULL, block = NULL, mode = "" ) )
    mode <- match.arg( mode, c( "auto", "fixed", "random" ) )
    n <- table( pairing )
    repeated <- names( n )[ n > 1L ]
    nested <- length( repeated ) > 0L &&
        all( vapply( repeated, function( p ) length( unique( group[ pairing == p ] ) ) == 1L, logical( 1 ) ) )
    if ( mode == "random" || ( mode == "auto" && nested ) )
        return( list( pairing = NULL, block = pairing, mode = "random" ) )
    list( pairing = pairing, block = NULL, mode = "fixed" )
}


#' Design matrix of a comparison
#' @title Design Matrix
#' @description
#' Intercept; one indicator per non-baseline group (in order of first
#' appearance); one indicator per pairing level except the first (sorted);
#' covariates, numeric ones centred, others as indicators of every level
#' except the first (sorted).  Covariates constant over the samples, or
#' within every pairing level, are dropped.  Stops naming the term that
#' makes the design rank-deficient.
#' @param group Group of each sample.
#' @param baseline Baseline group.
#' @param pairing Pairing level of each sample, or \code{NULL}.
#' @param covariates \code{data.frame} of covariates, or \code{NULL}.
#' @return Matrix with an attribute group.columns (group -> column).
#' @keywords internal
design.matrix <- function( group, baseline, pairing = NULL, covariates = NULL ) {
    if ( !( baseline %in% group ) )
        stop( sprintf( "Baseline group '%s' has no samples in this comparison.", baseline ) )
    full.rank <- function( m ) qr( m )$rank == ncol( m )
    X <- matrix( 1, length( group ), 1L, dimnames = list( NULL, "(Intercept)" ) )
    group.cols <- integer( 0 )
    for ( g in setdiff( unique( group ), baseline ) ) {
        X <- cbind( X, as.numeric( group == g ) )
        colnames( X )[ ncol( X ) ] <- sprintf( "group[%s]", g )
        group.cols[ g ] <- ncol( X )
    }
    if ( !is.null( pairing ) ) {
        for ( lvl in sort( unique( pairing ) )[ -1L ] ) {
            X <- cbind( X, as.numeric( pairing == lvl ) )
            colnames( X )[ ncol( X ) ] <- sprintf( "pair[%s]", lvl )
        }
        if ( !full.rank( X ) )
            stop( "The pairing is confounded with the groups: no individual has samples in more than one group." )
    }
    for ( name in names( covariates ) ) {
        v <- covariates[[ name ]]
        if ( length( unique( v ) ) < 2L ) next
        if ( !is.null( pairing ) &&
             all( vapply( split( v, pairing ), function( x ) length( unique( x ) ) == 1L, logical( 1 ) ) ) ) next
        if ( is.numeric( v ) ) {
            X <- cbind( X, v - mean( v ) )
            colnames( X )[ ncol( X ) ] <- name
        } else for ( lvl in sort( unique( as.character( v ) ) )[ -1L ] ) {
            X <- cbind( X, as.numeric( as.character( v ) == lvl ) )
            colnames( X )[ ncol( X ) ] <- sprintf( "%s[%s]", name, lvl )
        }
        if ( !full.rank( X ) )
            stop( sprintf( "Covariate '%s' is confounded with group or pairing in this comparison.", name ) )
    }
    if ( !full.rank( X ) )
        stop( "The design matrix is rank-deficient." )
    attr( X, "group.columns" ) <- group.cols
    X
}


#' Moderated fit of one set of samples, read out for each contrast
#' @title Moderated Fit
#' @param values Matrix samples x features.
#' @param group Group of each sample.
#' @param contrasts List of c(baseline, other), one baseline.
#' @param pairing,block Fixed pairing or random-effect individual, or
#'   \code{NULL}.
#' @param covariates \code{data.frame} or \code{NULL}.
#' @return \code{data.frame} with feature, comparison, logFC, SE, df.total,
#'   AveExpr, t, plus attribute note (random-effect correlation).
#' @keywords internal
moderated.fit <- function( values, group, contrasts, pairing = NULL, block = NULL, covariates = NULL ) {
    base <- contrasts[[ 1L ]][ 1L ]
    X <- design.matrix( group, base, pairing, covariates )
    expr <- t( values )
    if ( nrow( expr ) < 3L )
        stop( sprintf( "Only %d feature(s) can be tested in %s vs %s; at least 3 are needed.",
                       nrow( expr ), contrasts[[ 1L ]][ 2L ], base ) )
    note <- character( 0 )
    if ( is.null( block ) ) {
        fit <- limma::lmFit( expr, X )
    } else {
        dc <- limma::duplicateCorrelation( expr, X, block = block )
        rho <- if ( is.finite( dc$consensus.correlation ) ) dc$consensus.correlation else 0
        note <- sprintf( "Correlation between samples of one individual: %.3f.", rho )
        fit <- limma::lmFit( expr, X, block = block, correlation = rho )
    }
    fit <- limma::eBayes( fit )
    gcols <- attr( X, "group.columns" )
    out <- do.call( rbind, lapply( contrasts, function( cn ) {
        k <- gcols[[ cn[ 2L ] ]]
        se <- fit$stdev.unscaled[ , k ] * sqrt( fit$s2.post )
        data.frame( feature = colnames( values ), comparison = sprintf( "%s vs %s", cn[ 2L ], cn[ 1L ] ),
                    logFC = unname( fit$coefficients[ , k ] ), SE = unname( se ),
                    df.total = unname( fit$df.total ), AveExpr = unname( fit$Amean ),
                    t = unname( fit$t[ , k ] ), stringsAsFactors = FALSE )
    } ) )
    attr( out, "note" ) <- note
    out
}


#' Stage-wise FDR over families of hypotheses
#' @title Stage-wise Adjustment
#' @description
#' Families are screened by their Simes p-value with Benjamini-Hochberg
#' across families; within each of the R selected families (of C testable)
#' Benjamini-Hochberg values are multiplied by C/R (Benjamini and Bogomolov
#' 2014), so they are compared with the FDR level directly.  Hypotheses in
#' unselected families get 1.
#' @param p P-values.
#' @param family Family of each.
#' @param alpha FDR level.
#' @return List with family.p, family.adj and adjusted.
#' @keywords internal
stagewise.adjust <- function( p, family, alpha ) {
    fam <- factor( family, levels = unique( family ) )
    simes <- vapply( split( p, fam ), function( x ) {
        x <- sort( x[ is.finite( x ) ] )
        if ( length( x ) == 0L ) NA_real_ else min( 1, min( x * length( x ) / seq_along( x ) ) )
    }, numeric( 1 ) )
    fam.adj <- stats::p.adjust( simes, "BH" )
    n.fam <- sum( is.finite( simes ) )
    selected <- names( simes )[ is.finite( fam.adj ) & fam.adj <= alpha ]
    adjusted <- ifelse( is.finite( p ), 1, NA_real_ )
    for ( f in selected ) {
        m <- which( fam == f )
        adjusted[ m ] <- pmin( stats::p.adjust( p[ m ], "BH" ) * n.fam / length( selected ), 1 )
    }
    list( family.p = unname( simes[ as.character( fam ) ] ), family.adj = unname( fam.adj[ as.character( fam ) ] ),
          adjusted = adjusted )
}


#' P-values, FDR and significance of a results table
#' @title Apply Significance
#' @param res Results with logFC, SE, df.total, comparison.
#' @param pval.threshold FDR level.
#' @param fc.threshold Smallest |logFC| of interest (TREAT threshold).
#' @param fdr.scope "global" or "per.comparison".
#' @param use.treat TREAT p-values for |logFC| above the threshold; else
#'   ordinary p-values with the threshold as a filter.
#' @param population.first Stage-wise FDR with populations as families
#'   (tables with a population column).
#' @return \code{res} with threshold, P.Value, adj.P.Val (per comparison),
#'   adj.P.Val.global, the stage-wise columns when used, and significant.
#' @keywords internal
apply.significance <- function( res, pval.threshold, fc.threshold, fdr.scope = "global",
                                use.treat = TRUE, population.first = FALSE ) {
    tau <- if ( use.treat ) fc.threshold else 0
    b <- abs( res$logFC )
    p <- stats::pt( ( b - tau ) / res$SE, res$df.total, lower.tail = FALSE ) +
         stats::pt( ( b + tau ) / res$SE, res$df.total, lower.tail = FALSE )
    res$threshold <- tau
    res$P.Value <- pmin( p, 1 )
    res$adj.P.Val <- stats::ave( res$P.Value, res$comparison, FUN = function( x ) stats::p.adjust( x, "BH" ) )
    res$adj.P.Val.global <- stats::p.adjust( res$P.Value, "BH" )
    passes <- !is.na( res$logFC ) & abs( res$logFC ) >= fc.threshold
    if ( population.first && "population" %in% names( res ) ) {
        scope <- if ( fdr.scope == "global" ) rep( "all", nrow( res ) ) else res$comparison
        res$population.P.Value <- NA_real_
        res$population.adj.P.Val <- NA_real_
        res$stagewise.adj.P.Val <- NA_real_
        for ( s in unique( scope ) ) {
            i <- which( scope == s )
            st <- stagewise.adjust( res$P.Value[ i ], paste( res$comparison[ i ], res$population[ i ], sep = "\r" ),
                                    pval.threshold )
            res$population.P.Value[ i ] <- st$family.p
            res$population.adj.P.Val[ i ] <- st$family.adj
            res$stagewise.adj.P.Val[ i ] <- st$adjusted
        }
        res$significant <- !is.na( res$stagewise.adj.P.Val ) & res$stagewise.adj.P.Val <= pval.threshold & passes
    } else {
        q <- if ( fdr.scope == "global" ) res$adj.P.Val.global else res$adj.P.Val
        res$significant <- !is.na( q ) & q <= pval.threshold & passes
    }
    res
}


#' Moderated test of every feature for every contrast
#' @title Moderated Test
#' @param values Matrix samples x features.
#' @param group Group of each sample.
#' @param contrasts List of c(baseline, other).
#' @param mode "reference" (one fit) or "pairwise" (each pair on its own
#'   samples).
#' @param pairing,block,covariates See \code{moderated.fit}.
#' @return Stacked results (before significance) with attribute notes.
#' @keywords internal
moderated.test <- function( values, group, contrasts, mode, pairing = NULL, block = NULL,
                            covariates = NULL ) {
    sub <- function( v, keep ) if ( is.null( v ) ) NULL else v[ keep ]
    if ( mode == "reference" ) {
        res <- moderated.fit( values, group, contrasts, pairing, block, covariates )
        return( res )
    }
    parts <- lapply( contrasts, function( cn ) {
        keep <- group %in% cn
        moderated.fit( values[ keep, , drop = FALSE ], group[ keep ], list( cn ), sub( pairing, keep ),
                       sub( block, keep ), if ( is.null( covariates ) ) NULL else covariates[ keep, , drop = FALSE ] )
    } )
    res <- do.call( rbind, parts )
    attr( res, "note" ) <- unlist( lapply( parts, attr, "note" ) )
    res
}


#' Compare gated populations between groups of samples
#'
#' @title Differential Statistics
#'
#' @description
#' Tests the population frequencies (log2 odds within the parent) and the
#' marker medians of every population between groups of samples, as
#' Honeychrome Automated Gating does; the method is described in the file
#' header of \code{differential_statistics.R}.
#'
#' @param gating Result of \code{\link{run.autogating}}.
#' @param sample.table \code{data.frame}, or path to a CSV, with a
#'   \code{sample} column and the group, pairing and covariate columns.
#' @param group.column Column holding each sample's group; samples with an
#'   empty group are left out.
#' @param group.order Group order (the first qualifying group is the
#'   baseline); groups not listed follow in sorted order.
#' @param contrast.mode "reference" (every group against \code{reference},
#'   one fit) or "pairwise" (every pair, each on its own samples).
#' @param reference Baseline group of "reference" mode.
#' @param pairing.column Optional column naming each sample's individual.
#' @param pairing.mode "auto", "fixed" or "random" (see the file header).
#' @param covariates Covariate columns adjusted for.
#' @param pval.threshold FDR level.
#' @param fc.threshold Threshold on |log2 odds ratio| for frequencies.
#' @param marker.threshold Threshold on |difference| of marker medians.
#' @param use.treat TREAT p-values (else the thresholds filter estimates).
#' @param fdr.scope "global" or "per.comparison".
#' @param use.stats.parent Use statistics parents as denominators.
#' @param min.parent.events Smallest parent measured.
#' @param run.markers Test marker medians.
#' @param markers.population.first Stage-wise FDR with populations as
#'   families.
#' @param exclude.gating.markers Leave out channels that define each
#'   population.
#' @param min.marker.events Smallest population whose medians are used.
#' @param output.dir Optional directory; the result tables are written
#'   there as CSV files.
#' @param figures Also write the figures of the results to
#'   \code{output.dir} (\code{\link{write.differential.figures}}).
#'
#' @return List with frequencies (results), frequency.values (\% of the
#'   denominator, samples x populations), markers, marker.values,
#'   comparisons, samples, groups (sample -> group), group.levels,
#'   pairing.mode and notes.
#'
#' @seealso \code{\link{run.autogating}}
#'
#' @export
differential.statistics <- function(
    gating, sample.table, group.column = "group", group.order = NULL,
    contrast.mode = "reference", reference = NULL, pairing.column = NULL,
    pairing.mode = "auto", covariates = character( 0 ),
    pval.threshold = ds.defaults$pval.threshold, fc.threshold = ds.defaults$fc.threshold,
    marker.threshold = ds.defaults$marker.threshold, use.treat = TRUE, fdr.scope = "global",
    use.stats.parent = TRUE, min.parent.events = ds.defaults$min.parent.events,
    run.markers = TRUE, markers.population.first = TRUE, exclude.gating.markers = TRUE,
    min.marker.events = ds.defaults$min.marker.events, output.dir = NULL, figures = TRUE ) {

    contrast.mode <- match.arg( contrast.mode, c( "reference", "pairwise" ) )
    fdr.scope <- match.arg( fdr.scope, c( "global", "per.comparison" ) )
    st <- if ( is.character( sample.table ) )
        utils::read.csv( sample.table, stringsAsFactors = FALSE, check.names = FALSE ) else sample.table
    for ( col in c( "sample", group.column, pairing.column, covariates ) )
        if ( !( col %in% names( st ) ) ) stop( sprintf( "The sample table has no column '%s'.", col ) )
    notes <- character( 0 )

    # ==========================================================================
    # 1. Samples, groups and contrasts
    # ==========================================================================

    gated <- unique( gating$population.table$sample )
    st <- st[ st$sample %in% gated & !is.na( st[[ group.column ]] ) & nzchar( trimws( st[[ group.column ]] ) ), , drop = FALSE ]
    if ( length( covariates ) ) {
        incomplete <- !stats::complete.cases( st[ , covariates, drop = FALSE ] )
        if ( any( incomplete ) ) {
            notes <- c( notes, sprintf( "Samples left out for a missing covariate value: %s",
                                        paste( st$sample[ incomplete ], collapse = ", " ) ) )
            st <- st[ !incomplete, , drop = FALSE ]
        }
    }
    sizes <- table( as.character( st[[ group.column ]] ) )
    order <- c( group.order, sort( setdiff( names( sizes ), group.order ), method = "radix" ) )
    groups <- order[ order %in% names( sizes )[ sizes >= ds.defaults$min.samples.per.group ] ]
    if ( length( groups ) < 2L )
        stop( sprintf( "At least two groups need %d or more gated samples each.", ds.defaults$min.samples.per.group ) )
    if ( contrast.mode == "reference" ) {
        if ( is.null( reference ) || !( reference %in% groups ) ) reference <- groups[ 1L ]
        contrasts <- lapply( setdiff( groups, reference ), function( g ) c( reference, g ) )
    } else {
        contrasts <- utils::combn( groups, 2L, simplify = FALSE )
    }
    st <- st[ st[[ group.column ]] %in% groups, , drop = FALSE ]
    samples <- st$sample
    group <- as.character( st[[ group.column ]] )
    pr <- resolve.pairing( group, if ( is.null( pairing.column ) ) NULL else as.character( st[[ pairing.column ]] ),
                           pairing.mode )
    cov <- if ( length( covariates ) ) st[ , covariates, drop = FALSE ] else NULL
    test <- function( values ) moderated.test( values, group, contrasts, contrast.mode,
                                               pr$pairing, pr$block, cov )

    # ==========================================================================
    # 2. Frequencies
    # ==========================================================================

    ft <- population.frequency.tables( gating$population.table, samples, use.stats.parent, min.parent.events )
    if ( length( ft$dropped ) )
        notes <- c( notes, sprintf( paste( "Populations measurable in fewer than %d samples (parent below %d events",
                                           "in the others) or with events in fewer than %d samples were left out: %s" ),
                                    ds.defaults$min.measured.samples, min.parent.events,
                                    ds.defaults$min.nonzero.samples, paste( ft$dropped, collapse = ", " ) ) )
    if ( length( ft$partial ) )
        notes <- c( notes, sprintf( "%d population(s) were tested on the samples where they could be measured: %s",
                                    length( ft$partial ), paste( ft$partial, collapse = ", " ) ) )
    labels <- unique( gating$population.table[ , c( "key", "label" ) ] )
    freq <- test( population.log2.odds( ft$counts, ft$denominators ) )
    notes <- c( notes, attr( freq, "note" ) )
    freq$label <- labels$label[ match( freq$feature, labels$key ) ]
    freq <- apply.significance( freq, pval.threshold, fc.threshold, fdr.scope, use.treat )

    # ==========================================================================
    # 3. Marker medians
    # ==========================================================================

    markers <- NULL
    mt <- NULL
    if ( run.markers ) {
        mt <- marker.median.table( gating$medians, gating$population.table, samples,
                                   gating$model$gate.definitions, exclude.gating.markers, min.marker.events )
        if ( length( mt$dropped ) )
            notes <- c( notes, sprintf( "Populations with %d or more events in fewer than %d samples were left out of the marker tests: %s",
                                        min.marker.events, ds.defaults$min.measured.samples,
                                        paste( mt$dropped, collapse = ", " ) ) )
        if ( ncol( mt$values ) >= 3L ) {
            markers <- test( mt$values )
            markers <- merge( markers, mt$meta, by = "feature", sort = FALSE )
            markers <- apply.significance( markers, pval.threshold, marker.threshold, fdr.scope, use.treat,
                                           population.first = markers.population.first )
        } else {
            notes <- c( notes, "Fewer than three marker medians could be tested; the marker test was not run." )
        }
    }

    out <- list( frequencies = freq, frequency.values = 100 * ft$counts / ft$denominators,
                 markers = markers, marker.values = mt$values,
                 comparisons = vapply( contrasts, function( cn ) sprintf( "%s vs %s", cn[ 2L ], cn[ 1L ] ), "" ),
                 samples = samples, groups = stats::setNames( group, samples ), group.levels = groups,
                 pairing.mode = pr$mode, notes = notes )
    for ( msg in notes ) message( msg )
    if ( !is.null( output.dir ) ) {
        dir.create( output.dir, recursive = TRUE, showWarnings = FALSE )
        utils::write.csv( freq, file.path( output.dir, "frequency_tests.csv" ), row.names = FALSE )
        utils::write.csv( data.frame( sample = samples, out$frequency.values, check.names = FALSE ),
                          file.path( output.dir, "frequency_values.csv" ), row.names = FALSE )
        if ( !is.null( markers ) ) {
            utils::write.csv( markers, file.path( output.dir, "marker_tests.csv" ), row.names = FALSE )
            utils::write.csv( data.frame( sample = samples, mt$values, check.names = FALSE ),
                              file.path( output.dir, "marker_values.csv" ), row.names = FALSE )
        }
        if ( figures )
            write.differential.figures( out, output.dir )
    }
    out
}
