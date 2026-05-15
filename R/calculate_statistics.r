# calculate_statistics.r
#
# Population and activation statistics functions.
#
# Public API:
#   calculate.population.statistics()
#   calculate.activation.statistics()


#' Calculate per-sample population statistics
#'
#' @title Calculate Population Statistics
#'
#' @description
#' Iterates over all population gate definitions and, for each non-ignored
#' population with a valid \code{stats.parent} specification, calculates the
#' fraction of parent-population events inside the gate for every sample.
#' Results are returned as a wide-format data frame with one row per sample
#' and one column per \emph{population :in: parent} combination.
#'
#' Relative parent references of the form \code{"parent.N"} (meaning N levels
#' up in the gate hierarchy) are resolved automatically via
#' \code{resolve.stats.parent()}.
#'
#' @param fcs.gates.sample Named list of gate results, one element per sample.
#'   Each element is a gate results list as produced by
#'   \code{\link{do.gate.generic}}.
#' @param fcs.population.gates List of population gate definitions (same
#'   structure as the input to \code{\link{do.population.gates}}).
#' @param ap Autogating parameter object; uses \code{ap$ignore},
#'   \code{ap$stats.separator}, \code{ap$parent.regexp}, and
#'   \code{ap$parent.regexp.sub}.
#'
#' @return A \code{data.frame} with one row per sample (column
#'   \code{Sample}) and one column per \emph{population :in: parent}
#'   statistic.  Values are numeric fractions in \eqn{[0, 1]}.
#'
#' @seealso \code{\link{calculate.activation.statistics}},
#'   \code{\link{do.population.gates}}
#'
#' @export
calculate.population.statistics <- function( fcs.gates.sample,
    fcs.population.gates, ap ) {

    samples <- names( fcs.gates.sample )
    stats   <- data.frame( Sample = samples, stringsAsFactors = FALSE )

    for ( gd in fcs.population.gates ) {

        gate.name    <- gd$gate.name
        populations  <- gd$populations
        stats.parent <- gd$stats.parent

        sp.gate  <- stats.parent$gate
        sp.popul <- stats.parent$popul

        # Skip entirely if stats.parent is the ignore sentinel on both sides
        if ( identical( sp.gate,  ap$ignore ) ||
             identical( sp.popul, ap$ignore ) )
            next

        stopifnot( length( sp.gate ) == length( sp.popul ) )

        for ( sp.idx in seq_along( sp.gate ) ) {

            sg <- sp.gate[ sp.idx ]
            sp <- sp.popul[ sp.idx ]

            # Resolve relative parent references ("parent.N")
            .resolved    <- resolve.stats.parent(sg, sp, gate.name, fcs.gates.sample[[1L]], ap)
            sg           <- .resolved[[1]]
            sp           <- .resolved[[2]]

            parent.label <- fcs.gates.sample[[ 1L ]][[ sg ]][[ sp ]]$label

            for ( pop.name in names( populations ) ) {

                if ( pop.name == ap$ignore )
                    next

                pop.label <- fcs.gates.sample[[ 1L ]][[ gate.name ]][[ pop.name ]]$label
                col.label <- paste0( pop.label, ap$stats.separator, parent.label )

                for ( samp in samples ) {
                    pop.idx    <- fcs.gates.sample[[ samp ]][[ gate.name ]][[ pop.name ]]$idx
                    parent.idx <- fcs.gates.sample[[ samp ]][[ sg ]][[ sp ]]$idx

                    if ( !all( pop.idx %in% parent.idx ) ) {
                      missing.idx <- pop.idx[ !( pop.idx %in% parent.idx ) ]
                      message( "=== Parent check failed ===" )
                      message( "  Sample:     ", samp )
                      message( "  Gate:       ", gate.name )
                      message( "  Population: ", pop.name )
                      message( "  Parent gate/popul: ", sg, " / ", sp )
                      message( "  pop.idx length:    ", length( pop.idx ) )
                      message( "  parent.idx length: ", length( parent.idx ) )
                      message( "  Missing from parent (first 10): ",
                               paste( utils::head( missing.idx, 10 ), collapse = ", " ) )
                      stop( "Parent containment violated: see diagnostics above" )
                    }

                    stats[ stats$Sample == samp, col.label ] <-
                        mean( parent.idx %in% pop.idx )
                }
            }
        }
    }

    stats[ is.na(stats) ] <- 0

    return(stats)
}


#' Calculate per-sample activation statistics
#'
#' @title Calculate Activation Statistics
#'
#' @description
#' For each marker in each activation gate definition, computes three
#' per-sample statistics for events in the positive (activated) population:
#' the fraction of parent events that are positive, the median back-transformed
#' expression value among positive events, and the inter-quartile range (IQR)
#' of back-transformed expression among positive events.
#'
#' Display-only rows (those with \code{gate.type == ""} or an empty
#' \code{markers} vector) are silently skipped.
#'
#' @param fcs.gates.sample Named list of gate results, one element per sample.
#' @param fcs.activation.gates List of activation gate definitions (same
#'   structure as the input to \code{\link{do.activation.gates}}).
#' @param fcs.expr.data Numeric matrix.  Full expression matrix for all
#'   samples; rows are events, columns are antigens.  Rows are indexed by the
#'   \code{idx} vectors stored in gate results.
#' @param ap Autogating parameter object; uses \code{ap$ignore},
#'   \code{ap$stats.separator}, \code{ap$activation.population.label}, and
#'   \code{ap$transform.inv} (named list of inverse transform functions, one
#'   per antigen).
#'
#' @return A \code{data.frame} with one row per sample and three columns per
#'   marker-in-parent combination:
#'   \describe{
#'     \item{\code{<marker>+ fraction :in: <parent>}}{Numeric. Fraction of
#'       parent events positive for the marker.}
#'     \item{\code{<marker>+ median :in: <parent>}}{Numeric. Median
#'       back-transformed expression among positive events.}
#'     \item{\code{<marker>+ iqr :in: <parent>}}{Numeric. IQR of
#'       back-transformed expression among positive events.}
#'   }
#'
#' @seealso \code{\link{calculate.population.statistics}},
#'   \code{\link{do.activation.gates}}
#'
#' @export
calculate.activation.statistics <- function( fcs.gates.sample,
    fcs.activation.gates, fcs.expr.data, ap ) {

    samples <- names( fcs.gates.sample )
    stats   <- data.frame( Sample = samples, stringsAsFactors = FALSE )

    for ( ad in fcs.activation.gates ) {

        if ( ad$gate.type == "" || length( ad$markers ) == 0L )
            next

        parent.gate  <- ad$parent.gate
        parent.popul <- ad$parent.popul

        for ( marker in ad$markers ) {

            gate.name    <- sprintf( "%s.%s.%s", parent.gate, parent.popul, marker )
            parent.label <- fcs.gates.sample[[ 1L ]][[ parent.gate ]][[ parent.popul ]]$label

            lbl.frac   <- sprintf( "%s+ fraction%s%s", marker, ap$stats.separator, parent.label )
            lbl.median <- sprintf( "%s+ median%s%s",   marker, ap$stats.separator, parent.label )
            lbl.iqr    <- sprintf( "%s+ iqr%s%s",      marker, ap$stats.separator, parent.label )

            for ( samp in samples ) {
                pop.idx    <- fcs.gates.sample[[ samp ]][[ gate.name ]][[ ap$activation.population.label ]]$idx
                parent.idx <- fcs.gates.sample[[ samp ]][[ parent.gate ]][[ parent.popul ]]$idx

                stopifnot( all( pop.idx %in% parent.idx ) )

                marker.expr <- ap$transform.inv[[ marker ]](
                    fcs.expr.data[ pop.idx, marker ] )

                stats[ stats$Sample == samp, lbl.frac   ] <- mean( parent.idx %in% pop.idx )
                stats[ stats$Sample == samp, lbl.median ] <- stats::median( marker.expr )
                stats[ stats$Sample == samp, lbl.iqr    ] <- stats::IQR( marker.expr )
            }
        }
    }

    stats[ is.na(stats) ] <- 0

    return(stats)
}


#' Resolve a relative parent gate reference to a concrete gate/population pair
#'
#' @title Resolve Relative Statistics Parent Reference
#'
#' @description
#' Converts a \code{"parent.N"} relative reference (meaning N levels up in
#' the gate hierarchy) into a concrete gate name and population name by
#' walking backwards through the gate hierarchy stored in
#' \code{gates.one.sample}.  If \code{sg} and \code{sp} are already concrete
#' names (i.e. do not match \code{ap$parent.regexp}), they are returned
#' unchanged.
#'
#' @param sg Character scalar.  Gate name, or a relative reference of the
#'   form \code{"parent.N"}.
#' @param sp Character scalar.  Population name within \code{sg}; must equal
#'   \code{sg} when using a relative reference.
#' @param gate.name Character scalar.  Name of the current gate (the starting
#'   point for upward traversal).
#' @param gates.one.sample Named list.  Gate results for a single sample
#'   (used to look up \code{parent.gate} at each level).
#' @param ap Autogating parameter object; uses \code{ap$parent.regexp} and
#'   \code{ap$parent.regexp.sub}.
#'
#' @return A length-2 list \code{list(gate.name, population.name)} with the
#'   resolved concrete names.
#'
#' @keywords internal
resolve.stats.parent <- function( sg, sp, gate.name, gates.one.sample, ap ) {

  parent.regexp     <- ap$parent.regexp
  parent.regexp.sub <- ap$parent.regexp.sub

  if ( grepl( parent.regexp, sg ) && identical( sg, sp ) ) {
    n            <- as.integer( sub( parent.regexp.sub, "\\1", sg ) )
    current.gate <- gate.name

    while ( n > 0L ) {
      parent.vec <- gates.one.sample[[ current.gate ]][[ 1L ]]$parent.gate
      popul.vec  <- gates.one.sample[[ current.gate ]][[ 1L ]]$parent.popul

      stopifnot( length( parent.vec ) == 1L &&
                   length( popul.vec  ) == 1L )

      current.gate <- parent.vec
      n            <- n - 1L
    }

    sg <- current.gate
    sp <- popul.vec
  }

  return( list( sg, sp ) )
}
