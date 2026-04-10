# analyze_aux.r
#
# Utility functions for the FCS autogating pipeline.
#
# Contents:
#   trim.quantiles()          — shared data-trimming helper
#   make.rect.boundary()      — shared rectangular gate boundary constructor
#   get.density.palette()     — colour palette for density scatter plots
#   get.pop.label.pos()       — figure label positioning
#   read.compensation.csv()   — load and reorder compensation matrix
#   read.transformation()     — build flowjo_biexp transform lists
#   load.gate.definitions()   — load an R list or RDS gate definition file


#' Trim a two-column event matrix to interior quantile ranges
#'
#' @title Trim Events to Interior Quantile Ranges
#'
#' @description
#' Subsets rows of a two-column numeric matrix to those that fall within the
#' interior quantile range on both axes simultaneously.  Used by all gate
#' calculation functions to remove extreme outliers before fitting thresholds
#' or boundaries.
#'
#' @param data Numeric matrix with at least two columns.
#' @param qx Numeric scalar in \eqn{[0, 0.5)}.  Lower-tail quantile to trim
#'   on column 1 (both the lower \code{qx} and upper \code{1 - qx} tails are
#'   removed).  Use \code{0} to skip x-axis trimming.
#' @param qy Numeric scalar in \eqn{[0, 0.5)}.  Lower-tail quantile to trim
#'   on column 2.  Use \code{0} to skip y-axis trimming.
#'
#' @return Subset of \code{data} rows lying within
#'   \eqn{[Q_{qx}, Q_{1-qx}]} on column 1 and
#'   \eqn{[Q_{qy}, Q_{1-qy}]} on column 2.  Original matrix dimensions are
#'   preserved via \code{drop = FALSE}.
#'
#' @seealso \code{\link{calculate.gate.1dsep}},
#'   \code{\link{calculate.gate.2dsep}}, \code{\link{calculate.gate.free}}
#'
#' @keywords internal
trim.quantiles <- function( data, qx, qy ) {
    keep <- rep( TRUE, nrow( data ) )
    if ( qx > 0 )
        keep <- keep &
            data[ , 1 ] >= stats::quantile( data[ , 1 ], qx ) &
            data[ , 1 ] <= stats::quantile( data[ , 1 ], 1 - qx )
    if ( qy > 0 )
        keep <- keep &
            data[ , 2 ] >= stats::quantile( data[ , 2 ], qy ) &
            data[ , 2 ] <= stats::quantile( data[ , 2 ], 1 - qy )
    data[ keep, , drop = FALSE ]
}


#' Build a standard 4-vertex rectangular gate boundary
#'
#' @title Construct Rectangular Gate Boundary
#'
#' @description
#' Returns a 4 × 2 numeric matrix representing a closed clockwise rectangle
#' defined by its axis-aligned extents.  Optionally swaps the x and y columns
#' to handle gates on the y-axis (\code{sep.axis == 2}).  This is the
#' canonical boundary representation used throughout the pipeline.
#'
#' @param xmin Numeric scalar.  Left edge (or bottom if \code{swap = TRUE}).
#' @param xmax Numeric scalar.  Right edge (or top if \code{swap = TRUE}).
#' @param ymin Numeric scalar.  Bottom edge (or left if \code{swap = TRUE}).
#' @param ymax Numeric scalar.  Top edge (or right if \code{swap = TRUE}).
#' @param swap Logical.  If \code{TRUE}, swap x and y columns after
#'   construction (used when the separation axis is y).  Default \code{FALSE}.
#'
#' @return A 4 × 2 numeric matrix with no dimnames.  Rows represent the four
#'   corners: bottom-left, bottom-right, top-right, top-left.
#'
#' @seealso \code{\link{calculate.gate.1dsep}},
#'   \code{\link{calculate.gate.2dsep}},
#'   \code{\link{calculate.boundary.all}}
#'
#' @keywords internal
make.rect.boundary <- function( xmin, xmax, ymin, ymax, swap = FALSE ) {
    b <- rbind(
        c( xmin, ymin ),
        c( xmax, ymin ),
        c( xmax, ymax ),
        c( xmin, ymax )
    )
    if ( swap ) b[ , c( 1, 2 ) ] <- b[ , c( 2, 1 ) ]
    b
}


#' Generate a density-mapped colour palette for gate scatter plots
#'
#' @title Density Colour Palette for Gate Plots
#'
#' @description
#' Maps density values to a rainbow colour palette (blue → cyan → green →
#' yellow → red) using the empirical CDF of the density values, so that the
#' palette adapts to the actual density distribution rather than the linear
#' range.  Returns a vector of hex colour strings, one per density grid point.
#'
#' @param fcs.density Numeric vector.  Interpolated density values for each
#'   event (from \code{\link{compute.gate.density}}).
#' @param ap Autogating parameter object; uses \code{ap$palette.base.n} (base
#'   palette resolution) and \code{ap$palette.n} (number of output colours).
#'
#' @return Character vector of length \code{ap$palette.n} containing hex
#'   colour codes.
#'
#' @seealso \code{\link{build.gate.plot}}
#'
#' @keywords internal
get.density.palette <- function( fcs.density, ap ) {
    rainbow.palette <- grDevices::colorRampPalette(
        c( "blue", "cyan", "green", "yellow", "red" ) )( ap$palette.base.n )

    fcs.density <- stats::na.omit( fcs.density )

    fcs.density.grid <- seq( min( fcs.density ), max( fcs.density ),
        length.out = ap$palette.n )

    density.palette.idx <-
        round( stats::ecdf( fcs.density )( fcs.density.grid ) * ap$palette.base.n )

    rainbow.palette[ density.palette.idx ]
}


#' Compute label position coordinates for a gate population
#'
#' @title Population Label Position
#'
#' @description
#' Returns \code{(x, y)} coordinates for placing a population label on a gate
#' plot, based on a compass position code and the extent of the population
#' boundary polygon.  The label is offset from the nearest boundary corner or
#' edge by an amount proportional to the label text dimensions and the axis
#' ranges.
#'
#' @param label.pos Integer position code (1–8) or any other value for
#'   centre:
#'   1 = below bottom-left corner,
#'   2 = below bottom-right corner,
#'   3 = above top-left corner,
#'   4 = above top-right corner,
#'   5 = left of bottom-left corner,
#'   6 = right of bottom-right corner,
#'   7 = left of top-left corner,
#'   8 = right of top-right corner,
#'   other = plot centre.
#' @param pop.boundary Numeric matrix (n × 2).  Gate boundary vertices.
#' @param limit.x Length-2 numeric vector.  x-axis limits.
#' @param limit.y Length-2 numeric vector.  y-axis limits.
#' @param pop.label Character scalar.  Label text; newlines are counted for
#'   multi-line height estimation.
#' @param ap Autogating parameter object; uses
#'   \code{ap$figure.label.width.factor} and
#'   \code{ap$figure.label.height.factor}.
#'
#' @return Length-2 numeric vector \code{c(x, y)}.
#'
#' @seealso \code{\link{add.gate.boundary}}
#'
#' @keywords internal
get.pop.label.pos <- function( label.pos, pop.boundary, limit.x, limit.y,
    pop.label, ap ) {

    pop.xmax <- max( pop.boundary[ , 1 ] )
    pop.ymax <- max( pop.boundary[ , 2 ] )
    pop.xmin <- min( pop.boundary[ , 1 ] )
    pop.ymin <- min( pop.boundary[ , 2 ] )

    label.line   <- scan( text = pop.label, what = character(), sep = "\n",
        quiet = TRUE )
    label.width  <- max( sapply( label.line, nchar ) )
    label.height <- length( label.line )

    dx <- ap$figure.label.width.factor  * label.width  * diff( limit.x )
    dy <- ap$figure.label.height.factor * label.height * diff( limit.y )

    switch( as.character( label.pos ),
        "1" = c( pop.xmin, pop.ymin - dy ),
        "2" = c( pop.xmax, pop.ymin - dy ),
        "3" = c( pop.xmin, pop.ymax + dy ),
        "4" = c( pop.xmax, pop.ymax + dy ),
        "5" = c( pop.xmin - dx, pop.ymin ),
        "6" = c( pop.xmax + dx, pop.ymin ),
        "7" = c( pop.xmin - dx, pop.ymax ),
        "8" = c( pop.xmax + dx, pop.ymax ),
              c( mean( limit.x ), mean( limit.y ) )
    )
}





#' Build forward and inverse biexponential transform lists
#'
#' @title Build flowjo_biexp Transformation Lists
#'
#' @description
#' Constructs named lists of forward and inverse
#' \code{\link[flowWorkspace]{flowjo_biexp}} transform functions for all
#' fluorescence channels, optionally applying per-experiment parameter
#' overrides from a CSV file.  Scatter channels receive a linear
#' (identity) transform.  The resulting lists are keyed by antigen name.
#'
#' @param transformation.filepath Character or \code{NULL}.  Path to a CSV
#'   file with columns \code{dye} plus any of \code{length},
#'   \code{max.range}, \code{pos}, \code{neg}, \code{width} for per-channel
#'   override.  \code{NULL} uses cytometer defaults throughout.
#' @param fcs.panel \code{data.frame}.  Panel definition with columns
#'   \code{dye} and \code{antigen}.
#' @param fcs.panel.scatter.idx Integer vector.  Scatter-channel row indices.
#' @param fcs.panel.lineage.idx Integer vector.  Lineage-channel row indices.
#' @param fcs.panel.activation.idx Integer vector.  Activation-channel row
#'   indices.
#' @param ap Autogating parameter object; uses
#'   \code{ap$transformation.defaults}.
#'
#' @return A length-2 list:
#'   \describe{
#'     \item{\code{[[1]]}}{Named list of forward transform functions (one per
#'       antigen), keyed by antigen name.}
#'     \item{\code{[[2]]}}{Named list of inverse transform functions (one per
#'       antigen), keyed by antigen name.}
#'   }
#'
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom flowCore linearTransform
#'
#' @seealso \code{\link{get.autogating.param}}
#'
#' @keywords internal
read.transformation <- function( transformation.filepath, fcs.panel,
    fcs.panel.scatter.idx, fcs.panel.lineage.idx, fcs.panel.activation.idx,
    ap ) {

    fluor.idx    <- c( fcs.panel.lineage.idx, fcs.panel.activation.idx )
    fluor.dye    <- fcs.panel$dye[ fluor.idx ]
    fluor.antigen <- fcs.panel$antigen[ fluor.idx ]

    # Start from cytometer defaults: one row per fluorescence channel
    defaults <- ap$transformation.defaults
    param.df <- data.frame(
        dye       = fluor.dye,
        length    = defaults$length,
        max.range = defaults$max.range,
        pos       = defaults$pos,
        neg       = defaults$neg,
        width     = defaults$width,
        stringsAsFactors = FALSE
    )

    # Apply per-experiment overrides if a CSV was provided
    if ( !is.null( transformation.filepath ) ) {
        ovr <- utils::read.csv( transformation.filepath, stringsAsFactors = FALSE )
        stopifnot( "dye" %in% colnames( ovr ) )

        for ( col in setdiff( colnames( ovr ), "dye" ) ) {
            if ( col %in% colnames( param.df ) ) {
                idx <- match( ovr$dye, param.df$dye )
                param.df[ idx[ !is.na( idx ) ], col ] <-
                    ovr[ !is.na( idx ), col ]
            }
        }
    }

    make.transform <- function( inverse ) {
        lapply( seq_len( nrow( param.df ) ), function( i )
            flowjo_biexp(
                channelRange = param.df$length[ i ],
                maxValue     = param.df$max.range[ i ],
                pos          = param.df$pos[ i ],
                neg          = param.df$neg[ i ],
                widthBasis   = param.df$width[ i ],
                inverse      = inverse
            )
        )
    }

    fcs.transform     <- make.transform( FALSE )
    fcs.transform.inv <- make.transform( TRUE )

    reorder <- match( fluor.dye, param.df$dye )
    fcs.transform     <- fcs.transform[ reorder ]
    fcs.transform.inv <- fcs.transform.inv[ reorder ]
    names( fcs.transform )     <- fluor.antigen
    names( fcs.transform.inv ) <- fluor.antigen

    # Scatter channels: identity transform
    for ( ant in fcs.panel$antigen[ fcs.panel.scatter.idx ] ) {
        fcs.transform[[ ant ]]     <- linearTransform()
        fcs.transform.inv[[ ant ]] <- linearTransform()
    }

    list( fcs.transform, fcs.transform.inv )
}


#' Load a gate definition list from an R file or RDS file
#'
#' @title Load Gate Definitions from File
#'
#' @description
#' Loads a gate definition list from either an R source file or an RDS file.
#' For R source files, the file is sourced into a fresh environment and the
#' single variable it defines is returned.  If the file defines more than one
#' variable, an error is raised.  If \code{filepath} is \code{NULL}, an empty
#' list is returned.
#'
#' @param filepath Character scalar or \code{NULL}.  Path to the gate
#'   definition file.  Extension determines loading method: \code{".rds"}
#'   uses \code{\link[base]{readRDS}}; any other extension uses
#'   \code{\link[base]{source}}.
#'
#' @return A list containing the gate definitions (population or activation),
#'   or an empty list if \code{filepath} is \code{NULL}.
#'
#' @seealso \code{\link{get.autogating.param}},
#'   \code{\link{do.population.gates}}, \code{\link{do.activation.gates}}
#'
#' @keywords internal
load.gate.definitions <- function( filepath ) {
    if ( is.null( filepath ) ) return( list() )

    ext <- tolower( tools::file_ext( filepath ) )

    if ( ext == "rds" ) {
        readRDS( filepath )
    } else {
        env <- new.env( parent = baseenv() )
        source( filepath, local = env )
        objs <- ls( env )
        if ( length( objs ) != 1L )
            stop( sprintf(
                "Gate definition file '%s' must assign exactly one variable.",
                filepath ) )
        env[[ objs[ 1 ] ]]
    }
}
