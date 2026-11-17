# threshold_algorithms.R
#
# One-dimensional threshold algorithms on the transformed (display) scale,
# as in Honeychrome Automated Gating.
#
#   tail     the negative population is assumed symmetric about its mode
#            (the highest density peak); the events below the mode are
#            mirrored to model it, and the threshold leaves a fraction
#            tail.fraction of that population above it.
#   mixture  a Gaussian mixture (2 or 3 components, unequal variances) is
#            fitted; the threshold is where two components are equally
#            dense: the lowest and highest ("outer"), the lowest two
#            ("lower") or the highest two ("upper").
#   kde_min  the lowest point of the density between its two highest
#            peaks; bimodal is the same rule.
#   otsu     the histogram split that maximises the between-class variance
#            (Otsu 1979).
#
# Kernel densities use a Gaussian kernel with Scott's bandwidth on 512
# points spanning the data; above 20,000 events the events are linearly
# binned on a finer grid and convolved with the kernel (agrees with the
# exact sum to about 1e-5 of the peak).
#
# The mixture is fitted by expectation-maximisation, initialised from a
# one-dimensional k-means started at evenly spaced quantiles, stopping when
# the mean log-likelihood changes by less than 1e-3 (as scikit-learn's
# GaussianMixture). It does not reproduce scikit-learn's random k-means++
# start, so thresholds can differ from Honeychrome's in the last digits of
# a converged fit.


#' Drop events outside interior quantiles
#'
#' @title Trim Events
#'
#' @param data Numeric vector or matrix.
#' @param q Quantile; events outside \eqn{[Q_q, Q_{1-q}]} on any column are
#'   dropped.
#'
#' @return The trimmed vector or matrix.
#'
#' @keywords internal
trim.events <- function( data, q = 0.001 ) {
    if ( q <= 0 || NROW( data ) == 0L )
        return( data )
    if ( is.null( dim( data ) ) ) {
        lo <- stats::quantile( data, q, names = FALSE )
        hi <- stats::quantile( data, 1 - q, names = FALSE )
        return( data[ data >= lo & data <= hi ] )
    }
    keep <- rep( TRUE, nrow( data ) )
    for ( j in seq_len( ncol( data ) ) ) {
        lo <- stats::quantile( data[ , j ], q, names = FALSE )
        hi <- stats::quantile( data[ , j ], 1 - q, names = FALSE )
        keep <- keep & data[ , j ] >= lo & data[ , j ] <= hi
    }
    data[ keep, , drop = FALSE ]
}


#' Restrict values to a window of their range
#'
#' @title Crop to Region
#'
#' @param col Numeric vector.
#' @param lo.frac,hi.frac Fractions of the empirical range (not
#'   percentiles).
#'
#' @return The values inside the window.
#'
#' @keywords internal
crop.to.region <- function( col, lo.frac, hi.frac ) {
    if ( length( col ) == 0L ) return( col )
    c.min <- min( col )
    c.max <- max( col )
    col[ col >= ( 1 - lo.frac ) * c.min + lo.frac * c.max &
         col <= ( 1 - hi.frac ) * c.min + hi.frac * c.max ]
}


#' Linear convolution by FFT
#'
#' @title FFT Convolution
#'
#' @param a,b Numeric vectors.
#'
#' @return The full linear convolution (length \code{length(a) +
#'   length(b) - 1}).
#'
#' @keywords internal
fft.convolve <- function( a, b ) {
    n.out <- length( a ) + length( b ) - 1L
    n.fft <- stats::nextn( n.out )
    fa <- stats::fft( c( a, numeric( n.fft - length( a ) ) ) )
    fb <- stats::fft( c( b, numeric( n.fft - length( b ) ) ) )
    Re( stats::fft( fa * fb, inverse = TRUE ) )[ seq_len( n.out ) ] / n.fft
}


#' Gaussian kernel density on a grid
#'
#' @title Kernel Density Grid
#'
#' @param x Numeric vector.
#' @param n.grid Number of grid points spanning the range of \code{x}.
#'
#' @return List with \code{grid} and \code{density}.
#'
#' @keywords internal
kde.grid <- function( x, n.grid = 512L ) {
    x <- as.numeric( x )
    n <- length( x )
    grid <- seq( min( x ), max( x ), length.out = n.grid )
    bw <- stats::sd( x ) * n^( -0.2 )

    if ( n > ag.defaults$kde.exact.max.events && all( is.finite( x ) ) ) {
        d <- binned.kde( x, grid, bw )
        if ( !is.null( d ) )
            return( list( grid = grid, density = d ) )
    }

    dens  <- numeric( n.grid )
    chunk <- max( 1L, floor( 2e6 / n ) )
    for ( s in seq( 1L, n.grid, by = chunk ) ) {
        idx <- s:min( n.grid, s + chunk - 1L )
        z <- outer( x, grid[ idx ], "-" ) / bw
        dens[ idx ] <- colSums( exp( -0.5 * z * z ) )
    }
    list( grid = grid, density = dens / ( n * bw * sqrt( 2 * pi ) ) )
}


#' Binned Gaussian kernel density
#'
#' @title Binned Kernel Density
#'
#' @param x Numeric vector.
#' @param grid Evenly spaced output grid.
#' @param bw Kernel bandwidth.
#' @param max.refine Largest refinement of the grid.
#'
#' @return Density at \code{grid}, or \code{NULL} when the bandwidth is too
#'   small for the range.
#'
#' @keywords internal
binned.kde <- function( x, grid, bw, max.refine = 2048L ) {
    n <- length( x )
    lo <- grid[ 1L ]
    hi <- grid[ length( grid ) ]
    if ( !( bw > 0 ) || hi <= lo )
        return( NULL )
    n.grid <- length( grid )
    refine <- max( 8L, as.integer( ceiling( 4 * ( hi - lo ) / ( ( n.grid - 1L ) * bw ) ) ) )
    if ( refine > max.refine )
        return( NULL )
    m <- ( n.grid - 1L ) * refine + 1L
    step <- ( hi - lo ) / ( m - 1L )
    pos  <- ( x - lo ) / step
    left <- pmin( pmax( floor( pos ), 0 ), m - 2L )
    frac <- pos - left
    bin  <- as.integer( left ) + 1L
    weights <- numeric( m )
    w1 <- rowsum( 1 - frac, bin, reorder = FALSE )
    w2 <- rowsum( frac, bin + 1L, reorder = FALSE )
    i1 <- as.integer( rownames( w1 ) )
    i2 <- as.integer( rownames( w2 ) )
    weights[ i1 ] <- w1[ , 1 ]
    weights[ i2 ] <- weights[ i2 ] + w2[ , 1 ]
    half <- min( as.integer( ceiling( 37.5 * bw / step ) ), m - 1L )
    offsets <- ( -half:half ) * step
    kernel <- exp( -0.5 * ( offsets / bw )^2 ) / ( bw * sqrt( 2 * pi ) * n )
    density <- fft.convolve( weights, kernel )[ ( half + 1L ):( half + m ) ]
    density[ seq( 1L, m, by = refine ) ]
}


#' Location of the highest density peak
#'
#' @title KDE Mode
#'
#' @param x Numeric vector.
#'
#' @return Numeric.
#'
#' @keywords internal
kde.mode <- function( x ) {
    x <- as.numeric( x )
    if ( length( x ) < 3L || max( x ) == min( x ) )
        return( if ( length( x ) ) stats::median( x ) else 0 )
    k <- kde.grid( x )
    k$grid[ which.max( k$density ) ]
}


#' Local maxima of a sequence
#'
#' @title Find Peaks
#'
#' @description
#' Indices of local maxima (a flat top counts once, at its middle); the
#' first and last points are never peaks.
#'
#' @param d Numeric vector.
#'
#' @return Integer indices.
#'
#' @keywords internal
find.peaks <- function( d ) {
    n <- length( d )
    peaks <- integer( 0 )
    i <- 2L
    while ( i < n ) {
        if ( d[ i - 1L ] < d[ i ] ) {
            j <- i
            while ( j < n && d[ j + 1L ] == d[ i ] ) j <- j + 1L
            if ( j < n && d[ j + 1L ] < d[ i ] ) {
                peaks <- c( peaks, ( i + j ) %/% 2L )
                i <- j
            }
        }
        i <- i + 1L
    }
    peaks
}


#' Tail threshold
#'
#' @title Tail Threshold
#'
#' @description
#' Upper tail of the negative population, assumed symmetric about its mode:
#' \eqn{2 \cdot mode - Q_{2 f}(\mathrm{lower\ half})} for tail fraction
#' \eqn{f}.
#'
#' @param x Numeric vector.
#' @param tail.fraction Fraction of the negative population left above the
#'   threshold.
#'
#' @return Numeric threshold.
#'
#' @export
threshold.tail <- function( x, tail.fraction = ag.defaults$tail.fraction ) {
    x <- as.numeric( x )
    if ( length( x ) < ag.defaults$min.events )
        return( if ( length( x ) ) stats::median( x ) else 0 )
    mode <- kde.mode( x )
    lower <- x[ x <= mode ]
    if ( length( lower ) < 5L )
        return( stats::quantile( x, 1 - tail.fraction, names = FALSE ) )
    q <- min( max( 2 * tail.fraction, 0 ), 1 )
    2 * mode - stats::quantile( lower, q, names = FALSE )
}


#' Density valley threshold
#'
#' @title KDE Minimum Threshold
#'
#' @description
#' Lowest density between the two highest peaks; falls back to the tail
#' threshold at a tail fraction of 0.01 when the density has a single peak.
#'
#' @param x Numeric vector.
#'
#' @return Numeric threshold.
#'
#' @export
threshold.kde.min <- function( x ) {
    x <- as.numeric( x )
    if ( length( x ) < 20L || max( x ) == min( x ) )
        return( if ( length( x ) ) stats::median( x ) else 0 )
    k <- kde.grid( x )
    peaks <- find.peaks( k$density )
    if ( length( peaks ) < 2L )
        return( threshold.tail( x, ag.defaults$kde.min.fallback.tail ) )
    top <- sort( peaks[ order( k$density[ peaks ], decreasing = TRUE ) ][ 1:2 ] )
    k$grid[ top[ 1 ] - 1L + which.min( k$density[ top[ 1 ]:top[ 2 ] ] ) ]
}


#' Otsu threshold
#'
#' @title Otsu Threshold
#'
#' @description
#' The 256-bin histogram split that maximises the between-class variance;
#' the threshold is the upper edge of the split bin.
#'
#' @param x Numeric vector.
#' @param n.bins Number of histogram bins.
#'
#' @return Numeric threshold.
#'
#' @export
threshold.otsu <- function( x, n.bins = 256L ) {
    x <- as.numeric( x )
    if ( length( x ) < 20L || max( x ) == min( x ) )
        return( if ( length( x ) ) stats::median( x ) else 0 )
    lo <- min( x )
    hi <- max( x )
    edges <- seq( lo, hi, length.out = n.bins + 1L )
    idx <- floor( ( x - lo ) * ( n.bins / ( hi - lo ) ) )
    idx[ idx >= n.bins ] <- n.bins - 1L
    dec <- x < edges[ idx + 1L ]
    idx[ dec ] <- idx[ dec ] - 1L
    inc <- x >= edges[ idx + 2L ] & idx != n.bins - 1L
    idx[ inc ] <- idx[ inc ] + 1L
    h <- tabulate( idx + 1L, n.bins )
    centres <- ( edges[ -1L ] + edges[ -( n.bins + 1L ) ] ) / 2
    hc <- h * centres
    w0 <- cumsum( h )
    w1 <- rev( cumsum( rev( h ) ) )
    m0 <- cumsum( hc ) / pmax( w0, 1e-12 )
    m1 <- rev( cumsum( rev( hc ) ) / pmax( cumsum( rev( h ) ), 1e-12 ) )
    between <- w0[ -n.bins ] * w1[ -1L ] * ( m0[ -n.bins ] - m1[ -1L ] )^2
    edges[ which.max( between ) + 1L ]
}


#' Fit a one-dimensional Gaussian mixture
#'
#' @title Fit Gaussian Mixture
#'
#' @description
#' Expectation-maximisation with unequal variances, initialised from a
#' one-dimensional k-means started at evenly spaced quantiles; a variance
#' floor of 1e-6 is added as in scikit-learn.
#'
#' @param x Numeric vector.
#' @param n.components Number of components.
#' @param tol Convergence tolerance on the mean log-likelihood.
#' @param max.iter Maximum EM iterations.
#'
#' @return List with \code{mean}, \code{weight}, \code{sd}, \code{loglik}
#'   (total log-likelihood) and \code{bic}.
#'
#' @keywords internal
fit.mixture <- function( x, n.components, tol = 1e-3, max.iter = 200L ) {
    x <- as.numeric( x )
    n <- length( x )
    k <- as.integer( n.components )
    reg <- 1e-6

    centres <- stats::quantile( x, ( seq_len( k ) - 0.5 ) / k, names = FALSE )
    for ( it in seq_len( 300L ) ) {
        cl <- if ( k == 1L ) rep( 1L, n ) else
            max.col( -abs( outer( x, centres, "-" ) ), ties.method = "first" )
        new <- vapply( seq_len( k ), function( j ) {
            if ( any( cl == j ) ) mean( x[ cl == j ] ) else centres[ j ]
        }, numeric( 1 ) )
        if ( all( new == centres ) ) break
        centres <- new
    }
    resp <- matrix( 0, n, k )
    resp[ cbind( seq_len( n ), cl ) ] <- 1

    m.step <- function( resp ) {
        nk <- colSums( resp ) + 10 * .Machine$double.eps
        mu <- colSums( resp * x ) / nk
        v  <- colSums( resp * outer( x, mu, "-" )^2 ) / nk + reg
        list( weight = nk / n, mean = mu, var = v )
    }
    e.step <- function( p ) {
        lp <- vapply( seq_len( k ), function( j )
            stats::dnorm( x, p$mean[ j ], sqrt( p$var[ j ] ), log = TRUE ) + log( p$weight[ j ] ),
            numeric( n ) )
        lp <- matrix( lp, ncol = k )
        mx <- lp[ , 1L ]
        for ( j in seq_len( k )[ -1L ] ) mx <- pmax( mx, lp[ , j ] )
        norm <- mx + log( rowSums( exp( lp - mx ) ) )
        list( loglik = norm, resp = exp( lp - norm ) )
    }

    p <- m.step( resp )
    lower <- -Inf
    for ( it in seq_len( max.iter ) ) {
        e <- e.step( p )
        p <- m.step( e$resp )
        lb <- mean( e$loglik )
        if ( abs( lb - lower ) < tol ) break
        lower <- lb
    }
    ll <- sum( e.step( p )$loglik )
    list( mean = p$mean, weight = p$weight, sd = sqrt( p$var ), loglik = ll,
          bic = -2 * ll + ( 3 * k - 1 ) * log( n ) )
}


#' Gaussian mixture crossover threshold
#'
#' @title Mixture Threshold
#'
#' @description
#' Point between two components' means where their weighted densities are
#' equal: the lowest and highest component (\code{cut = "outer"}), the
#' lowest two (\code{"lower"}) or the highest two (\code{"upper"}); with two
#' components the three are the same.  Falls back to the midpoint of the
#' means when the densities do not cross between them.
#'
#' @param x Numeric vector.
#' @param n.components 2 or 3.
#' @param cut \code{"outer"}, \code{"lower"} or \code{"upper"}.
#'
#' @return Numeric threshold.
#'
#' @export
threshold.mixture <- function( x, n.components = 2L, cut = "outer" ) {
    x <- as.numeric( x )
    if ( length( x ) < 20L )
        return( if ( length( x ) ) stats::median( x ) else 0 )
    fit <- fit.mixture( x, n.components )
    o <- order( fit$mean )
    k <- length( o )
    ij <- switch( cut, lower = c( 1L, 2L ), upper = c( k - 1L, k ), c( 1L, k ) )
    a <- o[ ij[ 1 ] ]
    b <- o[ ij[ 2 ] ]
    f <- function( v ) fit$weight[ a ] * stats::dnorm( v, fit$mean[ a ], fit$sd[ a ] ) -
        fit$weight[ b ] * stats::dnorm( v, fit$mean[ b ], fit$sd[ b ] )
    tryCatch(
        stats::uniroot( f, c( fit$mean[ a ], fit$mean[ b ] ), tol = 1e-12 )$root,
        error = function( e ) ( fit$mean[ a ] + fit$mean[ b ] ) / 2 )
}


#' Run a threshold algorithm
#'
#' @title Compute Threshold
#'
#' @param x Numeric vector (transformed scale).
#' @param algorithm One of \code{ag.threshold.algorithms}; unknown names use
#'   \code{"tail"}.
#' @param params Named list: \code{tail.fraction},
#'   \code{mixture.components}, \code{mixture.cut}.
#'
#' @return Numeric threshold.
#'
#' @export
compute.threshold <- function( x, algorithm, params = list() ) {
    tf  <- as.numeric( params[[ "tail.fraction" ]] %||% ag.defaults$tail.fraction )
    nc  <- as.integer( params[[ "mixture.components" ]] %||% 2L )
    cut <- as.character( params[[ "mixture.cut" ]] %||% "outer" )
    switch( algorithm,
        kde_min = threshold.kde.min( x ),
        mixture = threshold.mixture( x, nc, cut ),
        otsu    = threshold.otsu( x ),
        bimodal = threshold.kde.min( x ),
        threshold.tail( x, tf ) )
}


#' Threshold factor of one axis
#'
#' @title Threshold Factor
#'
#' @param gate.param Named list of gate parameters.
#' @param suffix \code{""}, \code{".x"} or \code{".y"}.
#'
#' @return \code{threshold.factor<suffix>}, else \code{threshold.factor},
#'   else 1.
#'
#' @keywords internal
threshold.factor <- function( gate.param, suffix = "" ) {
    v <- gate.param[[ paste0( "threshold.factor", suffix ) ]] %||%
        gate.param[[ "threshold.factor" ]] %||% ag.defaults$threshold.factor
    v <- suppressWarnings( as.numeric( v ) )
    if ( length( v ) == 1L && !is.na( v ) ) v else ag.defaults$threshold.factor
}


#' Multiply a threshold in raw units
#'
#' @title Nudge Threshold
#'
#' @description
#' Converts a display-scale threshold to raw units, multiplies it by
#' \code{factor} (2 = twice the intensity, i.e. stricter), and converts it
#' back.  A threshold at or below zero in raw units is left unchanged.
#'
#' @param t Threshold on the display scale.
#' @param factor Multiplier.
#' @param transform.params Transform parameters of the channel (see
#'   \code{\link{make.transform}}); \code{NULL} for untransformed data.
#'
#' @return List with \code{threshold} and \code{note} (empty unless the
#'   factor could not be applied).
#'
#' @keywords internal
nudge.threshold <- function( t, factor, transform.params = NULL ) {
    if ( factor == 1 || !is.finite( factor ) || factor <= 0 )
        return( list( threshold = t, note = "" ) )
    tf  <- make.transform( transform.params )
    raw <- tf$inverse( t )
    if ( !( raw > 0 ) )
        return( list( threshold = t, note = sprintf(
            "threshold %.4g is not above zero in raw units; factor %g not applied", raw, factor ) ) )
    list( threshold = tf$forward( raw * factor ), note = "" )
}


#' Evidence for two populations rather than one
#'
#' @title Bimodality Delta BIC
#'
#' @description
#' BIC of a one-component minus a two-component Gaussian mixture on the
#' quantile-trimmed values (a seeded subsample of at most
#' \code{max.events}).  Large positive values favour two populations.
#'
#' @param x Numeric vector.
#' @param max.events Largest number of events fitted.
#' @param trim Quantile trim.
#' @param seed Seed of the subsample (the caller's random stream is
#'   restored).
#'
#' @return Numeric \eqn{\Delta}BIC.
#'
#' @keywords internal
bimodality.delta.bic <- function( x, max.events = 5000L, trim = 0.001, seed = 0L ) {
    x <- as.numeric( x )
    x <- trim.events( x[ is.finite( x ) ], trim )
    if ( length( x ) < 50L )
        return( 0 )
    if ( length( x ) > max.events )
        x <- x[ seeded.call( seed, function() sample.int( length( x ), max.events ) ) ]
    fit.mixture( x, 1L )$bic - fit.mixture( x, 2L )$bic
}


#' Suggest a threshold algorithm for a gate
#'
#' @title Recommend Algorithm
#'
#' @description
#' Threshold gates get \code{"mixture"} when any separation axis looks
#' bimodal (\eqn{\Delta}BIC above \code{min.delta.bic}), otherwise
#' \code{"tail"}.  Other gate types keep their algorithm.
#'
#' @param gate.def Gate definition.
#' @param data Matrix of the parent population's events (columns named by
#'   channel).
#' @param min.delta.bic \eqn{\Delta}BIC threshold.
#'
#' @return Character algorithm.
#'
#' @export
recommend.algorithm <- function( gate.def, data, min.delta.bic = 10 ) {
    gate.type <- gate.type.of( gate.def )
    if ( !( gate.type %in% ag.threshold.types ) )
        return( if ( has.value( gate.def$algorithm ) ) gate.def$algorithm else
                    if ( gate.type == "free" ) "all" else gate.type )
    axes <- gate.def$gate.marker.x
    if ( gate.type == "2dsep" ) axes <- c( axes, gate.def$gate.marker.y )
    for ( ch in axes )
        if ( ch %in% colnames( data ) && nrow( data ) >= 100L &&
             bimodality.delta.bic( data[ , ch ] ) > min.delta.bic )
            return( "mixture" )
    "tail"
}
