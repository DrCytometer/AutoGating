#' Threshold by Otsu's method
#'
#' @title Otsu's Method Threshold
#'
#' @description
#' Maximises the between-class variance across 255 candidate thresholds on a
#' histogram with 256 equal-width bins.  This is the classical Otsu (1979)
#' image segmentation criterion adapted for flow cytometry data.
#'
#' @param marker.data Numeric vector.  Values to threshold.
#' @param threshold.param Not used by this method; present for interface
#'   consistency with other threshold functions.
#'
#' @return Numeric scalar.  The bin boundary that maximises between-class
#'   variance.
#'
#' @references
#' Otsu, N. (1979). A threshold selection method from gray-level histograms.
#' \emph{IEEE Transactions on Systems, Man, and Cybernetics}, 9(1), 62–66.
#' \doi{10.1109/TSMC.1979.4310076}
#'
#' @seealso \code{\link{calculate.threshold.tail}},
#'   \code{\link{calculate.gate.1dsep}}
#'
#' @keywords internal

calculate.threshold.otsu <- function( marker.data, threshold.param ) {
  breaks <- seq( min( marker.data ), max( marker.data ), length.out = 256 )
  counts <- graphics::hist( marker.data, breaks = breaks, plot = FALSE )$counts
  total  <- sum( counts )
  best.var <- -Inf
  best.t   <- breaks[2]
  for ( i in seq_along( counts ) ) {
    w0 <- sum( counts[ seq_len(i) ] ) / total
    w1 <- 1 - w0
    if ( w0 == 0 || w1 == 0 ) next
    mu0 <- sum( counts[ seq_len(i) ] * breaks[ seq_len(i) ] ) / ( w0 * total )
    mu1 <- sum( counts[ -seq_len(i) ] * breaks[ -seq_len(i) ] ) / ( w1 * total )
    bv  <- w0 * w1 * ( mu0 - mu1 )^2
    if ( bv > best.var ) { best.var <- bv; best.t <- breaks[i] }
  }
  best.t
}
