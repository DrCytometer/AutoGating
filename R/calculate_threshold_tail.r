# calculate_threshold_tail.r

calculate.threshold.tail <- function( marker.data, threshold.param )
{
    x <- marker.data

    if ( is.null( threshold.param ) )
        tail.prob <- fcs.default.tail.probability
    else
        tail.prob <- as.numeric( threshold.param )

    x.mode <- mlv( marker.data, method = "shorth" )

    x.symm.low <- x[ x < x.mode ]
    x.symm.high <- 2 * x.mode - x.symm.low
    x.symm <- c( x.symm.low, x.symm.high )

    x.symm.rev.high <- x[ x > x.mode ]
    x.symm.rev.low <- 2 * x.mode - x.symm.rev.high
    x.symm.rev <- c( x.symm.rev.low, x.symm.rev.high )

    x.thr <- as.numeric( quantile( x.symm, tail.prob ) )
    x.thr.rev <- as.numeric( quantile( x.symm.rev, tail.prob ) )

    min( x.thr, x.thr.rev )
}

