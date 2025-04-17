# calculate_boundary_all.r

calculate.boundary.all <- function( marker.data, boundary.param )
{
    marker.data.xmin <- min( marker.data[ , 1 ] )
    marker.data.xmax <- max( marker.data[ , 1 ] )
    marker.data.ymin <- min( marker.data[ , 2 ] )
    marker.data.ymax <- max( marker.data[ , 2 ] )

    gate.xmin <- + 1.01 * marker.data.xmin - 0.01 * marker.data.xmax
    gate.xmax <- - 0.01 * marker.data.xmin + 1.01 * marker.data.xmax
    gate.ymin <- + 1.01 * marker.data.ymin - 0.01 * marker.data.ymax
    gate.ymax <- - 0.01 * marker.data.ymin + 1.01 * marker.data.ymax

    gate.boundary <- rbind(
        c( gate.xmin, gate.ymin ),
        c( gate.xmax, gate.ymin ),
        c( gate.xmax, gate.ymax ),
        c( gate.xmin, gate.ymax ) )

    gate.boundary
}
