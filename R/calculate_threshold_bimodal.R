# calculate_threshold_bimodal.r

calculate.threshold.bimodal <- function( marker.data, threshold.param, agp )
{
  flow.data <- as.matrix( marker.data )
  colnames( flow.data ) <- "marker"

  threshold.calc <- flowStats::rangeGate( flowCore::flowFrame( flow.data ), "marker",
                               borderQuant = 0, alpha = "min" )

  as.numeric( threshold.calc@min )
}
