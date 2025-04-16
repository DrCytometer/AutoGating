# get_density_palette.r

get.density.palette <- function( fcs.density )
{
  rainbow.palette <- colorRampPalette(
    c( "blue", "cyan", "green", "yellow", "red" ) )( fcs.palette.base.n )

  fcs.density <- na.omit( fcs.density )

  fcs.density.grid <- seq( min( fcs.density ), max( fcs.density ),
                           length.out = fcs.palette.n )

  density.palette.idx <-
    round( ecdf( fcs.density )( fcs.density.grid ) * fcs.palette.base.n )

  rainbow.palette[ density.palette.idx ]
}
