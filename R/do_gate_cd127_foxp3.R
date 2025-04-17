# do_gate_cd127_foxp3.r

do.gate.cd127.foxp3 <- function( fcs.exp.data, fcs.gates, fcs.calculate,
                                 figure.dir )
{
  fcs.gate.label <- "Treg"

  fcs.gate.marker <- c( "CD4", "FOXP3" )

  fcs.gate.parent <- c( "cd127", "population" )

  fcs.gate.parent.idx <- fcs.gates[[ fcs.gate.parent[ 1 ] ]][[
    fcs.gate.parent[ 2 ] ]]$idx

  fcs.gate.data <- fcs.exp.data[ fcs.gate.parent.idx, fcs.gate.marker ]

  fcs.gate.data.dens <- fcs.gate.data[
    fcs.gate.data[ , 2 ] > quantile( fcs.gate.data[ , 2 ], 1e-4 ) &
      fcs.gate.data[ , 2 ] < quantile( fcs.gate.data[ , 2 ], 1 - 1e-4 ), ]

  fcs.gate.boundary <- fcs.gates$foxp3.cd4$treg$boundary

  fcs.gates$cd127.foxp3 <- list(
    population = list( parent = fcs.gate.parent,
                       boundary = fcs.gate.boundary ) )

  fcs.gate.pip <- point.in.polygon(
    fcs.gate.data[ , 1 ], fcs.gate.data[ , 2 ],
    fcs.gate.boundary[ , 1 ], fcs.gate.boundary[ , 2 ] )

  fcs.gate.fraction <- mean( fcs.gate.pip == 1 )

  fcs.gates$cd127.foxp3$population$idx <-
    fcs.gate.parent.idx[ fcs.gate.pip == 1 ]

  fcs.gate.data.ggp <- data.frame(
    x = fcs.gate.data[ , fcs.gate.marker[ 1 ] ],
    y = fcs.gate.data[ , fcs.gate.marker[ 2 ] ],
    z = interp.surface( kde2d(
      fcs.gate.data.dens[ , fcs.gate.marker[ 1 ] ],
      fcs.gate.data.dens[ , fcs.gate.marker[ 2 ] ], n = 100 ),
      fcs.gate.data ) )

  fcs.gate.boundary.ggp <- data.frame(
    x = c( fcs.gate.boundary[ , 1 ], fcs.gate.boundary[ 1, 1 ] ),
    y = c( fcs.gate.boundary[ , 2 ], fcs.gate.boundary[ 1, 2 ] ) )

  x.transform <- fcs.transform[[ fcs.gate.marker[ 1 ] ]]
  y.transform <- fcs.transform[[ fcs.gate.marker[ 2 ] ]]

  fig.limit.x <- x.transform( c( -2e2, 265e3 ) )
  fig.limit.y <- y.transform( c( -2e3, 265e3 ) )

  fig.gate.label <- sprintf( "%s\n%.1f%%", fcs.gate.label,
                             100 * fcs.gate.fraction )
  fig.gate.label.pos <- get.gate.label.pos( 2, fcs.gate.boundary,
                                            fig.limit.x, fig.limit.y, fig.gate.label )

  density.palette <- get.density.palette( fcs.gate.data.ggp$z )

  ggplot( fcs.gate.data.ggp, aes( x, y, color = z ) ) +
    scale_x_continuous(
      name = fcs.marker.label.figure[ fcs.gate.marker[ 1 ] ],
      breaks = x.transform( c( -10^2, 0, 10^(2:5) ) ),
      labels = c( expression( -10^2 ), 0, expression( 10^2 ),
                  expression( 10^3 ), expression( 10^4 ), expression( 10^5 ) ),
      limits = fig.limit.x ) +
    scale_y_continuous(
      name = fcs.marker.label.figure[ fcs.gate.marker[ 2 ] ],
      breaks = y.transform( c( -10^(3:2), 0, 10^(2:5) ) ),
      labels = c( expression( -10^3 ), expression( -10^2 ), 0,
                  expression( 10^2 ), expression( 10^3 ), expression( 10^4 ),
                  expression( 10^5 ) ),
      limits = fig.limit.y ) +
    geom_point( size = 0.1 ) +
    scale_color_gradientn( "", labels = NULL,
                           colors = density.palette,
                           guide = guide_colorbar( barheight = 30 ) ) +
    geom_path( aes( x, y, color = NULL ), data = fcs.gate.boundary.ggp ) +
    theme_bw() +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank() ) +
    annotate( "text", label = fig.gate.label, size = 5,
              x = fig.gate.label.pos[ 1 ], y = fig.gate.label.pos[ 2 ] )

  ggsave( file.path( figure.dir, "cd127_foxp3.png" ),
          width = 9, height = 7.5 )

  fcs.gates
}
