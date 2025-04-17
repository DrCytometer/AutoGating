# do_gate_memory_cd8.r

do.gate.memory.cd8 <- function( fcs.exp.data, fcs.gates, fcs.calculate,
                                figure.dir )
{
  fcs.gate.a.label <- "Naive CD8+ T cells"
  fcs.gate.b.label <- "Central Memory CD8+ T cells"
  fcs.gate.c.label <- "Effector Memory CD8+ T cells"
  fcs.gate.d.label <- "Effector Memory RA CD8+ T cells"

  fcs.gate.marker <- c( "CD45RA", "CCR7" )

  fcs.gate.parent <- c( "cd4.cd8", "cd8" )

  fcs.gate.parent.idx <- fcs.gates[[ fcs.gate.parent[ 1 ] ]][[
    fcs.gate.parent[ 2 ] ]]$idx

  fcs.gate.data <- fcs.exp.data[ fcs.gate.parent.idx, fcs.gate.marker ]

  fcs.gate.data.dens <- fcs.gate.data[
    fcs.gate.data[ , 1 ] > quantile( fcs.gate.data[ , 1 ], 1e-5 ) &
      fcs.gate.data[ , 1 ] < quantile( fcs.gate.data[ , 1 ], 1 - 1e-5 ) &
      fcs.gate.data[ , 2 ] > quantile( fcs.gate.data[ , 2 ], 1e-5 ) &
      fcs.gate.data[ , 2 ] < quantile( fcs.gate.data[ , 2 ], 1 - 1e-5 ), ]

  fcs.gate.a.boundary <- fcs.gates$memory$naive$boundary
  fcs.gate.b.boundary <- fcs.gates$memory$cenmem$boundary
  fcs.gate.c.boundary <- fcs.gates$memory$effmem$boundary
  fcs.gate.d.boundary <- fcs.gates$memory$temra$boundary

  fcs.gate.a.pip <- point.in.polygon(
    fcs.gate.data[ , 1 ], fcs.gate.data[ , 2 ],
    fcs.gate.a.boundary[ , 1 ], fcs.gate.a.boundary[ , 2 ] )

  fcs.gate.b.pip <- point.in.polygon(
    fcs.gate.data[ , 1 ], fcs.gate.data[ , 2 ],
    fcs.gate.b.boundary[ , 1 ], fcs.gate.b.boundary[ , 2 ] )

  fcs.gate.c.pip <- point.in.polygon(
    fcs.gate.data[ , 1 ], fcs.gate.data[ , 2 ],
    fcs.gate.c.boundary[ , 1 ], fcs.gate.c.boundary[ , 2 ] )

  fcs.gate.d.pip <- point.in.polygon(
    fcs.gate.data[ , 1 ], fcs.gate.data[ , 2 ],
    fcs.gate.d.boundary[ , 1 ], fcs.gate.d.boundary[ , 2 ] )

  fcs.gate.a.fraction <- mean( fcs.gate.a.pip == 1 )
  fcs.gate.b.fraction <- mean( fcs.gate.b.pip == 1 )
  fcs.gate.c.fraction <- mean( fcs.gate.c.pip == 1 )
  fcs.gate.d.fraction <- mean( fcs.gate.d.pip == 1 )

  fcs.gates$memory.cd8 <- list(
    naive = list( parent = fcs.gate.parent,
                  boundary = fcs.gate.a.boundary,
                  idx = fcs.gate.parent.idx[ fcs.gate.a.pip == 1 ] ),
    cenmem = list( parent = fcs.gate.parent,
                   boundary = fcs.gate.b.boundary,
                   idx = fcs.gate.parent.idx[ fcs.gate.b.pip == 1 ] ),
    effmem = list( parent = fcs.gate.parent,
                   boundary = fcs.gate.c.boundary,
                   idx = fcs.gate.parent.idx[ fcs.gate.c.pip == 1 ] ),
    temra = list( parent = fcs.gate.parent,
                  boundary = fcs.gate.d.boundary,
                  idx = fcs.gate.parent.idx[ fcs.gate.d.pip == 1 ] )
  )

  if( nrow( fcs.gate.data ) == 0 )
    return( fcs.gates )

  if ( nrow( fcs.gate.data.dens ) > 0 )
    fcs.gate.data.ggp <- data.frame(
      x = fcs.gate.data[ , fcs.gate.marker[ 1 ] ],
      y = fcs.gate.data[ , fcs.gate.marker[ 2 ] ],
      z = interp.surface( kde2d(
        fcs.gate.data.dens[ , fcs.gate.marker[ 1 ] ],
        fcs.gate.data.dens[ , fcs.gate.marker[ 2 ] ], n = 100 ),
        fcs.gate.data ) )
  else
    fcs.gate.data.ggp <- data.frame(
      x = fcs.gate.data[ , fcs.gate.marker[ 1 ] ],
      y = fcs.gate.data[ , fcs.gate.marker[ 2 ] ],
      z = rep( 0, nrow( fcs.gate.data ) ) )

  fcs.gate.a.boundary.ggp <- data.frame(
    x = c( fcs.gate.a.boundary[ , 1 ], fcs.gate.a.boundary[ 1, 1 ] ),
    y = c( fcs.gate.a.boundary[ , 2 ], fcs.gate.a.boundary[ 1, 2 ] ) )

  fcs.gate.b.boundary.ggp <- data.frame(
    x = c( fcs.gate.b.boundary[ , 1 ], fcs.gate.b.boundary[ 1, 1 ] ),
    y = c( fcs.gate.b.boundary[ , 2 ], fcs.gate.b.boundary[ 1, 2 ] ) )

  fcs.gate.c.boundary.ggp <- data.frame(
    x = c( fcs.gate.c.boundary[ , 1 ], fcs.gate.c.boundary[ 1, 1 ] ),
    y = c( fcs.gate.c.boundary[ , 2 ], fcs.gate.c.boundary[ 1, 2 ] ) )

  fcs.gate.d.boundary.ggp <- data.frame(
    x = c( fcs.gate.d.boundary[ , 1 ], fcs.gate.d.boundary[ 1, 1 ] ),
    y = c( fcs.gate.d.boundary[ , 2 ], fcs.gate.d.boundary[ 1, 2 ] ) )

  x.transform <- fcs.transform[[ fcs.gate.marker[ 1 ] ]]
  y.transform <- fcs.transform[[ fcs.gate.marker[ 2 ] ]]

  fig.limit.x <- x.transform( c( -2e3, 5e4 ) )
  fig.limit.y <- y.transform( c( -1e3, 1e5 ) )

  fig.gate.a.label <- sprintf( "%s\n%.1f%%", fcs.gate.a.label,
                               100 * fcs.gate.a.fraction )
  fig.gate.b.label <- sprintf( "%s\n%.1f%%", fcs.gate.b.label,
                               100 * fcs.gate.b.fraction )
  fig.gate.c.label <- sprintf( "%s\n%.1f%%", fcs.gate.c.label,
                               100 * fcs.gate.c.fraction )
  fig.gate.d.label <- sprintf( "%s\n%.1f%%", fcs.gate.d.label,
                               100 * fcs.gate.d.fraction )

  fig.gate.a.label.pos <- get.gate.label.pos( 1, fcs.gate.a.boundary,
                                              fig.limit.x, fig.limit.y, fig.gate.a.label )
  fig.gate.b.label.pos <- get.gate.label.pos( 2, fcs.gate.b.boundary,
                                              fig.limit.x, fig.limit.y, fig.gate.b.label )
  fig.gate.c.label.pos <- get.gate.label.pos( 3, fcs.gate.c.boundary,
                                              fig.limit.x, fig.limit.y, fig.gate.c.label )
  fig.gate.d.label.pos <- get.gate.label.pos( 4, fcs.gate.d.boundary,
                                              fig.limit.x, fig.limit.y, fig.gate.d.label )

  density.palette <- get.density.palette( fcs.gate.data.ggp$z )

  ggplot( fcs.gate.data.ggp, aes( x, y, color = z ) ) +
    scale_x_continuous(
      name = fcs.marker.label.figure[ fcs.gate.marker[ 1 ] ],
      breaks = x.transform( c( -10^3, 0, 10^(3:4) ) ),
      labels = c( expression( -10^3 ), 0, expression( 10^3 ),
                  expression( 10^4 ) ),
      limits = fig.limit.x ) +
    scale_y_continuous(
      name = fcs.marker.label.figure[ fcs.gate.marker[ 2 ] ],
      breaks = y.transform( c( -10^3, 0, 10^(3:4) ) ),
      labels = c( expression( -10^3 ), 0, expression( 10^3 ),
                  expression( 10^4 ) ),
      limits = fig.limit.y ) +
    geom_point( size = 0.1 ) +
    scale_color_gradientn( "", labels = NULL,
                           colors = density.palette,
                           guide = guide_colorbar( barheight = 30 ) ) +
    geom_path( aes( x, y, color = NULL ), data = fcs.gate.a.boundary.ggp ) +
    geom_path( aes( x, y, color = NULL ), data = fcs.gate.b.boundary.ggp ) +
    geom_path( aes( x, y, color = NULL ), data = fcs.gate.c.boundary.ggp ) +
    geom_path( aes( x, y, color = NULL ), data = fcs.gate.d.boundary.ggp ) +
    theme_bw() +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank() ) +
    annotate( "text", label = fig.gate.a.label, size = 5,
              x = fig.gate.a.label.pos[ 1 ], y = fig.gate.a.label.pos[ 2 ] ) +
    annotate( "text", label = fig.gate.b.label, size = 5,
              x = fig.gate.b.label.pos[ 1 ], y = fig.gate.b.label.pos[ 2 ] ) +
    annotate( "text", label = fig.gate.c.label, size = 5,
              x = fig.gate.c.label.pos[ 1 ], y = fig.gate.c.label.pos[ 2 ] ) +
    annotate( "text", label = fig.gate.d.label, size = 5,
              x = fig.gate.d.label.pos[ 1 ], y = fig.gate.d.label.pos[ 2 ] )

  ggsave( file.path( figure.dir, "memory_cd8.png" ),
          width = 9, height = 7.5 )

  fcs.gates
}
