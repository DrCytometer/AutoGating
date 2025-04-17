# do_gate_foxp3_cd4.r

do.gate.foxp3.cd4 <- function( fcs.exp.data, fcs.gates, fcs.calculate,
                               figure.dir )
{
  fcs.gate.a.label <- "Treg"
  fcs.gate.b.label <- "Non-Treg"

  fcs.gate.marker <- c( "CD4", "FOXP3" )

  fcs.gate.parent <- c( "cd4.cd8", "cd4" )

  fcs.gate.parent.idx <- fcs.gates[[ fcs.gate.parent[ 1 ] ]][[
    fcs.gate.parent[ 2 ] ]]$idx

  fcs.gate.data <- fcs.exp.data[ fcs.gate.parent.idx, fcs.gate.marker ]

  fcs.gate.data.dens <- fcs.gate.data[
    fcs.gate.data[ , 1 ] > quantile( fcs.gate.data[ , 1 ], 1e-5 ) &
      fcs.gate.data[ , 1 ] < quantile( fcs.gate.data[ , 1 ], 1 - 1e-5 ) &
      fcs.gate.data[ , 2 ] > quantile( fcs.gate.data[ , 2 ], 1e-5 ) &
      fcs.gate.data[ , 2 ] < quantile( fcs.gate.data[ , 2 ], 1 - 1e-5 ), ]

  if ( fcs.calculate )
  {
    fcs.gate.data.calc <- fcs.gate.data[
      fcs.gate.data[ , 1 ] > quantile( fcs.gate.data[ , 1 ], 1e-5 ) &
        fcs.gate.data[ , 1 ] < quantile( fcs.gate.data[ , 1 ], 1 - 1e-5 ) &
        fcs.gate.data[ , 2 ] > quantile( fcs.gate.data[ , 2 ], 1e-5 ) &
        fcs.gate.data[ , 2 ] < quantile( fcs.gate.data[ , 2 ], 1 - 1e-5 ), ]

    # fcs.gate.calc.y <- rangeGate( flowFrame( - fcs.gate.data.calc ),
    #     fcs.gate.marker[ 2 ], borderQuant = 0, alpha = "min" )
    # fcs.gate.ythr <- - fcs.gate.calc.y@min
    fcs.gate.ythr <- gate.tail( fcs.gate.data.calc[ , 2 ], 0.99 )

    fcs.gate.a.xmin <- 1.05 * min( fcs.gate.data[ , 1 ] ) +
      -0.05 * max( fcs.gate.data[ , 1 ] )
    fcs.gate.a.xmax <- -0.05 * min( fcs.gate.data[ , 1 ] ) +
      1.05 * max( fcs.gate.data[ , 1 ] )
    fcs.gate.a.ymin <- fcs.gate.ythr
    fcs.gate.a.ymax <- -0.01 * min( fcs.gate.data[ , 2 ] ) +
      1.01 * max( fcs.gate.data[ , 2 ] )

    # note manual threshold for b.ymin
    fcs.gate.b.xmin <- 1.05 * min( fcs.gate.data[ , 1 ] ) +
      -0.05 * max( fcs.gate.data[ , 1 ] )
    fcs.gate.b.xmax <- -0.05 * min( fcs.gate.data[ , 1 ] ) +
      1.05 * max( fcs.gate.data[ , 1 ] )
    fcs.gate.b.ymin <- as.numeric( quantile( fcs.gate.data[ , 2 ], 2e-5 ) )
    fcs.gate.b.ymax <- fcs.gate.ythr

    fcs.gate.a.boundary <- rbind(
      c( fcs.gate.a.xmin, fcs.gate.a.ymin ),
      c( fcs.gate.a.xmax, fcs.gate.a.ymin ),
      c( fcs.gate.a.xmax, fcs.gate.a.ymax ),
      c( fcs.gate.a.xmin, fcs.gate.a.ymax ) )
    colnames( fcs.gate.a.boundary ) <- fcs.gate.marker

    fcs.gate.b.boundary <- rbind(
      c( fcs.gate.b.xmin, fcs.gate.b.ymin ),
      c( fcs.gate.b.xmax, fcs.gate.b.ymin ),
      c( fcs.gate.b.xmax, fcs.gate.b.ymax ),
      c( fcs.gate.b.xmin, fcs.gate.b.ymax ) )
    colnames( fcs.gate.b.boundary ) <- fcs.gate.marker

    fcs.gates$foxp3.cd4 <- list(
      treg = list( parent = fcs.gate.parent,
                   boundary = fcs.gate.a.boundary ),
      non.treg = list( parent = fcs.gate.parent,
                       boundary = fcs.gate.b.boundary ) )
  }

  fcs.gate.a.boundary <- fcs.gates$foxp3.cd4$treg$boundary
  fcs.gate.b.boundary <- fcs.gates$foxp3.cd4$non.treg$boundary

  fcs.gate.a.pip <- point.in.polygon(
    fcs.gate.data[ , 1 ], fcs.gate.data[ , 2 ],
    fcs.gate.a.boundary[ , 1 ], fcs.gate.a.boundary[ , 2 ] )

  fcs.gate.b.pip <- point.in.polygon(
    fcs.gate.data[ , 1 ], fcs.gate.data[ , 2 ],
    fcs.gate.b.boundary[ , 1 ], fcs.gate.b.boundary[ , 2 ] )

  fcs.gate.a.fraction <- mean( fcs.gate.a.pip == 1 )
  fcs.gate.b.fraction <- mean( fcs.gate.b.pip == 1 )

  fcs.gates$foxp3.cd4$treg$idx <- fcs.gate.parent.idx[ fcs.gate.a.pip == 1 ]
  fcs.gates$foxp3.cd4$non.treg$idx <- fcs.gate.parent.idx[
    fcs.gate.b.pip == 1 ]

  fcs.gate.data.ggp <- data.frame(
    x = fcs.gate.data[ , fcs.gate.marker[ 1 ] ],
    y = fcs.gate.data[ , fcs.gate.marker[ 2 ] ],
    z = interp.surface( kde2d(
      fcs.gate.data.dens[ , fcs.gate.marker[ 1 ] ],
      fcs.gate.data.dens[ , fcs.gate.marker[ 2 ] ], n = 100 ),
      fcs.gate.data ) )

  fcs.gate.a.boundary.ggp <- data.frame(
    x = c( fcs.gate.a.boundary[ , 1 ], fcs.gate.a.boundary[ 1, 1 ] ),
    y = c( fcs.gate.a.boundary[ , 2 ], fcs.gate.a.boundary[ 1, 2 ] ) )

  fcs.gate.b.boundary.ggp <- data.frame(
    x = c( fcs.gate.b.boundary[ , 1 ], fcs.gate.b.boundary[ 1, 1 ] ),
    y = c( fcs.gate.b.boundary[ , 2 ], fcs.gate.b.boundary[ 1, 2 ] ) )

  x.transform <- fcs.transform[[ fcs.gate.marker[ 1 ] ]]
  y.transform <- fcs.transform[[ fcs.gate.marker[ 2 ] ]]

  fig.limit.x <- x.transform( c( -2e2, 265e3 ) )
  fig.limit.y <- y.transform( c( -2e3, 265e3 ) )

  fig.gate.a.label <- sprintf( "%s\n%.1f%%", fcs.gate.a.label,
                               100 * fcs.gate.a.fraction )
  fig.gate.b.label <- sprintf( "%s\n%.1f%%", fcs.gate.b.label,
                               100 * fcs.gate.b.fraction )

  fig.gate.a.label.pos <- get.gate.label.pos( 5, fcs.gate.a.boundary,
                                              fig.limit.x, fig.limit.y, fig.gate.a.label )
  fig.gate.b.label.pos <- get.gate.label.pos( 8, fcs.gate.b.boundary,
                                              fig.limit.x, fig.limit.y, fig.gate.b.label )

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
    geom_path( aes( x, y, color = NULL ), data = fcs.gate.a.boundary.ggp ) +
    geom_path( aes( x, y, color = NULL ), data = fcs.gate.b.boundary.ggp ) +
    theme_bw() +
    theme( panel.grid.major = element_blank(),
           panel.grid.minor = element_blank() ) +
    annotate( "text", label = fig.gate.a.label, size = 5,
              x = fig.gate.a.label.pos[ 1 ], y = fig.gate.a.label.pos[ 2 ] ) +
    annotate( "text", label = fig.gate.b.label, size = 5,
              x = fig.gate.b.label.pos[ 1 ], y = fig.gate.b.label.pos[ 2 ] )

  ggsave( file.path( figure.dir, "foxp3_cd4.png" ),
          width = 9, height = 7.5 )

  fcs.gates
}
