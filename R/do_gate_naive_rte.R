# do_gate_naive_rte.r

do.gate.naive.rte <- function( fcs.exp.data, fcs.gates, fcs.calculate,
                               figure.dir )
{
  fcs.gate.label <- "T-cell Recent Thymic Emigrants"
  fcs.gate.marker <- c( "CCR7", "CD31" )

  fcs.gate.parent <- list( c( "memory.cd4", "naive" ),
                           c( "memory.cd8", "naive" ) )

  fcs.gate.parent.idx <- sort( unlist( sapply( fcs.gate.parent,
                                               function( fgp ) fcs.gates[[ fgp[ 1 ] ]][[ fgp[ 2 ] ]]$idx ) ) )

  fcs.gate.data <- fcs.exp.data[ fcs.gate.parent.idx, fcs.gate.marker,
                                 drop = FALSE ]

  fcs.gate.data.dens <- fcs.gate.data[
    fcs.gate.data[ , 1 ] > quantile( fcs.gate.data[ , 1 ], 1e-5 ) &
      fcs.gate.data[ , 1 ] < quantile( fcs.gate.data[ , 1 ], 1 - 1e-5 ) &
      fcs.gate.data[ , 2 ] > quantile( fcs.gate.data[ , 2 ], 1e-5 ) &
      fcs.gate.data[ , 2 ] < quantile( fcs.gate.data[ , 2 ], 1 - 1e-5 ), ,
    drop = FALSE ]

  if ( fcs.calculate )
  {
    fcs.gate.data.calc <- fcs.gate.data[
      fcs.gate.data[ , 1 ] > quantile( fcs.gate.data[ , 1 ], 1e-5 ) &
        fcs.gate.data[ , 1 ] < quantile( fcs.gate.data[ , 1 ], 1 - 1e-5 ) &
        fcs.gate.data[ , 2 ] > quantile( fcs.gate.data[ , 2 ], 1e-5 ) &
        fcs.gate.data[ , 2 ] < quantile( fcs.gate.data[ , 2 ], 1 - 1e-5 ), ]

    fcs.gate.calc.y <- rangeGate( flowFrame( fcs.gate.data.calc ),
                                  fcs.gate.marker[ 2 ], borderQuant = 0, alpha = "min" )
    fcs.gate.ythr <- fcs.gate.calc.y@min

    fcs.gate.xmin <- 1.01 * min( fcs.gate.data.calc[ , 1 ] ) +
      - 0.01 * max( fcs.gate.data.calc[ , 1 ] )
    fcs.gate.xmax <- -0.01 * min( fcs.gate.data.calc[ , 1 ] ) +
      1.01 * max( fcs.gate.data.calc[ , 1 ] )
    fcs.gate.ymin <- fcs.gate.ythr
    fcs.gate.ymax <- -0.01 * min( fcs.gate.data.calc[ , 2 ] ) +
      1.01 * max( fcs.gate.data.calc[ , 2 ] )

    fcs.gate.boundary <- rbind(
      c( fcs.gate.xmin, fcs.gate.ymin ),
      c( fcs.gate.xmax, fcs.gate.ymin ),
      c( fcs.gate.xmax, fcs.gate.ymax ),
      c( fcs.gate.xmin, fcs.gate.ymax ) )
    colnames( fcs.gate.boundary ) <- fcs.gate.marker

    fcs.gates$naive.rte <- list(
      population = list( parent = fcs.gate.parent,
                         boundary = fcs.gate.boundary ) )
  }

  fcs.gate.boundary <- fcs.gates$naive.rte$population$boundary

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

  fcs.gate.boundary.ggp <- data.frame(
    x = c( fcs.gate.boundary[ , 1 ], fcs.gate.boundary[ 1, 1 ] ),
    y = c( fcs.gate.boundary[ , 2 ], fcs.gate.boundary[ 1, 2 ] ) )

  x.transform <- fcs.transform[[ fcs.gate.marker[ 1 ] ]]
  y.transform <- fcs.transform[[ fcs.gate.marker[ 2 ] ]]

  fig.limit.x <- x.transform( c( -1e2, 3e5 ) )
  fig.limit.y <- y.transform( c( -6e2, 1e5 ) )
  fig.gate.label <- fcs.gate.label
  fig.gate.label.pos <- get.gate.label.pos( 1, fcs.gate.boundary,
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
      breaks = y.transform( c( -10^2, 0, 10^(2:5) ) ),
      labels = c( expression( -10^2 ), 0, expression( 10^2 ),
                  expression( 10^3 ), expression( 10^4 ), expression( 10^5 ) ),
      limits = fig.limit.y ) +
    ggtitle( "CD4+ and CD8+ Naive T cells - only for calculation of gate" ) +
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

  ggsave( file.path( figure.dir, "naive_rte.png" ),
          width = 9, height = 7.5 )

  fcs.gates
}
