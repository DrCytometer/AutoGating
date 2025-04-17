# do_gate_memory.r

do.gate.memory <- function( fcs.exp.data, fcs.gates, fcs.calculate,
                            figure.dir )
{
  fcs.gate.a.label <- "Naive T cells"
  fcs.gate.b.label <- "Central Memory T cells"
  fcs.gate.c.label <- "Effector Memory T cells"
  fcs.gate.d.label <- "Effector Memory RA T cells"

  fcs.gate.marker <- c( "CD45RA", "CCR7" )

  fcs.gate.parent <- list( c( "cd4.cd8", "cd4" ), c( "cd4.cd8", "cd8" ) )

  fcs.gate.parent.idx <- sort( unlist( sapply( fcs.gate.parent,
                                               function( fgp ) fcs.gates[[ fgp[ 1 ] ]][[ fgp[ 2 ] ]]$idx ) ) )

  fcs.gate.data <- fcs.exp.data[ fcs.gate.parent.idx, fcs.gate.marker ]

  fcs.gate.data.dens <- fcs.gate.data[
    fcs.gate.data[ , 1 ] > quantile( fcs.gate.data[ , 1 ], 1e-5 ) &
      fcs.gate.data[ , 1 ] < quantile( fcs.gate.data[ , 1 ], 1 - 1e-5 ) &
      fcs.gate.data[ , 2 ] > quantile( fcs.gate.data[ , 2 ], 1e-5 ) &
      fcs.gate.data[ , 2 ] < quantile( fcs.gate.data[ , 2 ], 1 - 1e-5 ), ]

  if ( fcs.calculate )
  {
    fcs.gate.data.calc <- fcs.gate.data[
      fcs.gate.data[ , 1 ] > quantile( fcs.gate.data[ , 1 ], 1e-4 ) &
        fcs.gate.data[ , 1 ] < quantile( fcs.gate.data[ , 1 ], 0.9 ) &
        fcs.gate.data[ , 2 ] > quantile( fcs.gate.data[ , 2 ], 1e-4 ) &
        fcs.gate.data[ , 2 ] < quantile( fcs.gate.data[ , 2 ], 0.9 ), ]

    fcs.gate.calc.x <- rangeGate( flowFrame( fcs.gate.data.calc ),
                                  fcs.gate.marker[ 1 ], borderQuant = 0, alpha = "min" )
    fcs.gate.calc.y <- rangeGate( flowFrame( fcs.gate.data.calc ),
                                  fcs.gate.marker[ 2 ], borderQuant = 0, alpha = "min" )

    fcs.gate.xthr <- fcs.gate.calc.x@min
    fcs.gate.ythr <- fcs.gate.calc.y@min

    fcs.gate.a.xmin <- fcs.gate.xthr
    fcs.gate.a.xmax <- -0.01 * min( fcs.gate.data.dens[ , 1 ] ) +
      1.01 * max( fcs.gate.data.dens[ , 1 ] )
    fcs.gate.a.ymin <- fcs.gate.ythr
    fcs.gate.a.ymax <- -0.01 * min( fcs.gate.data.dens[ , 2 ] ) +
      1.01 * max( fcs.gate.data.dens[ , 2 ] )

    fcs.gate.b.xmin <- 1.01 * min( fcs.gate.data.dens[ , 1 ] ) +
      -0.01 * max( fcs.gate.data.dens[ , 1 ] )
    fcs.gate.b.xmax <- fcs.gate.xthr
    fcs.gate.b.ymin <- fcs.gate.ythr
    fcs.gate.b.ymax <- -0.01 * min( fcs.gate.data.dens[ , 2 ] ) +
      1.01 * max( fcs.gate.data.dens[ , 2 ] )

    fcs.gate.c.xmin <- 1.01 * min( fcs.gate.data.dens[ , 1 ] ) +
      -0.01 * max( fcs.gate.data.dens[ , 1 ] )
    fcs.gate.c.xmax <- fcs.gate.xthr
    fcs.gate.c.ymin <- 1.01 * min( fcs.gate.data.calc[ , 2 ] ) + # note calc
      -0.01 * max( fcs.gate.data.dens[ , 2 ] )
    fcs.gate.c.ymax <- fcs.gate.ythr

    fcs.gate.d.xmin <- fcs.gate.xthr
    fcs.gate.d.xmax <- -0.01 * min( fcs.gate.data.dens[ , 1 ] ) +
      1.01 * max( fcs.gate.data.dens[ , 1 ] )
    fcs.gate.d.ymin <- 1.01 * min( fcs.gate.data.calc[ , 2 ] ) + # note calc
      -0.01 * max( fcs.gate.data.dens[ , 2 ] )
    fcs.gate.d.ymax <- fcs.gate.ythr

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

    fcs.gate.c.boundary <- rbind(
      c( fcs.gate.c.xmin, fcs.gate.c.ymin ),
      c( fcs.gate.c.xmax, fcs.gate.c.ymin ),
      c( fcs.gate.c.xmax, fcs.gate.c.ymax ),
      c( fcs.gate.c.xmin, fcs.gate.c.ymax ) )
    colnames( fcs.gate.c.boundary ) <- fcs.gate.marker

    fcs.gate.d.boundary <- rbind(
      c( fcs.gate.d.xmin, fcs.gate.d.ymin ),
      c( fcs.gate.d.xmax, fcs.gate.d.ymin ),
      c( fcs.gate.d.xmax, fcs.gate.d.ymax ),
      c( fcs.gate.d.xmin, fcs.gate.d.ymax ) )
    colnames( fcs.gate.d.boundary ) <- fcs.gate.marker

    fcs.gates$memory <- list(
      naive = list( parent = fcs.gate.parent,
                    boundary = fcs.gate.a.boundary ),
      cenmem = list( parent = fcs.gate.parent,
                     boundary = fcs.gate.b.boundary ),
      effmem = list( parent = fcs.gate.parent,
                     boundary = fcs.gate.c.boundary ),
      temra = list( parent = fcs.gate.parent,
                    boundary = fcs.gate.d.boundary ) )
  }

  fcs.gate.a.boundary <- fcs.gates$memory$naive$boundary
  fcs.gate.b.boundary <- fcs.gates$memory$cenmem$boundary
  fcs.gate.c.boundary <- fcs.gates$memory$effmem$boundary
  fcs.gate.d.boundary <- fcs.gates$memory$temra$boundary

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

  fig.gate.a.label <- fcs.gate.a.label
  fig.gate.b.label <- fcs.gate.b.label
  fig.gate.c.label <- fcs.gate.c.label
  fig.gate.d.label <- fcs.gate.d.label

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
    ggtitle( "CD4+ and CD8+ T cells - only for calculation of gates" ) +
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

  ggsave( file.path( figure.dir, "memory.png" ),
          width = 9, height = 7.5 )

  fcs.gates
}
