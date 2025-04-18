# get_pop_label_pos.r

get.pop.label.pos <- function( label.pos, pop.boundary, limit.x, limit.y,
                               pop.label, agp )
{
  pop.xmax <- max( pop.boundary[ , 1 ] )
  pop.ymax <- max( pop.boundary[ , 2 ] )
  pop.xmin <- min( pop.boundary[ , 1 ] )
  pop.ymin <- min( pop.boundary[ , 2 ] )

  label.line <- scan( text = pop.label, what = character(), sep = "\n",
                      quiet = TRUE )
  label.width <- max( sapply( label.line, nchar ) )
  label.height <- length( label.line )

  if ( label.pos == 1 ) {
    label.pos.x <- pop.xmin
    label.pos.y <- pop.ymin -
      agp$fcs.figure.label.height.factor * label.height *
      ( limit.y[ 2 ] - limit.y[ 1 ] )
  }
  else if ( label.pos == 2 ) {
    label.pos.x <- pop.xmax
    label.pos.y <- pop.ymin -
      agp$fcs.figure.label.height.factor * label.height *
      ( limit.y[ 2 ] - limit.y[ 1 ] )
  }
  else if ( label.pos == 3 ) {
    label.pos.x <- pop.xmin
    label.pos.y <- pop.ymax +
      agp$fcs.figure.label.height.factor * label.height *
      ( limit.y[ 2 ] - limit.y[ 1 ] )
  }
  else if ( label.pos == 4 ) {
    label.pos.x <- pop.xmax
    label.pos.y <- pop.ymax +
      agp$fcs.figure.label.height.factor * label.height *
      ( limit.y[ 2 ] - limit.y[ 1 ] )
  }
  else if ( label.pos == 5 ) {
    label.pos.x <- pop.xmin -
      agp$fcs.figure.label.width.factor * label.width *
      ( limit.x[ 2 ] - limit.x[ 1 ] )
    label.pos.y <- pop.ymin
  }
  else if ( label.pos == 6 ) {
    label.pos.x <- pop.xmax +
      agp$fcs.figure.label.width.factor * label.width *
      ( limit.x[ 2 ] - limit.x[ 1 ] )
    label.pos.y <- pop.ymin
  }
  else if ( label.pos == 7 ) {
    label.pos.x <- pop.xmin -
      agp$fcs.figure.label.width.factor * label.width *
      ( limit.x[ 2 ] - limit.x[ 1 ] )
    label.pos.y <- pop.ymax
  }
  else if ( label.pos == 8 ) {
    label.pos.x <- pop.xmax +
      agp$fcs.figure.label.width.factor * label.width *
      ( limit.x[ 2 ] - limit.x[ 1 ] )
    label.pos.y <- pop.ymax
  }
  else {
    label.pos.x <- ( limit.x[ 1 ] + limit.x[ 2 ] ) / 2
    label.pos.y <- ( limit.y[ 1 ] + limit.y[ 2 ] ) / 2
  }

  c( label.pos.x, label.pos.y )
}
