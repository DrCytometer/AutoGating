#' Compute density, build plot, overlay boundaries, and save as PNG
#'
#' @title Plot Gate and Save to File
#'
#' @description
#' Orchestrates the full figure-generation workflow for a single gate:
#' \enumerate{
#'   \item Optionally downsamples events for faster rendering.
#'   \item Calls \code{\link{compute.gate.density}} for density estimation.
#'   \item Calls \code{\link{build.gate.plot}} to create the base scatter plot.
#'   \item Iterates over non-ignored populations, calling
#'     \code{\link{add.gate.boundary}} for each.
#'   \item Saves the completed figure as a PNG file in \code{figure.dir}.
#' }
#'
#' @param fcs.gates Named list.  Gate results; boundary, label, and fraction
#'   are read from \code{fcs.gates[[gate.name]][[pop.name]]}.
#' @param gate.data Numeric matrix (n × 2).  Parent-population events used
#'   for density estimation.
#' @param gate.name Character.  Key into \code{fcs.gates}.
#' @param gate.number Integer.  Used to construct the output filename.
#' @param gate.marker.x Character.  Antigen/channel name for the x-axis.
#' @param gate.marker.y Character.  Antigen/channel name for the y-axis.
#' @param populations Named list of population metadata (must contain
#'   \code{label} and \code{label.pos} for each population).
#' @param p Named list.  Resolved gate parameters.
#' @param figure.dir Character.  Output directory for PNG files.
#' @param ap Autogating parameter object; uses \code{ap$seed},
#'   \code{ap$gate.number.width}, \code{ap$figure.width},
#'   \code{ap$figure.height}, and \code{ap$ignore}.
#' @param max.points Integer.  Maximum number of events to render; events are
#'   downsampled by random sampling if \code{nrow(gate.data) > max.points}.
#'   Default \code{1e5}.
#'
#' @return \code{invisible(NULL)}.  Called for its side-effect of writing a
#'   PNG file.
#'
#' @importFrom ggplot2 ggsave
#' @importFrom ragg agg_png
#'
#' @seealso \code{\link{compute.gate.density}}, \code{\link{build.gate.plot}},
#'   \code{\link{add.gate.boundary}}, \code{\link{do.gate.generic}}
#'
#' @keywords internal
gate.plot <- function( fcs.gates, gate.data, gate.name, gate.number,
                       gate.marker.x, gate.marker.y, populations, p, figure.dir,
                       ap, max.points = 1e5 ) {

  # downsample (faster plotting)
  n.points <- nrow( gate.data )
  if ( n.points > max.points ) {
    # random sampling
    set.seed( ap$seed )
    gate.data <- gate.data[ sample( seq_len( nrow( gate.data ) ), max.points ), ]
  }

  ggp.data <- compute.gate.density( gate.data, p, ap )

  ggp <- build.gate.plot( ggp.data, gate.marker.x, gate.marker.y, p, ap )

  # Reconstruct limits here so add.gate.boundary can use them.
  # Scatter channels use scatter limits; fluorescence channels use expr limits.
  tx <- ap$transform[[ gate.marker.x ]]
  ty <- ap$transform[[ gate.marker.y ]]

  marker.limit <- function( marker ) {
    if ( marker %in% ap$singlet.scatter.parameter ||
         marker %in% ap$leukocyte.scatter.parameter )
      c( ap$figure.limits$scatter.min, ap$figure.limits$scatter.max )
    else
      c( ap$figure.limits$expr.min, ap$figure.limits$expr.max )
  }

  lim.x <- tx( marker.limit( gate.marker.x ) )
  lim.y <- ty( marker.limit( gate.marker.y ) )

  for ( pop.name in names( populations ) ) {

    message( sprintf( "  pop: %s  boundary: %s  fraction: %s",
                      pop.name,
                      if ( is.null( fcs.gates[[ gate.name ]][[ pop.name ]]$boundary ) )
                        "NULL" else paste( dim( fcs.gates[[ gate.name ]][[ pop.name ]]$boundary ), collapse = "x" ),
                      if ( is.null( fcs.gates[[ gate.name ]][[ pop.name ]]$fraction ) )
                        "NULL" else round( fcs.gates[[ gate.name ]][[ pop.name ]]$fraction, 3 )
    ) )

    if ( pop.name == ap$ignore )
      next

    pop.info  <- populations[[ pop.name ]]
    label.pos <- pop.info$label.pos

    if ( identical( label.pos, ap$ignore ) )
      next

    ggp <- add.gate.boundary(
      ggp,
      boundary  = fcs.gates[[ gate.name ]][[ pop.name ]]$boundary,
      label     = pop.info$label,
      label.pos = label.pos,
      fraction  = fcs.gates[[ gate.name ]][[ pop.name ]]$fraction,
      lim.x     = lim.x,
      lim.y     = lim.y,
      p         = p,
      ap        = ap,
      gate.marker.x = gate.marker.x,
      gate.marker.y = gate.marker.y
    )
  }

  ggsave(
    file.path( figure.dir,
               sprintf( "%0*d - %s.png",
                        ap$gate.number.width, gate.number, gate.name ) ),
    ggp,
    device = ragg::agg_png,
    width  = ap$figure.width,
    height = ap$figure.height
  )

  return(invisible( NULL ))
}
