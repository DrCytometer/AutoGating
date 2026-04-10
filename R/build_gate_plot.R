#' Construct the base ggplot2 object for a gate scatter plot
#'
#' @title Build Base Gate Scatter Plot
#'
#' @description
#' Creates a density-coloured scatter plot of parent-population events using
#' \pkg{ggplot2} and \pkg{scattermore}.  Axis scales, breaks, labels, and
#' limits are set according to whether each channel is a scatter channel or a
#' fluorescence channel, using transformation functions stored in \code{ap}.
#' No gate boundaries are added here; those are overlaid separately by
#' \code{\link{add.gate.boundary}}.
#'
#' @param ggp.data \code{data.frame} with columns \code{x}, \code{y}, and
#'   \code{z} as returned by \code{\link{compute.gate.density}}.
#' @param gate.marker.x Character.  Antigen/channel name for the x-axis;
#'   used to look up the transform and break settings.
#' @param gate.marker.y Character.  Antigen/channel name for the y-axis.
#' @param p Named list.  Resolved gate parameters; uses
#'   \code{break.x}, \code{break.x.factor}, \code{break.y},
#'   \code{break.y.factor}, \code{break.scatter}, \code{break.scatter.factor},
#'   and \code{label.size}.
#' @param ap Autogating parameter object; uses
#'   \code{ap$transform}, \code{ap$singlet.scatter.parameter},
#'   \code{ap$leukocyte.scatter.parameter}, \code{ap$figure.limits},
#'   \code{ap$figure.point.size}, \code{ap$figure.axis.label.size},
#'   \code{ap$figure.axis.title.size}, and \code{ap$figure.margin}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object with a density-coloured
#'   scatter layer and formatted axes, ready for boundary overlays.
#'
#' @importFrom ggplot2 ggplot aes geom_path scale_x_continuous
#'   scale_y_continuous scale_color_gradientn annotate theme_bw theme
#'   element_text element_blank margin
#' @importFrom scattermore geom_scattermore
#'
#' @seealso \code{\link{compute.gate.density}},
#'   \code{\link{add.gate.boundary}}, \code{\link{gate.plot}}
#'
#' @keywords internal
build.gate.plot <- function( ggp.data, gate.marker.x, gate.marker.y, p, ap ) {

  density.palette <- if ( nrow( ggp.data ) > 0L )
    get.density.palette( ggp.data$z, ap )
  else
    "black"

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

  marker.breaks <- function( marker, axis ) {
    if ( marker %in% ap$singlet.scatter.parameter ||
         marker %in% ap$leukocyte.scatter.parameter )
      list( breaks = p$break.scatter, factor = p$break.scatter.factor )
    else if ( axis == "x" )
      list( breaks = p$break.x, factor = p$break.x.factor )
    else
      list( breaks = p$break.y, factor = p$break.y.factor )
  }

  bx <- marker.breaks( gate.marker.x, "x" )
  by <- marker.breaks( gate.marker.y, "y" )

  brk.x <- tx( bx$breaks * bx$factor )
  brk.y <- ty( by$breaks * by$factor )

  suffix <- function( factor )
    switch( as.character( log10( factor ) / 3 + 1 ), "1" = "", "2" = "K", "3" = "M", "" )

  lbl.x <- paste0( bx$breaks, suffix( bx$factor ) )
  lbl.y <- paste0( by$breaks, suffix( by$factor ) )

  #x.label <- ap$marker.label[[ gate.marker.x ]]
  #y.label <- ap$marker.label[[ gate.marker.y ]]

  ggplot( ggp.data, aes( x, y, color = z ) ) +
    scale_x_continuous(
      name   = gate.marker.x,
      breaks = brk.x,
      labels = lbl.x,
      limits = lim.x
    ) +
    scale_y_continuous(
      name   = gate.marker.y,
      breaks = brk.y,
      labels = lbl.y,
      limits = lim.y
    ) +
    geom_scattermore( size = ap$figure.point.size, na.rm = TRUE ) +
    scale_color_gradientn(
      "",
      labels = NULL,
      colors = density.palette
      #guide  = guide_colorbar( barheight = ap$density.barheight )
    ) +
    theme_bw() +
    theme(
      aspect.ratio     = 1,
      plot.margin      = margin(
        ap$figure.margin[ "up" ],
        ap$figure.margin[ "right" ],
        ap$figure.margin[ "down" ],
        ap$figure.margin[ "left" ]
      ),
      axis.text        = element_text( size = ap$figure.axis.label.size ),
      axis.title       = element_text( size = ap$figure.axis.title.size ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}
