#' Overlay a gate boundary polygon and fraction label on a ggplot
#'
#' @title Add Gate Boundary to Plot
#'
#' @description
#' Adds a closed boundary polygon path and an annotated text label (showing
#' the population name and percentage) to an existing \pkg{ggplot2} object.
#' Boundary vertices are clamped to the axis limits before plotting so that
#' out-of-range vertices appear at the plot edge rather than being dropped.
#' Label position is controlled by a compass code (1–8) via
#' \code{\link{get.pop.label.pos}}.
#'
#' @param ggp \code{\link[ggplot2]{ggplot}} object, typically produced by
#'   \code{\link{build.gate.plot}}.
#' @param boundary Numeric matrix (n × 2).  Boundary polygon vertices
#'   (columns = x, y in transformed space).
#' @param label Character scalar.  Population display label (e.g.
#'   \code{"CD4+ T Cells"}).
#' @param label.pos Integer in 1–8, or \code{ap$ignore} to suppress the
#'   label.  Position codes: 1 = bottom-left, 2 = bottom-right,
#'   3 = top-left, 4 = top-right, 5 = left-bottom, 6 = right-bottom,
#'   7 = left-top, 8 = right-top.
#' @param fraction Numeric scalar in \eqn{[0, 1]}.  Population fraction;
#'   displayed as a percentage alongside the label.
#' @param lim.x Length-2 numeric vector.  x-axis limits in transformed space.
#' @param lim.y Length-2 numeric vector.  y-axis limits in transformed space.
#' @param p Named list.  Resolved gate parameters; uses \code{label.size}.
#' @param ap Autogating parameter object; uses \code{ap$figure.label.alpha}
#'   and \code{ap$ignore}.
#' @param gate.marker.x Character.  x-axis channel name (currently reserved
#'   for future transformation use).
#' @param gate.marker.y Character.  y-axis channel name (currently reserved
#'   for future transformation use).
#'
#' @return Updated \code{\link[ggplot2]{ggplot}} object with boundary path
#'   and label annotation added.
#'
#' @importFrom ggplot2 geom_path annotate aes
#'
#' @seealso \code{\link{build.gate.plot}}, \code{\link{gate.plot}},
#'   \code{\link{get.pop.label.pos}}
#'
#' @keywords internal
add.gate.boundary <- function( ggp, boundary, label, label.pos, fraction,
                               lim.x, lim.y, p, ap,
                               gate.marker.x, gate.marker.y ) {

  #tx <- ap$transform[[ gate.marker.x ]]
  #ty <- ap$transform[[ gate.marker.y ]]

  #boundary.t <- boundary
  #boundary.t[ , 1L ] <- tx( boundary[ , 1L ] )
  #boundary.t[ , 2L ] <- ty( boundary[ , 2L ] )

  # Clamp to axis limits so out-of-range vertices plot at the edge
  # rather than being silently dropped by ggplot
  boundary[ , 1L ] <- pmax( lim.x[ 1L ], pmin( lim.x[ 2L ], boundary[ , 1L ] ) )
  boundary[ , 2L ] <- pmax( lim.y[ 1L ], pmin( lim.y[ 2L ], boundary[ , 2L ] ) )

  bnd.df <- data.frame(
    x = c( boundary[ , 1L ], boundary[ 1L, 1L ] ),
    y = c( boundary[ , 2L ], boundary[ 1L, 2L ] )
  )

  fig.label <- sprintf( "%s\n%.1f%%", label, 100 * fraction )

  fig.label.pos <- get.pop.label.pos(
    label.pos, boundary, lim.x, lim.y, fig.label, ap )

  ggp +
    geom_path( aes( x, y, color = NULL ), data = bnd.df ) +
    annotate( "label",
              label      = fig.label,
              size       = p$label.size,
              x          = fig.label.pos[ 1L ],
              y          = fig.label.pos[ 2L ],
              color      = "white",
              linewidth = NA,
              alpha      = ap$figure.label.alpha,
              na.rm = TRUE
    ) +
    annotate( "text",
              label = fig.label,
              size  = p$label.size,
              x     = fig.label.pos[ 1L ],
              y     = fig.label.pos[ 2L ],
              na.rm = TRUE
    )
}
