# calculate_boundary_singlets.r

calculate.boundary.singlets <- function( marker.data, boundary.param )
{
  gate.singlet.res <- gate.singlet( marker.data,
                                    maxit = 100, gate.grid.n = 20, gate.skip.low = 0.06,
                                    gate.skip.high = 0.10, spread.span = 1.0, spread.factor = 7,
                                    spread.spline.x =
                                      c( 0.00, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 0.90, 1.00 ),
                                    spread.spline.y =
                                      # c( 2.00, 2.00, 2.00, 1.30, 0.70, 0.45, 0.40, 0.35, 0.35 ) )
                                      c( 2.00, 2.10, 1.90, 1.30, 0.90, 0.45, 0.30, 0.25, 0.25 ) )

  gate.boundary <- gate.singlet.res@boundaries
  dimnames( gate.boundary ) <- NULL

  gate.boundary
}

