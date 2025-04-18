
gate.singlet <- function( x, area = "FSC-A", height = "FSC-H",
    sidescatter = NULL, filterId = "singlet",
    maxit = 100, gate.grid.n = 20,
    gate.skip.low = 0.015, gate.skip.high = 0.15,
    spread.span = 0.5, spread.factor = 10,
    spread.spline.x = c( 0, 0.2, 0.4, 1 ),
    spread.spline.y = c( 1, 0.5, 1, 0 ), ... )
{
    flowCore:::checkClass(x, "matrix")
    flowCore:::checkClass(area, "character")
    flowCore:::checkClass(height, "character")
    if (!is.null(sidescatter)) {
        flowCore:::checkClass(sidescatter, "character")
    }
    if (length(area) + length(height) != 2) {
        stop("Each of 'area' and 'height' must be 'character' vectors of length 1.")
    }

    # filter out extrema points

    x <- data.frame( x[
        x[ , area ] > gate.skip.low * max( x[ , area ] ) +
            ( 1 - gate.skip.low ) * min( x[ , area ] ) &
        x[ , area ] < gate.skip.high * min( x[ , area ] ) +
            ( 1 - gate.skip.high ) * max( x[ , area ] ) &
        x[ , height ] > gate.skip.low * max( x[ , height ] ) +
            ( 1 - gate.skip.low ) * min( x[ , height ] ) &
        x[ , height ] < gate.skip.high * min( x[ , height ] ) +
            ( 1 - gate.skip.high ) * max( x[ , height ] ),
        c( area, height, sidescatter ) ] )

    channel_names <- c(area, height)
    area <- make.names(area)
    height <- make.names(height)
    rlm_formula <- paste(make.names(height), make.names(area),
        sep = " ~ ")
    if (!is.null(sidescatter)) {
        sidescatter <- make.names(sidescatter)
        ssc_ratio <- paste0("I(", sidescatter, " / ", area, ")")
        rlm_formula <- paste(rlm_formula, sidescatter, ssc_ratio,
            sep = " + ")
    }
    rlm_formula <- as.formula(rlm_formula)
    rlm_fit <- MASS::rlm(rlm_formula, data = x, maxit = maxit, ...)
    if (!rlm_fit$converged) {
        warning("The IRLS algorithm employed in 'rlm' did not converge.")
    }

    x.delta <- ( max( x[[ area ]] ) - min( x[[ area ]] ) ) / ( gate.grid.n - 1 )
    x.grid <- seq( min( x[[ area ]] ), max( x[[ area ]] ), by = x.delta )
    x.point.idx <- sapply( x.grid, function( xg )
        which.min( abs( x[[ area ]] - xg ) ) )
    x.predict <- x[ x.point.idx, ]

    predictions <- predict(rlm_fit, x.predict, interval = "none" )

    x.spread <- vapply( 1 : gate.grid.n, function( idx ) {
        x.range.idx <- which(
            abs( x[[ area ]] - x.grid[ idx ] ) < x.delta * spread.span )

        # Check if x.range.idx is empty
        if (length(x.range.idx) > 20 ) {
          x.pred <- x[ x.range.idx, ]
          pred <- predict( rlm_fit, x.pred, interval = "none" )
          x.range.pos.bol <- x[[ height ]][ x.range.idx ] > pred
          x.range.pos.idx <- x.range.idx[ x.range.pos.bol ]
          pred.pos <- pred[ x.range.pos.bol ]

          spread_value <- quantile( x[[ height ]][ x.range.pos.idx ] - pred.pos,
                                    0.5, na.rm = TRUE )

          return( spread_value )
        } else {
          return( NA_real_ )
        }
    }, numeric(1) )

    x.spread[ is.na( x.spread )] <- min( x.spread[ !is.null(x.spread) ], na.rm = TRUE )
    x.spread[ is.null( x.spread ) ] <- min(x.spread[ !is.null(x.spread) ], na.rm = TRUE )

    x.spread <- x.spread *
        splinefun( spread.spline.x, spread.spline.y, method = "natural" )(
            ( 1 : gate.grid.n - 1 ) / ( gate.grid.n - 1 ) )

    gate_vertices <- rbind(
        t( sapply( 1 : gate.grid.n, function( idx )
            c( x.predict[[ area ]][ idx ],
                predictions[ idx ] - spread.factor * x.spread[ idx ] ) ) ),
        c( x.predict[[ area ]][ gate.grid.n ], predictions[ gate.grid.n ] ),
        t( sapply( gate.grid.n: 1, function( idx )
            c( x.predict[[ area ]][ idx ],
                predictions[ idx ] + spread.factor * x.spread[ idx ] ) ) ) )

    colnames(gate_vertices) <- channel_names
    flowCore::polygonGate(gate_vertices, filterId = filterId)
}

