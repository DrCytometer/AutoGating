
do.activation.gates <- function( fcs.expr.data, fcs.population.gates, 
    fcs.figure.dir, fcs.activation.thresholds, fcs.activation.gates, 
    fcs.activation.gates.calculate, fcs.activation.gates.param )
{
    for ( fcs.activation.marker in names( fcs.activation.gates ) )
        for ( fcs.activation.population in 
            names( fcs.activation.gates[[ fcs.activation.marker ]] ) )
            fcs.activation.gates <- do.act.gate( fcs.expr.data, 
                fcs.population.gates, fcs.figure.dir, 
                fcs.activation.thresholds[[ fcs.activation.marker ]]$threshold, 
                fcs.activation.marker, fcs.activation.population, 
                fcs.activation.gates, 
                fcs.activation.gates.calculate[ sprintf( "%s.%s", 
                    fcs.activation.marker, fcs.activation.population ) ], 
                fcs.activation.gates.param[[ fcs.activation.marker ]][[ 
                    fcs.activation.population ]] )
    
    fcs.activation.gates
}


do.activation.thresholds <- function( fcs.expr.data, 
    fcs.population.gates, fcs.figure.dir, fcs.activation.thresholds, 
    fcs.activation.thresholds.param )
{
    for ( activation.marker in names( fcs.activation.thresholds ) )
        if ( fcs.activation.thresholds[[ activation.marker ]] == 0 ) {
            # general threshold function
            fcs.activation.thresholds <- do.act.threshold( 
                fcs.expr.data, fcs.population.gates, fcs.figure.dir, 
                fcs.activation.thresholds, activation.marker, 
                fcs.activation.thresholds.param[[ activation.marker ]] )
        }
        else {
            # special threshold function
            threshold.function <- 
                get( paste0( "do.act.threshold.", activation.marker ) )
            fcs.activation.thresholds <- threshold.function( fcs.expr.data, 
                fcs.population.gates, fcs.figure.dir, fcs.activation.thresholds )
            
        }
    
    fcs.activation.thresholds
}


do.population.gates <- function( fcs.gates, fcs.gates.data, 
    fcs.popul.gates.definition, calculate.gates, figure.popgate.dir )
{
    for ( popul.gate.def in fcs.popul.gates.definition )
    {
        gate.number <- popul.gate.def$gate.number
        gate.name <- popul.gate.def$gate.name
        gate.marker <- c( popul.gate.def$gate.marker.x, 
            popul.gate.def$gate.marker.y )
        
        cat( sprintf( "%0*d - %s\n", fcs.gate.number.width, gate.number, 
            gate.name ) )
        
        fcs.gates <- do.gate.generic( fcs.gates, 
            fcs.gates.data[ , gate.marker ], popul.gate.def, calculate.gates, 
            figure.popgate.dir )
        
        str( fcs.gates[[ gate.name ]] )
        cat( "\n" )
    }
    
    fcs.gates
}


do.activation.gates <- function( fcs.gates, fcs.gates.data, 
    fcs.activ.gates.definition, calculate.gates, figure.actgate.dir )
{
    for ( activ.gate.def in fcs.activ.gates.definition )
    {
        marker.activation <- fcs.marker.activation[ 
            ! is.na( activ.gate.def[ fcs.marker.activation ] ) &
            activ.gate.def[ fcs.marker.activation ] != "" ]
        
        if ( length( marker.activation ) > 0 )
        {
            popul.number <- activ.gate.def$popul.number
            
            parent.gate <- activ.gate.def$parent.gate
            parent.popul <- activ.gate.def$parent.popul
            
            popul.label.base <- activ.gate.def$popul.label
            popul.label.pos <- activ.gate.def$popul.label.pos
            
            gate.algorithm.base <- activ.gate.def$gate.algorithm
            
            stopifnot( xor(
                activ.gate.def$gate.marker.x == fcs.gate.parameter.ignore, 
                activ.gate.def$gate.marker.y == fcs.gate.parameter.ignore 
            ) )
            
            marker.act.axis <- ifelse( 
                activ.gate.def$gate.marker.x == fcs.gate.parameter.ignore, 
                1, 2 )
            
            activ.gate.def$gate.number <- popul.number
            
            activ.gate.def$popul.label.pos <- paste( fcs.gate.parameter.ignore, 
                popul.label.pos, sep = fcs.gate.parameter.delimiter )
            
            for ( marker.act in marker.activation )
            {
                gate.name <- sprintf( "%s.%s.%s", parent.gate, parent.popul, 
                    marker.act )
                
                popul.label <- sprintf( "%s+ %s", marker.act, 
                    popul.label.base )
                
                activ.gate.def$gate.name <- gate.name
                
                activ.gate.def$popul.name <- paste( fcs.gate.parameter.ignore, 
                    fcs.activation.population.label, 
                    sep = fcs.gate.parameter.delimiter )
                
                activ.gate.def$popul.label <- paste( fcs.gate.parameter.ignore, 
                    popul.label, 
                    sep = fcs.gate.parameter.delimiter )
                
                if ( marker.act.axis == 1 )
                    activ.gate.def$gate.marker.x <- marker.act
                else
                    activ.gate.def$gate.marker.y <- marker.act
                
                activ.gate.def$gate.algorithm <- ifelse( 
                    activ.gate.def[[ marker.act ]] == 
                        fcs.activation.marker.selected, 
                    gate.algorithm.base, activ.gate.def[[ marker.act ]] )
                
                gate.marker <- c( activ.gate.def$gate.marker.x, 
                    activ.gate.def$gate.marker.y )
                
                cat( sprintf( "%0*d - %s\n", fcs.gate.number.width, 
                    popul.number, gate.name ) )
                
                fcs.gates <- do.gate.generic( fcs.gates, 
                    fcs.gates.data[ , gate.marker ], activ.gate.def, 
                    calculate.gates, figure.actgate.dir )
                
                str( fcs.gates[[ gate.name ]] )
                cat( "\n" )
            }
        }
    }
    
    fcs.gates
}


do.gate.generic <- function( fcs.gates, fcs.gate.data, fcs.gate.def, 
    calculate.gate, figure.dir )
{
    # get gate definition
    
    gate.number <- fcs.gate.def$gate.number
    
    gate.name <- fcs.gate.def$gate.name
    popul.name <- strsplit( fcs.gate.def$popul.name, 
        fcs.gate.parameter.delimiter )[[ 1 ]]
    
    gate.marker.x <- fcs.gate.def$gate.marker.x
    gate.marker.y <- fcs.gate.def$gate.marker.y
    
    parent.gate <- strsplit( fcs.gate.def$parent.gate, 
        fcs.gate.parameter.delimiter )[[ 1 ]]
    parent.popul <- strsplit( fcs.gate.def$parent.popul, 
        fcs.gate.parameter.delimiter )[[ 1 ]]
    
    popul.label <- strsplit( fcs.gate.def$popul.label, 
        fcs.gate.parameter.delimiter )[[ 1 ]]
    popul.label.pos <- strsplit( fcs.gate.def$popul.label.pos, 
        fcs.gate.parameter.delimiter )[[ 1 ]]
    
    gate.type <- fcs.gate.def$gate.type
    
    gate.algorithm <- strsplit( fcs.gate.def$gate.algorithm, 
        fcs.gate.parameter.delimiter )[[ 1 ]]
    
    gate.param.raw <- strsplit( fcs.gate.def$gate.param, 
        fcs.gate.parameter.delimiter )[[ 1 ]]
    gate.param <- vector( "list", length( gate.param.raw ) )
    
    if ( length( gate.param.raw ) > 0 )
        for ( gpr.idx in 1 : length( gate.param.raw ) )
        {
            gpar.raw <- gate.param.raw[ gpr.idx ]
            
            gpar.pos <- regexpr( fcs.gate.parameter.operator, gpar.raw )[[ 1 ]]
            gpar.name <- substr( gpar.raw, 1, gpar.pos - 1 )
            gpar.value <- substr( gpar.raw, gpar.pos + 1, nchar( gpar.raw ) )
            
            names( gate.param )[[ gpr.idx ]] <- gpar.name
            gate.param[[ gpr.idx ]] <- gpar.value
        }
    
    density.quantile.x <- fcs.gate.def$density.quantile.x
    density.quantile.y <- fcs.gate.def$density.quantile.y
    
    density.trimneg.x <- fcs.gate.def$density.trimneg.x
    if ( is.na( density.trimneg.x ) )
        density.trimneg.x <- fcs.default.density.trimneg.x

    density.trimneg.y <- fcs.gate.def$density.trimneg.y
    if ( is.na( density.trimneg.y ) )
        density.trimneg.y <- fcs.default.density.trimneg.y
    
    limit.xmin <- fcs.gate.def$limit.xmin
    limit.xmax <- fcs.gate.def$limit.xmax
    
    limit.ymin <- fcs.gate.def$limit.ymin
    limit.ymax <- fcs.gate.def$limit.ymax
    
    break.x <- as.numeric( strsplit( fcs.gate.def$break.x, 
        fcs.gate.parameter.delimiter )[[ 1 ]] )
    
    break.x.factor <- fcs.gate.def$break.x.factor
    
    break.y <- as.numeric( strsplit( fcs.gate.def$break.y, 
        fcs.gate.parameter.delimiter )[[ 1 ]] )
    
    break.y.factor <- fcs.gate.def$break.y.factor
    
    label.size <- fcs.gate.def$label.size
    if ( is.null( label.size ) )
        label.size <- fcs.figure.label.size
    
    # get data from parent population
    
    if ( gate.number > 0 )
    {
        stopifnot( length( parent.gate ) == length( parent.popul ) )
        
        gate.parent.idx <- sort( unlist( lapply( 1 : length( parent.gate ), 
            function ( pg.idx ) fcs.gates[[ parent.gate[ pg.idx ] ]][[ 
                parent.popul[ pg.idx ] ]]$idx ) ) )
        
        stopifnot( anyDuplicated( gate.parent.idx ) == 0 )
    }
    else
        gate.parent.idx <- 1 : nrow( fcs.gate.data )
    
    stopifnot( colnames( fcs.gate.data ) == c( gate.marker.x, gate.marker.y ) )
    
    gate.data <- fcs.gate.data[ gate.parent.idx, , drop = FALSE ]
    
    # calculate gate boundary if needed
    
    if ( calculate.gate )
    {
        calculate.gate.function <- get( paste0( "calculate.gate.", gate.type ) )
        
        gate.boundary <- calculate.gate.function( gate.data, popul.name, 
            gate.algorithm, gate.param, gate.number, gate.name, 
            fcs.gates, figure.dir )
        
        fcs.gates[[ gate.name ]] <- list()
        
        for ( pop.name in names( gate.boundary ) )
            fcs.gates[[ gate.name ]][[ pop.name ]] <- list( 
                parent.gate = parent.gate, parent.popul = parent.popul, 
                boundary = gate.boundary[[ pop.name ]] )
    }
    
    # obtain populations defined by gate boundaries
    
    for ( pop.name in names( fcs.gates[[ gate.name ]] ) )
    {
        pop.label <- popul.label[ pop.name == popul.name ]
        
        fcs.gates[[ gate.name ]][[ pop.name ]]$label <- pop.label
        
        pop.boundary <- fcs.gates[[ gate.name ]][[ pop.name ]]$boundary
        
        pop.pip <- point.in.polygon( gate.data[ , 1 ], gate.data[ , 2 ], 
            pop.boundary[ , 1 ], pop.boundary[ , 2 ] )
        
        fcs.gates[[ gate.name ]][[ pop.name ]]$idx <- 
            gate.parent.idx[ pop.pip == 1 ]
        
        fcs.gates[[ gate.name ]][[ pop.name ]]$fraction <- 
            mean( pop.pip == 1 )
    }
    
    # plot gate
    
    if ( gate.number > 0 )
    {
        # trim data for density calculation
        
        if ( nrow( gate.data ) > 0 )
        {
            gate.data.xmin <- min( gate.data[ , 1 ] )
            gate.data.dens.xmin <- quantile( gate.data[ , 1 ], 
                density.quantile.x )
            if ( gate.data.dens.xmin <= gate.data.xmin )
                gate.data.dens.xmin <- ( 1 + .Machine$double.eps^0.5 ) * 
                    gate.data.xmin
            if ( gate.data.dens.xmin < 0 && density.trimneg.x )
                gate.data.dens.xmin <- 0
            
            gate.data.xmax <- max( gate.data[ , 1 ] )
            gate.data.dens.xmax <- quantile( gate.data[ , 1 ], 
                1 - density.quantile.x )
            if ( gate.data.dens.xmax >= gate.data.xmax )
                gate.data.dens.xmax <- ( 1 - .Machine$double.eps^0.5 ) * 
                    gate.data.xmax
            
            gate.data.ymin <- min( gate.data[ , 2 ] )
            gate.data.dens.ymin <- quantile( gate.data[ , 2 ], 
                density.quantile.y )
            if ( gate.data.dens.ymin <= gate.data.ymin )
                gate.data.dens.ymin <- ( 1 + .Machine$double.eps^0.5 ) * 
                    gate.data.ymin
            if ( gate.data.dens.ymin < 0 && density.trimneg.y )
                gate.data.dens.ymin <- 0
            
            gate.data.ymax <- max( gate.data[ , 2 ] )
            gate.data.dens.ymax <- quantile( gate.data[ , 2 ], 
                1 - density.quantile.y )
            if ( gate.data.dens.ymax >= gate.data.ymax )
                gate.data.dens.ymax <- ( 1 - .Machine$double.eps^0.5 ) * 
                    gate.data.ymax
            
            gate.data.dens <- gate.data[ 
                gate.data[ , 1 ] >= gate.data.dens.xmin & 
                gate.data[ , 1 ] <= gate.data.dens.xmax & 
                gate.data[ , 2 ] >= gate.data.dens.ymin & 
                gate.data[ , 2 ] <= gate.data.dens.ymax, , drop = FALSE ]
        }
        else
            gate.data.dens <- gate.data
        
        if ( nrow( gate.data.dens ) > 1 )
        {
          bdw.x <- 3 * dpik( gate.data[ , 1 ] )
          bdw.y <- 3 * dpik( gate.data[ , 2 ] )
          
          kde.density <- bkde2D( gate.data.dens, 
                                 bandwidth = c( bdw.x, bdw.y ), 
                         gridsize = c( fcs.density.grid.n, fcs.density.grid.n ) )
          
          names( kde.density ) <- c( "x", "y", "z" )
          
          gate.data.ggp <- data.frame( 
                x = gate.data[ , 1 ], 
                y = gate.data[ , 2 ], 
                z = interp.surface(  kde.density, gate.data ) 
            )
        }
        else 
        {
            gate.data.n <- nrow( gate.data )
            
            gate.data.ggp <- data.frame( 
                x = gate.data[ , 1 ], 
                y = gate.data[ , 2 ], 
                z = rep( 1 / gate.data.n, times = gate.data.n )
            )
        }
        
        if ( length( gate.data.ggp$z ) > 0 )
            density.palette <- get.density.palette( gate.data.ggp$z )
        else
            density.palette <- "black"
        
        transform.x <- fcs.transform[[ gate.marker.x ]]
        transform.y <- fcs.transform[[ gate.marker.y ]]
        
        limit.x <- transform.x( c( limit.xmin, limit.xmax ) )
        limit.y <- transform.y( c( limit.ymin, limit.ymax ) )
        
        ggp <- ggplot( gate.data.ggp, aes( x, y, color = z ) ) + 
            scale_x_continuous( 
                name = fcs.marker.label.figure[ gate.marker.x ], 
                breaks = transform.x( break.x * break.x.factor ), 
                labels = paste0( break.x, 
                    switch( log10( break.x.factor ) / 3 + 1, "", "K", "M" ) ), 
                limits = limit.x ) + 
            scale_y_continuous( 
                name = fcs.marker.label.figure[ gate.marker.y ], 
                breaks = transform.y( break.y * break.y.factor ), 
                labels = paste0( break.y, 
                    switch( log10( break.y.factor ) / 3 + 1, "", "K", "M" ) ), 
                limits = limit.y ) + 
            geom_scattermore( size = fcs.figure.point.size * 3.3 ) + 
            scale_color_gradientn( "", labels = NULL, colors = density.palette, 
                guide = guide_colorbar( barheight = fcs.density.barheight ) ) + 
            theme_bw() + 
            theme( aspect.ratio = 1, 
                plot.margin = margin( fcs.figure.margin.up, 
                    fcs.figure.margin.ri, fcs.figure.margin.do, 
                    fcs.figure.margin.le ), 
                axis.text = element_text( size = fcs.figure.axis.label.size ), 
                axis.title = element_text( size = fcs.figure.axis.title.size ), 
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank() )
        
        # add populations to plot
        
        for ( pn.idx in 1 : length( popul.name ) )
        {
            pop.name <- popul.name[ pn.idx ]
            
            if ( pop.name != fcs.gate.parameter.ignore )
            {
                pop.label <- popul.label[ pn.idx ]
                pop.label.pos <- popul.label.pos[ pn.idx ]
                pop.boundary <- fcs.gates[[ gate.name ]][[ 
                    pop.name ]]$boundary
                pop.fraction <- fcs.gates[[ gate.name ]][[ 
                    pop.name ]]$fraction
                
                pop.boundary.ggp <- data.frame( 
                    x = c( pop.boundary[ , 1 ], pop.boundary[ 1, 1 ] ), 
                    y = c( pop.boundary[ , 2 ], pop.boundary[ 1, 2 ] ) )
                
                fig.label <- sprintf( "%s\n%.1f%%", pop.label, 
                    100 * pop.fraction )
                fig.label.pos <- get.pop.label.pos( pop.label.pos, 
                    pop.boundary, limit.x, limit.y, fig.label )
                
                ggp <- ggp + 
                    geom_path( aes( x, y, color = NULL ), 
                        data = pop.boundary.ggp ) + 
                    annotate( "label", label = fig.label, size = label.size, 
                        x = fig.label.pos[ 1 ], y = fig.label.pos[ 2 ], 
                        color = "white", label.size = NA, 
                        alpha = fcs.figure.label.alpha ) + 
                    annotate( "text", label = fig.label, size = label.size, 
                        x = fig.label.pos[ 1 ], y = fig.label.pos[ 2 ] )
            }
        }
        
        ggsave( file.path( figure.dir, 
                sprintf( "%0*d - %s.jpg", fcs.gate.number.width, gate.number, 
                    gate.name ) ), 
            ggp, width = fcs.figure.width, height = fcs.figure.height )
    }
    
    fcs.gates
}

