# analyze.r
#
# Pipeline entry point for FCS autogating.
#
# Usage: edit the get.autogating.param() call below to match your experiment,
# then source this file.  Everything else is driven by the parameter object
# and the gate definition files.


library( fields )
library( flowClust )
library( flowStats )
library( flowWorkspace )
library( flowCore )
library( ggplot2 )
library( MASS )
library( modeest )
library( moments )
library( parallel )
library( RColorBrewer )
library( sp )
library( scattermore )
library( mclust )
library( ragg )


# ---- source pipeline functions ------------------------------------------------

fcs.src.dir <- "C:/Users/ob240/OneDrive - University of Cambridge/Bioinformatics/AutoGating/AutoGating/R"
#devtools::load_all(fcs.src.dir)
invisible(lapply(list.files(fcs.src.dir, full.names = TRUE), source))

working.dir <- "C:/Users/ob240/OneDrive - University of Cambridge/Bioinformatics/AutoGating/20260408_test"
setwd(working.dir)

# ---- build parameter object ---------------------------------------------------

ap <- get.autogating.param(
    cytometer              = "symphony",
    data.dir               = "C:/Users/ob240/OneDrive - University of Cambridge/Bioinformatics/AutoGating/10.01_fcs_data_devel_st1",
    panel.filename         = "fcs_panel_st1.csv",
    experiment.filename    = "fcs_experiment_st1.csv",
    compensation.filename  = "fcs_compensation_st1.csv",
    transformation.filename = "fcs_transformation_st1.csv",          # NULL = use cytometer defaults
    population.gates.filename  = "fcs_population_gates_st1.r",
    activation.gates.filename  = "fcs_activation_gates_st1.r",
    figure.dir.basename    = "figure",
    statistics.dir         = "statistics",
    population.statistics.filename = "population_statistics.csv",
    activation.statistics.filename = "activation_statistics.csv",
    seed                   = 20180101L,
    worker.n               = 10L
)


# ---- set random seed ----------------------------------------------------------

set.seed( ap$seed )


# ==============================================================================
# 1. Read experiment manifest
# ==============================================================================

fcs.experiment <- read.csv( ap$experiment.filename, stringsAsFactors = FALSE )

stopifnot( anyDuplicated( fcs.experiment$filename ) == 0 )
stopifnot( anyDuplicated( fcs.experiment$sample )   == 0 )

# gate.calculation defaults to TRUE for all samples unless explicitly set FALSE
fcs.experiment$gate.calculation[ is.na( fcs.experiment$gate.calculation ) |
    fcs.experiment$gate.calculation == "" ] <- FALSE
fcs.experiment$gate.calculation <- as.logical( fcs.experiment$gate.calculation )

fcs.sample              <- fcs.experiment$sample
fcs.sample.n            <- length( fcs.sample )
fcs.sample.filename     <- fcs.experiment$filename
fcs.sample.gate.calc    <- fcs.sample[ fcs.experiment$gate.calculation ]

message( sprintf( "%d samples total, %d used for gate calculation.",
    fcs.sample.n, length( fcs.sample.gate.calc ) ) )


# ==============================================================================
# 2. Read panel and define markers
# ==============================================================================

fcs.panel <- read.csv( ap$panel.filename, stringsAsFactors = FALSE )

# Drop the AutoF channel — it is not used in any gate
fcs.panel <- fcs.panel[ fcs.panel$antigen != "AutoF", ]

fcs.panel$type <- factor( fcs.panel$type,
    levels = c( "scatter", "time", "lineage", "activation" ) )

fcs.panel.scatter.idx    <- which( fcs.panel$type == "scatter" )
fcs.panel.time.idx       <- which( fcs.panel$type == "time" )
fcs.panel.lineage.idx    <- which( fcs.panel$type == "lineage" )
fcs.panel.activation.idx <- which( fcs.panel$type == "activation" )

stopifnot( sort( c( fcs.panel.scatter.idx, fcs.panel.time.idx,
    fcs.panel.lineage.idx, fcs.panel.activation.idx ) ) ==
    seq_len( nrow( fcs.panel ) ) )

fcs.marker.scatter    <- fcs.panel$antigen[ fcs.panel.scatter.idx ]
fcs.marker.time       <- fcs.panel$antigen[ fcs.panel.time.idx ]
fcs.marker.lineage    <- fcs.panel$antigen[ fcs.panel.lineage.idx ]
fcs.marker.activation <- fcs.panel$antigen[ fcs.panel.activation.idx ]

fcs.marker <- c( fcs.marker.scatter, fcs.marker.time,
                 fcs.marker.lineage, fcs.marker.activation )

fcs.marker.label <- c(
    fcs.marker.scatter,
    fcs.marker.time,
    paste0( fcs.panel$dye[ fcs.panel.lineage.idx ],    " - ", fcs.marker.lineage ),
    paste0( fcs.panel$dye[ fcs.panel.activation.idx ], " - ", fcs.marker.activation )
)

fcs.marker.label.figure        <- fcs.marker.label
names( fcs.marker.label.figure ) <- fcs.marker

fcs.marker.label.table         <- fcs.marker
names( fcs.marker.label.table )  <- fcs.marker


# ==============================================================================
# 3. Read compensation and transformation
# ==============================================================================

fcs.compensation <- read.compensation.csv( ap$compensation.filename, fcs.panel,
    fcs.panel.lineage.idx, fcs.panel.activation.idx )

fcs.transform.both <- read.transformation( ap$transformation.filename, fcs.panel,
    fcs.panel.scatter.idx, fcs.panel.lineage.idx, fcs.panel.activation.idx, ap )

fcs.transform     <- fcs.transform.both[[ 1 ]]
fcs.transform.inv <- fcs.transform.both[[ 2 ]]

ap$transform <- fcs.transform
ap$transform.inv <- fcs.transform.inv


# ==============================================================================
# 4. Read FCS files; apply compensation and transformation
# ==============================================================================

fcs.flow.set <- suppressWarnings( read.flowSet( files = fcs.sample.filename,
    path = ap$data.dir, transformation = NULL, truncate_max_range = TRUE ) )

stopifnot( fcs.panel$dye %in% colnames( fcs.flow.set ) )

fcs.flow.set <- fcs.flow.set[ , fcs.panel$dye ]
colnames( fcs.flow.set ) <- fcs.panel$antigen

fcs.flow.set <- compensate( fcs.flow.set, compensation( fcs.compensation ) )
fcs.flow.set <- transform( fcs.flow.set,
    transformList( names( fcs.transform ), fcs.transform ) )


# ==============================================================================
# 5. Extract and label expression data
# ==============================================================================

# Determine zero-padded event ID width from the largest sample
fcs.sample.event.n.max <- max( sapply( seq_len( fcs.sample.n ),
    function( i ) nrow( exprs( fcs.flow.set[[ i ]] ) ) ) )
fcs.event.number.width <- floor( log10( fcs.sample.event.n.max ) ) + 1L
fcs.event.regexp       <- sprintf( "\\.[0-9]{%d}$", fcs.event.number.width )

fcs.expr.data <- do.call( rbind, lapply( seq_len( fcs.sample.n ), function( i ) {
    ed   <- exprs( fcs.flow.set[[ i ]] )
    samp <- fcs.sample[ i ]

    rownames( ed ) <- sprintf( "%s.%0*d", samp, fcs.event.number.width,
        seq_len( nrow( ed ) ) )

    stopifnot( sort( colnames( ed ) ) == sort( fcs.panel$antigen ) )

    ed <- ed[ , fcs.marker ]
    attr( ed, "ranges" )            <- NULL
    attr( colnames( ed ), "names" ) <- NULL
    ed
} ) )

fcs.expr.data[ , "Time" ] <- fcs.expr.data[ , "Time" ] / 100

fcs.event        <- rownames( fcs.expr.data )
fcs.event.n      <- length( fcs.event )
fcs.event.sample <- factor( sub( fcs.event.regexp, "", fcs.event ),
    levels = fcs.sample )

message( sprintf( "%d total events across %d samples.", fcs.event.n, fcs.sample.n ) )


# ==============================================================================
# 6. Load gate definitions
# ==============================================================================

fcs.population.gates <- load.gate.definitions( ap$population.gates.filename )
fcs.activation.gates <- load.gate.definitions( ap$activation.gates.filename )


# ==============================================================================
# 7. Create output directories
# ==============================================================================

fcs.figure.popgate.dir <- c(
    all    = paste0( ap$figure.dir.basename, "_popgate_all" ),
    setNames(
        sprintf( "%s_popgate_%s", ap$figure.dir.basename, fcs.sample ),
        fcs.sample
    )
)

fcs.figure.actgate.dir <- c(
    all    = paste0( ap$figure.dir.basename, "_actgate_all" ),
    setNames(
        sprintf( "%s_actgate_%s", ap$figure.dir.basename, fcs.sample ),
        fcs.sample
    )
)

for ( d in c( fcs.figure.popgate.dir, fcs.figure.actgate.dir,
              ap$statistics.dir ) )
    if ( !file.exists( d ) )
        dir.create( d, recursive = TRUE )


# ==============================================================================
# 8. Calculate gate boundaries on pooled gate-calculation samples
# ==============================================================================

fcs.gates.all      <- list()
fcs.gates.all.data <- fcs.expr.data[ fcs.event.sample %in% fcs.sample.gate.calc, ]

fcs.gates.all <- do.population.gates(
    fcs.gates.all, fcs.gates.all.data,
    fcs.population.gates, ap,
    calculate.gates = TRUE,
    figure.dir      = fcs.figure.popgate.dir[ "all" ]
)

fcs.gates.all <- do.activation.gates(
    fcs.gates.all, fcs.gates.all.data,
    fcs.activation.gates, ap,
    calculate.gates = TRUE,
    figure.dir      = fcs.figure.actgate.dir[ "all" ]
)


# ==============================================================================
# 9. Apply gates to each individual sample (parallel)
# ==============================================================================

# check fcs.figure.popgate.dir
# turn off plotting by argument
# change to ragg::png
fcs.gates.sample <- lapply( fcs.sample, function( samp ) {

    samp.data  <- fcs.expr.data[ fcs.event.sample == samp, ]
    samp.gates <- fcs.gates.all   # carry boundaries calculated on pooled data

    samp.gates <- do.population.gates(
        samp.gates, samp.data,
        fcs.population.gates, ap,
        calculate.gates = FALSE,
        figure.dir      = fcs.figure.popgate.dir[ samp ]
    )

    samp.gates <- do.activation.gates(
        samp.gates, samp.data,
        fcs.activation.gates, ap,
        calculate.gates = FALSE,
        figure.dir      = fcs.figure.actgate.dir[ samp ]
    )

    samp.gates

} )

names( fcs.gates.sample ) <- fcs.sample


# ==============================================================================
# 10. Calculate and write statistics
# ==============================================================================

population.stats <- calculate.population.statistics(
    fcs.gates.sample, fcs.population.gates, ap )

activation.stats <- calculate.activation.statistics(
    fcs.gates.sample, fcs.activation.gates, fcs.expr.data, ap )

write.table( population.stats,
    file  = file.path( ap$statistics.dir, ap$population.statistics.filename ),
    sep   = ",", quote = FALSE, row.names = FALSE )

write.table( activation.stats,
    file  = file.path( ap$statistics.dir, ap$activation.statistics.filename ),
    sep   = ",", quote = FALSE, row.names = FALSE )

message( "Done." )

## differential abundance analysis
library( limma )
library( pheatmap )

diff.res <- calculate.differential.abundance(
  population.stats,
  activation.stats,
  ap,
  experiment.filename = "fcs_experiment_st1.csv"
)








