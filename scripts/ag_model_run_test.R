# ag_model_run_test.R
#
# A trained Honeychrome model applied by AutoGating to the samples exported
# unmixed from Honeychrome, against Honeychrome's own run of the same model
# on the same samples (hc_harness scenario ag_model_run). Both use the
# experiment's transforms (Honeychrome's, written by the scenario), and both
# convert the model's stored boundaries to them.
#
# Compared per sample and population: event count and fraction of the
# parent. Outputs in file.path( hc.dir, "r" ): counts.csv (every population
# in either run), summary.txt.

library( AutoGating )

hc.dir  <- "~/dev/claude-sandboxes/hc_harness/results/ag_model_run/20261117T132714Z"
model   <- "~/Experiments/UCL IMMPROVE/AG_MODEL/IMMPROVE_stacked_panel2.agmodel"
fcs.dir <- "~/Experiments/UCL IMMPROVE/Unmixed/IC controls"
hc.samples <- c( BuffyCoat_Donor_E = "Raw/IC controls/C4 IC_BuffyCoat_Donor_E_Plate 4.fcs",
                 PBMC_Donor_1      = "Raw/IC controls/C5 IC_PBMC_Donor_1_Plate 4.fcs" )

out.dir <- file.path( hc.dir, "r" )
dir.create( out.dir, showWarnings = FALSE, recursive = TRUE )
summary.lines <- character( 0 )
say <- function( ... ) {
    line <- sprintf( ... )
    summary.lines <<- c( summary.lines, line )
    cat( line, "\n" )
}


# ---- Experiment with Honeychrome's transforms ---------------------------------

fcs.files <- list.files( fcs.dir, pattern = "\\.fcs$" )
panel <- create.panel.file( file.path( fcs.dir, fcs.files[ 1 ] ), file.path( out.dir, "panel" ) )
panel <- panel[ panel$dye != "event_id", ]
utils::write.csv( panel, file.path( out.dir, "panel.csv" ), row.names = FALSE )

hc.tf <- jsonlite::read_json( file.path( hc.dir, "transforms.json" ) )
tf.dyes <- panel$dye[ panel$type != "time" ]
tf.tbl <- data.frame( dye = tf.dyes, stringsAsFactors = FALSE )
for ( col in c( "id", "scale.t", "linear.a", "logicle.w", "logicle.m", "logicle.a", "log.m" ) )
    tf.tbl[[ col ]] <- vapply( tf.dyes, function( d ) as.numeric( hc.tf[[ d ]][[ gsub( ".", "_", col, fixed = TRUE ) ]] ),
                               numeric( 1 ) )
utils::write.csv( tf.tbl, file.path( out.dir, "transformation.csv" ), row.names = FALSE )

create.experiment.file( data.frame( Sample = names( hc.samples ) ), fcs.dir,
                        filename = file.path( out.dir, "experiment.csv" ) )

ap <- get.autogating.param(
    cytometer               = "aurora",
    data.dir                = fcs.dir,
    panel.filename          = file.path( out.dir, "panel.csv" ),
    experiment.filename     = file.path( out.dir, "experiment.csv" ),
    transformation.filename = file.path( out.dir, "transformation.csv" ),
    model.filename          = model,
    statistics.dir          = file.path( out.dir, "statistics" ),
    figure.training         = FALSE )
fcs.setup <- setup.autogating( ap )


# ---- Run the model --------------------------------------------------------------

gating <- run.autogating( fcs.setup, read.agmodel( model, fcs.setup$fcs.panel ), ap )
say( "R: %d gates; flagged %d gate results", length( gating$model$gate.definitions ),
     sum( !is.na( gating$status$flag ) ) )


# ---- Compare ------------------------------------------------------------------------

hc.pop <- utils::read.csv( file.path( hc.dir, "populations.csv" ), stringsAsFactors = FALSE,
                           check.names = FALSE )
hc.pop$sample <- names( hc.samples )[ match( hc.pop$sample, hc.samples ) ]
r.pop <- gating$population.table
m <- merge( r.pop[ , c( "sample", "key", "n.events", "n.parent", "fraction.of.parent" ) ],
            hc.pop[ , c( "sample", "key", "n_events", "n_parent", "fraction_of_parent" ) ],
            by = c( "sample", "key" ), all = TRUE )
m$diff.events <- m$n.events - m$n_events
m$diff.fraction <- m$fraction.of.parent - m$fraction_of_parent
utils::write.csv( m, file.path( out.dir, "counts.csv" ), row.names = FALSE )

both <- !is.na( m$n.events ) & !is.na( m$n_events )
say( "Populations: R %d, Honeychrome %d, in both %d, only in R %d, only in Honeychrome %d",
     nrow( r.pop ), nrow( hc.pop ), sum( both ), sum( is.na( m$n_events ) ), sum( is.na( m$n.events ) ) )
say( "Identical counts: %d of %d; differing: %d", sum( m$diff.events[ both ] == 0 ), sum( both ),
     sum( m$diff.events[ both ] != 0 ) )
d <- m[ both & m$diff.events != 0, ]
if ( nrow( d ) ) {
    say( "Largest |count difference| %d events; largest |fraction difference| %.3g",
         max( abs( d$diff.events ) ), max( abs( d$diff.fraction ), na.rm = TRUE ) )
    d <- d[ order( -abs( d$diff.events ) ), ]
    for ( i in seq_len( min( 20L, nrow( d ) ) ) )
        say( "   %s  %s: R %d, Honeychrome %d (parent %d / %d)", d$sample[ i ], d$key[ i ],
             d$n.events[ i ], d$n_events[ i ], d$n.parent[ i ], d$n_parent[ i ] )
}

writeLines( summary.lines, file.path( out.dir, "summary.txt" ) )
