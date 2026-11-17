# ag_r_match_test.R
#
# AutoGating on the same inputs as Honeychrome's Automated Gating run of the
# hc_harness scenario ag_r_match (UCL IMMPROVE IC controls, IMMPROVE FlowJo
# template), and the comparison of the two:
#
#   0. events     Honeychrome's review events against the same rows of the
#                 exported unmixed FCS files, transformed here with
#                 Honeychrome's transforms (largest difference)
#   1. conversion convert.imported.gates() on Honeychrome's review events
#                 against Honeychrome's conversion evidence: gate type,
#                 algorithm of each axis, share drawn positive, share below
#                 zero and agreement (the best fit to the shapes as drawn)
#   2. engine     Honeychrome's trained model applied here to every event
#                 of both samples, against Honeychrome's population counts
#   3. pipeline   the gates as converted here, trained on Honeychrome's
#                 training events and applied to both samples, against
#                 Honeychrome's population counts
#   4. training   Honeychrome's gates trained here on Honeychrome's training
#                 events and applied to both samples, against Honeychrome's
#                 population counts
#   5. thresholds each 1dsep tail and mixture gate calculated here on the
#                 parent events Honeychrome's model gives it, against
#                 Honeychrome's trained threshold
#
# Inputs: hc.dir, a run folder of scenarios/ag_r_match.py; fcs.dir, the two
# samples exported unmixed from Honeychrome; wspt, the FlowJo template.
# Outputs: CSV tables and summary.txt in file.path( hc.dir, "r" ).

library( AutoGating )

hc.dir  <- "~/dev/claude-sandboxes/hc_harness/results/ag_r_match/20261117T102742Z"
fcs.dir <- "~/Experiments/UCL IMMPROVE/Unmixed/IC controls"
wspt    <- "~/dev/claude-sandboxes/FlowJo_workspace/20251212_IMMPROVE panel template_new_panel_v2.wspt"

out.dir <- file.path( hc.dir, "r" )
dir.create( out.dir, showWarnings = FALSE, recursive = TRUE )
summary.lines <- character( 0 )
say <- function( ... ) {
    line <- sprintf( ... )
    summary.lines <<- c( summary.lines, line )
    cat( line, "\n" )
}


# ---- Experiment with Honeychrome's transforms --------------------------------

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

hc.samples <- c( BuffyCoat_Donor_E = "Raw/IC controls/C4 IC_BuffyCoat_Donor_E_Plate 4.fcs",
                 PBMC_Donor_1      = "Raw/IC controls/C5 IC_PBMC_Donor_1_Plate 4.fcs" )
create.experiment.file( data.frame( Sample = names( hc.samples ) ), fcs.dir,
                        filename = file.path( out.dir, "experiment.csv" ) )

ap <- get.autogating.param(
    cytometer               = "aurora",
    data.dir                = fcs.dir,
    panel.filename          = file.path( out.dir, "panel.csv" ),
    experiment.filename     = file.path( out.dir, "experiment.csv" ),
    transformation.filename = file.path( out.dir, "transformation.csv" ),
    figure.dir              = file.path( out.dir, "figures" ),
    statistics.dir          = file.path( out.dir, "statistics" ),
    figure.training         = FALSE )
fcs.setup <- setup.autogating( ap )
ap <- fcs.setup$ap
transforms <- ap$transform.params


# ---- 0. Honeychrome's events ---------------------------------------------------

events <- jsonlite::read_json( file.path( hc.dir, "events.json" ), simplifyVector = TRUE )
read.events <- function( file, shape ) {
    m <- matrix( readBin( file, "numeric", n = prod( shape ), size = 4L ),
                 nrow = shape[ 1 ], byrow = TRUE )
    colnames( m ) <- events$channels
    m <- m[ , panel$dye, drop = FALSE ]
    colnames( m ) <- panel$antigen
    m[ , colnames( fcs.setup$fcs.expr.data ), drop = FALSE ]
}
review <- read.events( file.path( hc.dir, "review_events.f32" ), events$review$shape )
training <- read.events( file.path( hc.dir, "training_events.f32" ), events$training_shape )

review.sample <- names( hc.samples )[ hc.samples == events$review$sample ]
review.rows <- which( fcs.setup$fcs.event.sample == review.sample )[ events$review$event_index + 1L ]
fl <- setdiff( colnames( review ), c( fcs.setup$fcs.marker.scatter, "Time" ) )
sc <- fcs.setup$fcs.marker.scatter
say( "0. events: %d review events of %s; largest difference from the exported FCS: fluorescence %.2g, scatter %.2g (display units)",
     nrow( review ), review.sample,
     max( abs( review[ , fl ] - fcs.setup$fcs.expr.data[ review.rows, fl ] ) ),
     max( abs( review[ , sc ] - fcs.setup$fcs.expr.data[ review.rows, sc ] ) ) )


# ---- 1. Conversion on the review events ---------------------------------------------

imported <- import.flowjo.workspace( wspt, fcs.setup, ap )
converted <- convert.imported.gates( imported$gate.definitions, review, transforms = transforms,
                                     panel = fcs.setup$fcs.panel )
r.prop <- converted$proposals

hc.ev <- jsonlite::read_json( file.path( hc.dir, "evidence.json" ) )
hc.model <- jsonlite::read_json( file.path( hc.dir, "model.agmodel" ) )
hc.defs <- hc.model$gate_definitions
names( hc.defs ) <- vapply( hc.defs, `[[`, character( 1 ), "gate_name" )
num <- function( v ) if ( is.null( v ) ) NA_real_ else as.numeric( v )
hc.axis.algorithm <- function( g, axis ) {
    if ( is.null( g ) ) return( NA_character_ )
    gp <- g$gate_param
    if ( identical( g$gate_type, "2dsep" ) )
        return( if ( is.null( gp[[ paste0( "algorithm_", axis ) ]] ) ) g$algorithm
                else gp[[ paste0( "algorithm_", axis ) ]] )
    if ( axis == "x" ) g$algorithm else NA_character_
}
cmp <- do.call( rbind, lapply( r.prop$gate, function( gate ) {
    ev <- if ( is.null( hc.ev[[ gate ]] ) ) list() else hc.ev[[ gate ]]
    fit <- unlist( ev$fit )
    g <- hc.defs[[ gate ]]
    data.frame( gate = gate,
                hc.type = if ( is.null( g ) ) NA_character_ else g$gate_type,
                hc.algorithm.x = hc.axis.algorithm( g, "x" ),
                hc.algorithm.y = hc.axis.algorithm( g, "y" ),
                hc.split = !is.null( ev$split ),
                hc.positive.x = num( ev$positive$x ), hc.positive.y = num( ev$positive$y ),
                hc.below.zero.x = num( ev$below_zero$x ), hc.below.zero.y = num( ev$below_zero$y ),
                hc.agreement = if ( length( fit ) ) max( fit ) else NA_real_,
                stringsAsFactors = FALSE )
} ) )
cmp <- cbind( r.prop[ , c( "gate", "type", "algorithm.x", "algorithm.y", "positive.x", "positive.y",
                           "below.zero.x", "below.zero.y", "agreement" ) ], cmp[ , -1 ] )
names( cmp )[ 2:9 ] <- paste0( "r.", names( cmp )[ 2:9 ] )
utils::write.csv( cmp, file.path( out.dir, "conversion_comparison.csv" ), row.names = FALSE )

same <- function( a, b ) ( is.na( a ) & is.na( b ) ) | ( !is.na( a ) & !is.na( b ) & a == b )
unsplit <- !cmp$hc.split & !is.na( cmp$hc.type )
say( "1. conversion: %d gates; %d split into 1-D gates by Honeychrome (not compared); of the other %d:",
     nrow( cmp ), sum( cmp$hc.split ), sum( unsplit ) )
say( "   same type %d, same algorithm on x %d, on y %d",
     sum( same( cmp$r.type, cmp$hc.type )[ unsplit ] ),
     sum( same( cmp$r.algorithm.x, cmp$hc.algorithm.x )[ unsplit ] ),
     sum( same( cmp$r.algorithm.y, cmp$hc.algorithm.y )[ unsplit ] ) )
say( "   largest difference: positive share x %.2g, y %.2g; below zero x %.2g, y %.2g",
     max( abs( cmp$r.positive.x - cmp$hc.positive.x ), na.rm = TRUE ),
     max( abs( cmp$r.positive.y - cmp$hc.positive.y ), na.rm = TRUE ),
     max( abs( cmp$r.below.zero.x - cmp$hc.below.zero.x ), na.rm = TRUE ),
     max( abs( cmp$r.below.zero.y - cmp$hc.below.zero.y ), na.rm = TRUE ) )
scored <- !is.na( cmp$r.agreement ) & !is.na( cmp$hc.agreement )
d <- abs( cmp$r.agreement - cmp$hc.agreement )[ scored ]
say( "   agreement: scored in R %d, in Honeychrome %d, both %d; |difference| median %.2g, 95th percentile %.2g, max %.2g, > 0.05 in %d",
     sum( !is.na( cmp$r.agreement ) ), sum( !is.na( cmp$hc.agreement ) ), sum( scored ),
     stats::median( d ), stats::quantile( d, 0.95 ), max( d ), sum( d > 0.05 ) )


# ---- Population counts ---------------------------------------------------------------

hc.pop <- utils::read.csv( file.path( hc.dir, "populations.csv" ), stringsAsFactors = FALSE,
                           check.names = FALSE )
hc.pop$sample <- names( hc.samples )[ match( hc.pop$sample, hc.samples ) ]
compare.counts <- function( gating, tag ) {
    r <- gating$population.table
    m <- merge( r[ , c( "sample", "key", "n.events", "fraction.of.parent" ) ],
                hc.pop[ , c( "sample", "key", "n_events", "fraction_of_parent" ) ],
                by = c( "sample", "key" ) )
    m$diff.events <- m$n.events - m$n_events
    m$diff.fraction <- m$fraction.of.parent - m$fraction_of_parent
    utils::write.csv( m, file.path( out.dir, paste0( tag, "_counts.csv" ) ), row.names = FALSE )
    big <- abs( m$diff.fraction ) > 0.01 & !is.na( m$diff.fraction )
    say( "   %d populations in both (R %d, Honeychrome %d rows); identical counts %d; |fraction difference| > 0.01 in %d; largest %.3g",
         nrow( m ), nrow( r ), nrow( hc.pop ), sum( m$diff.events == 0, na.rm = TRUE ), sum( big ),
         max( abs( m$diff.fraction ), na.rm = TRUE ) )
    m
}


# ---- 2. Honeychrome's model applied here ------------------------------------------

model.hc <- read.agmodel( file.path( hc.dir, "model.agmodel" ), fcs.setup$fcs.panel )
ap.engine <- ap
ap.engine$statistics.dir <- file.path( out.dir, "statistics_engine" )
gating.engine <- run.autogating( fcs.setup, model.hc, ap.engine )
say( "2. engine (Honeychrome's model):" )
engine <- compare.counts( gating.engine, "engine" )


# ---- 3. Converted here, trained on Honeychrome's training events --------------------

model.r <- new.agmodel( converted$gate.definitions )
trained <- train.boundaries( model.r$gate.definitions, training, transforms = transforms )
model.r$trained.boundaries <- trained$boundaries
model.r$transforms <- transforms
ap.pipeline <- ap
ap.pipeline$statistics.dir <- file.path( out.dir, "statistics_pipeline" )
gating.pipeline <- run.autogating( fcs.setup, model.r, ap.pipeline )
say( "3. pipeline (converted and trained here on Honeychrome's training events):" )
pipeline <- compare.counts( gating.pipeline, "pipeline" )



# ---- 4. Honeychrome's gates trained here ------------------------------------------

model.t <- model.hc
model.t$trained.boundaries <- train.boundaries( model.hc$gate.definitions, training,
                                                transforms = transforms )$boundaries
ap.training <- ap
ap.training$statistics.dir <- file.path( out.dir, "statistics_training" )
gating.training <- run.autogating( fcs.setup, model.t, ap.training )
say( "4. training (Honeychrome's gates trained here on Honeychrome's training events):" )
training.counts <- compare.counts( gating.training, "training" )


# ---- 5. Thresholds on Honeychrome's parent events --------------------------------

fixed <- run.gating( model.hc$gate.definitions, training, reference = model.hc$trained.boundaries,
                     transforms = transforms )
rows <- list()
for ( g in model.hc$gate.definitions ) {
    if ( !identical( g$gate.type, "1dsep" ) || !( g$algorithm %in% c( "tail", "mixture" ) ) )
        next
    parent <- fixed$parent.events[[ g$gate.name ]]
    hc.b <- model.hc$trained.boundaries[[ g$gate.name ]]
    if ( is.null( parent ) || is.null( hc.b ) )
        next
    b <- tryCatch( calculate.gate( training[ parent, , drop = FALSE ], g ), error = function( e ) NULL )
    if ( length( b ) == 0L )
        next
    rows[[ length( rows ) + 1L ]] <- data.frame( gate = g$gate.name, algorithm = g$algorithm,
                                                 n.parent = length( parent ),
                                                 r = b[[ 1 ]]$threshold.x, hc = hc.b[[ 1 ]]$threshold.x,
                                                 stringsAsFactors = FALSE )
}
thresholds <- do.call( rbind, rows )
thresholds$difference <- thresholds$r - thresholds$hc
utils::write.csv( thresholds, file.path( out.dir, "threshold_comparison.csv" ), row.names = FALSE )
say( "5. thresholds on Honeychrome's parent events (display units):" )
for ( a in c( "tail", "mixture" ) ) {
    x <- abs( thresholds$difference[ thresholds$algorithm == a ] )
    say( "   %s: %d gates; |difference| < 1e-6 in %d, > 1e-3 in %d, > 0.01 in %d; largest %.3g",
         a, length( x ), sum( x < 1e-6 ), sum( x > 1e-3 ), sum( x > 0.01 ), max( x ) )
}

writeLines( summary.lines, file.path( out.dir, "summary.txt" ) )
