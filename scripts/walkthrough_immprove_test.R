# walkthrough_immprove_test.R
#
# AutoGating walkthrough on the UCL IMMPROVE batch controls: the two buffy
# coat / PBMC "IC" samples, unmixed in Honeychrome, gated with the IMMPROVE
# FlowJo template (20251212_IMMPROVE panel template_new_panel_v2.wspt).
#
# Steps: panel and experiment files -> parameters -> read the data ->
# import the FlowJo template -> convert drawn gates to calculated gates ->
# train -> save the .agmodel -> gate every sample -> figures.
#
# Run from any working directory; set the three paths below.  Outputs go to
# out.dir.

library( AutoGating )

data.dir <- "~/Experiments/UCL IMMPROVE/Unmixed/IC controls"
wspt     <- "~/dev/claude-sandboxes/FlowJo_workspace/20251212_IMMPROVE panel template_new_panel_v2.wspt"
out.dir  <- "~/Experiments/UCL IMMPROVE/AutoGating_walkthrough"

dir.create( out.dir, showWarnings = FALSE, recursive = TRUE )
fcs.files <- list.files( data.dir, pattern = "\\.fcs$" )


# ---- 1. Panel file ----------------------------------------------------------
# Drafted from the first FCS file (autofluorescence channels are left out);
# Honeychrome's event_id, when present, is not gated either.
panel <- create.panel.file( file.path( data.dir, fcs.files[ 1 ] ),
                            filename = file.path( out.dir, "panel" ) )
panel <- panel[ panel$dye != "event_id", ]

# The IMMPROVE panel stacks two markers on some fluorophores (one on T
# cells, one on B or myeloid cells).  The exported FCS names only the first,
# so the antigen names are completed here; they are the axis titles of
# every figure and the column names of every table.
stacked <- c(
    "Spark UV 387"   = "CD8_CD11b",     "BUV395"       = "CD4_IgM",
    "BUV563"         = "CD127_FCRL5",   "BUV661"       = "ICOS_CD80",
    "BUV737"         = "CD161_CD21",    "SBUV770"      = "CD28_CD20",
    "BV421"          = "Perforin_CD40", "BV605"        = "TCRVa72_CD1c",
    "BV650"          = "IgG_TCRgd",     "BV711"        = "CTLA4_BAFFR",
    "BV750"          = "LAG3_CD10",     "SBV515"       = "CD103_CD163",
    "BB700"          = "CD25_CD141",    "RB705"        = "TCF1_CD138",
    "RB545"          = "TCRVa24_TCRVg9_CD24",
    "PerCP-Fire 806" = "KLRG1_IgD",     "PE"           = "CD134_IgE",
    "PE-eFluor 610"  = "GATA3_CD23",    "PE-Cy5"       = "NKp46_CD86",
    "APC"            = "CD49a_IgA",     "R718"         = "Tbet_CD123" )
hit <- panel$dye %in% names( stacked )
panel$antigen[ hit ] <- stacked[ panel$dye[ hit ] ]
utils::write.csv( panel, file.path( out.dir, "panel.csv" ), row.names = FALSE )


# ---- 2. Experiment file -----------------------------------------------------
# One row per FCS file; both samples are used to train the gates.
samples <- data.frame( Sample = c( "BuffyCoat_Donor_E", "PBMC_Donor_1" ),
                       source = c( "buffy coat", "PBMC" ) )
create.experiment.file( samples, data.dir,
                        filename      = file.path( out.dir, "experiment.csv" ),
                        sample.col    = "Sample",
                        condition.col = "source" )


# ---- 3. Parameters ----------------------------------------------------------
ap <- get.autogating.param(
    cytometer           = "aurora",
    data.dir            = data.dir,
    panel.filename      = file.path( out.dir, "panel.csv" ),
    experiment.filename = file.path( out.dir, "experiment.csv" ),
    figure.dir          = file.path( out.dir, "figures" ),
    statistics.dir      = file.path( out.dir, "statistics" ),
    events.per.sample   = 50000L
)


# ---- 4. Read the data -------------------------------------------------------
fcs.setup <- setup.autogating( ap )
ap <- fcs.setup$ap        # carries the transforms set from the data


# ---- 5. Import the FlowJo template ------------------------------------------
# Matches the template's parameters to the panel (by fluorophore here) and
# reads every gate as drawn ("imported").
imported <- import.flowjo.workspace( wspt, fcs.setup, ap )
table( imported$mapping$method )


# ---- 6. Convert drawn gates to calculated gates ------------------------------
# On the events the model is trained on: 50,000 per training sample.
ts <- training.set( fcs.setup, ap )
converted <- convert.imported.gates( imported$gate.definitions,
                                     ts$data,
                                     transforms = ap$transform.params,
                                     panel      = fcs.setup$fcs.panel )
table( converted$proposals$type, converted$proposals$algorithm.x )
utils::write.csv( converted$proposals, file.path( out.dir, "conversion_proposals.csv" ),
                  row.names = FALSE )


# ---- 7. Train and save the model --------------------------------------------
model <- define.gates( fcs.setup, ap, model = new.agmodel( converted$gate.definitions ) )
write.agmodel( model, file.path( out.dir, "IMMPROVE_IC.agmodel" ), fcs.setup$fcs.panel )


# ---- 8. Gate every sample ----------------------------------------------------
gating <- run.autogating( fcs.setup, model, ap )
head( gating$population.table )
table( gating$status$flag, useNA = "ifany" )


# ---- 9. Selected gates on each sample -----------------------------------------
# (figure.samples = TRUE in the parameters draws every gate on every sample.)
show <- c( "Cells", "Single Cells (Cells)", "FSC-A, viability subset", "CD14+CD11b+ & CD14-",
           "B cells (CD14-)", "CD3",
           "non-Trans non-Plasma & Plasmablasts & Transitional B cells",
           "Activated & Atypical & Mature Naive & Memory" )
gate.names <- vapply( model$gate.definitions, `[[`, character( 1 ), "gate.name" )
for ( s in fcs.setup$fcs.sample ) {
    rows <- fcs.setup$fcs.event.sample == s
    write.gate.figures( model$gate.definitions[ gate.names %in% show ], gating$results[[ s ]],
                        fcs.setup$fcs.expr.data[ rows, ], fcs.setup, ap,
                        dir = file.path( ap$figure.dir, "samples", s ), title = s )
}


# ---- 10. (Optional) review and edit in the app -------------------------------
# edited <- autogating.app( fcs.setup, ap, model )
# model  <- edited$model
