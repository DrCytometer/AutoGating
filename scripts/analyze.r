# analyze.r
#
# AutoGating pipeline on the ST1 example (BD instrument, 25 samples).
#
# Set fcs.data.dir to the folder holding the ST1 FCS files, then source this
# file.  The gate definitions in scripts/ST1 are in the earlier R format;
# define.gates() converts them (convert.gate.definitions) and trains the
# boundaries on the samples marked gate.calculation in the experiment file.
# The trained model is saved as an .agmodel file that Honeychrome opens.
# Each gate is drawn on the training events in output/figures/training; set
# figure.samples = TRUE to draw every gate on every sample as well.

library( AutoGating )

fcs.data.dir <- "~/path/to/10.01_fcs_data_devel_st1"
st1.dir      <- "~/dev/claude-sandboxes/AutoGating-sandbox/scripts/ST1"
output.dir   <- file.path( st1.dir, "output" )
dir.create( output.dir, showWarnings = FALSE )

ap <- get.autogating.param(
    cytometer                 = "symphony",
    data.dir                  = fcs.data.dir,
    panel.filename            = file.path( st1.dir, "fcs_panel_st1.csv" ),
    experiment.filename       = file.path( st1.dir, "fcs_experiment_st1.csv" ),
    compensation.filename     = file.path( st1.dir, "fcs_compensation_st1.csv" ),
    transformation.filename   = file.path( st1.dir, "fcs_transformation_st1.csv" ),
    population.gates.filename = file.path( st1.dir, "fcs_population_gates_st1.r" ),
    activation.gates.filename = file.path( st1.dir, "fcs_activation_gates_st1.r" ),
    statistics.dir            = file.path( output.dir, "statistics" ),
    figure.dir                = file.path( output.dir, "figures" ),
    worker.n                  = 4L
)

fcs.setup <- setup.autogating( ap )
model     <- define.gates( fcs.setup, ap )
write.agmodel( model, file.path( output.dir, "st1.agmodel" ), fcs.setup$fcs.panel )
gating    <- run.autogating( fcs.setup, model, ap )
