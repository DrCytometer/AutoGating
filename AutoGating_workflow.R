# possible AutoGating workflow

library( AutoGating ) # or whatever we want to call it

# need checks
experiment <- read.experiment( parameter.folder, "fcs_experiment_st1.csv" )

autogate.param <- get.autogating.param( "symphony", experiment )

parameter.folder <- "./10.02_analysis_devel_st1_original/"
data.folder <- "./10.01_fcs_data_devel_st1"

create.directories( autogate.param, experiment )

# need checks
panel <- read.panel.design( parameter.folder, "fcs_panel_st1.csv" )

# define the gates
flow.gates <- set.up.autogating( parameter.folder, data.folder,
                                 "fcs_transformation_st1.csv",
                                 experiment, panel, autogate.param,
                                 "fcs_population_gates_st1.csv", "fcs_activation_gates_st1.csv",
                                 compensate = TRUE, "fcs_compensation_st1.csv",
                                 setup = TRUE )
# use the gates
autogate()
