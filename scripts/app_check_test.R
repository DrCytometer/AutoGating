# app_check_test.R
#
# Opens the AutoGating app on the simulated Aurora experiment of
# figures_check_test.R (sourced first: it builds fcs.setup, ap and model
# and writes its figures), for a check of the app by hand:
#   - Gates: select gates, edit a population name (children follow), click
#     the CD4 x CD8 figure to set manual thresholds, add and delete a gate,
#     Undo;
#   - Transforms: change W of a fluorescence channel and Apply;
#   - Save and run: save the model, run on every sample, Finish.
# The model and experiment returned by Finish are in `result`.

source( file.path( "scripts", "figures_check_test.R" ) )

result <- autogating.app( fcs.setup, ap, model )
if ( !is.null( result ) )
    cat( sprintf( "Returned a model with %d gate(s).\n", length( result$model$gate.definitions ) ) )
