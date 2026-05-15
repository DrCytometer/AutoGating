# globals.R
#
# Package-level global variable declarations.
#
# Declares variables that appear as unbound symbols in ggplot2 aesthetics or
# NSE calls so that R CMD check does not raise "no visible binding for global
# variable" NOTEs.  No runtime state is stored here; all pipeline state flows
# through the 'ap' parameter object and function arguments.

utils::globalVariables( c(
    # ---- ggplot2 aesthetic variables used in build_gate_plot.R ---------------
    "x",        # mapped in aes(x = x, ...)
    "y",        # mapped in aes(... y = y, ...)
    "z",         # mapped in aes(colour = z) for density colouring
    "logFC",
    "neglog10p",
    "Significance",
    "Feature"
) )
