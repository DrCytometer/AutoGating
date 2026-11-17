# globals.R
#
# Package-level declarations.
#
# Names used as unbound symbols in ggplot2 aesthetics are declared so that
# R CMD check does not report "no visible binding for global variable".  No
# runtime state is stored here; pipeline state flows through the 'ap'
# parameter object and function arguments.

utils::globalVariables( c(
    # gate figures
    "x", "y", "z", "at", "group", "text", "hjust", "vjust",
    # statistics figures
    "logFC", "neglog10p", "significance", "display", "comparison",
    "value", "group.position", "star"
) )

# Null default: a unless it is NULL, else b (the same as base R 4.4.0 and
# later, defined here for earlier versions of R).
`%||%` <- function( a, b ) if ( is.null( a ) ) b else a
