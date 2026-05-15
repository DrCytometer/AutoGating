# create_experiment_file.R
#
# Utility to scaffold the experiment manifest CSV for the AutoGating pipeline.


#' Create an experiment manifest CSV from a metadata table and FCS directory
#'
#' @title Create Experiment Manifest CSV
#'
#' @description
#' Scaffolds the \code{fcs_experiment.csv} file required by the AutoGating
#' pipeline.  For each row of \code{metadata}, the function searches
#' \code{fcs.dir} for an FCS filename that contains the sample identifier
#' (case-insensitive partial match), then writes a CSV with columns
#' \code{filename}, \code{sample}, \code{gate.calculation}, and, if a
#' \code{condition.col} is specified, \code{Group}.
#'
#' The resulting file is written to \code{filename} and returned invisibly as
#' a data frame so the caller can inspect it before proceeding.
#'
#' @details
#' The match between sample identifiers and FCS filenames is performed by
#' \code{\link[base]{grep}} with \code{ignore.case = TRUE}.  If a sample
#' identifier matches more than one file, the first match is used and a
#' warning is issued.  If no match is found, \code{NA} is stored in the
#' \code{filename} column and a warning is issued; these rows should be
#' corrected before running the pipeline.
#'
#' All samples are initialised with \code{gate.calculation = TRUE}.  Set
#' individual rows to \code{FALSE} in the saved CSV to reuse boundaries
#' calculated from other samples (see \code{CLAUDE.md} for details).
#'
#' @param metadata data.frame.  Sample metadata table.  Must contain at least
#'   a column of sample identifiers (named by \code{sample.col}).
#' @param fcs.dir Character scalar.  Path to the directory containing FCS
#'   files.  Passed to \code{\link[base]{list.files}}.
#' @param filename Character scalar.  Output CSV filename (including path if
#'   needed).  Default \code{"fcs_experiment.csv"}.
#' @param sample.col Character scalar.  Name of the column in \code{metadata}
#'   that holds sample identifiers used for filename matching.
#'   Default \code{"Sample"}.
#' @param condition.col Character scalar or \code{NULL}.  Name of the column
#'   in \code{metadata} to write as the \code{Group} column in the output.
#'   \code{NULL} omits this column.  Default \code{NULL}.
#'
#' @return The experiment manifest data frame (invisibly), with columns:
#'   \describe{
#'     \item{\code{filename}}{Character.  FCS filename (relative to
#'       \code{fcs.dir}), or \code{NA} if no match was found.}
#'     \item{\code{sample}}{Character.  Sample identifier from
#'       \code{metadata[[sample.col]]}.}
#'     \item{\code{gate.calculation}}{Logical.  \code{TRUE} for all rows.}
#'     \item{\code{Group}}{Character (only present when \code{condition.col}
#'       is not \code{NULL}).  Condition / group label from
#'       \code{metadata[[condition.col]]}.}
#'   }
#'
#' @seealso \code{\link{create.panel.file}}, \code{\link{get.autogating.param}}
#'
#' @importFrom utils read.csv write.csv
#'
#' @export
create.experiment.file <- function(
    metadata,
    fcs.dir,
    filename      = "fcs_experiment.csv",
    sample.col    = "Sample",
    condition.col = NULL
) {
    stopifnot( is.data.frame( metadata ) )
    stopifnot( sample.col %in% colnames( metadata ) )

    if ( !is.null( condition.col ) )
        stopifnot( condition.col %in% colnames( metadata ) )

    fcs.files   <- list.files( fcs.dir )
    sample.ids  <- metadata[[ sample.col ]]

    # Match each sample identifier to an FCS filename
    matched.filenames <- vapply( sample.ids, function( s ) {
        hits <- grep( s, fcs.files, ignore.case = TRUE, value = TRUE )
        if ( length( hits ) == 0L ) {
            warning( sprintf(
                "No FCS file found matching sample identifier '%s'.", s ) )
            return( NA_character_ )
        }
        if ( length( hits ) > 1L ) {
            warning( sprintf(
                "Multiple FCS files match sample identifier '%s'; using first: %s",
                s, hits[ 1L ] ) )
        }
        hits[ 1L ]
    }, character( 1L ) )

    experiment.df <- data.frame(
        filename         = matched.filenames,
        sample           = sample.ids,
        gate.calculation = TRUE,
        stringsAsFactors = FALSE
    )

    if ( !is.null( condition.col ) )
        experiment.df$Group <- metadata[[ condition.col ]]

    n.missing <- sum( is.na( experiment.df$filename ) )
    if ( n.missing > 0L )
        warning( sprintf(
            "%d sample(s) could not be matched to an FCS file. Edit '%s' before running the pipeline.",
            n.missing, filename ) )

    utils::write.csv( experiment.df, file = filename, row.names = FALSE )
    message( sprintf( "Experiment manifest written to: %s  (%d samples)",
                      filename, nrow( experiment.df ) ) )

    invisible( experiment.df )
}
