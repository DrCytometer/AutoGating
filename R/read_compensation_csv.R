#' Load and reorder a compensation matrix from a CSV file
#'
#' @title Read Compensation Matrix CSV
#'
#' @description
#' Reads a square compensation matrix from a CSV file, subsets and reorders
#' its rows and columns to match the lineage and activation antigens in the
#' panel (ordered by dye name), and assigns antigen names to both row and
#' column names.  The first column of the CSV (typically row labels) is
#' discarded.
#'
#' The CSV holds the \emph{spillover} matrix, as exported by FlowJo or BD
#' FACSDiva: rows are fluorochromes, columns are detectors, the diagonal is
#' 1, and entry \eqn{(i, j)} is the fraction of fluorochrome \eqn{i}'s signal
#' detected in detector \eqn{j}.  \code{\link{setup.autogating}} compensates
#' with \eqn{X S^{-1}}.
#'
#' @param compensation.filepath Character.  Path to the compensation CSV file.
#'   The CSV must contain a column for each dye listed in the panel for the
#'   lineage and activation channels (column names are matched after applying
#'   \code{\link[base]{make.names}}).
#' @param fcs.panel \code{data.frame}.  Panel definition with at least columns
#'   \code{dye} and \code{antigen}.
#' @param fcs.panel.lineage.idx Integer vector.  Row indices into
#'   \code{fcs.panel} identifying lineage channels.
#' @param fcs.panel.activation.idx Integer vector.  Row indices into
#'   \code{fcs.panel} identifying activation channels.
#'
#' @return Numeric matrix.  Square compensation matrix with antigen names as
#'   both row and column names, ordered to match the panel.
#'
#' @seealso \code{\link{get.autogating.param}}
#'
#' @export
read.compensation.csv <- function( compensation.filepath, fcs.panel,
                                   fcs.panel.lineage.idx, fcs.panel.activation.idx ) {

  fcs.compensation <- utils::read.csv(
    compensation.filepath,
    stringsAsFactors = FALSE
  )

  fcs.compensation <- data.matrix( fcs.compensation[ , -1 ] )

  fcs.compensation.dye <- fcs.panel$dye[
    c( fcs.panel.lineage.idx, fcs.panel.activation.idx ) ]

  fcs.compensation.antigen <- fcs.panel$antigen[
    c( fcs.panel.lineage.idx, fcs.panel.activation.idx ) ]

  stopifnot(
    sort( make.names( fcs.compensation.dye ) %in% sort( colnames( fcs.compensation ) ) )
  )

  reorder <- match( make.names( fcs.compensation.dye ),
                    colnames( fcs.compensation ) )

  fcs.compensation <- fcs.compensation[ reorder, reorder ]
  rownames( fcs.compensation ) <- fcs.compensation.antigen
  colnames( fcs.compensation ) <- fcs.compensation.antigen

  return(fcs.compensation)
}
