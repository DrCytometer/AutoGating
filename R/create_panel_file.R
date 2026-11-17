# create_panel_file.R
#
# Utility to scaffold the panel CSV for the AutoGating pipeline.


#' Create a panel definition CSV from a single FCS file
#'
#' @title Create Panel Definition CSV
#'
#' @description
#' Reads the channel metadata from a single FCS file and writes a three-column
#' panel CSV (\code{dye}, \code{antigen}, \code{type}) to disk.  The resulting
#' file is intended as a starting scaffold: review and edit the \code{antigen}
#' and \code{type} columns before using it with the AutoGating pipeline.
#'
#' The function:
#' \enumerate{
#'   \item Reads the FCS keywords with
#'     \code{\link[AutoSpectral]{readFCS}}.
#'   \item Removes autofluorescence channels: descriptions \code{AF} or
#'     \code{AF Index}, or containing \code{Autofluorescence} (e.g.
#'     Honeychrome's \code{Autofluorescence Abundance} and
#'     \code{Autofluorescence Index}).
#'   \item Fills missing marker descriptions (\code{NA}) with the channel
#'     name so that every row has a non-empty \code{antigen} value.
#'   \item Assigns a default \code{type} based on channel name:
#'     \code{"time"} for channels matching \code{time} (case-insensitive),
#'     \code{"scatter"} for FSC/SSC channels, and \code{"lineage"} for all
#'     others.  Activation channels must be relabelled manually.
#' }
#'
#' @param fcs.file Character scalar.  Path to the FCS file to read.
#' @param filename Character scalar.  Base name (without extension) of the
#'   output CSV file.  The \code{".csv"} extension is appended automatically.
#'   Default \code{"panel_file"}.
#'
#' @return The panel data frame (invisibly), with columns:
#'   \describe{
#'     \item{\code{dye}}{Character.  FCS channel name (e.g. \code{"BV421-A"}).}
#'     \item{\code{antigen}}{Character.  Marker description from the FCS
#'       parameter block, or the channel name when the description is absent.}
#'     \item{\code{type}}{Character.  Default channel type: one of
#'       \code{"scatter"}, \code{"time"}, or \code{"lineage"}.  Channels
#'       intended as activation markers must be changed to
#'       \code{"activation"} manually.}
#'   }
#'
#' @seealso \code{\link{create.experiment.file}},
#'   \code{\link{get.autogating.param}}
#'
#' @importFrom AutoSpectral readFCS
#' @importFrom utils write.csv
#'
#' @export
create.panel.file <- function( fcs.file, filename = "panel_file" ) {

    # Only the keywords are needed; read a single event
    fcs.kw   <- readFCS( fcs.file, return.keywords = TRUE,
                         start.row = 1, end.row = 1 )$keywords
    n.par    <- as.integer( fcs.kw[[ "$PAR" ]] )
    channels <- vapply( seq_len( n.par ), function( i )
        fcs.kw[[ paste0( "$P", i, "N" ) ]], character( 1L ) )
    markers  <- vapply( seq_len( n.par ), function( i ) {
        desc <- fcs.kw[[ paste0( "$P", i, "S" ) ]]
        if ( is.null( desc ) || trimws( desc ) == "" ) NA_character_ else desc
    }, character( 1L ) )

    # Remove autofluorescence channels - not used in any gate
    non.af.idx <- which( !( markers %in% c( "AF", "AF Index" ) |
                            grepl( "Autofluorescence", markers, fixed = TRUE ) ) )
    channels <- channels[ non.af.idx ]
    markers  <- markers[  non.af.idx ]

    # Fill missing descriptions with the channel name
    missing.marker.idx <- which( is.na( markers ) )
    markers[ missing.marker.idx ] <- channels[ missing.marker.idx ]

    # Assign default channel type using base R conditionals
    # (ordered: time -> scatter -> lineage, so scatter doesn't subsume time)
    channel.type <- ifelse(
        grepl( "time", channels, ignore.case = TRUE ),
        "time",
        ifelse(
            grepl( "^(FSC|SSC)", channels, ignore.case = FALSE ),
            "scatter",
            "lineage"
        )
    )

    panel.df <- data.frame(
        dye      = channels,
        antigen  = markers,
        type     = channel.type,
        stringsAsFactors = FALSE
    )

    out.path <- paste0( filename, ".csv" )
    utils::write.csv( panel.df, file = out.path, row.names = FALSE )
    message( sprintf( "Panel file written to: %s  (%d channels)",
                      out.path, nrow( panel.df ) ) )

    invisible( panel.df )
}
