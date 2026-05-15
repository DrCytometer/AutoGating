# calculate_differential_abundance.R
#
# Differential abundance analysis for the FCS autogating pipeline.
#
# Reads Group (and optional covariate) columns from the experiment CSV,
# runs limma on the combined population + activation statistics, and
# produces a pheatmap (top-k features) and volcano plot(s).
#
# Public API:
#   calculate.differential.abundance()
#
# Experiment CSV conventions
# --------------------------
# The pipeline's original experiment CSV has (at minimum):
#   filename         - FCS file path, one row per sample
#   sample / Sample  - human-readable sample name (case-insensitive match)
#   gate.calculation - logical flag controlling gate recalculation
#
# For differential abundance a *separate* (or modified) experiment CSV may be
# supplied via the `experiment.filename` argument.  That file is allowed to:
#   - contain a subset of the samples present in the statistics data.frames
#   - add or modify the following columns:
#       Group          mandatory factor defining comparison groups
#       Block          optional blocking variable (donor, time-point, etc.);
#                      triggers limma::duplicateCorrelation()
#       Covariate.*    any column whose name starts with "Covariate." is
#                      added as a nuisance covariate in the design matrix
#
# The reload file must NOT change the structural identity columns
# (filename, sample / Sample, gate.calculation) for any row it contains;
# those are cross-checked against the original ap$experiment.filename.


# ==============================================================================
# calculate.differential.abundance
# ==============================================================================

#' Differential abundance analysis of population and activation statistics
#'
#' @title Differential Abundance Analysis
#'
#' @description
#' Performs group-level differential abundance analysis on the tidy statistics
#' tables produced by \code{\link{calculate.population.statistics}} and
#' \code{\link{calculate.activation.statistics}}.
#'
#' The function:
#' \enumerate{
#'   \item Merges the population and activation statistics tables by
#'     \code{Sample}.
#'   \item Reads group assignments (and optional blocking or covariate columns)
#'     from an experiment CSV — either \code{ap$experiment.filename} or a
#'     separate file supplied via \code{experiment.filename}.
#'   \item Fits a \code{\link[limma]{lmFit}} linear model with a no-intercept
#'     group design.  If a \code{Block} column is present,
#'     \code{\link[limma]{duplicateCorrelation}} is used to account for
#'     repeated-measures or paired structure.  Columns named
#'     \code{Covariate.*} are appended to the design matrix as nuisance
#'     regressors.
#'   \item Runs \code{\link[limma]{eBayes}} and extracts pairwise contrasts
#'     via \code{\link[limma]{topTable}} (Benjamini-Hochberg FDR correction).
#'   \item Saves a \code{\link[pheatmap]{pheatmap}} of the top-\eqn{k} features
#'     ranked by F-statistic, and one volcano plot per contrast, to
#'     \code{figure.dir}.
#' }
#'
#' When \code{experiment.filename} points to a different file from
#' \code{ap$experiment.filename}, the reloaded file is validated: every sample
#' it names must appear in both the original experiment file and the statistics
#' tables, and the structural columns (\code{filename},
#' \code{gate.calculation}) must be unchanged.  Only \code{Group},
#' \code{Block}, and \code{Covariate.*} columns may differ.  Supplying a
#' subset of samples is permitted; the analysis is restricted to those rows.
#'
#' @param population.stats data.frame returned by
#'   \code{\link{calculate.population.statistics}}.  First column must be
#'   \code{"Sample"}.
#' @param activation.stats data.frame returned by
#'   \code{\link{calculate.activation.statistics}}, or \code{NULL} to perform
#'   population-level analysis only.  First column must be \code{"Sample"}.
#'   Default \code{NULL}.
#' @param ap Autogating parameter object (see
#'   \code{\link{get.autogating.param}}).  Used for
#'   \code{ap$experiment.filename} and \code{ap$statistics.dir}.
#' @param experiment.filename Character scalar or \code{NULL}.  Path to a
#'   modified experiment CSV that adds \code{Group}, \code{Block}, and/or
#'   \code{Covariate.*} columns.  When \code{NULL} (default),
#'   \code{ap$experiment.filename} is used directly.
#' @param figure.dir Character scalar.  Directory for output figures.
#'   Created recursively if absent.  Default
#'   \code{file.path(ap$statistics.dir, "figures")}.
#' @param top.k Integer.  Number of features to display in the heatmap,
#'   ranked by F-statistic.  Default \code{30L}.
#' @param fdr.threshold Numeric.  Adjusted p-value significance threshold for
#'   colouring points on the volcano plot.  Default \code{0.05}.
#' @param lfc.threshold Numeric.  Absolute log\eqn{_2} fold-change threshold
#'   for the volcano significance call.  Default \code{0.5}.
#' @param volcano.groups Character vector of length \eqn{\ge 2} specifying
#'   which group levels to compare.  Required when there are more than 4
#'   groups (to avoid generating an unmanageable number of plots); ignored
#'   otherwise.  When \code{NULL} and \eqn{\le 4} groups, all pairwise
#'   combinations are used.  Default \code{NULL}.
#' @param heatmap.filename Character scalar.  Base filename (without
#'   extension) for the heatmap PNG.
#'   Default \code{"differential_abundance_heatmap"}.
#' @param volcano.filename.prefix Character scalar.  Prefix for volcano PNG
#'   filenames; the contrast name is appended.
#'   Default \code{"differential_abundance_volcano"}.
#' @param scale.rows Logical.  If \code{TRUE} (default), rows are z-scored
#'   before hierarchical clustering in the heatmap.
#' @param width.heatmap Numeric.  Heatmap figure width in inches.
#'   Default \code{10}.
#' @param height.heatmap Numeric.  Heatmap figure height in inches.
#'   Default \code{10}.
#' @param width.volcano Numeric.  Volcano figure width in inches.
#'   Default \code{7}.
#' @param height.volcano Numeric.  Volcano figure height in inches.
#'   Default \code{6}.
#'
#' @return A named list (returned invisibly) with elements:
#'   \describe{
#'     \item{\code{fit}}{The \code{MArrayLM} object returned by
#'       \code{\link[limma]{eBayes}} after contrast fitting.}
#'     \item{\code{results}}{Named list of \code{topTable} data frames, one
#'       per pairwise contrast (e.g. \code{"A_vs_B"}).  Rows are features,
#'       columns are the standard \code{limma} output columns
#'       (\code{logFC}, \code{AveExpr}, \code{t}, \code{P.Value},
#'       \code{adj.P.Val}, \code{B}).}
#'     \item{\code{expr.mat}}{Numeric matrix, features \eqn{\times} samples,
#'       used as input to \code{lmFit}.}
#'     \item{\code{meta}}{data.frame of sample metadata after validation and
#'       subsetting (samples with blank \code{Group} removed).}
#'   }
#'
#' @seealso \code{\link{calculate.population.statistics}},
#'   \code{\link{calculate.activation.statistics}},
#'   \code{\link{run.autogating}}
#'
#' @importFrom limma lmFit duplicateCorrelation contrasts.fit eBayes
#'   makeContrasts topTable
#' @importFrom pheatmap pheatmap
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_manual geom_vline
#'   geom_hline xlim labs theme_bw theme element_text ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices colorRampPalette png dev.off
#' @importFrom stats model.matrix sd
#' @importFrom utils combn read.csv
#'
#' @export
calculate.differential.abundance <- function(
    population.stats,
    activation.stats        = NULL,
    ap,
    experiment.filename     = NULL,
    figure.dir              = file.path( ap$statistics.dir, "figures" ),
    top.k                   = 30L,
    fdr.threshold           = 0.05,
    lfc.threshold           = 0.5,
    volcano.groups          = NULL,
    heatmap.filename        = "differential_abundance_heatmap",
    volcano.filename.prefix = "differential_abundance_volcano",
    scale.rows              = TRUE,
    width.heatmap           = 10,
    height.heatmap          = 10,
    width.volcano           = 7,
    height.volcano          = 6
) {

  # ---- package checks -------------------------------------------------------

  for ( pkg in c( "limma", "pheatmap", "ggplot2", "ggrepel" ) ) {
    if ( !requireNamespace( pkg, quietly = TRUE ) )
      stop( sprintf( "Package '%s' is required. Install it with install.packages('%s').",
                     pkg, pkg ) )
  }

  # ---- merge statistics tables ---------------------------------------------

  stats <- .da.merge.stats( population.stats, activation.stats )

  # ---- resolve and validate experiment file --------------------------------

  # Determine which file to use for metadata
  reload.path <- if ( !is.null( experiment.filename ) ) experiment.filename
  else                                   ap$experiment.filename

  meta <- .da.load.meta( reload.path )

  # When the user has supplied a different file, cross-check it against the
  # original experiment CSV and against the statistics data.frames.
  if ( !is.null( experiment.filename ) &&
       !identical( normalizePath( experiment.filename, mustWork = FALSE ),
                   normalizePath( ap$experiment.filename, mustWork = FALSE ) ) ) {

    message( sprintf( "Re-loading experiment file: %s", experiment.filename ) )
    original.meta <- .da.load.meta( ap$experiment.filename )
    .da.validate.reload( meta, original.meta, stats )

  } else {
    # Using the original file - standard alignment check only
    missing.from.meta <- setdiff( stats$Sample, meta$Sample )
    if ( length( missing.from.meta ) > 0L )
      warning( sprintf(
        "%d sample(s) in statistics tables have no match in the experiment CSV and will be dropped: %s",
        length( missing.from.meta ),
        paste( missing.from.meta, collapse = ", " ) ) )
  }

  # ---- align samples (meta drives subsetting after validation) -------------

  # All samples in meta must be in stats (enforced above in reload path;
  # applied here for both paths).
  missing.from.stats <- setdiff( meta$Sample, stats$Sample )
  if ( length( missing.from.stats ) > 0L )
    stop( sprintf(
      "The following sample(s) appear in the experiment file but are absent from the statistics tables: %s",
      paste( missing.from.stats, collapse = ", " ) ) )

  # Restrict stats to the samples present in meta (allows subsetting)
  stats <- stats[ stats$Sample %in% meta$Sample, , drop = FALSE ]

  # Enforce identical row order (meta order is canonical)
  stats <- stats[ match( meta$Sample, stats$Sample ), , drop = FALSE ]
  stopifnot( identical( meta$Sample, stats$Sample ) )

  # ---- build feature matrix (samples x features) ---------------------------

  # Select feature columns explicitly by name (exclude Sample; guard against
  # any other non-numeric columns that may have survived the merge).
  feat.cols <- setdiff( colnames( stats ), "Sample" )

  non.numeric <- feat.cols[ !vapply( stats[ feat.cols ],
                                     is.numeric, logical( 1L ) ) ]
  if ( length( non.numeric ) > 0L ) {
    warning( sprintf(
      "Dropping %d non-numeric column(s) from the feature matrix: %s",
      length( non.numeric ),
      paste( non.numeric, collapse = ", " ) ) )
    feat.cols <- setdiff( feat.cols, non.numeric )
  }

  if ( length( feat.cols ) == 0L )
    stop( "No numeric feature columns remain after filtering. Check the statistics data.frames." )

  # Guard against duplicate feature names (e.g. identical column names across
  # the population and activation stats tables after merging).  limma's eBayes
  # calls cut() on fit$Amean which becomes NULL when rownames are duplicated,
  # producing the cryptic "x must be numeric" error from cut.default.
  if ( anyDuplicated( feat.cols ) ) {
    dup.names <- unique( feat.cols[ duplicated( feat.cols ) ] )
    warning( sprintf(
      "Duplicate feature name(s) detected and disambiguated with a numeric suffix: %s",
      paste( dup.names, collapse = ", " ) ) )
    feat.cols.unique <- make.unique( feat.cols, sep = "." )
  } else {
    feat.cols.unique <- feat.cols
  }

  # Select columns by position (feat.cols holds the original names, which may
  # be duplicated in stats); then apply the deduplicated names to the matrix.
  # We use a positional index built from the full colnames vector so that both
  # copies of a duplicated name are correctly selected (match() would only find
  # the first occurrence of each).
  all.cols     <- colnames( stats )
  feat.col.idx <- which( all.cols != "Sample" )   # same order as feat.cols
  feat.mat <- as.matrix( stats[ , feat.col.idx, drop = FALSE ] )
  colnames( feat.mat ) <- feat.cols.unique
  storage.mode( feat.mat ) <- "double"
  rownames( feat.mat ) <- stats$Sample

  # Replace NA / NaN / Inf with 0
  feat.mat[ !is.finite( feat.mat ) ] <- 0

  # Transpose to features x samples (limma convention)
  expr.mat <- t( feat.mat )

  # Sanity check: limma requires non-NULL, non-empty, unique rownames
  if ( is.null( rownames( expr.mat ) ) || anyDuplicated( rownames( expr.mat ) ) )
    stop( "Internal error: expr.mat rownames are NULL or duplicated after feature construction." )

  # ---- group factor --------------------------------------------------------

  if ( !"Group" %in% colnames( meta ) )
    stop( "Experiment CSV has no 'Group' column. Add one to enable differential abundance analysis." )

  # Treat empty strings and whitespace-only values as missing
  group.raw <- trimws( as.character( meta$Group ) )
  blank.idx  <- group.raw == "" | is.na( group.raw ) | group.raw == "NA"

  if ( all( blank.idx ) )
    stop( "All samples have a blank Group assignment. Please fill in the 'Group' column." )

  if ( any( blank.idx ) ) {
    dropped <- meta$Sample[ blank.idx ]
    warning( sprintf(
      "%d sample(s) have no Group assignment and will be excluded from the analysis: %s",
      length( dropped ),
      paste( dropped, collapse = ", " ) ) )
    meta  <- meta[  !blank.idx, , drop = FALSE ]
    stats <- stats[ !blank.idx, , drop = FALSE ]
    expr.mat <- expr.mat[ , !blank.idx, drop = FALSE ]
  }

  meta$Group <- factor( group.raw[ !blank.idx ] )
  group.levels <- levels( meta$Group )
  n.groups <- length( group.levels )

  if ( n.groups < 2L )
    stop( paste0(
      "At least two groups are required for differential abundance analysis. ",
      "Only one group found after removing blank assignments: '",
      group.levels, "'." ) )

  # Ensure every group has at least 2 samples (required for limma variance
  # estimation - a group of 1 contributes no within-group information).
  group.sizes <- table( meta$Group )
  singleton.groups <- names( group.sizes )[ group.sizes < 2L ]
  if ( length( singleton.groups ) > 0L )
    stop( sprintf(
      "The following group(s) contain fewer than 2 samples, which is insufficient for variance estimation: %s.\nEither assign more samples to these groups or exclude them.",
      paste( sprintf( "'%s' (n=%d)", singleton.groups,
                      group.sizes[ singleton.groups ] ),
             collapse = ", " ) ) )

  # ---- decide which pairwise contrasts to run ------------------------------

  contrast.groups <- .da.resolve.contrasts(
    n.groups     = n.groups,
    group.levels = group.levels,
    volcano.groups = volcano.groups
  )

  # ---- build design matrix + covariates ------------------------------------

  design.result <- .da.build.design( meta )
  design  <- design.result$design
  block   <- design.result$block    # NULL or a vector for duplicateCorrelation

  # ---- optional duplicate correlation (paired / repeated measures) ---------

  corfit <- NULL
  if ( !is.null( block ) ) {
    message( "Estimating intra-block correlation for repeated-measures / paired design ..." )
    corfit <- limma::duplicateCorrelation( expr.mat, design, block = block )
    message( sprintf( "  consensus correlation = %.3f", corfit$consensus.correlation ) )
  }

  # ---- limma fit -----------------------------------------------------------

  if ( is.null( corfit ) ) {
    fit <- limma::lmFit( expr.mat, design )
  } else {
    fit <- limma::lmFit( expr.mat, design,
                         block       = block,
                         correlation = corfit$consensus.correlation )
  }

  # Build contrast matrix for all pairwise comparisons among selected groups
  contrast.strings <- vapply( contrast.groups, function( pair ) {
    sprintf( "Group%s-Group%s", pair[ 2L ], pair[ 1L ] )
  }, character( 1L ) )
  names( contrast.strings ) <- vapply( contrast.groups, function( pair ) {
    paste( pair, collapse = "_vs_" )
  }, character( 1L ) )

  contrast.mat <- limma::makeContrasts(
    contrasts = contrast.strings,
    levels    = design
  )
  colnames( contrast.mat ) <- names( contrast.strings )

  fit2 <- limma::contrasts.fit( fit, contrast.mat )

  # contrasts.fit() does not propagate fit$Amean, which eBayes needs when
  # trend=TRUE or robust=TRUE.  Copy it forward explicitly so downstream
  # calls with those options do not fail with "x must be numeric" in cut().
  fit2$Amean <- fit$Amean

  fit2 <- limma::eBayes( fit2 )

  # ---- extract results for each contrast -----------------------------------

  results <- lapply( names( contrast.strings ), function( cname ) {
    limma::topTable( fit2, coef = cname,
                     number  = Inf,
                     sort.by = "P",
                     adjust.method = "BH" )
  } )
  names( results ) <- names( contrast.strings )

  # ---- create output directory ---------------------------------------------

  dir.create( figure.dir, recursive = TRUE, showWarnings = FALSE )

  # ---- heatmap (top-k features across all contrasts) -----------------------

  .da.plot.heatmap(
    fit2          = fit2,
    expr.mat      = expr.mat,
    meta          = meta,
    top.k         = top.k,
    scale.rows    = scale.rows,
    figure.dir    = figure.dir,
    filename      = heatmap.filename,
    width         = width.heatmap,
    height        = height.heatmap
  )

  # ---- volcano plots -------------------------------------------------------

  for ( cname in names( results ) ) {
    .da.plot.volcano(
      toptable      = results[[ cname ]],
      contrast.name = cname,
      fdr.threshold = fdr.threshold,
      lfc.threshold = lfc.threshold,
      figure.dir    = figure.dir,
      filename.prefix = volcano.filename.prefix,
      width         = width.volcano,
      height        = height.volcano
    )
  }

  message( sprintf( "Differential abundance complete. Figures written to: %s", figure.dir ) )

  invisible( list(
    fit     = fit2,
    results = results,
    expr.mat = expr.mat,
    meta    = meta
  ) )
}


# ==============================================================================
# Internal helpers
# ==============================================================================

# ------------------------------------------------------------------------------
# .da.merge.stats
# ------------------------------------------------------------------------------

#' Merge population and activation statistics tables
#'
#' @description
#' Outer-joins \code{population.stats} and \code{activation.stats} by the
#' \code{Sample} column, producing a single wide data frame.  If
#' \code{activation.stats} is \code{NULL}, \code{population.stats} is returned
#' unchanged.
#'
#' @param population.stats data.frame with a \code{Sample} column.
#' @param activation.stats data.frame with a \code{Sample} column, or
#'   \code{NULL}.
#'
#' @return data.frame of merged statistics.
#'
#' @keywords internal
.da.merge.stats <- function( population.stats, activation.stats ) {

  pop <- population.stats

  if ( !"Sample" %in% colnames( pop ) )
    stop( "population.stats must have a 'Sample' column." )

  if ( !is.null( activation.stats ) ) {
    if ( !"Sample" %in% colnames( activation.stats ) )
      stop( "activation.stats must have a 'Sample' column." )

    # Merge by Sample; outer join so neither side loses rows
    merged <- merge( pop, activation.stats, by = "Sample", all = TRUE )
    return( merged )
  }

  pop
}


# ------------------------------------------------------------------------------
# .da.load.meta
# ------------------------------------------------------------------------------

#' Read an experiment CSV and normalise the Sample column
#'
#' @description
#' Reads a CSV file and ensures the sample-identity column is named
#' \code{"Sample"} (accepting \code{"sample"} and \code{"SAMPLE"} as
#' case-insensitive aliases).
#'
#' @param path Character scalar.  Path to the experiment CSV.
#'
#' @return data.frame with at minimum a \code{Sample} column, plus whatever
#'   additional columns the file contains (\code{filename},
#'   \code{gate.calculation}, \code{Group}, \code{Block},
#'   \code{Covariate.*}, etc.).
#'
#' @keywords internal
.da.load.meta <- function( path ) {

  if ( is.null( path ) || !file.exists( path ) )
    stop( sprintf( "Experiment file not found: '%s'", path ) )

  meta <- utils::read.csv( path, stringsAsFactors = FALSE )

  # Normalise the sample-identity column to "Sample" (title-case)
  col.map <- c( "Sample", "sample", "SAMPLE" )
  hit <- intersect( col.map, colnames( meta ) )[ 1L ]
  if ( is.na( hit ) )
    stop( sprintf(
      "Experiment CSV '%s' must have a 'Sample' column (case-insensitive).",
      path ) )
  if ( hit != "Sample" )
    colnames( meta )[ colnames( meta ) == hit ] <- "Sample"

  meta
}


# ------------------------------------------------------------------------------
# .da.validate.reload
# ------------------------------------------------------------------------------

#' Validate a reloaded experiment CSV against the original and statistics tables
#'
#' @description
#' Enforces three rules when the user supplies a different experiment CSV for
#' differential abundance analysis:
#' \enumerate{
#'   \item Every sample in \code{reload} must appear in \code{stats} (subset
#'     is allowed; extras in \code{stats} are fine).
#'   \item Every sample in \code{reload} must appear in \code{original}
#'     (samples cannot be invented that were never gated).
#'   \item For samples shared between \code{reload} and \code{original}, the
#'     structural identity columns (\code{filename}, \code{gate.calculation})
#'     must be unchanged.  Only \code{Group}, \code{Block}, and
#'     \code{Covariate.*} columns may differ.
#' }
#'
#' @param reload data.frame.  The reloaded experiment metadata.
#' @param original data.frame.  The original experiment metadata from
#'   \code{ap$experiment.filename}.
#' @param stats data.frame.  Merged statistics table (output of
#'   \code{.da.merge.stats}).
#'
#' @return \code{NULL} invisibly.  Stops with an informative error on any
#'   violation.
#'
#' @keywords internal
.da.validate.reload <- function( reload, original, stats ) {

  # ---- 1. All reload samples must appear in stats --------------------------

  missing.from.stats <- setdiff( reload$Sample, stats$Sample )
  if ( length( missing.from.stats ) > 0L )
    stop( sprintf( paste0(
      "The following sample(s) in the reloaded experiment file are absent ",
      "from the statistics data.frames:\n  %s\n",
      "All samples in the reloaded file must have been processed by the ",
      "gating pipeline." ),
      paste( missing.from.stats, collapse = ", " ) ) )

  # ---- 2. All reload samples must appear in the original experiment file ---

  missing.from.original <- setdiff( reload$Sample, original$Sample )
  if ( length( missing.from.original ) > 0L )
    stop( sprintf( paste0(
      "The following sample(s) in the reloaded experiment file were not ",
      "present in the original experiment file (ap$experiment.filename):\n  %s\n",
      "Only samples that were part of the original gating run are permitted." ),
      paste( missing.from.original, collapse = ", " ) ) )

  # ---- 3. Structural columns must be unchanged for shared samples ----------

  # Columns that must not change between original and reload
  structural.cols <- c( "filename", "gate.calculation" )

  # Only check columns that actually exist in both files
  cols.to.check <- intersect( structural.cols,
                              intersect( colnames( reload ), colnames( original ) ) )

  shared.samples <- intersect( reload$Sample, original$Sample )

  reload.sub   <- reload[   match( shared.samples, reload$Sample   ), , drop = FALSE ]
  original.sub <- original[ match( shared.samples, original$Sample ), , drop = FALSE ]

  problems <- character( 0L )

  for ( col in cols.to.check ) {

    reload.vals   <- as.character( reload.sub[[ col ]] )
    original.vals <- as.character( original.sub[[ col ]] )

    mismatch.idx <- which( reload.vals != original.vals )

    if ( length( mismatch.idx ) > 0L ) {
      details <- paste( sprintf(
        "  Sample '%s': original = '%s', reloaded = '%s'",
        shared.samples[ mismatch.idx ],
        original.vals[ mismatch.idx ],
        reload.vals[   mismatch.idx ] ), collapse = "\n" )
      problems <- c( problems,
                     sprintf( "Column '%s' differs for %d sample(s):\n%s",
                              col, length( mismatch.idx ), details ) )
    }
  }

  if ( length( problems ) > 0L )
    stop( paste0(
      "Reloaded experiment file has inconsistencies in structural columns ",
      "(only Group, Block, and Covariate.* columns may change):\n\n",
      paste( problems, collapse = "\n\n" ) ) )

  # ---- informational messages ----------------------------------------------

  n.original <- nrow( original )
  n.reload   <- nrow( reload )
  if ( n.reload < n.original )
    message( sprintf(
      "  Using %d of %d sample(s) from original experiment (subset mode).",
      n.reload, n.original ) )

  invisible( NULL )
}


# ------------------------------------------------------------------------------
# .da.resolve.contrasts
# ------------------------------------------------------------------------------

#' Resolve pairwise group contrasts for differential abundance
#'
#' @description
#' Returns a list of \code{character(2)} pairs (reference, treatment) for all
#' pairwise comparisons among the chosen groups.  When there are more than 4
#' groups, the user must supply \code{volcano.groups} explicitly; all pairwise
#' combinations would otherwise produce an unmanageable number of plots.
#'
#' @param n.groups Integer.  Number of unique groups.
#' @param group.levels Character vector of group level names.
#' @param volcano.groups Character vector of selected group levels, or
#'   \code{NULL} (use all levels when \eqn{n \le 4}).
#'
#' @return A list of \code{character(2)} vectors, each giving the reference
#'   and treatment levels for one contrast.
#'
#' @keywords internal
.da.resolve.contrasts <- function( n.groups, group.levels, volcano.groups ) {

  if ( n.groups > 4L ) {
    if ( is.null( volcano.groups ) )
      stop( paste0(
        "There are ", n.groups, " groups, which would produce ",
        choose( n.groups, 2L ), " pairwise volcano plots.\n",
        "Please pass a character vector to the 'volcano.groups' argument ",
        "specifying the 2-4 group levels you want to compare.\n",
        "Available groups: ", paste( group.levels, collapse = ", " ) ) )

    # Validate user-supplied groups
    bad <- setdiff( volcano.groups, group.levels )
    if ( length( bad ) > 0L )
      stop( sprintf( "volcano.groups contains unknown group(s): %s.\nAvailable: %s",
                     paste( bad, collapse = ", " ),
                     paste( group.levels, collapse = ", " ) ) )
    if ( length( volcano.groups ) < 2L )
      stop( "volcano.groups must contain at least 2 group names." )

    active.levels <- volcano.groups
  } else {
    active.levels <- group.levels
  }

  # All pairwise combinations (reference = alphabetically first)
  pairs <- utils::combn( active.levels, 2L, simplify = FALSE )

  pairs
}


# ------------------------------------------------------------------------------
# .da.build.design
# ------------------------------------------------------------------------------

#' Build a limma design matrix from sample metadata
#'
#' @description
#' Constructs a no-intercept group design matrix from the \code{Group} factor
#' in \code{meta}, then appends any \code{Covariate.*} columns as additional
#' regressors (dummy-encoded for character/factor covariates; numeric
#' otherwise).  Detects a \code{Block} column for use with
#' \code{\link[limma]{duplicateCorrelation}} in paired or
#' repeated-measures designs.
#'
#' Covariate column conventions:
#' \describe{
#'   \item{\code{Group}}{Mandatory.  Factor defining comparison groups.
#'     Encoded without an intercept so \code{\link[limma]{makeContrasts}}
#'     works cleanly.}
#'   \item{\code{Block}}{Optional.  Paired or repeated-measures blocking
#'     variable.  Triggers \code{\link[limma]{duplicateCorrelation}} in the
#'     caller.  Ignored if entirely blank.}
#'   \item{\code{Covariate.*}}{Optional.  Any column whose name starts with
#'     \code{"Covariate."} is added as a nuisance regressor.  Character or
#'     factor columns are dummy-encoded (one level dropped to avoid rank
#'     deficiency); numeric columns are used as-is.  Entirely blank columns
#'     are silently skipped with a warning.}
#' }
#'
#' @param meta data.frame.  Sample metadata with at minimum a \code{Group}
#'   factor column and, optionally, \code{Block} and \code{Covariate.*}
#'   columns.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{design}}{Numeric matrix.  The design matrix for
#'       \code{\link[limma]{lmFit}}.}
#'     \item{\code{block}}{Character vector of block labels, or \code{NULL}
#'       if no \code{Block} column is present or it is entirely blank.}
#'   }
#'
#' @keywords internal
.da.build.design <- function( meta ) {

  group <- meta$Group   # already a factor

  # Detect covariate columns (prefix "Covariate.")
  cov.cols <- grep( "^Covariate\\.", colnames( meta ), value = TRUE )

  # Start with a no-intercept group design so Group coefficients are
  # interpretable and makeContrasts works cleanly.
  design <- model.matrix( ~ 0 + Group, data = meta )
  colnames( design ) <- gsub( "^Group", "Group", colnames( design ) )

  # Append covariates - skip any column that is entirely blank
  for ( cc in cov.cols ) {
    cov.vec   <- meta[[ cc ]]
    cov.label <- sub( "^Covariate\\.", "", cc )

    # Treat empty strings and whitespace-only values as missing
    if ( is.character( cov.vec ) || is.factor( cov.vec ) ) {
      cov.chr <- trimws( as.character( cov.vec ) )
      if ( all( cov.chr == "" | is.na( cov.chr ) ) ) {
        warning( sprintf(
          "Covariate column '%s' is entirely blank and will be ignored.", cc ) )
        next
      }
      # Dummy-encode; drop one level to avoid rank deficiency
      cov.factor <- factor( cov.chr )
      cov.dm     <- model.matrix( ~ cov.factor )[ , -1L, drop = FALSE ]
      colnames( cov.dm ) <- paste0( cov.label, "_",
                                    levels( cov.factor )[ -1L ] )
    } else {
      cov.num <- suppressWarnings( as.numeric( cov.vec ) )
      if ( all( is.na( cov.num ) ) ) {
        warning( sprintf(
          "Covariate column '%s' is entirely blank/NA and will be ignored.", cc ) )
        next
      }
      cov.dm <- matrix( cov.num, ncol = 1L )
      colnames( cov.dm ) <- cov.label
    }

    design <- cbind( design, cov.dm )
  }

  # Detect optional Block column for paired / repeated-measures.
  # Treat as absent if the column exists but is entirely blank.
  block <- NULL
  if ( "Block" %in% colnames( meta ) ) {
    block.raw <- trimws( as.character( meta$Block ) )
    if ( all( block.raw == "" | is.na( block.raw ) ) ) {
      warning( "Block column is present but entirely blank and will be ignored." )
    } else {
      block <- block.raw
    }
  }

  list( design = design, block = block )
}


# ------------------------------------------------------------------------------
# .da.plot.heatmap
# ------------------------------------------------------------------------------

#' Save a pheatmap of the top-k differentially abundant features
#'
#' @description
#' Ranks features by their F-statistic across all contrasts, subsets the
#' expression matrix to the top \code{top.k}, optionally z-scores rows, and
#' saves a \code{\link[pheatmap]{pheatmap}} PNG with \code{Group} and any
#' \code{Covariate.*} annotation columns.
#'
#' @param fit2 \code{MArrayLM} object after \code{\link[limma]{eBayes}}.
#' @param expr.mat Numeric matrix, features x samples.
#' @param meta data.frame of sample metadata.
#' @param top.k Integer.  Number of features to display.
#' @param scale.rows Logical.  Z-score rows before clustering.
#' @param figure.dir Character.  Output directory.
#' @param filename Character.  Base filename (no extension).
#' @param width Numeric.  Figure width in inches.
#' @param height Numeric.  Figure height in inches.
#'
#' @return \code{NULL} invisibly.
#'
#' @keywords internal
.da.plot.heatmap <- function( fit2, expr.mat, meta, top.k, scale.rows,
                              figure.dir, filename, width, height ) {

  # Rank features by F-statistic (summary across all contrasts)
  f.stat <- fit2$F
  names( f.stat ) <- rownames( fit2 )
  top.features <- names( sort( f.stat, decreasing = TRUE ) )[ seq_len(
    min( top.k, length( f.stat ) ) ) ]

  # mat is features x samples; columns are sample names
  mat <- expr.mat[ top.features, , drop = FALSE ]

  # ---- build annotation data.frame -----------------------------------------
  # pheatmap matches annotation_col rows to heatmap columns by name.
  # rownames(annotation.col) must equal colnames(mat), in the same order.

  ann.cols <- intersect( c( "Group",
                            grep( "^Covariate\\.", colnames( meta ),
                                  value = TRUE ) ),
                         colnames( meta ) )
  annotation.col <- as.data.frame( meta[ , ann.cols, drop = FALSE ],
                                   stringsAsFactors = FALSE )
  rownames( annotation.col ) <- meta$Sample
  colnames( annotation.col ) <- sub( "^Covariate\\.", "", colnames( annotation.col ) )

  # Align annotation rows to exactly the sample columns present in mat
  mat.samples <- colnames( mat )
  missing.ann <- setdiff( mat.samples, rownames( annotation.col ) )
  if ( length( missing.ann ) > 0L )
    stop( sprintf(
      "Internal error: sample(s) in expr.mat have no annotation row: %s",
      paste( missing.ann, collapse = ", " ) ) )
  annotation.col <- annotation.col[ mat.samples, , drop = FALSE ]

  # Drop any annotation column that is entirely NA or blank - these arise when
  # the experiment CSV has a Covariate.* or Block column that was never filled
  # in.  Passing all-NA columns to pheatmap causes errors in its colour-scale
  # code.
  all.blank.ann <- vapply( annotation.col, function( col ) {
    all( is.na( col ) | trimws( as.character( col ) ) == "" )
  }, logical( 1L ) )
  if ( any( all.blank.ann ) ) {
    warning( sprintf(
      "Dropping entirely blank annotation column(s) from heatmap: %s",
      paste( colnames( annotation.col )[ all.blank.ann ], collapse = ", " ) ) )
    annotation.col <- annotation.col[ , !all.blank.ann, drop = FALSE ]
  }

  # If no annotation columns remain, pass NULL so pheatmap omits the sidebar
  if ( ncol( annotation.col ) == 0L )
    annotation.col <- NULL

  # Coerce Group to character so pheatmap treats it as categorical, not numeric
  if ( !is.null( annotation.col ) && "Group" %in% colnames( annotation.col ) )
    annotation.col$Group <- as.character( annotation.col$Group )

  # ---- z-score rows --------------------------------------------------------
  if ( scale.rows ) {
    row.sd <- apply( mat, 1L, stats::sd, na.rm = TRUE )
    row.sd[ row.sd == 0 | is.na( row.sd ) ] <- 1
    mat <- sweep( mat - rowMeans( mat, na.rm = TRUE ), 1L, row.sd, "/" )
  }

  # Ensure no non-finite values survive into pheatmap (would break dist())
  mat[ !is.finite( mat ) ] <- 0

  # ---- draw ----------------------------------------------------------------
  heat.cols <- grDevices::colorRampPalette(
    c( "#2166AC", "#F7F7F7", "#D6604D" ) )( 100L )

  out.path <- file.path( figure.dir, paste0( filename, ".png" ) )

  grDevices::png( out.path, width = width, height = height,
                  units = "in", res = 150L )
  on.exit( grDevices::dev.off(), add = TRUE )   # always close device

  pheatmap::pheatmap(
    mat,
    color              = heat.cols,
    annotation_col     = annotation.col,
    cluster_rows       = TRUE,
    cluster_cols       = TRUE,
    show_colnames      = TRUE,
    show_rownames      = TRUE,
    fontsize_row       = max( 5L, 10L - nrow( mat ) %/% 10L ),
    main               = sprintf( "Top %d differentially abundant features",
                                  nrow( mat ) ),
    border_color       = NA
  )

  message( sprintf( "  Heatmap saved: %s", out.path ) )
}


# ------------------------------------------------------------------------------
# .da.plot.volcano
# ------------------------------------------------------------------------------

#' Save a volcano plot for one pairwise contrast
#'
#' @description
#' Draws a \code{ggplot2} volcano plot (log\eqn{_2} fold change vs
#' \eqn{-\log_{10}} adjusted p-value) for a single pairwise contrast, labels
#' the top 15 significant hits with \code{\link[ggrepel]{geom_text_repel}},
#' and saves it as a PNG.
#'
#' @param toptable data.frame returned by \code{\link[limma]{topTable}}.
#' @param contrast.name Character scalar.  Name of the contrast
#'   (e.g. \code{"A_vs_B"}), used in the title and output filename.
#' @param fdr.threshold Numeric.  Adjusted p-value threshold.
#' @param lfc.threshold Numeric.  Absolute log\eqn{_2} fold-change threshold.
#' @param figure.dir Character.  Output directory.
#' @param filename.prefix Character.  Filename prefix; contrast name appended.
#' @param width Numeric.  Figure width in inches.
#' @param height Numeric.  Figure height in inches.
#'
#' @return \code{NULL} invisibly.
#'
#' @keywords internal
.da.plot.volcano <- function( toptable, contrast.name, fdr.threshold,
                              lfc.threshold, figure.dir, filename.prefix,
                              width, height ) {

  # topTable returns log2FC in "logFC", adj.P.Val
  df <- toptable
  df$Feature  <- rownames( df )
  df$neglog10p <- -log10( df$adj.P.Val + .Machine$double.eps )

  # Significance classification
  df$Significance <- "NS"
  df$Significance[ abs( df$logFC ) >= lfc.threshold &
                     df$adj.P.Val   <  fdr.threshold  ] <- "Significant"
  df$Significance[ abs( df$logFC ) >= lfc.threshold &
                     df$logFC > 0   &
                     df$adj.P.Val   <  fdr.threshold  ] <- "Up"
  df$Significance[ abs( df$logFC ) >= lfc.threshold &
                     df$logFC < 0   &
                     df$adj.P.Val   <  fdr.threshold  ] <- "Down"
  df$Significance <- factor( df$Significance,
                             levels = c( "NS", "Up", "Down", "Significant" ) )

  # Label the most significant hits (top 15 by adjusted p-value)
  n.label <- min( 15L, sum( df$Significance != "NS" ) )
  label.df <- df[ order( df$adj.P.Val )[ seq_len( n.label ) ], ]

  # Build x-axis limits symmetrically
  lfc.max <- max( abs( df$logFC ), na.rm = TRUE )
  x.lim   <- c( -lfc.max, lfc.max ) * 1.05

  col.map <- c(
    NS           = "grey70",
    Up           = "#D6604D",
    Down         = "#4393C3",
    Significant  = "#74C476"
  )

  # Contrast label for title: "A_vs_B" -> "A vs B"
  title.str <- sprintf( "Volcano: %s\n(FDR < %.2f, |log2FC| >= %.2f)",
                        gsub( "_", " ", contrast.name ),
                        fdr.threshold, lfc.threshold )

  ggp <- ggplot2::ggplot(
    df,
    ggplot2::aes( x = logFC, y = neglog10p,
                  colour = Significance ) ) +
    ggplot2::geom_point( alpha = 0.6, size = 1.8 ) +
    ggplot2::scale_colour_manual( values = col.map, drop = FALSE ) +
    ggrepel::geom_text_repel(
      data            = label.df,
      ggplot2::aes( label = Feature ),
      size            = 3,
      max.overlaps    = 20,
      segment.colour  = "grey50",
      show.legend     = FALSE
    ) +
    ggplot2::geom_vline( xintercept = c( -lfc.threshold, lfc.threshold ),
                         linetype = "dashed", colour = "grey40" ) +
    ggplot2::geom_hline( yintercept = -log10( fdr.threshold ),
                         linetype = "dashed", colour = "grey40" ) +
    ggplot2::xlim( x.lim ) +
    ggplot2::labs(
      title  = title.str,
      x      = expression( log[2]~"fold change" ),
      y      = expression( -log[10]~"(adjusted p-value)" ),
      colour = NULL
    ) +
    ggplot2::theme_bw( base_size = 12 ) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text( size = 11, face = "bold" ),
      legend.position = "bottom"
    )

  # Safe filename: replace spaces and special chars
  safe.contrast <- gsub( "[^A-Za-z0-9_.-]", "_", contrast.name )
  out.path <- file.path( figure.dir,
                         paste0( filename.prefix, "_", safe.contrast, ".png" ) )

  ggplot2::ggsave( out.path, ggp,
                   width = width, height = height, dpi = 150L )

  message( sprintf( "  Volcano saved:  %s", out.path ) )
}
