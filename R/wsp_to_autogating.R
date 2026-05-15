# wsp_to_autogating.R
#
# Convert a FlowJo .wsp workspace file to the population gate and activation
# gate definition lists used by the AutoGating pipeline.
#
# Public API:
#   wsp.to.autogating()       — main entry point
#   get.cytometer.zero.params() — cytometer scatter / biexp parameters
#   write.gate.definitions()  — write gate lists to R source files
#
# Internal helpers:
#   parse.wsp.xml()              — extract gate tree from WSP XML
#   extract.gates.recursive()    — recursive XML walker
#   extract.gate.definition()    — parse a single Gate XML element
#   resolve.dye.to.antigen()     — map dye name to antigen name
#   classify.gates()             — assign role, gate.type, algorithm per gate
#   clean.gate.name()            — sanitise a gate name for use as a list key
#   build.population.gates()     — assemble fcs.population.gates list
#   build.activation.gates()     — assemble fcs.activation.gates list
#   dirname_gate() / basename_gate() / get.sibling.gates() — path utilities
#   write.pop.file() / write.act.file() — R-source serialisers


# =============================================================================
# Package-level constants
# =============================================================================

.wsp.ns <- c(
  gating     = "http://www.isac-net.org/std/Gating-ML/v2.0/gating",
  datatypes  = "http://www.isac-net.org/std/Gating-ML/v2.0/datatypes",
  transforms = "http://www.isac-net.org/std/Gating-ML/v2.0/transformations"
)

# Known viability marker names (case-insensitive matching)
.viability.patterns <- c(
  "viability", "livedead", "live.dead", "live/dead", "zombie", "viakrome",
  "fixable", "aqua", "7-aad", "7aad", "dapi", "sytox", "propidium",
  "live_dead", "l/d", "ld"
)

# Known scatter channel prefixes
.scatter.prefixes <- c("FSC", "SSC")


# =============================================================================
# wsp.to.autogating — main entry point
# =============================================================================

#' Convert a FlowJo workspace to AutoGating gate definition lists
#'
#' @title Convert FlowJo Workspace to AutoGating Format
#'
#' @description
#' Reads a FlowJo \code{.wsp} (or \code{.wspt}) workspace file, extracts the
#' gate hierarchy, maps FlowJo channel names to antigen names using the panel
#' CSV, classifies each gate as a population or activation gate, and returns
#' the \code{fcs.population.gates} and \code{fcs.activation.gates} lists
#' compatible with the AutoGating pipeline.
#'
#' The function proceeds in five stages:
#' \enumerate{
#'   \item Parse the WSP XML into a flat list of gate records
#'     (\code{parse.wsp.xml}).
#'   \item Map dye channel names to antigen names via the panel CSV
#'     (\code{resolve.dye.to.antigen}).
#'   \item Classify each gate as singlet, leukocyte/scatter, viability,
#'     population, activation, or NOT gate across six classification passes
#'     (\code{classify.gates}).
#'   \item Build the population gate list (\code{build.population.gates}).
#'   \item Build the activation gate list (\code{build.activation.gates}).
#' }
#'
#' Sample-level gates are preferred over group-level gates when both exist;
#' fall-back to group-level is automatic.
#'
#' @param wsp.filename Character scalar.  Path to the \code{.wsp} or
#'   \code{.wspt} FlowJo workspace file.
#' @param panel.filename Character scalar.  Path to the panel CSV with columns
#'   \code{dye}, \code{antigen}, and \code{type}.
#' @param cytometer Character scalar.  Instrument identifier used to look up
#'   scatter channel names and the biexponential near-zero point.  One of
#'   \code{"aurora"}, \code{"id7000"}, \code{"symphony"}, or any string that
#'   falls through to the built-in fallback defaults.
#'   Default \code{"aurora"}.
#' @param group.name Character scalar.  Name of the FlowJo group to read
#'   gates from when sample-level gates are unavailable.
#'   Default \code{"All Samples"}.
#' @param sample.index Integer scalar.  Index of the sample within the
#'   workspace to read gates from (1-based).  Set to \code{NULL} to force
#'   group-level gate reading.  Default \code{1L}.
#' @param prefer.sample Logical.  If \code{TRUE} (default), sample-level
#'   gates are used when available; if \code{FALSE}, group-level gates are
#'   always used.
#' @param two.marker.strategy Character scalar.  Strategy for 2-D lineage
#'   gates on two fluorescence markers: \code{"2dsep"} (default) creates a
#'   single \code{2dsep} gate; \code{"sequential"} converts it to a
#'   \code{1dsep} gate on each axis separately.
#' @param verbose Logical.  Print diagnostic progress messages.
#'   Default \code{TRUE}.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{population.gates}}{List of population gate definitions
#'       compatible with \code{\link{do.population.gates}}.}
#'     \item{\code{activation.gates}}{List of activation gate definitions
#'       compatible with \code{\link{do.activation.gates}}.}
#'     \item{\code{gate.tree}}{The flat list of classified gate records for
#'       inspection and debugging.}
#'     \item{\code{panel}}{The panel data frame as read and normalised from
#'       \code{panel.filename}.}
#'   }
#'
#' @seealso \code{\link{write.gate.definitions}},
#'   \code{\link{do.population.gates}}, \code{\link{do.activation.gates}},
#'   \code{\link{get.autogating.param}}
#'
#' @importFrom utils read.csv
#' @importFrom xml2 read_xml xml_find_all xml_find_first xml_children
#'   xml_name xml_attr
#'
#' @export
wsp.to.autogating <- function(
    wsp.filename,
    panel.filename,
    cytometer             = "aurora",
    group.name            = "All Samples",
    sample.index          = 1L,
    prefer.sample         = TRUE,
    two.marker.strategy   = c( "2dsep", "sequential" ),
    verbose               = TRUE
) {
  two.marker.strategy <- match.arg(two.marker.strategy)

  # ---- read panel -----------------------------------------------------------
  panel <- utils::read.csv(panel.filename, stringsAsFactors = FALSE)
  panel$dye     <- trimws(panel$dye)
  panel$antigen <- trimws(panel$antigen)
  panel$type    <- trimws(tolower(panel$type))

  # Build dye -> antigen lookup
  dye.to.antigen <- stats::setNames(panel$antigen, panel$dye)

  # Activation markers from panel
  activation.antigens <- panel$antigen[panel$type == "activation"]
  lineage.antigens    <- panel$antigen[panel$type == "lineage"]
  scatter.dyes        <- panel$dye[panel$type == "scatter"]

  if (verbose)
    message(sprintf("Panel: %d channels (%d lineage, %d activation, %d scatter)",
                    nrow(panel), length(lineage.antigens),
                    length(activation.antigens), length(scatter.dyes)))

  # ---- cytometer-specific parameters ----------------------------------------
  cyto.params <- get.cytometer.zero.params(cytometer)

  # ---- parse WSP XML --------------------------------------------------------
  gate.tree <- parse.wsp.xml(
    wsp.filename,
    group.name    = group.name,
    sample.index  = sample.index,
    prefer.sample = prefer.sample,
    verbose       = verbose
  )

  if (verbose)
    message(sprintf("Parsed %d gates from workspace", length(gate.tree)))

  # ---- map dye names to antigen names ---------------------------------------
  gate.tree <- lapply(gate.tree, function(g) {
    g$antigen.x <- resolve.dye.to.antigen(g$dim.x, dye.to.antigen, scatter.dyes)
    g$antigen.y <- resolve.dye.to.antigen(g$dim.y, dye.to.antigen, scatter.dyes)
    g
  })

  # ---- classify gates -------------------------------------------------------
  gate.tree <- classify.gates(
    gate.tree,
    activation.antigens = activation.antigens,
    scatter.dyes        = scatter.dyes,
    dye.to.antigen      = dye.to.antigen,
    cyto.params         = cyto.params,
    two.marker.strategy = two.marker.strategy,
    verbose             = verbose
  )

  # ---- build output lists ---------------------------------------------------
  pop.gates <- build.population.gates(gate.tree, cyto.params, verbose)
  act.gates <- build.activation.gates(gate.tree, pop.gates, verbose)

  list(
    population.gates = pop.gates,
    activation.gates = act.gates,
    gate.tree        = gate.tree,
    panel            = panel
  )
}


# =============================================================================
# get.cytometer.zero.params
# =============================================================================

#' Return cytometer-specific scatter and biexponential parameters
#'
#' @title Cytometer-Specific Parameters for Gate Classification
#'
#' @description
#' Returns a named list of instrument-specific constants used during gate
#' classification: the approximate biexponential near-zero raw value (used to
#' determine whether a gate spans the negative region), scatter channel names
#' split by area/height/width, and the default FSC/SSC pair for the root gate.
#'
#' Known cytometers: \code{"aurora"}, \code{"id7000"}, \code{"symphony"}.
#' Unrecognised strings fall back to Symphony-like defaults.
#'
#' @param cytometer Character scalar.  Cytometer identifier (case-insensitive).
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{\code{biexp.zero}}{Numeric.  Raw value corresponding to biexp
#'       zero (always \code{0} for current instruments).}
#'     \item{\code{biexp.near.zero}}{Numeric.  Raw value below which a gate
#'       boundary is considered to span the negative region.}
#'     \item{\code{scatter.max}}{Numeric.  Maximum raw scatter value for the
#'       instrument.}
#'     \item{\code{scatter.area}}{Character vector of area scatter channels.}
#'     \item{\code{scatter.height}}{Character vector of height scatter
#'       channels.}
#'     \item{\code{scatter.width}}{Character vector of width scatter channels.}
#'     \item{\code{default.scatter}}{Length-2 character vector giving the
#'       default x and y scatter channels for the root gate.}
#'     \item{\code{all.scatter}}{Character vector of all scatter channels
#'       (union of area, height, width).}
#'   }
#'
#' @seealso \code{\link{wsp.to.autogating}}
#'
#' @export
get.cytometer.zero.params <- function(cytometer) {

  cytometer <- tolower(cytometer)

  params <- switch(cytometer,
    aurora = list(
      biexp.zero       = 0,
      biexp.near.zero  = 5000,
      scatter.max      = 4194304,
      scatter.area     = c("FSC-A", "SSC-A", "SSC-B-A"),
      scatter.height   = c("FSC-H", "SSC-H", "SSC-B-H"),
      scatter.width    = c("FSC-W", "SSC-W", "SSC-B-W"),
      default.scatter  = c("FSC-A", "SSC-B-A")
    ),
    id7000 = list(
      biexp.zero       = 0,
      biexp.near.zero  = 2000,
      scatter.max      = 1048576,
      scatter.area     = c("FSC-A", "SSC-A"),
      scatter.height   = c("FSC-H", "SSC-H"),
      scatter.width    = c("FSC-W", "SSC-W"),
      default.scatter  = c("FSC-A", "SSC-A")
    ),
    symphony = list(
      biexp.zero       = 0,
      biexp.near.zero  = 1500,
      scatter.max      = 262144,
      scatter.area     = c("FSC-A", "SSC-A"),
      scatter.height   = c("FSC-H", "SSC-H"),
      scatter.width    = c("FSC-W", "SSC-W"),
      default.scatter  = c("FSC-A", "SSC-A")
    ),
    # Default fallback
    list(
      biexp.zero       = 0,
      biexp.near.zero  = 3000,
      scatter.max      = 262144,
      scatter.area     = c("FSC-A", "SSC-A"),
      scatter.height   = c("FSC-H", "SSC-H"),
      scatter.width    = c("FSC-W", "SSC-W"),
      default.scatter  = c("FSC-A", "SSC-A")
    )
  )

  params$all.scatter <- c(params$scatter.area, params$scatter.height,
                          params$scatter.width)
  params
}


# =============================================================================
# parse.wsp.xml
# =============================================================================

#' Parse a FlowJo workspace XML and return a flat gate tree
#'
#' @title Parse FlowJo Workspace XML
#'
#' @description
#' Reads a FlowJo \code{.wsp} file, attempts to extract gates from the
#' specified sample (sample-level gates are more detailed), and falls back to
#' group-level gates if no sample gates are found.  Returns a flat list of
#' gate records produced by \code{extract.gates.recursive}.
#'
#' @param wsp.filename Character scalar.  Path to the workspace file.
#' @param group.name Character scalar.  Name of the FlowJo group for fallback
#'   extraction.  Default \code{"All Samples"}.
#' @param sample.index Integer or \code{NULL}.  Index of the sample to read
#'   from (1-based).  \code{NULL} forces group-level extraction.
#' @param prefer.sample Logical.  Whether to prefer sample-level gates.
#' @param verbose Logical.  Print extraction source messages.
#'
#' @return A flat list of gate records.  Each record is a named list with
#'   fields: \code{name}, \code{path}, \code{parent.path}, \code{depth},
#'   \code{gate.class}, \code{ndims}, \code{dim.x}, \code{dim.y},
#'   \code{vertices}, \code{min.x}, \code{max.x}, \code{min.y},
#'   \code{max.y}, \code{is.not.gate}, \code{not.ref.id},
#'   \code{gate.id}, \code{events.inside}, \code{children}.
#'
#' @importFrom xml2 read_xml xml_find_all xml_find_first xml_children
#'   xml_name xml_attr
#'
#' @keywords internal
parse.wsp.xml <- function(
    wsp.filename,
    group.name    = "All Samples",
    sample.index  = 1L,
    prefer.sample = TRUE,
    verbose       = TRUE
) {
  doc <- xml2::read_xml(wsp.filename)

  gate.tree <- list()
  source.label <- ""

  # Helper: check if xml node is valid (not NULL and not xml_missing)
  is.valid.node <- function(x) {
    !is.null(x) && !inherits(x, "xml_missing") && length(x) > 0
  }

  # ---- try sample-level gates -----------------------------------------------
  if (prefer.sample && !is.null(sample.index)) {
    samples <- xml2::xml_find_all(doc, "//SampleList/Sample")
    if (length(samples) >= sample.index) {
      sn <- xml2::xml_find_first(samples[[sample.index]], ".//SampleNode")
      if (is.valid.node(sn)) {
        sn.name <- xml2::xml_attr(sn, "name")
        subpops <- xml2::xml_find_first(sn, ".//Subpopulations")
        if (is.valid.node(subpops) && length(xml2::xml_children(subpops)) > 0) {
          source.label <- sprintf("sample '%s'", sn.name)
          gate.tree <- extract.gates.recursive(subpops, parent.path = "", depth = 0)
        }
      }
    }
  }

  # ---- fall back to group-level gates ---------------------------------------
  if (length(gate.tree) == 0) {
    group.nodes <- xml2::xml_find_all(doc, "//Groups/GroupNode")
    for (gn in group.nodes) {
      gn.name <- xml2::xml_attr(gn, "name")
      if (gn.name == group.name) {
        subpops <- xml2::xml_find_first(gn, ".//Subpopulations")
        if (is.valid.node(subpops) && length(xml2::xml_children(subpops)) > 0) {
          source.label <- sprintf("group '%s'", gn.name)
          gate.tree <- extract.gates.recursive(subpops, parent.path = "", depth = 0)
        }
        break
      }
    }
  }

  if (length(gate.tree) == 0)
    warning("No gates found in workspace. Check group.name and sample.index.")

  if (verbose && nchar(source.label) > 0)
    message(sprintf("Gates extracted from %s", source.label))

  gate.tree
}


# =============================================================================
# extract.gates.recursive
# =============================================================================

#' Recursively walk an XML Subpopulations node and collect gate records
#'
#' @description
#' Iterates over direct children of a \code{<Subpopulations>} XML node,
#' processes each \code{<Population>} or \code{<NotNode>} element by calling
#' \code{extract.gate.definition}, then recurses into its own
#' \code{<Subpopulations>} child.  Returns a flat list of gate records
#' ordered depth-first (parent before children).
#'
#' @param subpops.node An \code{xml_node} pointing to a
#'   \code{<Subpopulations>} element.
#' @param parent.path Character scalar.  Slash-separated path of the parent
#'   gate (empty string for root-level gates).
#' @param depth Integer.  Depth of the current node in the tree (0 = root).
#'
#' @return A flat list of gate records.
#'
#' @keywords internal
extract.gates.recursive <- function(subpops.node, parent.path, depth) {
  gates <- list()

  children <- xml2::xml_children(subpops.node)
  for (child in children) {
    tag <- xml2::xml_name(child)
    if (!(tag %in% c("Population", "NotNode")))
      next

    gate.name <- xml2::xml_attr(child, "name")
    is.not    <- (tag == "NotNode")

    # Build path
    gate.path <- if (nchar(parent.path) == 0) gate.name
                 else paste(parent.path, gate.name, sep = "/")

    # ---- extract gate definition ------------------------------------------
    gate.info <- extract.gate.definition(child, is.not)
    gate.info$name        <- gate.name
    gate.info$path        <- gate.path
    gate.info$parent.path <- parent.path
    gate.info$depth       <- depth
    gate.info$is.not.gate <- is.not

    # ---- recurse into children --------------------------------------------
    child.subpops <- xml2::xml_find_first(child, ".//Subpopulations[1]")
    child.paths   <- character(0)

    if (!inherits(child.subpops, "xml_missing") && !is.null(child.subpops) &&
        length(xml2::xml_children(child.subpops)) > 0) {
      child.gates <- extract.gates.recursive(child.subpops, gate.path, depth + 1)
      child.paths <- vapply(child.gates, `[[`, character(1), "path")
      gates <- c(gates, child.gates)
    }

    gate.info$children <- child.paths
    gates <- c(list(gate.info), gates)
  }

  gates
}


# =============================================================================
# extract.gate.definition
# =============================================================================

#' Parse the gate geometry from a single FlowJo Population node
#'
#' @description
#' Reads a \code{<Gate>} child of a \code{<Population>} or \code{<NotNode>}
#' XML element and extracts the gate class (Polygon, Rectangle, Ellipsoid, or
#' NotGate), the channel dimension names, bounding-box coordinates, and vertex
#' matrix.  Boolean NOT gates are detected from both the \code{<NotNode>} tag
#' and the \code{<gating:BooleanGate>} element.
#'
#' @param pop.node An \code{xml_node} for the \code{<Population>} or
#'   \code{<NotNode>} element.
#' @param is.not Logical.  Whether this node was encountered as a
#'   \code{<NotNode>} tag.  Default \code{FALSE}.
#'
#' @return A named list with fields: \code{gate.class}, \code{ndims},
#'   \code{dim.x}, \code{dim.y}, \code{vertices}, \code{min.x},
#'   \code{max.x}, \code{min.y}, \code{max.y}, \code{not.ref.id},
#'   \code{gate.id}, \code{events.inside}.
#'
#' @keywords internal
extract.gate.definition <- function(pop.node, is.not = FALSE) {

  info <- list(
    gate.class    = NA_character_,
    ndims         = 0L,
    dim.x         = NA_character_,
    dim.y         = NA_character_,
    vertices      = NULL,
    min.x         = NA_real_,
    max.x         = NA_real_,
    min.y         = NA_real_,
    max.y         = NA_real_,
    not.ref.id    = NA_character_,
    gate.id       = NA_character_,
    events.inside = TRUE
  )

  gate.node <- xml2::xml_find_first(pop.node, ".//Gate[1]")

  if (inherits(gate.node, "xml_missing") ||
      is.null(gate.node) || length(gate.node) == 0) {
    # NOT gates (NotNode) don't have a Gate element with geometry
    if (is.not) {
      info$gate.class <- "NotGate"
    }
    return(info)
  }

  info$gate.id <- tryCatch({
    id <- xml2::xml_attr(gate.node, paste0("{", .wsp.ns[["gating"]], "}id"))
    if (is.na(id))
      id <- xml2::xml_attr(gate.node, "id")
    id
  }, error = function(e) NA_character_)

  # ---- look for Polygon, Rectangle, Ellipsoid gates ----------------------
  gate.types <- c("PolygonGate", "RectangleGate", "EllipsoidGate")
  for (gt in gate.types) {
    xpath <- sprintf(".//gating:%s", gt)
    g <- tryCatch(
      xml2::xml_find_first(gate.node, xpath, ns = .wsp.ns),
      error = function(e) NULL
    )
    if (is.null(g) || inherits(g, "xml_missing")) next

    info$gate.class <- gt
    info$events.inside <- as.logical(
      ifelse(is.na(xml2::xml_attr(g, "eventsInside")), TRUE,
             xml2::xml_attr(g, "eventsInside") != "0"))

    dims <- xml2::xml_find_all(g, ".//gating:dimension", ns = .wsp.ns)
    info$ndims <- length(dims)

    for (i in seq_along(dims)) {
      fcs.dim <- xml2::xml_find_first(dims[[i]], ".//datatypes:fcs-dimension",
                                      ns = .wsp.ns)
      dname <- NA_character_
      if (!inherits(fcs.dim, "xml_missing") && length(fcs.dim) > 0) {
        dname <- xml2::xml_attr(fcs.dim,
                                paste0("{", .wsp.ns[["datatypes"]], "}name"))
        if (is.na(dname))
          dname <- xml2::xml_attr(fcs.dim, "name")
      }

      gmin.str <- xml2::xml_attr(dims[[i]],
                                 paste0("{", .wsp.ns[["gating"]], "}min"))
      gmax.str <- xml2::xml_attr(dims[[i]],
                                 paste0("{", .wsp.ns[["gating"]], "}max"))
      if (is.na(gmin.str)) gmin.str <- xml2::xml_attr(dims[[i]], "min")
      if (is.na(gmax.str)) gmax.str <- xml2::xml_attr(dims[[i]], "max")
      gmin <- suppressWarnings(as.numeric(gmin.str))
      gmax <- suppressWarnings(as.numeric(gmax.str))

      if (i == 1) {
        info$dim.x <- dname; info$min.x <- gmin; info$max.x <- gmax
      } else if (i == 2) {
        info$dim.y <- dname; info$min.y <- gmin; info$max.y <- gmax
      }
    }

    # Vertices (Polygon gates)
    verts <- xml2::xml_find_all(g, ".//gating:vertex", ns = .wsp.ns)
    if (length(verts) > 0) {
      coords <- t(sapply(verts, function(v) {
        cc <- xml2::xml_find_all(v, ".//gating:coordinate", ns = .wsp.ns)
        vals <- sapply(cc, function(c) {
          val <- xml2::xml_attr(c, paste0("{", .wsp.ns[["datatypes"]], "}value"))
          if (is.na(val)) val <- xml2::xml_attr(c, "value")
          as.numeric(val)
        })
        if (length(vals) < 2) vals <- c(vals, NA_real_)
        vals[1:2]
      }))
      colnames(coords) <- c("x", "y")
      info$vertices <- coords
      if (is.na(info$min.x)) info$min.x <- min(coords[, "x"], na.rm = TRUE)
      if (is.na(info$max.x)) info$max.x <- max(coords[, "x"], na.rm = TRUE)
      if (is.na(info$min.y)) info$min.y <- min(coords[, "y"], na.rm = TRUE)
      if (is.na(info$max.y)) info$max.y <- max(coords[, "y"], na.rm = TRUE)
    }

    break
  }

  # ---- Boolean gate (NOT) ------------------------------------------------
  bg <- tryCatch(
    xml2::xml_find_first(gate.node, ".//gating:BooleanGate", ns = .wsp.ns),
    error = function(e) NULL
  )
  if (!is.null(bg) && !inherits(bg, "xml_missing")) {
    not.elem <- tryCatch(
      xml2::xml_find_first(bg, ".//gating:not", ns = .wsp.ns),
      error = function(e) NULL
    )
    if (!is.null(not.elem) && !inherits(not.elem, "xml_missing")) {
      info$gate.class <- "NotGate"
      refs <- xml2::xml_find_all(not.elem, ".//gating:gateReference",
                                 ns = .wsp.ns)
      if (length(refs) > 0) {
        ref.id <- xml2::xml_attr(refs[[1]],
                                 paste0("{", .wsp.ns[["gating"]], "}ref"))
        if (is.na(ref.id)) ref.id <- xml2::xml_attr(refs[[1]], "ref")
        info$not.ref.id <- ref.id
      }
    }
  }

  info
}


# =============================================================================
# resolve.dye.to.antigen
# =============================================================================

#' Map a FlowJo channel (dye) name to its antigen name
#'
#' @description
#' Looks up a dye name in the panel-derived \code{dye.to.antigen} map.
#' Scatter channels (FSC/SSC) are returned unchanged.  If the dye is not
#' found in the panel, the original dye name is returned as a fallback with
#' no error.
#'
#' @param dye.name Character scalar.  The channel name from the FlowJo XML.
#' @param dye.to.antigen Named character vector of antigen names keyed by dye.
#' @param scatter.dyes Character vector of scatter dye names from the panel.
#'
#' @return Character scalar.  The corresponding antigen name, the scatter
#'   channel name unchanged, or \code{NA_character_} when \code{dye.name} is
#'   \code{NA}.
#'
#' @keywords internal
resolve.dye.to.antigen <- function(dye.name, dye.to.antigen, scatter.dyes) {
  if (is.na(dye.name)) return(NA_character_)

  # Scatter channels map to themselves
  if (dye.name %in% scatter.dyes || grepl("^(FSC|SSC)", dye.name))
    return(dye.name)

  # Look up in panel
  if (dye.name %in% names(dye.to.antigen))
    return(dye.to.antigen[[dye.name]])

  # If not found, return as-is
  dye.name
}


# =============================================================================
# classify.gates
# =============================================================================

#' Classify gates from the parsed WSP gate tree
#'
#' @description
#' Annotates every gate record in \code{gate.tree} with classification fields
#' (\code{role}, \code{gate.type}, \code{algorithm}, \code{is.singlet},
#' \code{is.scatter}, \code{is.viability}, \code{is.activation},
#' \code{is.replicate}, \code{polarity}) by applying a six-pass rule set:
#'
#' \describe{
#'   \item{Pass 1}{Assigns a default role and tests each gate against rules for
#'     singlet (area vs height), leukocyte/scatter (two different scatter
#'     channels), viability (name/marker matching), activation (1-D gates,
#'     activation-panel markers, name pattern), and NOT gates.}
#'   \item{Pass 2}{Re-classifies activation gates that have non-activation
#'     children as population gates.}
#'   \item{Pass 3}{Detects \code{2dsep} gates: multiple gates at the same
#'     hierarchy level on the same marker pair.}
#'   \item{Pass 4}{Resolves 2-D lineage gates by checking whether each axis
#'     spans the biexponential zero (Rule 19/20).}
#'   \item{Pass 5}{Detects replicate population gates (identical boundary
#'     coordinates on the same parameters).}
#'   \item{Pass 6}{Detects replicate activation gates.}
#' }
#'
#' @param gate.tree List of gate records from \code{parse.wsp.xml}.
#' @param activation.antigens Character vector of activation marker names from
#'   the panel.
#' @param scatter.dyes Character vector of scatter dye names from the panel.
#' @param dye.to.antigen Named character vector mapping dye to antigen.
#' @param cyto.params List from \code{get.cytometer.zero.params}.
#' @param two.marker.strategy Character scalar.  \code{"2dsep"} or
#'   \code{"sequential"}.
#' @param verbose Logical.  Print classification summary.
#'
#' @return The annotated \code{gate.tree} list.
#'
#' @keywords internal
classify.gates <- function(
    gate.tree,
    activation.antigens,
    scatter.dyes,
    dye.to.antigen,
    cyto.params,
    two.marker.strategy = "2dsep",
    verbose = TRUE
) {
  all.scatter <- cyto.params$all.scatter

  # ---- index gates by path for quick lookup ---------------------------------
  path.idx <- stats::setNames(seq_along(gate.tree),
                        vapply(gate.tree, `[[`, character(1), "path"))

  # ---- pass 1: basic classification ----------------------------------------
  for (i in seq_along(gate.tree)) {
    g <- gate.tree[[i]]

    g$role             <- "population"
    g$gate.type        <- "1dsep"
    g$algorithm        <- "mixture"
    g$is.singlet       <- FALSE
    g$is.scatter       <- FALSE
    g$is.viability     <- FALSE
    g$is.activation    <- FALSE
    g$is.replicate     <- FALSE
    g$replicate.origin <- NA_character_
    g$polarity         <- "positive"

    dim.x <- g$dim.x
    dim.y <- g$dim.y
    ant.x <- g$antigen.x
    ant.y <- g$antigen.y

    is.scatter.x <- !is.na(dim.x) && dim.x %in% all.scatter
    is.scatter.y <- !is.na(dim.y) && dim.y %in% all.scatter

    # Rule 10: singlet gate
    if (is.scatter.x && is.scatter.y) {
      base.x   <- sub("-(A|H|W)$", "", dim.x)
      base.y   <- sub("-(A|H|W)$", "", dim.y)
      suffix.x <- sub(".*-", "", dim.x)
      suffix.y <- sub(".*-", "", dim.y)
      if (base.x == base.y && suffix.x != suffix.y &&
          any(c(suffix.x, suffix.y) %in% c("A", "H"))) {
        g$is.singlet <- TRUE
        g$role       <- "population"
        g$gate.type  <- "free"
        g$algorithm  <- "singlets"
      }
    }

    # Rule 11: leukocyte/scatter gate
    if (is.scatter.x && is.scatter.y && !g$is.singlet) {
      base.x <- sub("-(A|H|W)$", "", dim.x)
      base.y <- sub("-(A|H|W)$", "", dim.y)
      if (base.x != base.y) {
        g$is.scatter <- TRUE
        g$role       <- "population"
        g$gate.type  <- "1dsep"
        g$algorithm  <- "mixture"
      }
    }

    # Rule 12: viability gate
    if (!is.na(ant.x) || !is.na(ant.y)) {
      viab.match.x    <- !is.na(ant.x) && any(sapply(.viability.patterns, function(p) grepl(p, tolower(ant.x))))
      viab.match.y    <- !is.na(ant.y) && any(sapply(.viability.patterns, function(p) grepl(p, tolower(ant.y))))
      viab.match.dimx <- !is.na(dim.x) && any(sapply(.viability.patterns, function(p) grepl(p, tolower(dim.x))))
      viab.match.dimy <- !is.na(dim.y) && any(sapply(.viability.patterns, function(p) grepl(p, tolower(dim.y))))
      if (viab.match.x || viab.match.y || viab.match.dimx || viab.match.dimy) {
        g$is.viability <- TRUE
        g$role         <- "population"
        g$gate.type    <- "1dsep"
        g$algorithm    <- "mixture"
      }
    }

    # Rule 5: activation gate recognition
    is.activation <- FALSE
    if (g$ndims == 1) is.activation <- TRUE
    if (g$ndims == 1) {
      if (!is.na(ant.x) && ant.x %in% activation.antigens) is.activation <- TRUE
      if (!is.na(ant.y) && ant.y %in% activation.antigens) is.activation <- TRUE
    } else if (g$ndims == 2) {
      both.act <- (!is.na(ant.x) && ant.x %in% activation.antigens) &&
                  (!is.na(ant.y) && ant.y %in% activation.antigens)
      if (both.act) is.activation <- TRUE
    }
    if (grepl("^[A-Za-z0-9-]+\\+$", g$name)) is.activation <- TRUE
    if (g$is.singlet || g$is.scatter || g$is.viability) is.activation <- FALSE
    if (is.activation) {
      g$is.activation <- TRUE
      g$role          <- "activation"
      g$gate.type     <- "1dsep"
      g$algorithm     <- "tail"
    }

    # Rule 21: NOT gate
    if (g$is.not.gate || g$gate.class == "NotGate") {
      g$role     <- "not.gate"
      g$polarity <- "negative"
    }

    # Rule 16: positive/negative determination
    if (!g$is.activation && !g$is.not.gate && g$ndims >= 1) {
      if (grepl("-$", g$name) || grepl("(?i)(neg|^no[^n]|^non-|^CD\\d+-$)", g$name))
        g$polarity <- "negative"
      if (!is.na(g$min.y) && g$min.y < -cyto.params$biexp.near.zero)
        g$polarity <- "negative"
    }

    gate.tree[[i]] <- g
  }

  # ---- pass 2: activation gates with non-activation children ---------------
  for (i in seq_along(gate.tree)) {
    g <- gate.tree[[i]]
    if (g$is.activation && length(g$children) > 0) {
      has.pop.children <- FALSE
      for (cp in g$children) {
        if (cp %in% names(path.idx)) {
          child.g <- gate.tree[[ path.idx[[cp]] ]]
          if (!child.g$is.activation) { has.pop.children <- TRUE; break }
        }
      }
      if (has.pop.children) {
        g$is.activation <- FALSE
        g$role          <- "population"
        g$gate.type     <- "1dsep"
        g$algorithm     <- "mixture"
        gate.tree[[i]]  <- g
      }
    }
  }

  # ---- pass 3: detect 2dsep ------------------------------------------------
  parent.groups <- list()
  for (i in seq_along(gate.tree)) {
    g <- gate.tree[[i]]
    if (g$role != "population" || g$is.singlet || g$is.scatter) next
    if (is.na(g$antigen.x) && is.na(g$antigen.y)) next
    key <- paste(g$parent.path,
                 paste(sort(c(g$dim.x, g$dim.y)), collapse = "|"),
                 sep = "###")
    parent.groups[[key]] <- c(parent.groups[[key]], i)
  }
  for (key in names(parent.groups)) {
    indices <- parent.groups[[key]]
    if (length(indices) >= 2) {
      for (idx in indices) {
        if (gate.tree[[idx]]$ndims >= 2)
          gate.tree[[idx]]$gate.type <- "2dsep"
      }
    }
  }

  # ---- pass 4: two-marker gate resolution ----------------------------------
  for (i in seq_along(gate.tree)) {
    g <- gate.tree[[i]]
    if (g$role != "population" || g$ndims != 2) next
    if (g$is.singlet || g$is.scatter || g$is.viability) next
    ant.x <- g$antigen.x; ant.y <- g$antigen.y
    if (is.na(ant.x) || is.na(ant.y)) next
    if (ant.x %in% cyto.params$all.scatter && ant.y %in% cyto.params$all.scatter) next
    near.zero    <- cyto.params$biexp.near.zero
    x.spans.zero <- (!is.na(g$min.x) && !is.na(g$max.x) && g$min.x < near.zero && g$max.x > near.zero)
    y.spans.zero <- (!is.na(g$min.y) && !is.na(g$max.y) && g$min.y < near.zero && g$max.y > near.zero)
    if (x.spans.zero && !y.spans.zero) {
      g$gate.type <- "1dsep"
    } else if (!x.spans.zero && y.spans.zero) {
      g$gate.type <- "1dsep"
    } else if (!x.spans.zero && !y.spans.zero) {
      g$gate.type <- if (two.marker.strategy == "2dsep") "2dsep" else "1dsep"
    }
    gate.tree[[i]] <- g
  }

  # ---- pass 5: replicate population gates ----------------------------------
  gate.signatures <- list()
  for (i in seq_along(gate.tree)) {
    g <- gate.tree[[i]]
    if (g$role == "activation" || g$is.not.gate) next
    sig <- paste(g$dim.x, g$dim.y,
                 round(g$min.x, 2), round(g$max.x, 2),
                 round(g$min.y, 2), round(g$max.y, 2), sep = "|")
    if (sig %in% names(gate.signatures)) {
      g$is.replicate     <- TRUE
      g$replicate.origin <- gate.signatures[[sig]]
      g$algorithm        <- "replicate"
      gate.tree[[i]]     <- g
    } else {
      gate.signatures[[sig]] <- g$path
    }
  }

  # ---- pass 6: replicate activation gates ----------------------------------
  act.signatures <- list()
  for (i in seq_along(gate.tree)) {
    g <- gate.tree[[i]]
    if (g$role != "activation") next
    sig <- paste(g$dim.x, round(g$min.x, 2), round(g$max.x, 2), sep = "|")
    if (sig %in% names(act.signatures)) {
      g$is.replicate     <- TRUE
      g$replicate.origin <- act.signatures[[sig]]
      gate.tree[[i]]     <- g
    } else {
      act.signatures[[sig]] <- g$path
    }
  }

  if (verbose) {
    roles <- vapply(gate.tree, `[[`, character(1), "role")
    message(sprintf("Classification: %d population, %d activation, %d NOT gates",
                    sum(roles == "population"), sum(roles == "activation"),
                    sum(roles == "not.gate")))
    n.rep <- sum(vapply(gate.tree, `[[`, logical(1), "is.replicate"))
    if (n.rep > 0) message(sprintf("  %d replicate gates detected", n.rep))
  }

  gate.tree
}


# =============================================================================
# clean.gate.name
# =============================================================================

#' Sanitise a FlowJo gate name for use as an AutoGating list key
#'
#' @description
#' Converts a raw FlowJo gate name (which may contain spaces, slashes,
#' parentheses, and trailing \code{+}/\code{-} markers) into a lowercase,
#' dot-separated key suitable for use as an R list name.
#'
#' @param name Character scalar.  Raw gate name from the FlowJo workspace.
#'
#' @return Character scalar.  Cleaned, lowercase gate name.
#'
#' @keywords internal
clean.gate.name <- function(name) {
  x <- name
  x <- gsub("[+]$", "", x)
  x <- gsub("[-]$", "", x)
  x <- gsub("[/ ,;:]+", ".", x)
  x <- gsub("[()\\[\\]{}]", "", x, perl = TRUE)
  x <- gsub("[^A-Za-z0-9._-]", "", x)
  x <- gsub("\\.{2,}", ".", x)
  x <- gsub("^\\.|\\.$", "", x)
  x <- tolower(x)
  x <- gsub("[+]", ".pos", x)
  x <- gsub("[-](?=[A-Z])", ".neg.", x, perl = TRUE, ignore.case = TRUE)
  x <- gsub("\\.{2,}", ".", x)
  x <- gsub("^\\.|\\.$", "", x)
  x
}


# =============================================================================
# build.population.gates
# =============================================================================

#' Build the fcs.population.gates list from a classified gate tree
#'
#' @description
#' Iterates breadth-first over the classified gate tree, creating one
#' AutoGating population gate definition per population-role gate.  A root
#' \code{"all"} gate is always prepended.  NOT gates are processed in a second
#' pass after the main population gates.
#'
#' @param gate.tree List of classified gate records (output of
#'   \code{classify.gates}).
#' @param cyto.params List from \code{get.cytometer.zero.params}.
#' @param verbose Logical.  Print a count of gates built.
#'
#' @return A named list of population gate definitions compatible with
#'   \code{\link{do.population.gates}}.
#'
#' @keywords internal
build.population.gates <- function(gate.tree, cyto.params, verbose = TRUE) {

  pop.gates <- list()
  pop.idx   <- list()
  gate.num  <- 0L

  # ---- Gate 0: root --------------------------------------------------------
  all.gate <- list(
    gate.number   = 0L,
    gate.name     = "all",
    gate.type     = "free",
    gate.marker.x = cyto.params$default.scatter[1],
    gate.marker.y = cyto.params$default.scatter[2],
    parent.gate   = NULL,
    parent.popul  = NULL,
    populations   = list(
      all = list(label = "All", label.pos = 0L)
    ),
    algorithm     = "all",
    gate.param    = list(quantile.x = 0, quantile.y = 0),
    stats.parent  = list(gate = "*", popul = "*")
  )
  pop.gates[["all"]] <- all.gate
  pop.idx[["__ROOT__"]] <- "all"

  depths    <- vapply(gate.tree, `[[`, numeric(1), "depth")
  order.idx <- order(depths)
  viable.gate.name <- "all"

  for (oi in order.idx) {
    g <- gate.tree[[oi]]
    if (g$role != "population" && g$role != "not.gate") next
    if (g$is.not.gate && g$role == "not.gate") next

    gate.num  <- gate.num + 1L
    gate.name <- clean.gate.name(g$name)

    if (gate.name %in% names(pop.gates)) {
      parent.clean <- clean.gate.name(basename(g$parent.path))
      gate.name <- paste(parent.clean, gate.name, sep = ".")
    }
    if (gate.name %in% names(pop.gates))
      gate.name <- paste0(gate.name, ".", gate.num)

    pop.idx[[g$path]] <- gate.name

    # ---- parent -----------------------------------------------------------
    parent.gate.name <- NULL; parent.popul.name <- NULL
    if (nchar(g$parent.path) > 0 && g$parent.path %in% names(pop.idx)) {
      parent.gate.name  <- pop.idx[[g$parent.path]]
      parent.popul.name <- parent.gate.name
    } else if (nchar(g$parent.path) > 0) {
      pp <- g$parent.path
      while (nchar(pp) > 0 && !(pp %in% names(pop.idx))) pp <- dirname_gate(pp)
      if (pp %in% names(pop.idx)) {
        parent.gate.name  <- pop.idx[[pp]]
        parent.popul.name <- parent.gate.name
      } else {
        parent.gate.name  <- "all"
        parent.popul.name <- "all"
      }
    } else {
      parent.gate.name  <- "all"
      parent.popul.name <- "all"
    }

    # ---- markers ----------------------------------------------------------
    marker.x <- if (!is.na(g$antigen.x)) g$antigen.x else g$dim.x
    marker.y <- if (!is.na(g$antigen.y)) g$antigen.y else g$dim.y
    if (is.na(marker.x)) marker.x <- cyto.params$default.scatter[1]
    if (is.na(marker.y)) marker.y <- cyto.params$default.scatter[2]
    if (g$ndims == 1) { marker.y <- marker.x; marker.x <- cyto.params$default.scatter[1] }

    # ---- gate.type / algorithm --------------------------------------------
    gate.type <- g$gate.type
    algorithm <- g$algorithm
    if (g$is.replicate) {
      algorithm <- "replicate"
      gate.type <- gate.tree[[ which(vapply(gate.tree, function(x)
        x$path == g$replicate.origin, logical(1)))[1] ]]$gate.type
    }
    if (g$is.singlet) { marker.x <- g$dim.x; marker.y <- g$dim.y }

    # ---- populations list -------------------------------------------------
    popul.list <- list()
    siblings <- get.sibling.gates(g, gate.tree)
    sibling.on.same.markers <- Filter(function(s) {
      same.dims <- (identical(s$dim.x, g$dim.x) && identical(s$dim.y, g$dim.y)) ||
                   (identical(s$dim.x, g$dim.y) && identical(s$dim.y, g$dim.x))
      same.dims && s$role == "population"
    }, siblings)

    if (length(sibling.on.same.markers) > 0 && gate.type == "2dsep") {
      popul.name <- clean.gate.name(g$name)
      popul.list[[popul.name]] <- list(label = g$name, label.pos = 4L)
      popul.list[["*"]] <- list(label = "*", label.pos = "*")
    } else {
      popul.name <- gate.name
      popul.list[[popul.name]] <- list(label = g$name, label.pos = 4L)
      if (g$polarity == "negative") {
        popul.list <- list()
        popul.list[[gate.name]] <- list(label = g$name, label.pos = 4L)
        popul.list[["*"]] <- list(label = "*", label.pos = "*")
      } else if (!g$is.singlet && !g$is.scatter) {
        popul.list[["*"]] <- list(label = "*", label.pos = "*")
      }
    }

    # ---- gate.param -------------------------------------------------------
    gate.param <- list()
    if (g$is.singlet) {
      gate.param <- list(quantile.x = 0.0001, quantile.y = 0.0001,
                         density.quantile.x = 0.0001, density.quantile.y = 0.0001,
                         density.trimneg.x = TRUE, density.trimneg.y = TRUE)
    } else if (g$is.scatter) {
      gate.param <- list(quantile.x = 0, density.quantile.x = 0,
                         density.trimneg.x = TRUE,
                         sep.axis = 1L, region.min = 0.01, region.max = 0.30)
    } else if (gate.type == "1dsep") {
      gate.param <- list(quantile.x = 0, sep.axis = 2L,
                         region.min = 0.01, region.max = 0.99)
    } else if (gate.type == "2dsep") {
      gate.param <- list(quantile.x = 0.0001, quantile.y = 0.0001,
                         density.quantile.x = 0.0001, density.quantile.y = 0.0001,
                         region.xmin = 0.01, region.xmax = 0.99,
                         region.ymin = 0.01, region.ymax = 0.99)
    }
    if (g$is.replicate && !is.na(g$replicate.origin)) {
      origin.clean <- pop.idx[[g$replicate.origin]]
      if (!is.null(origin.clean)) gate.param <- list(origin.gate = origin.clean)
    }

    if (g$is.viability || g$is.scatter || g$is.singlet)
      viable.gate.name <- gate.name

    gd <- list(
      gate.number   = gate.num,
      gate.name     = gate.name,
      gate.type     = gate.type,
      gate.marker.x = marker.x,
      gate.marker.y = marker.y,
      parent.gate   = parent.gate.name,
      parent.popul  = parent.popul.name,
      populations   = popul.list,
      algorithm     = algorithm,
      gate.param    = gate.param,
      stats.parent  = list(gate = "*", popul = "*")
    )
    attr(gd, "flowjo.label") <- g$name
    pop.gates[[gate.name]] <- gd
  }

  # ---- NOT gates -----------------------------------------------------------
  for (oi in order.idx) {
    g <- gate.tree[[oi]]
    if (!g$is.not.gate) next
    ref.name <- sub("-$", "", g$name)
    ref.path <- if (nchar(g$parent.path) > 0) paste(g$parent.path, ref.name, sep = "/")
                else ref.name
    if (ref.path %in% names(pop.idx)) {
      gate.num      <- gate.num + 1L
      gate.name     <- clean.gate.name(g$name)
      if (gate.name %in% names(pop.gates)) gate.name <- paste0(gate.name, ".not")
      ref.gate.name <- pop.idx[[ref.path]]
      ref.gd        <- pop.gates[[ref.gate.name]]
      gd <- list(
        gate.number   = gate.num,
        gate.name     = gate.name,
        gate.type     = ref.gd$gate.type,
        gate.marker.x = ref.gd$gate.marker.x,
        gate.marker.y = ref.gd$gate.marker.y,
        parent.gate   = ref.gd$parent.gate,
        parent.popul  = ref.gd$parent.popul,
        populations   = list(),
        algorithm     = "not",
        gate.param    = list(reference.gate = ref.gate.name),
        stats.parent  = list(gate = "*", popul = "*")
      )
      gd$populations[[gate.name]] <- list(label = g$name, label.pos = 2L)
      pop.gates[[gate.name]] <- gd
      pop.idx[[g$path]]      <- gate.name
    } else {
      if (verbose)
        message(sprintf("NOT gate '%s': could not find reference gate '%s'",
                        g$name, ref.name))
    }
  }

  if (verbose)
    message(sprintf("Built %d population gate definitions", length(pop.gates)))
  pop.gates
}


# =============================================================================
# build.activation.gates
# =============================================================================

#' Build the fcs.activation.gates list from a classified gate tree
#'
#' @description
#' For every population gate, collects the activation marker gates that sit
#' directly beneath it in the gate hierarchy and assembles one activation gate
#' definition entry (with a \code{markers} character vector).  Population
#' gates with no activation children produce display-only entries
#' (\code{gate.type = ""}, \code{markers = character(0)}).
#'
#' @param gate.tree List of classified gate records (output of
#'   \code{classify.gates}).
#' @param pop.gates Named list of population gate definitions (output of
#'   \code{build.population.gates}).
#' @param verbose Logical.  Print a summary count.
#'
#' @return A list of activation gate definitions compatible with
#'   \code{\link{do.activation.gates}}.
#'
#' @keywords internal
build.activation.gates <- function(gate.tree, pop.gates, verbose = TRUE) {

  act.gates <- list()
  path.to.idx <- stats::setNames(seq_along(gate.tree),
                           vapply(gate.tree, `[[`, character(1), "path"))

  # Group activation gates by parent path
  act.by.parent <- list()
  for (i in seq_along(gate.tree)) {
    g <- gate.tree[[i]]
    if (g$role != "activation") next
    key <- g$parent.path
    if (!(key %in% names(act.by.parent))) act.by.parent[[key]] <- list()
    act.by.parent[[key]] <- c(act.by.parent[[key]], list(g))
  }

  popul.number <- 0L
  for (pg.name in names(pop.gates)) {
    pg <- pop.gates[[pg.name]]

    main.popul  <- setdiff(names(pg$populations), "*")
    if (length(main.popul) == 0) main.popul <- pg.name
    popul.label <- pg$populations[[ main.popul[1] ]]$label
    if (is.null(popul.label)) popul.label <- pg.name

    # Match activation gates to this population by path basename
    matching.path <- NULL
    for (pp in names(act.by.parent)) {
      pp.basename <- basename_gate(pp)
      if (clean.gate.name(pp.basename) == pg.name ||
          clean.gate.name(pp) == pg.name) {
        matching.path <- pp; break
      }
    }

    act.markers         <- character(0)
    marker.alg.override <- list()

    if (!is.null(matching.path)) {
      for (ag in act.by.parent[[matching.path]]) {
        marker.name <- if (!is.na(ag$antigen.x) && ag$ndims == 1) ag$antigen.x
                       else if (!is.na(ag$antigen.y))             ag$antigen.y
                       else if (!is.na(ag$dim.x))                 ag$dim.x
                       else NA_character_
        if (!is.na(marker.name)) {
          act.markers <- c(act.markers, marker.name)
          if (ag$is.replicate) marker.alg.override[[marker.name]] <- "replicate"
        }
      }
    }
    act.markers <- unique(act.markers)

    parent.marker <- pg$gate.marker.y
    if (is.null(parent.marker) || is.na(parent.marker)) parent.marker <- "*"

    ad <- list(
      popul.number  = popul.number,
      parent.gate   = pg.name,
      parent.popul  = main.popul[1],
      popul.label   = popul.label,
      label.pos     = if (length(act.markers) > 0) 3L else 0L,
      gate.marker.x = if (length(act.markers) > 0) parent.marker else "*",
      gate.marker.y = "*",
      gate.type     = if (length(act.markers) > 0) "1dsep" else "",
      algorithm     = if (length(act.markers) > 0) "tail" else NULL,
      gate.param    = if (length(act.markers) > 0)
        list(quantile.x = 0, sep.axis = 2L, region.min = 0.01, region.max = 0.99)
        else list(),
      markers                   = act.markers,
      marker.algorithm.override = marker.alg.override
    )

    act.gates    <- c(act.gates, list(ad))
    popul.number <- popul.number + 1L
  }

  if (verbose) {
    n.with  <- sum(vapply(act.gates, function(a) length(a$markers) > 0, logical(1)))
    n.total <- sum(vapply(act.gates, function(a) length(a$markers),    integer(1)))
    message(sprintf("Built %d activation gate entries (%d with markers, %d total markers)",
                    length(act.gates), n.with, n.total))
  }
  act.gates
}


# =============================================================================
# Utility functions
# =============================================================================

#' Gate path dirname (everything before the last slash)
#' @keywords internal
dirname_gate <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  if (length(parts) <= 1) return("")
  paste(parts[-length(parts)], collapse = "/")
}

#' Gate path basename (last component after the last slash)
#' @keywords internal
basename_gate <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  parts[length(parts)]
}

#' Return sibling gates (same parent, different path)
#' @keywords internal
get.sibling.gates <- function(g, gate.tree) {
  Filter(function(x) {
    x$parent.path == g$parent.path && x$path != g$path
  }, gate.tree)
}


# =============================================================================
# write.gate.definitions
# =============================================================================

#' Write population and activation gate definitions to R source files
#'
#' @title Save Gate Definitions to R Source Files
#'
#' @description
#' Serialises the \code{population.gates} and \code{activation.gates} lists
#' returned by \code{\link{wsp.to.autogating}} to \code{.r} (or \code{.rds})
#' files that can be sourced or loaded directly by the AutoGating pipeline.
#' Delegates the actual formatting to \code{write.pop.file} and
#' \code{write.act.file}.
#'
#' @param result Named list returned by \code{\link{wsp.to.autogating}},
#'   containing at least \code{population.gates} and \code{activation.gates}.
#' @param pop.filename Character scalar.  Output path for the population gate
#'   definitions.  Default \code{"fcs_population_gates.r"}.
#' @param act.filename Character scalar.  Output path for the activation gate
#'   definitions.  Default \code{"fcs_activation_gates.r"}.
#'
#' @return \code{NULL} invisibly.
#'
#' @seealso \code{\link{wsp.to.autogating}}, \code{\link{load.gate.definitions}}
#'
#' @export
write.gate.definitions <- function(
    result,
    pop.filename = "fcs_population_gates.r",
    act.filename = "fcs_activation_gates.r"
) {
  write.pop.file(result$population.gates, pop.filename)
  message(sprintf("Population gates written to: %s", pop.filename))
  write.act.file(result$activation.gates, act.filename)
  message(sprintf("Activation gates written to: %s", act.filename))
}


# =============================================================================
# write.pop.file / write.act.file — R-source serialisers
# =============================================================================

#' Serialise a population gate list to an R source file
#' @keywords internal
write.pop.file <- function(pop.gates, filename) {
  lines <- c(
    "# Auto-generated population gate definitions",
    sprintf("# Generated: %s", Sys.time()),
    "#",
    "# Each element of fcs.population.gates is a named list describing one gate.",
    "# See fcs_population_gates_st1.r for field documentation.",
    "",
    "fcs.population.gates <- list(",
    ""
  )

  for (i in seq_along(pop.gates)) {
    gd      <- pop.gates[[i]]
    is.last <- (i == length(pop.gates))

    lines <- c(lines, sprintf("    # Gate %d -- %s", gd$gate.number,
                               attr(gd, "flowjo.label") %||% gd$gate.name))
    lines <- c(lines, "    list(")
    lines <- c(lines, sprintf("        gate.number  = %dL,",  gd$gate.number))
    lines <- c(lines, sprintf('        gate.name    = "%s",', gd$gate.name))
    lines <- c(lines, sprintf('        gate.type    = "%s",', gd$gate.type))
    lines <- c(lines, sprintf('        gate.marker.x = "%s",', gd$gate.marker.x))
    lines <- c(lines, sprintf('        gate.marker.y = "%s",', gd$gate.marker.y))

    if (is.null(gd$parent.gate)) {
      lines <- c(lines, "        parent.gate  = NULL,")
    } else if (length(gd$parent.gate) == 1) {
      lines <- c(lines, sprintf('        parent.gate  = "%s",', gd$parent.gate))
    } else {
      lines <- c(lines, sprintf('        parent.gate  = c( %s ),',
                                 paste(sprintf('"%s"', gd$parent.gate), collapse = ", ")))
    }
    if (is.null(gd$parent.popul)) {
      lines <- c(lines, "        parent.popul = NULL,")
    } else if (length(gd$parent.popul) == 1) {
      lines <- c(lines, sprintf('        parent.popul = "%s",', gd$parent.popul))
    } else {
      lines <- c(lines, sprintf('        parent.popul = c( %s ),',
                                 paste(sprintf('"%s"', gd$parent.popul), collapse = ", ")))
    }

    lines <- c(lines, "        populations  = list(")
    pop.names <- names(gd$populations)
    for (j in seq_along(pop.names)) {
      pn    <- pop.names[j]
      pi    <- gd$populations[[pn]]
      lpos  <- pi$label.pos
      lpos.str <- if (identical(lpos, "*")) '"*"' else sprintf("%dL", as.integer(lpos))
      pn.str   <- if (pn == "*") "`*`" else pn
      comma    <- if (j < length(pop.names)) "," else ""
      lines <- c(lines, sprintf('            %s = list( label = "%s", label.pos = %s )%s',
                                 pn.str, pi$label, lpos.str, comma))
    }
    lines <- c(lines, "        ),")

    if (length(gd$algorithm) == 1) {
      lines <- c(lines, sprintf('        algorithm    = "%s",', gd$algorithm))
    } else {
      lines <- c(lines, sprintf('        algorithm    = c( %s ),',
                                 paste(sprintf('"%s"', gd$algorithm), collapse = ", ")))
    }

    if (length(gd$gate.param) == 0) {
      lines <- c(lines, "        gate.param   = list(),")
    } else {
      lines <- c(lines, "        gate.param   = list(")
      gp.names <- names(gd$gate.param)
      for (k in seq_along(gp.names)) {
        val   <- gd$gate.param[[ gp.names[k] ]]
        comma <- if (k < length(gp.names)) "," else ""
        if (is.character(val)) {
          lines <- c(lines, sprintf('            %s = "%s"%s', gp.names[k], val, comma))
        } else if (is.logical(val)) {
          lines <- c(lines, sprintf('            %s = %s%s', gp.names[k],
                                    toupper(as.character(val)), comma))
        } else {
          lines <- c(lines, sprintf('            %s = %s%s', gp.names[k],
                                    format(val, digits = 4), comma))
        }
      }
      lines <- c(lines, "        ),")
    }

    sg <- gd$stats.parent$gate;  if (is.null(sg)) sg <- "*"
    sp <- gd$stats.parent$popul; if (is.null(sp)) sp <- "*"
    if (length(sg) == 1) {
      lines <- c(lines, sprintf('        stats.parent = list( gate = "%s", popul = "%s" )',
                                 sg, sp))
    } else {
      lines <- c(lines, sprintf('        stats.parent = list( gate = c( %s ), popul = c( %s ) )',
                                 paste(sprintf('"%s"', sg), collapse = ", "),
                                 paste(sprintf('"%s"', sp), collapse = ", ")))
    }

    lines <- c(lines, sprintf("    )%s", if (!is.last) "," else ""), "")
  }

  lines <- c(lines, ")", "",
    "# Name the list by gate.name for convenient lookup",
    'names( fcs.population.gates ) <- sapply( fcs.population.gates, `[[`, "gate.name" )')
  writeLines(lines, filename)
}


#' Serialise an activation gate list to an R source file
#' @keywords internal
write.act.file <- function(act.gates, filename) {
  lines <- c(
    "# Auto-generated activation gate definitions",
    sprintf("# Generated: %s", Sys.time()),
    "#",
    "# Each element defines the activation gating performed within one cell",
    "# population.  See fcs_activation_gates_st1.r for field documentation.",
    "",
    "fcs.activation.gates <- list(",
    ""
  )

  for (i in seq_along(act.gates)) {
    ad      <- act.gates[[i]]
    is.last <- (i == length(act.gates))

    if (ad$gate.type == "" || length(ad$markers) == 0) {
      lines <- c(lines,
        sprintf('    list( popul.number = %dL, parent.gate = "%s", parent.popul = "%s",',
                ad$popul.number, ad$parent.gate, ad$parent.popul),
        sprintf('          popul.label = "%s", label.pos = %dL,',
                ad$popul.label, as.integer(ad$label.pos)),
        '          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",',
        '          algorithm = NULL, gate.param = list(), markers = character(0),',
        sprintf('          marker.algorithm.override = list() )%s',
                if (!is.last) "," else ""))
    } else {
      lines <- c(lines, "    list(")
      lines <- c(lines, sprintf("        popul.number = %dL,",   ad$popul.number))
      lines <- c(lines, sprintf('        parent.gate  = "%s",',  ad$parent.gate))
      lines <- c(lines, sprintf('        parent.popul = "%s",',  ad$parent.popul))
      lines <- c(lines, sprintf('        popul.label  = "%s",',  ad$popul.label))
      lines <- c(lines, sprintf("        label.pos    = %dL,",   as.integer(ad$label.pos)))
      lines <- c(lines, sprintf('        gate.marker.x = "%s",', ad$gate.marker.x))
      lines <- c(lines, sprintf('        gate.marker.y = "%s",', ad$gate.marker.y))
      lines <- c(lines, sprintf('        gate.type    = "%s",',  ad$gate.type))
      lines <- c(lines, sprintf('        algorithm    = "%s",',  ad$algorithm))

      if (length(ad$gate.param) == 0) {
        lines <- c(lines, "        gate.param   = list(),")
      } else {
        lines <- c(lines, "        gate.param   = list(")
        gp.names <- names(ad$gate.param)
        for (k in seq_along(gp.names)) {
          comma <- if (k < length(gp.names)) "," else ""
          lines <- c(lines, sprintf("            %s = %s%s",
                                     gp.names[k],
                                     format(ad$gate.param[[ gp.names[k] ]], digits = 4),
                                     comma))
        }
        lines <- c(lines, "        ),")
      }

      if (length(ad$markers) == 0) {
        lines <- c(lines, "        markers = character(0),")
      } else {
        lines <- c(lines, sprintf("        markers = c( %s ),",
                                   paste(sprintf('"%s"', ad$markers), collapse = ", ")))
      }

      if (length(ad$marker.algorithm.override) == 0) {
        lines <- c(lines, "        marker.algorithm.override = list()")
      } else {
        lines <- c(lines, "        marker.algorithm.override = list(")
        mao.names <- names(ad$marker.algorithm.override)
        for (k in seq_along(mao.names)) {
          comma <- if (k < length(mao.names)) "," else ""
          lines <- c(lines, sprintf('            %s = "%s"%s',
                                     mao.names[k],
                                     ad$marker.algorithm.override[[ mao.names[k] ]],
                                     comma))
        }
        lines <- c(lines, "        )")
      }

      lines <- c(lines, sprintf("    )%s", if (!is.last) "," else ""))
    }
    lines <- c(lines, "")
  }

  lines <- c(lines, ")")
  writeLines(lines, filename)
}


# NULL-coalescing operator (package-private)
`%||%` <- function(a, b) if (is.null(a)) b else a
