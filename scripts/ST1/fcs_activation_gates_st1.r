# fcs_activation_gates_st1.r
#
# Activation gate definitions for the ST1 panel (BD LSRFortessa).
#
# Each element defines the activation gating performed within one cell
# population.  For populations that are shown in figures but have no
# activation markers, only the display metadata is needed.
#
# Fields:
#   popul.number  Integer.  Matches the row number in the original CSV for
#                 traceability.
#   parent.gate   Gate name (character).
#   parent.popul  Population name within that gate (character).
#   popul.label   Display label for the parent population on figures.
#   label.pos     Integer 1-8 or "*" (suppress).
#   gate.marker.x Activation axis marker, or "*" if axis is y.
#   gate.marker.y Usually "*" (activation is 1-D on x, or the column marker).
#   gate.type     "1dsep" or "" (no gating, display-only row).
#   algorithm     Default algorithm for this population's activation gates.
#                 Individual markers may override via marker.algorithm.override.
#   gate.param    Named list of per-population overrides to cytometer defaults.
#   markers       Character vector of activation markers assayed in this
#                 population (using "X" convention from the original CSV).
#   marker.algorithm.override
#                 Named list: marker -> algorithm string, for any marker that
#                 should use a different algorithm than the population default.
#                 Omit entirely if all markers use the default.
#
# Display-only rows (no activation gating) have gate.type = "" and empty
# markers / gate.param.  CTLA-4 is the canonical spelling.


fcs.activation.gates <- list(

    # ---- display-only rows (no activation gates) --------------------------

    list( popul.number = 0L,  parent.gate = "all",          parent.popul = "all",
          popul.label = "All",                    label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list( popul.number = 1L,  parent.gate = "singlets",     parent.popul = "singlets",
          popul.label = "Singlets",               label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list( popul.number = 2L,  parent.gate = "leukocytes",   parent.popul = "leukocytes",
          popul.label = "Leukocytes",             label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list( popul.number = 3L,  parent.gate = "monocytes",    parent.popul = "no.monocytes",
          popul.label = "No Monocytes",           label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list( popul.number = 4L,  parent.gate = "monocytes",    parent.popul = "monocytes",
          popul.label = "Monocytes",              label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list( popul.number = 5L,  parent.gate = "live.no.monocytes", parent.popul = "live.no.monocytes",
          popul.label = "Live No Monocytes",      label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    # ---- populations with activation gates --------------------------------

    list(
        popul.number = 6L,
        parent.gate  = "live.monocytes",
        parent.popul = "live.monocytes",
        popul.label  = "Live Monocytes",
        label.pos    = 3L,
        gate.marker.x = "CD14",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD95", "CXCR3", "HLA-DR",
                     "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CCR2 = "mixture",
          CD95 = "mixture",
          `HLA-DR` = "mixture"
        )
    ),

    list(
        popul.number = 7L,
        parent.gate  = "t.cells",
        parent.popul = "no.t.cells",
        popul.label  = "No T Cells",
        label.pos    = 4L,
        gate.marker.x = "CD3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD95", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          `HLA-DR` = "mixture"
        )
    ),

    list( popul.number = 8L,  parent.gate = "t.cells",      parent.popul = "t.cells",
          popul.label = "T Cells",                label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list( popul.number = 9L,  parent.gate = "cd4.cd8",      parent.popul = "double.neg",
          popul.label = "CD4-/CD8- T Cells",      label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list( popul.number = 10L, parent.gate = "cd4.cd8",      parent.popul = "cd8.pos",
          popul.label = "CD8+ T Cells",           label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list(
        popul.number = 11L,
        parent.gate  = "cd4.cd8",
        parent.popul = "cd4.pos",
        popul.label  = "CD4+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD4",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25" ),
        marker.algorithm.override = list()
    ),

    list( popul.number = 12L, parent.gate = "cd4.cd8",      parent.popul = "double.pos",
          popul.label = "CD4+/CD8+ T Cells",      label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list(
        popul.number = 13L,
        parent.gate  = "cd4.memory",
        parent.popul = "tem",
        popul.label  = "Effector Memory CD4+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD4",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD95", "CD127", "CTLA-4", "CXCR5",
                     "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 14L,
        parent.gate  = "cd4.memory",
        parent.popul = "temra",
        popul.label  = "Effector Memory RA CD4+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD4",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD95", "CD127", "CTLA-4", "CXCR5",
                     "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 15L,
        parent.gate  = "cd4.memory",
        parent.popul = "tcm",
        popul.label  = "Central Memory CD4+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD4",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD95", "CD127", "CTLA-4", "CXCR5",
                     "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 16L,
        parent.gate  = "cd4.memory",
        parent.popul = "naive",
        popul.label  = "Naive CD4+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD4",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD127", "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    list( popul.number = 17L, parent.gate = "cd4.no.naive", parent.popul = "no.naive",
          popul.label = "No Naive CD4+ T Cells",  label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list(
        popul.number = 18L,
        parent.gate  = "cd8.memory",
        parent.popul = "tem",
        popul.label  = "Effector Memory CD8+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD8",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "CXCR3", "CXCR5", "HLA-DR", "ICOS",
                     "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
            CCR2  = "mixture",
            CCR4  = "mixture",
            CD28  = "mixture",
            CD95  = "mixture",
            CD127 = "mixture"
        )
    ),

    list(
        popul.number = 19L,
        parent.gate  = "cd8.memory",
        parent.popul = "temra",
        popul.label  = "Effector Memory RA CD8+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD8",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "CXCR3", "CXCR5", "HLA-DR", "ICOS",
                     "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list( CD28 = "mixture" )
    ),

    list(
        popul.number = 20L,
        parent.gate  = "cd8.memory",
        parent.popul = "tcm",
        popul.label  = "Central Memory CD8+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD8",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "CXCR3", "CXCR5", "HLA-DR", "ICOS",
                     "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CCR2  = "mixture",
          CCR4  = "mixture",
          CD28  = "mixture",
          CD95  = "mixture",
          CD127 = "mixture",
          CXCR3 = "mixture"
        )
    ),

    list(
        popul.number = 21L,
        parent.gate  = "cd8.memory",
        parent.popul = "naive",
        popul.label  = "Naive CD8+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD8",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD127", "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CD28  = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 22L,
        parent.gate  = "cd4.naive.cd31",
        parent.popul = "rte",
        popul.label  = "Recent Thymic Emigrant CD4+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD4",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD127", "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CD28  = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 23L,
        parent.gate  = "cd4.naive.cd95",
        parent.popul = "scm",
        popul.label  = "Stem Memory CD4+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD4",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD127", "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CD28  = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 24L,
        parent.gate  = "cd8.naive.cd31",
        parent.popul = "rte",
        popul.label  = "Recent Thymic Emigrant CD8+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD8",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD127", "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CD28  = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 25L,
        parent.gate  = "cd8.naive.cd95",
        parent.popul = "scm",
        popul.label  = "Stem Memory CD8+ T Cells",
        label.pos    = 3L,
        gate.marker.x = "CD8",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD127", "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CD28  = "mixture",
          CD127 = "mixture"
        )
    ),

    list( popul.number = 26L, parent.gate = "cd4.treg", parent.popul = "no.treg",
          popul.label = "No Treg CD4+ T Cells",   label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list(
        popul.number = 27L,
        parent.gate  = "cd4.treg",
        parent.popul = "treg",
        popul.label  = "Regulatory T Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CTLA-4",
                     "CXCR3", "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD25 = "mixture",
          CD28  = "mixture",
          CD95  = "mixture",
          CD127 = "mixture"
        )
    ),

    # Display-only rows for no-Treg memory subsets
    list( popul.number = 28L, parent.gate = "cd4.no.treg.memory", parent.popul = "tem",
          popul.label = "Effector Memory No Treg CD4+ T Cells",   label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 29L, parent.gate = "cd4.no.treg.memory", parent.popul = "temra",
          popul.label = "Effector Memory RA No Treg CD4+ T Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 30L, parent.gate = "cd4.no.treg.memory", parent.popul = "tcm",
          popul.label = "Central Memory No Treg CD4+ T Cells",     label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 31L, parent.gate = "cd4.no.treg.memory", parent.popul = "naive",
          popul.label = "Naive No Treg CD4+ T Cells",              label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 32L, parent.gate = "cd4.no.treg.naive.cd31", parent.popul = "rte",
          popul.label = "Recent Thymic Emigrant No Treg CD4+ T Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 33L, parent.gate = "cd4.no.treg.naive.cd95", parent.popul = "scm",
          popul.label = "Stem Memory No Treg CD4+ T Cells",         label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list(
        popul.number = 34L,
        parent.gate  = "cd4.treg.memory",
        parent.popul = "tem",
        popul.label  = "Effector Memory Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CTLA-4",
                     "CXCR3", "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD25 = "mixture",
          CD28  = "mixture",
          CD95  = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 35L,
        parent.gate  = "cd4.treg.memory",
        parent.popul = "temra",
        popul.label  = "Effector Memory RA Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CTLA-4",
                     "CXCR3", "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD25 = "mixture",
          CD28  = "mixture",
          CD95  = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 36L,
        parent.gate  = "cd4.treg.memory",
        parent.popul = "tcm",
        popul.label  = "Central Memory Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CTLA-4",
                     "CXCR3", "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD25 = "mixture",
          CD28  = "mixture",
          CD95  = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 37L,
        parent.gate  = "cd4.treg.memory",
        parent.popul = "naive",
        popul.label  = "Naive Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CD25 = "mixture",
          CD28  = "mixture"
        )
    ),

    list( popul.number = 38L, parent.gate = "cd4.treg.no.naive", parent.popul = "no.naive",
          popul.label = "No Naive Treg Cells",               label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list(
        popul.number = 39L,
        parent.gate  = "cd4.treg.naive.cd31",
        parent.popul = "rte",
        popul.label  = "Recent Thymic Emigrant Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CD25 = "mixture",
          CD28  = "mixture"
        )
    ),

    list(
        popul.number = 40L,
        parent.gate  = "cd4.treg.naive.cd95",
        parent.popul = "scm",
        popul.label  = "Stem Memory Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CD25 = "mixture",
          CD28  = "mixture"
        )
    ),

    list(
        popul.number = 41L,
        parent.gate  = "treg.follicular",
        parent.popul = "follicular",
        popul.label  = "Follicular Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "mixture",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CTLA-4",
                     "CXCR3", "HLA-DR", "ICOS", "KI67", "RORgt" ),
        marker.algorithm.override = list(
          CCR2 = "tail"
        )
    ),

    # Display-only rows for follicular Treg memory subsets
    list( popul.number = 42L, parent.gate = "treg.follicular.memory", parent.popul = "tem",
          popul.label = "Effector Memory Follicular Treg Cells",    label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 43L, parent.gate = "treg.follicular.memory", parent.popul = "temra",
          popul.label = "Effector Memory RA Follicular Treg Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 44L, parent.gate = "treg.follicular.memory", parent.popul = "tcm",
          popul.label = "Central Memory Follicular Treg Cells",     label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    # ---- treg.cd127 gated populations (rows 45–53) -----------------------

    list(
        popul.number = 45L,
        parent.gate  = "treg.cd127",
        parent.popul = "cd127.neg",
        popul.label  = "CD127- Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "CXCR3", "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD25 = "mixture",
          CD28  = "mixture",
          CD95 = "mixture"
        )
    ),

    list(
        popul.number = 46L,
        parent.gate  = "treg.cd127.memory",
        parent.popul = "tem",
        popul.label  = "Effector Memory CD127- Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "CXCR3", "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD25 = "mixture",
          CD28  = "mixture",
          CD95 = "mixture",
          `CTLA-4` = "mixture"
        )
    ),

    list(
        popul.number = 47L,
        parent.gate  = "treg.cd127.memory",
        parent.popul = "temra",
        popul.label  = "Effector Memory RA CD127- Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "CXCR3", "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD25 = "mixture",
          CD28  = "mixture",
          CD95 = "mixture",
          `CTLA-4` = "mixture"
        )
    ),

    list(
        popul.number = 48L,
        parent.gate  = "treg.cd127.memory",
        parent.popul = "tcm",
        popul.label  = "Central Memory CD127- Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "CXCR3", "HLA-DR", "ICOS", "KI67", "PD1", "RORgt" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD25 = "mixture",
          CD28  = "mixture",
          CD95 = "mixture",
          `CTLA-4` = "mixture"
        )
    ),

    list(
        popul.number = 49L,
        parent.gate  = "treg.cd127.memory",
        parent.popul = "naive",
        popul.label  = "Naive CD127- Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD127", "CTLA-4", "ICOS", "KI67" ),
        marker.algorithm.override = list(
          CD25 = "mixture",
          CD28  = "mixture"
        )
    ),

    list( popul.number = 50L, parent.gate = "treg.cd127.no.naive", parent.popul = "no.naive",
          popul.label = "No Naive CD127- Treg Cells",              label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list(
        popul.number = 51L,
        parent.gate  = "treg.cd127.naive.cd31",
        parent.popul = "rte",
        popul.label  = "Recent Thymic Emigrant CD127- Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD127", "CTLA-4", "ICOS", "KI67" ),
        marker.algorithm.override = list(
          CD25 = "mixture",
          CD28  = "mixture"
        )
    ),

    list(
        popul.number = 52L,
        parent.gate  = "treg.cd127.naive.cd95",
        parent.popul = "scm",
        popul.label  = "Stem Memory CD127- Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD127", "CTLA-4", "ICOS", "KI67" ),
        marker.algorithm.override = list(
          CD25 = "mixture",
          CD28  = "mixture"
        )
    ),

    list(
        popul.number = 53L,
        parent.gate  = "treg.cd127.follicular",
        parent.popul = "follicular",
        popul.label  = "Follicular CD127- Treg Cells",
        label.pos    = 3L,
        gate.marker.x = "FOXP3",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "mixture",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "CXCR3", "HLA-DR", "ICOS", "KI67", "RORgt" ),
        marker.algorithm.override = list(
          CCR2 = "tail",
          CD127 = "tail"
        )
    ),

    # ---- display-only rows 54–59 ------------------------------------------

    list( popul.number = 54L, parent.gate = "treg.cd127.follicular.memory", parent.popul = "tem",
          popul.label = "Effector Memory Follicular CD127- Treg Cells",    label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 55L, parent.gate = "treg.cd127.follicular.memory", parent.popul = "temra",
          popul.label = "Effector Memory RA Follicular CD127- Treg Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 56L, parent.gate = "treg.cd127.follicular.memory", parent.popul = "tcm",
          popul.label = "Central Memory Follicular CD127- Treg Cells",     label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 57L, parent.gate = "cd4.thelper",          parent.popul = "thelper",
          popul.label = "T Helper Cells",                            label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 58L, parent.gate = "thelper.cxcr3.ccr4",   parent.popul = "cxcr3.pos",
          popul.label = "CCR4- CXCR3+ T Helper Cells",              label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 59L, parent.gate = "thelper.cxcr3.ccr4",   parent.popul = "ccr4.pos",
          popul.label = "CCR4+ CXCR3- T Helper Cells",              label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    # ---- T helper gated populations (rows 60–77) --------------------------

    list(
        popul.number = 60L,
        parent.gate  = "thelper.cxcr3.rorgt",
        parent.popul = "th1",
        popul.label  = "T Helper 1 Cells",
        label.pos    = 4L,
        gate.marker.x = "CCR7",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CD25", "CD28", "CD95", "CD127", "CTLA-4",
                     "HLA-DR", "ICOS", "KI67", "PD1" ),
        marker.algorithm.override = list(
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 61L,
        parent.gate  = "thelper.cxcr3.rorgt",
        parent.popul = "th1.ast",
        popul.label  = "T Helper 1* Cells",
        label.pos    = 4L,
        gate.marker.x = "CCR7",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CD25", "CD28", "CD95", "CD127", "CTLA-4",
                     "HLA-DR", "ICOS", "KI67", "PD1" ),
        marker.algorithm.override = list(
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    # display-only rows 62–67
    list( popul.number = 62L, parent.gate = "thelper.th1.memory",     parent.popul = "tem",
          popul.label = "Effector Memory T Helper 1 Cells",           label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 63L, parent.gate = "thelper.th1.memory",     parent.popul = "temra",
          popul.label = "Effector Memory RA T Helper 1 Cells",        label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 64L, parent.gate = "thelper.th1.memory",     parent.popul = "tcm",
          popul.label = "Central Memory T Helper 1 Cells",            label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 65L, parent.gate = "thelper.th1.ast.memory", parent.popul = "tem",
          popul.label = "Effector Memory T Helper 1* Cells",          label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 66L, parent.gate = "thelper.th1.ast.memory", parent.popul = "temra",
          popul.label = "Effector Memory RA T Helper 1* Cells",       label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 67L, parent.gate = "thelper.th1.ast.memory", parent.popul = "tcm",
          popul.label = "Central Memory T Helper 1* Cells",           label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list(
        popul.number = 68L,
        parent.gate  = "thelper.ccr4.rorgt",
        parent.popul = "th2",
        popul.label  = "T Helper 2 Cells",
        label.pos    = 4L,
        gate.marker.x = "CCR7",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CD25", "CD28", "CD95", "CD127", "CTLA-4",
                     "HLA-DR", "ICOS", "KI67", "PD1" ),
        marker.algorithm.override = list(
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 69L,
        parent.gate  = "thelper.ccr4.rorgt",
        parent.popul = "th17",
        popul.label  = "T Helper 17 Cells",
        label.pos    = 4L,
        gate.marker.x = "CCR7",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD95", "CD127", "CTLA-4",
                     "HLA-DR", "ICOS", "KI67", "PD1" ),
        marker.algorithm.override = list(
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    # display-only rows 70–75
    list( popul.number = 70L, parent.gate = "thelper.th2.memory",  parent.popul = "tem",
          popul.label = "Effector Memory T Helper 2 Cells",        label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 71L, parent.gate = "thelper.th2.memory",  parent.popul = "temra",
          popul.label = "Effector Memory RA T Helper 2 Cells",     label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 72L, parent.gate = "thelper.th2.memory",  parent.popul = "tcm",
          popul.label = "Central Memory T Helper 2 Cells",         label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 73L, parent.gate = "thelper.th17.memory", parent.popul = "tem",
          popul.label = "Effector Memory T Helper 17 Cells",       label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 74L, parent.gate = "thelper.th17.memory", parent.popul = "temra",
          popul.label = "Effector Memory RA T Helper 17 Cells",    label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 75L, parent.gate = "thelper.th17.memory", parent.popul = "tcm",
          popul.label = "Central Memory T Helper 17 Cells",        label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    list(
        popul.number = 76L,
        parent.gate  = "thelper.th17.ccr2",
        parent.popul = "classic",
        popul.label  = "Classic T Helper 17 Cells",
        label.pos    = 4L,
        gate.marker.x = "CCR7",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD95", "CD127", "CTLA-4",
                     "HLA-DR", "ICOS", "KI67", "PD1" ),
        marker.algorithm.override = list(
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 77L,
        parent.gate  = "thelper.th17.ccr2",
        parent.popul = "dual",
        popul.label  = "Dual T Helper 17 Cells",
        label.pos    = 4L,
        gate.marker.x = "CCR7",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CD25", "CD28", "CD95", "CD127", "CTLA-4",
                     "HLA-DR", "ICOS", "KI67", "PD1" ),
        marker.algorithm.override = list(
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture",
          `CTLA-4` = "mixture",
          `HLA-DR` = "mixture"
        )
    ),

    # display-only rows 78–87
    list( popul.number = 78L, parent.gate = "thelper.th17.classic.memory", parent.popul = "tem",
          popul.label = "Effector Memory Classic T Helper 17 Cells",    label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 79L, parent.gate = "thelper.th17.classic.memory", parent.popul = "temra",
          popul.label = "Effector Memory RA Classic T Helper 17 Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 80L, parent.gate = "thelper.th17.classic.memory", parent.popul = "tcm",
          popul.label = "Central Memory Classic T Helper 17 Cells",     label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 81L, parent.gate = "thelper.th17.dual.memory",    parent.popul = "tem",
          popul.label = "Effector Memory Dual T Helper 17 Cells",        label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 82L, parent.gate = "thelper.th17.dual.memory",    parent.popul = "temra",
          popul.label = "Effector Memory RA Dual T Helper 17 Cells",     label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 83L, parent.gate = "thelper.th17.dual.memory",    parent.popul = "tcm",
          popul.label = "Central Memory Dual T Helper 17 Cells",         label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 84L, parent.gate = "thelper.follicular",          parent.popul = "follicular",
          popul.label = "Follicular T helper Cells",                      label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 85L, parent.gate = "thelper.follicular.memory",   parent.popul = "tem",
          popul.label = "Effector Memory Follicular T helper Cells",       label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 86L, parent.gate = "thelper.follicular.memory",   parent.popul = "temra",
          popul.label = "Effector Memory RA Follicular T helper Cells",    label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 87L, parent.gate = "thelper.follicular.memory",   parent.popul = "tcm",
          popul.label = "Central Memory Follicular T helper Cells",        label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),

    # ---- Tfh CXCR3/RORgT gated populations (rows 88–90) ------------------

    list(
        popul.number = 88L,
        parent.gate  = "thelper.follicular.cxcr3.rorgt",
        parent.popul = "th2",
        popul.label  = "Follicular T helper 2-like Cells",
        label.pos    = 4L,
        gate.marker.x = "PD1",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 89L,
        parent.gate  = "thelper.follicular.cxcr3.rorgt",
        parent.popul = "th1",
        popul.label  = "Follicular T helper 1-like Cells",
        label.pos    = 4L,
        gate.marker.x = "PD1",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    list(
        popul.number = 90L,
        parent.gate  = "thelper.follicular.cxcr3.rorgt",
        parent.popul = "th17",
        popul.label  = "Follicular T helper 17-like Cells",
        label.pos    = 4L,
        gate.marker.x = "PD1",
        gate.marker.y = "*",
        gate.type    = "1dsep",
        algorithm    = "tail",
        gate.param   = list(
            quantile.x = 0,
            sep.axis   = 2L,
            region.min = 0.01,
            region.max = 0.99
        ),
        markers = c( "CCR2", "CCR4", "CD25", "CD28", "CD95", "CD127",
                     "CTLA-4", "HLA-DR", "KI67" ),
        marker.algorithm.override = list(
          CCR4 = "mixture",
          CD28 = "mixture",
          CD95 = "mixture",
          CD127 = "mixture"
        )
    ),

    # display-only rows 91–111
    list( popul.number = 91L,  parent.gate = "thelper.follicular.th1.memory",  parent.popul = "tem",
          popul.label = "Effector Memory Follicular T helper 1-like Cells",    label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 92L,  parent.gate = "thelper.follicular.th1.memory",  parent.popul = "temra",
          popul.label = "Effector Memory RA Follicular T helper 1-like Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 93L,  parent.gate = "thelper.follicular.th1.memory",  parent.popul = "tcm",
          popul.label = "Central Memory Follicular T helper 1-like Cells",     label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 94L,  parent.gate = "thelper.follicular.th2.memory",  parent.popul = "tem",
          popul.label = "Effector Memory Follicular T helper 2-like Cells",    label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 95L,  parent.gate = "thelper.follicular.th2.memory",  parent.popul = "temra",
          popul.label = "Effector Memory RA Follicular T helper 2-like Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 96L,  parent.gate = "thelper.follicular.th2.memory",  parent.popul = "tcm",
          popul.label = "Central Memory Follicular T helper 2-like Cells",     label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 97L,  parent.gate = "thelper.follicular.th17.memory", parent.popul = "tem",
          popul.label = "Effector Memory Follicular T helper 17-like Cells",   label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 98L,  parent.gate = "thelper.follicular.th17.memory", parent.popul = "temra",
          popul.label = "Effector Memory RA Follicular T helper 17-like Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 99L,  parent.gate = "thelper.follicular.th17.memory", parent.popul = "tcm",
          popul.label = "Central Memory Follicular T helper 17-like Cells",    label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 100L, parent.gate = "thelper.follicular.tem.icos",     parent.popul = "quiesc1",
          popul.label = "Quiescent 1 Follicular T helper Cells",               label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 101L, parent.gate = "thelper.follicular.tem.icos",     parent.popul = "quiesc2",
          popul.label = "Quiescent 2 Follicular T helper Cells",               label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 102L, parent.gate = "thelper.follicular.tem.icos",     parent.popul = "rec.act",
          popul.label = "Recently Activated Follicular T helper Cells",        label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 103L, parent.gate = "thelper.follicular.th1.tem.icos", parent.popul = "quiesc1",
          popul.label = "Quiescent 1 Follicular T helper 1-like Cells",        label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 104L, parent.gate = "thelper.follicular.th1.tem.icos", parent.popul = "quiesc2",
          popul.label = "Quiescent 2 Follicular T helper 1-like Cells",        label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 105L, parent.gate = "thelper.follicular.th1.tem.icos", parent.popul = "rec.act",
          popul.label = "Recently Activated Follicular T helper 1-like Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 106L, parent.gate = "thelper.follicular.th2.tem.icos", parent.popul = "quiesc1",
          popul.label = "Quiescent 1 Follicular T helper 2-like Cells",        label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 107L, parent.gate = "thelper.follicular.th2.tem.icos", parent.popul = "quiesc2",
          popul.label = "Quiescent 2 Follicular T helper 2-like Cells",        label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 108L, parent.gate = "thelper.follicular.th2.tem.icos", parent.popul = "rec.act",
          popul.label = "Recently Activated Follicular T helper 2-like Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 109L, parent.gate = "thelper.follicular.th17.tem.icos", parent.popul = "quiesc1",
          popul.label = "Quiescent 1 Follicular T helper 17-like Cells",       label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 110L, parent.gate = "thelper.follicular.th17.tem.icos", parent.popul = "quiesc2",
          popul.label = "Quiescent 2 Follicular T helper 17-like Cells",       label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() ),
    list( popul.number = 111L, parent.gate = "thelper.follicular.th17.tem.icos", parent.popul = "rec.act",
          popul.label = "Recently Activated Follicular T helper 17-like Cells", label.pos = 0L,
          gate.marker.x = "*", gate.marker.y = "*", gate.type = "",
          algorithm = NULL, gate.param = list(), markers = character(0),
          marker.algorithm.override = list() )

)

# Name by parent.gate + "." + parent.popul for convenient lookup
names( fcs.activation.gates ) <- sprintf( "%s.%s",
    sapply( fcs.activation.gates, `[[`, "parent.gate" ),
    sapply( fcs.activation.gates, `[[`, "parent.popul" ) )
