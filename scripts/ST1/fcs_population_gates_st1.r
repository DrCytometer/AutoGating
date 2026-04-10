# fcs_population_gates_st1.r
#
# Population gate definitions for the ST1 panel (BD LSRFortessa).
#
# Each element of fcs.population.gates is a named list describing one gate.
# Fields:
#   gate.number     Integer.  0 = root (no figure produced).
#   gate.name       Character.  Unique key used throughout the pipeline.
#   gate.type       "free" | "1dsep" | "2dsep"
#   gate.marker.x   Antigen name for the x-axis.
#   gate.marker.y   Antigen name for the y-axis.
#   parent.gate     Character (or character vector for multi-parent gates).
#   parent.popul    Character (or character vector, parallel to parent.gate).
#   populations     Named list.  Each element is a sub-population:
#                     label     Display label shown on the figure.
#                     label.pos Integer 1-8 (corner/edge) or 0 (centre).
#                               Use ap$ignore ("*") to suppress.
#   algorithm       Character (1dsep) or length-2 character (2dsep: x then y).
#                   "replicate" copies boundaries from another gate.
#   gate.param      Named list of per-gate overrides to cytometer gate.defaults.
#   stats.parent    Named list(gate, popul).
#
# Conventions:
#   - Gates 3-6 have mixed scatter (x: FSC-A) / fluorescence (y) axes.
#     break.x is overridden to scatter breaks; break.y uses fluorescence default.
#   - density.trimneg.x/y default is TRUE/FALSE.  Blank in CSV = FALSE.
#   - label.size defaults to 4; replicate memory gates use 3 (per CSV).


fcs.population.gates <- list(
  
  # ---- Gate 0: root -------------------------------------------------------
  list(
    gate.number   = 0L,
    gate.name     = "all",
    gate.type     = "free",
    gate.marker.x = "FSC-A",
    gate.marker.y = "FSC-H",
    parent.gate   = NULL,
    parent.popul  = NULL,
    populations   = list(
      all = list( label = "All", label.pos = 0L )
    ),
    algorithm     = "all",
    gate.param    = list( quantile.x = 0, quantile.y = 0 ),
    stats.parent  = list( gate = "*", popul = "*" )
  ),
  
  # ---- Gate 1: singlets ---------------------------------------------------
  list(
    gate.number   = 1L,
    gate.name     = "singlets",
    gate.type     = "free",
    gate.marker.x = "FSC-A",
    gate.marker.y = "FSC-H",
    parent.gate   = "all",
    parent.popul  = "all",
    populations   = list(
      singlets = list( label = "Singlets", label.pos = 4L )
    ),
    algorithm     = "singlets",
    gate.param    = list(
      quantile.x         = 0.0001,
      quantile.y         = 0.0001,
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      density.trimneg.x  = TRUE,
      density.trimneg.y  = TRUE,
      break.x            = c( 0, 50, 100, 150, 200, 250 ),
      break.x.factor     = 1000,
      break.y            = c( 0, 50, 100, 150, 200, 250 ),
      break.y.factor     = 1000,
      limit.xmin         = 0,
      limit.xmax         = 265000,
      limit.ymin         = 0,
      limit.ymax         = 265000
    ),
    stats.parent  = list( gate = "*", popul = "*" )
  ),
  
  # ---- Gate 2: leukocytes -------------------------------------------------
  # sep.axis=1 (FSC-A). density.trimneg.y=TRUE. limit.ymax=0.85.
  list(
    gate.number   = 2L,
    gate.name     = "leukocytes",
    gate.type     = "1dsep",
    gate.marker.x = "FSC-A",
    gate.marker.y = "SSC-A",
    parent.gate   = "singlets",
    parent.popul  = "singlets",
    populations   = list(
      `*`        = list( label = "*",          label.pos = "*" ),
      leukocytes = list( label = "Leukocytes", label.pos = 4L )
    ),
    algorithm     = "mixture",
    gate.param    = list(
      quantile.x         = 0,
      quantile.y         = 0.0001,
      density.quantile.x = 0,
      density.quantile.y = 0.0001,
      density.trimneg.x  = TRUE,
      density.trimneg.y  = TRUE,
      sep.axis           = 1L,
      region.min         = 0.01,
      region.max         = 0.10,
      limit.xmin         = -0.01,
      limit.xmax         = 1.01,
      limit.ymin         = -0.01,
      limit.ymax         = 0.85,
      break.x            = c( 0, 50, 100, 150, 200, 250 ),
      break.x.factor     = 1000,
      break.y            = c( 0, 50, 100, 150, 200, 250 ),
      break.y.factor     = 1000
    ),
    stats.parent  = list( gate = "*", popul = "*" )
  ),
  
  # ---- Gate 3: monocytes --------------------------------------------------
  # x: FSC-A (scatter breaks), y: CD14 (fluorescence default breaks).
  list(
    gate.number   = 3L,
    gate.name     = "monocytes",
    gate.type     = "1dsep",
    gate.marker.x = "FSC-A",
    gate.marker.y = "CD14",
    parent.gate   = "leukocytes",
    parent.popul  = "leukocytes",
    populations   = list(
      no.monocytes = list( label = "No Monocytes", label.pos = 2L ),
      monocytes    = list( label = "Monocytes",    label.pos = 4L )
    ),
    algorithm     = "mixture",
    gate.param    = list(
      quantile.x         = 0,
      density.quantile.x = 0,
      density.trimneg.x  = TRUE,
      sep.axis           = 2L,
      region.min         = 0.01,
      region.max         = 0.99,
      break.x            = c( 0, 50, 100, 150, 200, 250 ),
      break.x.factor     = 1000
    ),
    stats.parent  = list( gate = "*", popul = "*" )
  ),
  
  # ---- Gate 4: live non-monocytes -----------------------------------------
  # x: FSC-A (scatter breaks), y: livedead (fluorescence default breaks).
  # region.min = 0.01.
  list(
    gate.number   = 4L,
    gate.name     = "live.no.monocytes",
    gate.type     = "1dsep",
    gate.marker.x = "FSC-A",
    gate.marker.y = "livedead",
    parent.gate   = "monocytes",
    parent.popul  = "no.monocytes",
    populations   = list(
      live.no.monocytes = list( label = "Live No Monocytes", label.pos = 4L ),
      `*`               = list( label = "*",                  label.pos = "*" )
    ),
    algorithm     = "mixture",
    gate.param    = list(
      quantile.x         = 0,
      density.quantile.x = 0,
      density.trimneg.x  = TRUE,
      sep.axis           = 2L,
      region.min         = 0.10,
      region.max         = 0.80,
      break.x            = c( 0, 50, 100, 150, 200, 250 ),
      break.x.factor     = 1000
    ),
    stats.parent  = list( gate = "*", popul = "*" )
  ),
  
  # ---- Gate 5: live monocytes ---------------------------------------------
  # x: FSC-A (scatter breaks), y: livedead. region.min = 0.30.
  list(
    gate.number   = 5L,
    gate.name     = "live.monocytes",
    gate.type     = "1dsep",
    gate.marker.x = "FSC-A",
    gate.marker.y = "livedead",
    parent.gate   = "monocytes",
    parent.popul  = "monocytes",
    populations   = list(
      live.monocytes = list( label = "Live Monocytes", label.pos = 4L ),
      `*`            = list( label = "*",               label.pos = "*" )
    ),
    algorithm     = "mixture",
    gate.param    = list(
      quantile.x         = 0,
      density.quantile.x = 0,
      density.trimneg.x  = TRUE,
      sep.axis           = 2L,
      region.min         = 0.20,
      region.max         = 0.99,
      break.x            = c( 0, 50, 100, 150, 200, 250 ),
      break.x.factor     = 1000
    ),
    stats.parent  = list( gate = "*", popul = "*" )
  ),
  
  # ---- Gate 6: T cells ----------------------------------------------------
  # x: FSC-A (scatter breaks), y: CD3.
  list(
    gate.number   = 6L,
    gate.name     = "t.cells",
    gate.type     = "1dsep",
    gate.marker.x = "FSC-A",
    gate.marker.y = "CD3",
    parent.gate   = "live.no.monocytes",
    parent.popul  = "live.no.monocytes",
    populations   = list(
      no.t.cells = list( label = "No T Cells", label.pos = 2L ),
      t.cells    = list( label = "T Cells",    label.pos = 4L )
    ),
    algorithm     = "mixture",
    gate.param    = list(
      quantile.x         = 0,
      density.quantile.x = 0,
      density.trimneg.x  = TRUE,
      sep.axis           = 2L,
      region.min         = 0.01,
      region.max         = 0.99,
      break.x            = c( 0, 50, 100, 150, 200, 250 ),
      break.x.factor     = 1000
    ),
    stats.parent  = list(
      gate  = c( "parent.1", "leukocytes", "singlets" ),
      popul = c( "parent.1", "leukocytes", "singlets" )
    )
  ),
  
  # ---- Gate 7: CD4 / CD8 --------------------------------------------------
  # quantile.y=0.0005. density.quantile.x/y=0.0001.
  list(
    gate.number   = 7L,
    gate.name     = "cd4.cd8",
    gate.type     = "2dsep",
    gate.marker.x = "CD8",
    gate.marker.y = "CD4",
    parent.gate   = "t.cells",
    parent.popul  = "t.cells",
    populations   = list(
      double.neg = list( label = "CD4-/CD8- T Cells",  label.pos = 1L ),
      cd8.pos    = list( label = "CD8+ T Cells",        label.pos = 2L ),
      cd4.pos    = list( label = "CD4+ T Cells",        label.pos = 3L ),
      double.pos = list( label = "CD4+/CD8+ T Cells",  label.pos = 4L )
    ),
    algorithm     = c( "mixture", "mixture" ),
    gate.param    = list(
      quantile.x         = 0.0001,
      quantile.y         = 0.0005,
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      region.xmin        = 0.01,
      region.xmax        = 0.99,
      region.ymin        = 0.01,
      region.ymax        = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 8: CD4 memory subsets -----------------------------------------
  # quantile.x/y=0.0001. density.quantile.x/y=0.0001. label.size=4.
  # stats.parent = cd4.cd8.
  list(
    gate.number   = 8L,
    gate.name     = "cd4.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "cd4.cd8",
    parent.popul  = "cd4.pos",
    populations   = list(
      tem   = list( label = "Effector Memory CD4+ T Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA CD4+ T Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory CD4+ T Cells",     label.pos = 3L ),
      naive = list( label = "Naive CD4+ T Cells",              label.pos = 4L )
    ),
    algorithm     = c( "mixture", "mixture" ),
    gate.param    = list(
      quantile.x         = 0.0001,
      quantile.y         = 0.0001,
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      region.max         = 0.3
    ),
    stats.parent  = list( gate = "cd4.cd8", popul = "*" )
  ),
  
  # ---- Gate 9: CD4 non-naive (joined L-shaped boundary) -------------------
  # quantile.x/y=0.0001. density.quantile.x/y=0.0001. label.size=4.
  list(
    gate.number   = 9L,
    gate.name     = "cd4.no.naive",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "cd4.cd8",
    parent.popul  = "cd4.pos",
    populations   = list(
      no.naive = list( label = "No Naive CD4+ T Cells", label.pos = 3L )
    ),
    algorithm     = c( "mixture", "mixture" ),
    gate.param    = list(
      join               = "1:2:3",
      quantile.x         = 0.0001,
      quantile.y         = 0.0001,
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      region.xmin        = 0.01,
      region.xmax        = 0.99,
      region.ymin        = 0.01,
      region.ymax        = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 10: CD8 memory subsets ----------------------------------------
  # algorithm = mixture;tail. threshold.param=0.999.
  # density.quantile.x/y=0.0001. label.size=4. stats.parent=cd4.cd8.
  list(
    gate.number   = 10L,
    gate.name     = "cd8.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "cd4.cd8",
    parent.popul  = "cd8.pos",
    populations   = list(
      tem   = list( label = "Effector Memory CD8+ T Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA CD8+ T Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory CD8+ T Cells",     label.pos = 3L ),
      naive = list( label = "Naive CD8+ T Cells",              label.pos = 4L )
    ),
    algorithm     = c( "mixture", "tail" ),
    gate.param    = list(
      quantile.x         = 0.0001,
      quantile.y         = 0.0001,
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      threshold.param    = 0.999
    ),
    stats.parent  = list( gate = "cd4.cd8", popul = "*" )
  ),
  
  # ---- Gate 11: CD4 naive CD31 (RTE) --------------------------------------
  # density.trimneg.x/y=FALSE (CSV blank). region.max=0.30. label.size=4.
  # stats.parent=cd4.cd8.
  list(
    gate.number   = 11L,
    gate.name     = "cd4.naive.cd31",
    gate.type     = "1dsep",
    gate.marker.x = "CD4",
    gate.marker.y = "CD31",
    parent.gate   = "cd4.memory",
    parent.popul  = "naive",
    populations   = list(
      `*`  = list( label = "*",                                   label.pos = "*" ),
      rte  = list( label = "Recent Thymic Emigrant CD4+ T Cells", label.pos = 3L )
    ),
    algorithm     = "mixture",
    gate.param    = list(
      quantile.x         = 0,
      density.quantile.x = 0,
      density.trimneg.x  = FALSE,
      sep.axis           = 2L,
      region.min         = 0.01,
      region.max         = 0.30
    ),
    stats.parent  = list( gate = "cd4.cd8", popul = "*" )
  ),
  
  # ---- Gate 12: CD4 naive CD95 (stem cell memory) -------------------------
  # density.trimneg.x/y=FALSE. label.size=4. stats.parent=cd4.cd8.
  list(
    gate.number   = 12L,
    gate.name     = "cd4.naive.cd95",
    gate.type     = "1dsep",
    gate.marker.x = "CD4",
    gate.marker.y = "CD95",
    parent.gate   = "cd4.memory",
    parent.popul  = "naive",
    populations   = list(
      `*`  = list( label = "*",                        label.pos = "*" ),
      scm  = list( label = "Stem Memory CD4+ T Cells", label.pos = 3L )
    ),
    algorithm     = "tail",
    gate.param    = list(
      quantile.x         = 0,
      density.quantile.x = 0,
      density.trimneg.x  = FALSE,
      sep.axis           = 2L,
      region.min         = 0.01,
      region.max         = 0.99
    ),
    stats.parent  = list( gate = "cd4.cd8", popul = "*" )
  ),
  
  # ---- Gate 13: CD8 naive CD31 (RTE) --------------------------------------
  # density.trimneg.x/y=FALSE. label.size=4.
  list(
    gate.number   = 13L,
    gate.name     = "cd8.naive.cd31",
    gate.type     = "1dsep",
    gate.marker.x = "CD8",
    gate.marker.y = "CD31",
    parent.gate   = "cd8.memory",
    parent.popul  = "naive",
    populations   = list(
      `*`  = list( label = "*",                                   label.pos = "*" ),
      rte  = list( label = "Recent Thymic Emigrant CD8+ T Cells", label.pos = 3L )
    ),
    algorithm     = "mixture",
    gate.param    = list(
      quantile.x         = 0,
      density.quantile.x = 0,
      density.trimneg.x  = FALSE,
      sep.axis           = 2L,
      region.min         = 0.01,
      region.max         = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 14: CD8 naive CD95 (stem cell memory) -------------------------
  # density.trimneg.x/y=FALSE. label.size=4.
  list(
    gate.number   = 14L,
    gate.name     = "cd8.naive.cd95",
    gate.type     = "1dsep",
    gate.marker.x = "CD8",
    gate.marker.y = "CD95",
    parent.gate   = "cd8.memory",
    parent.popul  = "naive",
    populations   = list(
      `*`  = list( label = "*",                        label.pos = "*" ),
      scm  = list( label = "Stem Memory CD8+ T Cells", label.pos = 3L )
    ),
    algorithm     = "tail",
    gate.param    = list(
      quantile.x         = 0,
      density.quantile.x = 0,
      density.trimneg.x  = FALSE,
      sep.axis           = 2L,
      region.min         = 0.01,
      region.max         = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 15: CD4 Treg --------------------------------------------------
  # sep.axis=1. quantile.x/y=0.0001. density.quantile.x/y=0.0001.
  # region.min=0.25.
  list(
    gate.number   = 15L,
    gate.name     = "cd4.treg",
    gate.type     = "1dsep",
    gate.marker.x = "FOXP3",
    gate.marker.y = "CD25",
    parent.gate   = "cd4.cd8",
    parent.popul  = "cd4.pos",
    populations   = list(
      no.treg = list( label = "No Treg CD4+ T Cells",  label.pos = 3L ),
      treg    = list( label = "Regulatory T Cells",     label.pos = 4L )
    ),
    algorithm     = "mixture",
    gate.param    = list(
      quantile.x         = 0.0001,
      quantile.y         = 0.0001,
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      sep.axis           = 1L,
      region.min         = 0.25,
      region.max         = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 16: CD4 non-Treg memory (replicate) ---------------------------
  # label.size=3. density.quantile.x/y=0.0001. stats.parent=cd4.treg.
  list(
    gate.number   = 16L,
    gate.name     = "cd4.no.treg.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "cd4.treg",
    parent.popul  = "no.treg",
    populations   = list(
      tem   = list( label = "Effector Memory No Treg CD4+ T Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA No Treg CD4+ T Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory No Treg CD4+ T Cells",     label.pos = 3L ),
      naive = list( label = "Naive No Treg CD4+ T Cells",              label.pos = 4L )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = "cd4.treg", popul = "*" )
  ),
  
  # ---- Gate 17: CD4 non-Treg naive CD31 (replicate) -----------------------
  # label.size=3.
  list(
    gate.number   = 17L,
    gate.name     = "cd4.no.treg.naive.cd31",
    gate.type     = "1dsep",
    gate.marker.x = "CD4",
    gate.marker.y = "CD31",
    parent.gate   = "cd4.no.treg.memory",
    parent.popul  = "naive",
    populations   = list(
      rte = list( label = "Recent Thymic Emigrant No Treg CD4+ T Cells", label.pos = 3L )
    ),
    algorithm     = "replicate",
    gate.param    = list( origin.gate = "cd4.naive.cd31", label.size = 3 ),
    stats.parent  = list( gate = c( "parent.1", "parent.2" ), popul = c( "parent.1", "parent.2" ) )
  ),
  
  # ---- Gate 18: CD4 non-Treg naive CD95 (replicate) -----------------------
  # label.size=3.
  list(
    gate.number   = 18L,
    gate.name     = "cd4.no.treg.naive.cd95",
    gate.type     = "1dsep",
    gate.marker.x = "CD4",
    gate.marker.y = "CD95",
    parent.gate   = "cd4.no.treg.memory",
    parent.popul  = "naive",
    populations   = list(
      scm = list( label = "Stem Memory No Treg CD4+ T Cells", label.pos = 3L )
    ),
    algorithm     = "replicate",
    gate.param    = list( origin.gate = "cd4.naive.cd95", label.size = 3 ),
    stats.parent  = list( gate = c( "parent.1", "parent.2" ), popul = c( "parent.1", "parent.2" ) )
  ),
  
  # ---- Gate 19: CD4 Treg memory (replicate) -------------------------------
  # label.size=3. density.quantile.x/y=0.0001. stats.parent=cd4.treg.
  list(
    gate.number   = 19L,
    gate.name     = "cd4.treg.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "cd4.treg",
    parent.popul  = "treg",
    populations   = list(
      tem   = list( label = "Effector Memory Treg Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA Treg Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory Treg Cells",     label.pos = 3L ),
      naive = list( label = "Naive Treg Cells",              label.pos = 4L )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = "cd4.treg", popul = "*" )
  ),
  
  # ---- Gate 20: CD4 Treg non-naive (replicate) ----------------------------
  # label.size=3. density.quantile.x/y=0.0001.
  list(
    gate.number   = 20L,
    gate.name     = "cd4.treg.no.naive",
    gate.type     = "free",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "cd4.treg",
    parent.popul  = "treg",
    populations   = list(
      no.naive = list( label = "No Naive Treg Cells", label.pos = 3L )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.no.naive",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 21: CD4 Treg naive CD31 (replicate) ---------------------------
  # label.size=3.
  list(
    gate.number   = 21L,
    gate.name     = "cd4.treg.naive.cd31",
    gate.type     = "1dsep",
    gate.marker.x = "CD4",
    gate.marker.y = "CD31",
    parent.gate   = "cd4.treg.memory",
    parent.popul  = "naive",
    populations   = list(
      rte = list( label = "Recent Thymic Emigrant Treg Cells", label.pos = 3L )
    ),
    algorithm     = "replicate",
    gate.param    = list( origin.gate = "cd4.naive.cd31", label.size = 3 ),
    stats.parent  = list( gate = c( "parent.1", "parent.2" ), popul = c( "parent.1", "parent.2" ) )
  ),
  
  # ---- Gate 22: CD4 Treg naive CD95 (replicate) ---------------------------
  # label.size=3.
  list(
    gate.number   = 22L,
    gate.name     = "cd4.treg.naive.cd95",
    gate.type     = "1dsep",
    gate.marker.x = "CD4",
    gate.marker.y = "CD95",
    parent.gate   = "cd4.treg.memory",
    parent.popul  = "naive",
    populations   = list(
      scm = list( label = "Stem Memory Treg Cells", label.pos = 3L )
    ),
    algorithm     = "replicate",
    gate.param    = list( origin.gate = "cd4.naive.cd95", label.size = 3 ),
    stats.parent  = list( gate = c( "parent.1", "parent.2" ), popul = c( "parent.1", "parent.2" ) )
  ),
  
  # ---- Gate 23: Treg follicular -------------------------------------------
  # Parent: cd4.treg.no.naive / no.naive. 2dsep tail;tail.
  # threshold.param=0.95. 4 pops: *;*;*;follicular.
  # quantile.x/y=0.0001. density.quantile.x/y=0.0001.
  list(
    gate.number   = 23L,
    gate.name     = "treg.follicular",
    gate.type     = "2dsep",
    gate.marker.x = "CXCR5",
    gate.marker.y = "PD1",
    parent.gate   = "cd4.treg.no.naive",
    parent.popul  = "no.naive",
    populations   = list(
      `*`        = list( label = "*",                     label.pos = "*" ),
      `*`        = list( label = "*",                     label.pos = "*" ),
      `*`        = list( label = "*",                     label.pos = "*" ),
      follicular = list( label = "Follicular Treg Cells", label.pos = 3L )
    ),
    algorithm     = c( "tail", "tail" ),
    gate.param    = list(
      quantile.x         = 0.0001,
      quantile.y         = 0.0001,
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      threshold.param    = 0.95,
      region.xmin        = 0.01,
      region.xmax        = 0.99,
      region.ymin        = 0.01,
      region.ymax        = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 24: Treg follicular memory (replicate) ------------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 24L,
    gate.name     = "treg.follicular.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "treg.follicular",
    parent.popul  = "follicular",
    populations   = list(
      tem   = list( label = "Effector Memory Follicular Treg Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA Follicular Treg Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory Follicular Treg Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                         label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 25: Treg CD127 ------------------------------------------------
  # Parent: cd4.treg / treg. Markers: CD25 (x), CD127 (y). sep.axis=2.
  # region.min=0.20. density.quantile.x=0. density.trimneg.x/y=FALSE.
  # stats.parent = parent.1 + parent.2.
  list(
    gate.number   = 25L,
    gate.name     = "treg.cd127",
    gate.type     = "1dsep",
    gate.marker.x = "CD25",
    gate.marker.y = "CD127",
    parent.gate   = "cd4.treg",
    parent.popul  = "treg",
    populations   = list(
      cd127.neg = list( label = "CD127- Treg Cells", label.pos = 2L ),
      `*`       = list( label = "*",                  label.pos = "*" )
    ),
    algorithm     = "mixture",
    gate.param    = list(
      quantile.x         = 0,
      density.quantile.x = 0,
      density.trimneg.x  = FALSE,
      sep.axis           = 2L,
      region.min         = 0.20,
      region.max         = 0.99
    ),
    stats.parent  = list( gate = c( "parent.1", "parent.2" ), popul = c( "parent.1", "parent.2" ) )
  ),
  
  # ---- Gate 26: Treg CD127 memory (replicate) -----------------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 26L,
    gate.name     = "treg.cd127.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "treg.cd127",
    parent.popul  = "cd127.neg",
    populations   = list(
      tem   = list( label = "Effector Memory CD127- Treg Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA CD127- Treg Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory CD127- Treg Cells",     label.pos = 3L ),
      naive = list( label = "Naive CD127- Treg Cells",              label.pos = 4L )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 27: Treg CD127 non-naive (replicate) --------------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 27L,
    gate.name     = "treg.cd127.no.naive",
    gate.type     = "free",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "treg.cd127",
    parent.popul  = "cd127.neg",
    populations   = list(
      no.naive = list( label = "No Naive CD127- Treg Cells", label.pos = 3L )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.no.naive",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 28: Treg CD127 naive CD31 (replicate) -------------------------
  # label.size=3.
  list(
    gate.number   = 28L,
    gate.name     = "treg.cd127.naive.cd31",
    gate.type     = "1dsep",
    gate.marker.x = "CD4",
    gate.marker.y = "CD31",
    parent.gate   = "treg.cd127.memory",
    parent.popul  = "naive",
    populations   = list(
      rte = list( label = "Recent Thymic Emigrant CD127- Treg Cells", label.pos = 3L )
    ),
    algorithm     = "replicate",
    gate.param    = list( origin.gate = "cd4.naive.cd31", label.size = 3 ),
    stats.parent  = list( gate = c( "parent.1", "parent.2" ), popul = c( "parent.1", "parent.2" ) )
  ),
  
  # ---- Gate 29: Treg CD127 naive CD95 (replicate) -------------------------
  # label.size=3.
  list(
    gate.number   = 29L,
    gate.name     = "treg.cd127.naive.cd95",
    gate.type     = "1dsep",
    gate.marker.x = "CD4",
    gate.marker.y = "CD95",
    parent.gate   = "treg.cd127.memory",
    parent.popul  = "naive",
    populations   = list(
      scm = list( label = "Stem Memory CD127- Treg Cells", label.pos = 3L )
    ),
    algorithm     = "replicate",
    gate.param    = list( origin.gate = "cd4.naive.cd95", label.size = 3 ),
    stats.parent  = list( gate = c( "parent.1", "parent.2" ), popul = c( "parent.1", "parent.2" ) )
  ),
  
  # ---- Gate 30: Treg CD127 follicular (replicate from treg.follicular) ----
  # density.quantile.x/y=0.0001. label.size=4.
  list(
    gate.number   = 30L,
    gate.name     = "treg.cd127.follicular",
    gate.type     = "2dsep",
    gate.marker.x = "CXCR5",
    gate.marker.y = "PD1",
    parent.gate   = "treg.cd127.no.naive",
    parent.popul  = "no.naive",
    populations   = list(
      follicular = list( label = "Follicular CD127- Treg Cells", label.pos = 3L )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "treg.follicular",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 31: Treg CD127 follicular memory (replicate) ------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 31L,
    gate.name     = "treg.cd127.follicular.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "treg.cd127.follicular",
    parent.popul  = "follicular",
    populations   = list(
      tem   = list( label = "Effector Memory Follicular CD127- Treg Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA Follicular CD127- Treg Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory Follicular CD127- Treg Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                                label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 32: CD4 T helper (non-Treg, non-naive; replicate) -------------
  # density.quantile.x/y=0.0001. label.size=4.
  list(
    gate.number   = 32L,
    gate.name     = "cd4.thelper",
    gate.type     = "free",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "cd4.treg",
    parent.popul  = "no.treg",
    populations   = list(
      thelper = list( label = "T Helper Cells", label.pos = 3L )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.no.naive",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 33: T helper CXCR3 / CCR4 ------------------------------------
  # quantile.x=0.0001, quantile.y=0.002. density.quantile.x/y=0.0001.
  # 4 pops: *;cxcr3.pos;ccr4.pos;*. algorithm=tail;mixture.
  list(
    gate.number   = 33L,
    gate.name     = "thelper.cxcr3.ccr4",
    gate.type     = "2dsep",
    gate.marker.x = "CXCR3",
    gate.marker.y = "CCR4",
    parent.gate   = "cd4.thelper",
    parent.popul  = "thelper",
    populations   = list(
      `*`       = list( label = "*",                            label.pos = "*" ),
      cxcr3.pos = list( label = "CCR4- CXCR3+ T Helper Cells", label.pos = 2L ),
      ccr4.pos  = list( label = "CCR4+ CXCR3- T Helper Cells", label.pos = 4L ),
      `*`       = list( label = "*",                            label.pos = "*" )
    ),
    algorithm     = c( "mixture", "mixture" ),
    gate.param    = list(
      quantile.x         = 0.0001,
      quantile.y         = 0.002,
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      region.xmin        = 0.01,
      region.xmax        = 0.99,
      region.ymin        = 0.01,
      region.ymax        = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 34: T helper 1 / 1* (CXCR3+, RORgT axis) ---------------------
  # quantile.y=0.0005. density.quantile.x=0. density.trimneg.x/y=FALSE.
  # sep.axis=2.
  list(
    gate.number   = 34L,
    gate.name     = "thelper.cxcr3.rorgt",
    gate.type     = "1dsep",
    gate.marker.x = "CXCR3",
    gate.marker.y = "RORgt",
    parent.gate   = "thelper.cxcr3.ccr4",
    parent.popul  = "cxcr3.pos",
    populations   = list(
      th1     = list( label = "T Helper 1 Cells",  label.pos = 1L ),
      th1.ast = list( label = "T Helper 1* Cells", label.pos = 3L )
    ),
    algorithm     = "tail",
    gate.param    = list(
      quantile.x         = 0,
      quantile.y         = 0.0005,
      density.quantile.x = 0,
      density.trimneg.x  = FALSE,
      sep.axis           = 2L,
      region.min         = 0.01,
      region.max         = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 35: Th1 memory (replicate) ------------------------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 35L,
    gate.name     = "thelper.th1.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "thelper.cxcr3.rorgt",
    parent.popul  = "th1",
    populations   = list(
      tem   = list( label = "Effector Memory T Helper 1 Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA T Helper 1 Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory T Helper 1 Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                    label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 36: Th1* memory (replicate) -----------------------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 36L,
    gate.name     = "thelper.th1.ast.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "thelper.cxcr3.rorgt",
    parent.popul  = "th1.ast",
    populations   = list(
      tem   = list( label = "Effector Memory T Helper 1* Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA T Helper 1* Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory T Helper 1* Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                     label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 37: T helper 2 / 17 (CCR4+, RORgT axis) ----------------------
  # quantile.y=0.0005. density.quantile.x=0. density.trimneg.x/y=FALSE.
  # sep.axis=2.
  list(
    gate.number   = 37L,
    gate.name     = "thelper.ccr4.rorgt",
    gate.type     = "1dsep",
    gate.marker.x = "CCR4",
    gate.marker.y = "RORgt",
    parent.gate   = "thelper.cxcr3.ccr4",
    parent.popul  = "ccr4.pos",
    populations   = list(
      th2  = list( label = "T Helper 2 Cells",  label.pos = 1L ),
      th17 = list( label = "T Helper 17 Cells", label.pos = 3L )
    ),
    algorithm     = "tail",
    gate.param    = list(
      quantile.x         = 0,
      quantile.y         = 0.0005,
      density.quantile.x = 0,
      density.trimneg.x  = FALSE,
      sep.axis           = 2L,
      region.min         = 0.01,
      region.max         = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 38: Th2 memory (replicate) ------------------------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 38L,
    gate.name     = "thelper.th2.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "thelper.ccr4.rorgt",
    parent.popul  = "th2",
    populations   = list(
      tem   = list( label = "Effector Memory T Helper 2 Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA T Helper 2 Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory T Helper 2 Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                    label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 39: Th17 memory (replicate) -----------------------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 39L,
    gate.name     = "thelper.th17.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "thelper.ccr4.rorgt",
    parent.popul  = "th17",
    populations   = list(
      tem   = list( label = "Effector Memory T Helper 17 Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA T Helper 17 Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory T Helper 17 Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                     label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 40: Th17 CCR2 subsets (classic / dual) ------------------------
  # quantile.y=0.0001. density.quantile.x=0. density.trimneg.x/y=FALSE.
  # sep.axis=2.
  list(
    gate.number   = 40L,
    gate.name     = "thelper.th17.ccr2",
    gate.type     = "1dsep",
    gate.marker.x = "CCR4",
    gate.marker.y = "CCR2",
    parent.gate   = "thelper.ccr4.rorgt",
    parent.popul  = "th17",
    populations   = list(
      classic = list( label = "Classic T Helper 17 Cells", label.pos = 1L ),
      dual    = list( label = "Dual T Helper 17 Cells",    label.pos = 3L )
    ),
    algorithm     = "tail",
    gate.param    = list(
      quantile.x         = 0,
      quantile.y         = 0.0001,
      density.quantile.x = 0,
      density.trimneg.x  = FALSE,
      sep.axis           = 2L,
      region.min         = 0.01,
      region.max         = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 41: Classic Th17 memory (replicate) ---------------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 41L,
    gate.name     = "thelper.th17.classic.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "thelper.th17.ccr2",
    parent.popul  = "classic",
    populations   = list(
      tem   = list( label = "Effector Memory Classic T Helper 17 Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA Classic T Helper 17 Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory Classic T Helper 17 Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                              label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 42: Dual Th17 memory (replicate) ------------------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 42L,
    gate.name     = "thelper.th17.dual.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "thelper.th17.ccr2",
    parent.popul  = "dual",
    populations   = list(
      tem   = list( label = "Effector Memory Dual T Helper 17 Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA Dual T Helper 17 Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory Dual T Helper 17 Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                          label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 43: T follicular helper (Tfh) ---------------------------------
  # sep.axis=1. quantile.x/y=0.0001. density.quantile.x/y=0.0001.
  list(
    gate.number   = 43L,
    gate.name     = "thelper.follicular",
    gate.type     = "1dsep",
    gate.marker.x = "CXCR5",
    gate.marker.y = "PD1",
    parent.gate   = "cd4.thelper",
    parent.popul  = "thelper",
    populations   = list(
      `*`        = list( label = "*",                         label.pos = "*" ),
      follicular = list( label = "Follicular T helper Cells", label.pos = 3L )
    ),
    algorithm     = "tail",
    gate.param    = list(
      quantile.x         = 0.0001,
      quantile.y         = 0.0001,
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      sep.axis           = 1L,
      region.min         = 0.01,
      region.max         = 0.99
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 44: Tfh memory (replicate) ------------------------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 44L,
    gate.name     = "thelper.follicular.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "thelper.follicular",
    parent.popul  = "follicular",
    populations   = list(
      tem   = list( label = "Effector Memory Follicular T helper Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA Follicular T helper Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory Follicular T helper Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                             label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 45: Tfh CXCR3 / RORgT subsets --------------------------------
  # quantile.x/y=0.0001. density.quantile.x/y=0.0001. label.size=3.
  # algorithm=bimodal;tail. 4 pops: th2;th1;th17;*.
  list(
    gate.number   = 45L,
    gate.name     = "thelper.follicular.cxcr3.rorgt",
    gate.type     = "2dsep",
    gate.marker.x = "CXCR3",
    gate.marker.y = "RORgt",
    parent.gate   = "thelper.follicular",
    parent.popul  = "follicular",
    populations   = list(
      th2  = list( label = "Follicular T helper 2-like Cells",  label.pos = 1L ),
      th1  = list( label = "Follicular T helper 1-like Cells",  label.pos = 2L ),
      th17 = list( label = "Follicular T helper 17-like Cells", label.pos = 3L ),
      `*`  = list( label = "*",                                  label.pos = "*" )
    ),
    algorithm     = c( "mixture", "tail" ),
    gate.param    = list(
      quantile.x         = 0.0001,
      quantile.y         = 0.0001,
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      region.xmin        = 0.01,
      region.xmax        = 0.99,
      region.ymin        = 0.01,
      region.ymax        = 0.99,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gates 46–48: Tfh subset memories (replicate) -----------------------
  # density.quantile.x/y=0.0001. label.size=3.
  list(
    gate.number   = 46L,
    gate.name     = "thelper.follicular.th1.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "thelper.follicular.cxcr3.rorgt",
    parent.popul  = "th1",
    populations   = list(
      tem   = list( label = "Effector Memory Follicular T helper 1-like Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA Follicular T helper 1-like Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory Follicular T helper 1-like Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                                    label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  list(
    gate.number   = 47L,
    gate.name     = "thelper.follicular.th2.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "thelper.follicular.cxcr3.rorgt",
    parent.popul  = "th2",
    populations   = list(
      tem   = list( label = "Effector Memory Follicular T helper 2-like Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA Follicular T helper 2-like Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory Follicular T helper 2-like Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                                    label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  list(
    gate.number   = 48L,
    gate.name     = "thelper.follicular.th17.memory",
    gate.type     = "2dsep",
    gate.marker.x = "CD45RA",
    gate.marker.y = "CCR7",
    parent.gate   = "thelper.follicular.cxcr3.rorgt",
    parent.popul  = "th17",
    populations   = list(
      tem   = list( label = "Effector Memory Follicular T helper 17-like Cells",    label.pos = 1L ),
      temra = list( label = "Effector Memory RA Follicular T helper 17-like Cells", label.pos = 2L ),
      tcm   = list( label = "Central Memory Follicular T helper 17-like Cells",     label.pos = 3L ),
      `*`   = list( label = "*",                                                     label.pos = "*" )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "cd4.memory",
      density.quantile.x = 0.0001,
      density.quantile.y = 0.0001,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gate 49: Tfh EM ICOS subsets ---------------------------------------
  # quantile.x=0.0001, quantile.y=0. density.quantile.x=0.0001,
  # density.quantile.y=0. label.size=3. algorithm=tail;tail.
  # 4 pops: quiesc1;*;quiesc2;rec.act. Multi-parent: tem+temra.
  list(
    gate.number   = 49L,
    gate.name     = "thelper.follicular.tem.icos",
    gate.type     = "2dsep",
    gate.marker.x = "ICOS",
    gate.marker.y = "PD1",
    parent.gate   = c( "thelper.follicular.memory", "thelper.follicular.memory" ),
    parent.popul  = c( "tem", "temra" ),
    populations   = list(
      quiesc1 = list( label = "Quiescent 1 Follicular T helper Cells",        label.pos = 1L ),
      `*`     = list( label = "*",                                             label.pos = "*" ),
      quiesc2 = list( label = "Quiescent 2 Follicular T helper Cells",        label.pos = 3L ),
      rec.act = list( label = "Recently Activated Follicular T helper Cells", label.pos = 4L )
    ),
    algorithm     = c( "tail", "tail" ),
    gate.param    = list(
      quantile.x         = 0.0001,
      quantile.y         = 0,
      density.quantile.x = 0.0001,
      density.quantile.y = 0,
      region.xmin        = 0.01,
      region.xmax        = 0.99,
      region.ymin        = 0.01,
      region.ymax        = 0.99,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  # ---- Gates 50–52: Tfh subset EM ICOS (replicate) ------------------------
  # density.quantile.x=0.0001, density.quantile.y=0. label.size=3.
  list(
    gate.number   = 50L,
    gate.name     = "thelper.follicular.th1.tem.icos",
    gate.type     = "2dsep",
    gate.marker.x = "ICOS",
    gate.marker.y = "PD1",
    parent.gate   = c( "thelper.follicular.th1.memory", "thelper.follicular.th1.memory" ),
    parent.popul  = c( "tem", "temra" ),
    populations   = list(
      quiesc1 = list( label = "Quiescent 1 Follicular T helper 1-like Cells",        label.pos = 1L ),
      quiesc2 = list( label = "Quiescent 2 Follicular T helper 1-like Cells",        label.pos = 3L ),
      rec.act = list( label = "Recently Activated Follicular T helper 1-like Cells", label.pos = 4L )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "thelper.follicular.tem.icos",
      density.quantile.x = 0.0001,
      density.quantile.y = 0,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  list(
    gate.number   = 51L,
    gate.name     = "thelper.follicular.th2.tem.icos",
    gate.type     = "2dsep",
    gate.marker.x = "ICOS",
    gate.marker.y = "PD1",
    parent.gate   = c( "thelper.follicular.th2.memory", "thelper.follicular.th2.memory" ),
    parent.popul  = c( "tem", "temra" ),
    populations   = list(
      quiesc1 = list( label = "Quiescent 1 Follicular T helper 2-like Cells",        label.pos = 1L ),
      quiesc2 = list( label = "Quiescent 2 Follicular T helper 2-like Cells",        label.pos = 3L ),
      rec.act = list( label = "Recently Activated Follicular T helper 2-like Cells", label.pos = 4L )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "thelper.follicular.tem.icos",
      density.quantile.x = 0.0001,
      density.quantile.y = 0,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  ),
  
  list(
    gate.number   = 52L,
    gate.name     = "thelper.follicular.th17.tem.icos",
    gate.type     = "2dsep",
    gate.marker.x = "ICOS",
    gate.marker.y = "PD1",
    parent.gate   = c( "thelper.follicular.th17.memory", "thelper.follicular.th17.memory" ),
    parent.popul  = c( "tem", "temra" ),
    populations   = list(
      quiesc1 = list( label = "Quiescent 1 Follicular T helper 17-like Cells",        label.pos = 1L ),
      quiesc2 = list( label = "Quiescent 2 Follicular T helper 17-like Cells",        label.pos = 3L ),
      rec.act = list( label = "Recently Activated Follicular T helper 17-like Cells", label.pos = 4L )
    ),
    algorithm     = "replicate",
    gate.param    = list(
      origin.gate        = "thelper.follicular.tem.icos",
      density.quantile.x = 0.0001,
      density.quantile.y = 0,
      label.size         = 3
    ),
    stats.parent  = list( gate = NULL, popul = NULL )
  )
  
)

# Name the list by gate.name for convenient lookup
names( fcs.population.gates ) <- sapply( fcs.population.gates, `[[`, "gate.name" )