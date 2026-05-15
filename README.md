# AutoGating GUI — Shiny App

Interactive gate definition editor for the AutoGating R package.

## Quick start

```r
# Install dependencies (once)
install.packages(c(
  "shiny", "bslib", "shinyjs",
  "flowCore", "ggplot2", "ggridges",
  "dplyr", "DT"
))
BiocManager::install("flowCore")   # if not already installed

# Source your AutoGating helpers, then run the app
source("get_autogating_param.r")
source("get_autogating_param_aurora.r")
# ... etc.
shiny::runApp("autogating_gui/")
```

## Tab overview

### 1 · Setup
- Upload one or more FCS files
- Upload the panel CSV (`dye` / `antigen` / `type` columns)
- Select the cytometer (initialises `ap` defaults)
- Optionally load existing population-gate or activation-gate definitions (`.r` or `.rds`)
- Click **Load & Initialise** to build the combined expression matrix and initialise
  the transformation parameter table

### 2 · Transformations
- Ridgeline (ggridges) plot of all fluorescence channels, downsampled to ≤ 30k events
- The highlighted channel (orange) matches the current selector on the left panel
- Adjust **width basis**, **pos**, **neg** per-channel using the slider (draggable)
  and the numeric input field (they stay in sync)
- The population shown in the ridgeline can be switched to any already-defined gate
  population (e.g. to assess transformation on live cells only)
- **Apply to channel** stores the updated parameters without recalculating gates
- **Save transformation CSV** exports `dye,pos,neg,width,length,max.range` — the
  format expected by `read.transformation()` in `analyze_aux.r`

### 3 · Gate Editor
Left panel (gate definition fields):
- Gate name, type (`free` / `1dsep` / `2dsep` / `replicate`)
- Parent gate + parent population (cascade selector)
- X and Y axis markers
- Algorithm — for 1dsep: `tail`, `bimodal`, `kde.min`, `otsu`, `mixture`;
  for 2dsep: x:y pair; for replicate: source gate
- Sub-population names (1dsep: neg/pos; 2dsep: DN/Xpos/Ypos/DP; use `*` to ignore)
- Reference/negative-control gate + population (used to prime the search region)
- Lower / upper search bounds (`region.min` / `region.max`) with slider
- Stats parent (immediate or explicit gate/popul)

Centre panel (plot):
- Scatter plot (hex-bin) for 2D or histogram for 1D / same-marker, downsampled to
  the user-controlled maximum
- **Single click** sets a gate-position hint (drawn as orange dashed line) that is
  passed as `user.hint.x` / `user.hint.y` into the gate parameter list
- **Double-click** clears the hint
- **Calculate gate** calls `calculate.gate.<type>()` and overlays the resulting
  boundary in green — the population is also indexed in `fcs_gates_temp` so that
  downstream gates can use it immediately
- **Add gate to definition** commits the gate to the `fcs.population.gates` list

Right panel: live preview of the R list that will be written on export.

### 4 · Gate Tree
- Table of all defined gates (number, name, type, parent, markers, populations,
  algorithm)
- Select a row to see the full R list representation
- **Delete selected gate** removes it (root gate `all` is protected)

### 5 · Export
- **Population gates** — download as `.r` (R source, `fcs.population.gates <- list(…)`)
  or `.rds`; variable name is configurable
- **Activation gates** — same, for `fcs.activation.gates`
- **Transformation CSV** — table preview + download; matches the format read by
  `read.transformation()` so it can be passed directly as `transformation.filename`
  in `get.autogating.param()`

## Gate definition format

The exported R source files produce a list that is directly compatible with
`load.gate.definitions()` in `analyze_aux.r`. Each element has the same fields
documented in `fcs_population_gates_st1.r`:

```
gate.number, gate.name, gate.type, gate.marker.x, gate.marker.y,
parent.gate, parent.popul, populations, algorithm, gate.param, stats.parent
```

## Integration notes

- Source all AutoGating helpers **before** `shiny::runApp()`. The app will call
  `calculate.gate.*()`, `apply.gate()`, etc. at runtime.
- The app stores a running `fcs_gates_temp` index as gates are calculated, so that
  child-gate plots immediately show only the parent population.
- The `user.hint.x` / `user.hint.y` gate-param fields can optionally be used by a
  custom `calculate.threshold.*()` function to bias the search toward the clicked
  position. The built-in algorithms ignore them unless you modify them to use this
  hint.
  
## Funding
This work has received funding from the KU Leuven C1 program, the European Union’s
Horizon 2020 research and innovation programme under grant agreement No 874707 
(EXIMIOUS), Wellcome Investigator Award, 222442/A/21/Z and UKRI Proactive 
Vaccinology Award, MR/Y004450/1 (IMMPROVE).
