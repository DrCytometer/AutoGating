# Plotting rewrite (gate figures and statistics figures)

## What the figures show

**Gate figures.** Each gate is drawn on the events of its parent population, exactly as it was applied: the thresholds (lines across the plot) or polygons (outlines), and each population's label with its percentage of the parent. Events are coloured by how crowded that part of the plot is (blue sparse, red dense), as in FlowJo's pseudocolour. A threshold gate with no second channel is a histogram. The subtitle gives the parent, its event count, where the boundary came from (trained, fixed, recalculated, reference) and any flag (drift, failed) with its message, so a flagged sample can be checked by eye.

**Axes.** Gates are calculated on the display scale (0 to 1 for the top of the instrument range); the figures use the same scale, but the tick labels are in raw units, worked out from each channel's own transform:
- linear channels (scatter, time) get evenly spaced ticks labelled with K / M / G for whatever range the instrument has (Aurora 0, 1M, 2M, 3M, 4M; Symphony 0, 50K to 250K);
- logicle and log channels get decade ticks (0, 10^3, 10^4, and negative decades); a decade too close to its neighbour to read (for example 10^3 next to 0 on an Aurora logicle with a wide linear region) stays as an unlabelled tick.

Whether a channel is scatter or fluorescence now comes from the panel `type` column, which also sets the axis range (`ap$figure.limits`: `scatter.min/max` for scatter, `expr.min/max` for fluorescence, raw units; the time channel uses its data range). Events beyond the range are drawn on the edge.

**Colour by density.** The colour is the rank of each event's local density: events are counted on a 128 x 128 grid over the plot, the counts are smoothed with a Gaussian of the usual normal-reference width (MASS `bandwidth.nrd` / 4, the same width `kde2d` uses), and each event takes the smoothed density interpolated from the four nearest grid points. Ranking makes the colours use the whole palette on every plot rather than being dominated by the single densest spot. This binned estimate is several times faster than `kde2d` and looks the same at plot resolution.

**Statistics figures** (`write.differential.figures`, called by `differential.statistics` when `output.dir` is set):
- volcano plots, one panel per comparison: difference across (log2 odds ratio for frequencies, difference of medians for markers), evidence up (-log10 p); significant results in red (higher in the second group) or blue (lower), labelled for the eight smallest p-values; dashed lines at the effect-size threshold;
- effect maps: every result significant in at least one comparison (up to 60), one column per comparison, starred where significant;
- group plots: each sample's value for the 24 most significant results by group (frequencies as % of parent, marker medians on the display scale), group median as a bar. Points are spread sideways by rank, not at random, so the figure is the same every time.

## When figures are drawn

- `get.autogating.param(figure.dir = "figures", figure.training = TRUE, figure.samples = FALSE)`.
- `define.gates`: with `figure.training`, one PNG per gate on the pooled training events in `figure.dir/training`.
- `run.autogating`: with `figure.samples`, one PNG per gate and sample in `figure.dir/<sample>` (drawn in the parallel workers). Off by default: ST1 has 415 gates.
- `write.gate.figures(gate.defs, result, data, fcs.setup, ap, dir)` and `gate.figure(...)` (a ggplot) for any gating result.
- Files are `NNN_<gate name>.png` (gate number, zero-padded), written with ragg; at most `ap$figure.max.events` (100,000) events per figure.
- Time per figure about 0.5 s for a scatter plot and 0.15 s for a histogram (6 x 6 in, 150 dpi), so ST1 training figures take about 3 minutes.

## Code changes

- New `R/gate_figures.R`: `axis.ticks`, `figure.axes`, `event.density`, `figure.channels`, `gate.annotations`, `gate.figure` (exported), `write.gate.figures` (exported).
- New `R/differential_figures.R`: `result.labels`, `volcano.figure`, `effect.map`, `group.figure`, `write.differential.figures` (exported). `differential.statistics` gains `figures = TRUE` and returns `group.levels`.
- Removed `R/plot_gate.R`, `R/build_gate_plot.R`, `R/add_gate_boundary.R`, `R/compute_gate_density.R`, and `get.density.palette` / `get.pop.label.pos` from `R/analyze_aux.r` (they read the old gate structure).
- Parameters: removed `ap$gate.defaults` (old break vectors, label sizes, algorithm defaults; nothing reads them) from the minimal and the four cytometer files, and the old palette, label-factor, margin and density-bar settings. Added `figure.max.events`, `figure.resolution`, `figure.label.size`, `figure.text.size`; `figure.dir.basename` replaced by `figure.dir`.
- Seeds: `seeded.call(seed, f)` (`R/analyze_aux.r`) sets the seed, calls `f` and restores the caller's random number stream; used by the training subsample, the bimodality subsample, the FlowJo review subsample and the figure subsample. `set.seed(ap$seed)` removed from `setup.autogating` (nothing there is random).
- S3-looking names renamed after the R CMD check note: `log2.odds` -> `population.log2.odds` (looked like a method of `log2`), `frequency.tables` -> `population.frequency.tables` (`stats::frequency` is a generic). Every function name in `R/` was checked against the S3 generics of base, stats, utils, graphics, grDevices, methods, tools, ggplot2, limma, MASS, jsonlite and xml2; none left.
- `%||%` defined in `R/globals.R` (base R has it only from 4.4.0; the package declares R >= 4.1.0).
- DESCRIPTION Imports gain grDevices, stats and utils (used with `::`).
- `create_experiment_file.R` docs: `gate.calculation` now described as training-sample selection.

## Checks

- Cloud (R 4.3.3, ggplot2 3.4.4): tests in `tests/testthat/test-figures.R` (axis labels for Aurora, Symphony and small ranges; logicle decades and spacing; axis ranges from the panel type; label corners and text for a three-quadrant population; subsample leaves the random stream untouched; figure from a gating result). All test files pass.
- `scripts/figures_check_test.R` draws every gate type on simulated Aurora data (training and per-sample figures) and the statistics figures on three simulated groups; checked by eye in the cloud. Minor ticks need ggplot2 >= 3.5.0 (not drawn on 3.4.4).

Fixes A2 #17 (scatter axes not recognised on Aurora/ID7000 and SSC on Symphony/Fortessa: axes now from the panel type and the channel transform) and #18 (Aurora labels off by 10,000: labels now computed from raw values for any range).
