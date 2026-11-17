# Gate engine and algorithms (2026-11-11)

Plan: project doc `claude/AUTOGATING_R_REVIEW_PLAN.md`, Part E Phase 2. Assumes the data-layer changes (`CHANGES_DATA_LAYER.md`) landed.

Ported from a snapshot of Honeychrome `ag_core.py` (working tree of honeychrome-sandbox at 2026-11-11 12:01 UTC, uncommitted edits included: MAD-based singlet band, mixture debris cut), saved at `dev/reference/ag_core.py` (sha1 740779fd1a2d1736b3ed576fefe2d65b66776e67). The validation script imports this copy.

New files
- `R/gate_model.R`: vocabulary, regions (`region.quadrants`, `region.sign`, `quadrant.mask`, `region.polygon`, `rect.polygon`), gate-definition helpers (`axis.algorithm`, `replicate.sources`, `gate.dependencies`, `ordered.gate.defs`, `effective.gate.def`, `gate.channels`, `separation.axis`, `population.regions`, shape helpers), `viability.channels` (AutoSpectral fluorophore database via `match.fluorophores`, optional panel dye lookup), `anchor.population`.
- `R/threshold_algorithms.R`: `trim.events`, `crop.to.region`, `kde.grid` (exact up to 20,000 events, binned FFT convolution above), `kde.mode`, `find.peaks`, `threshold.tail`, `threshold.kde.min`, `threshold.otsu`, `fit.mixture` (EM, unequal variances, k-means from quantile starts, tol 1e-3 on mean log-likelihood, variance floor 1e-6, BIC), `threshold.mixture` (2-3 components, outer/lower/upper cut), `compute.threshold`, `threshold.factor`, `nudge.threshold`, `bimodality.delta.bic`, `recommend.algorithm`.
- `R/gate_calculation.R`: `gate.context`, `calculate.gate` (singlets, 1dsep, 2dsep, free: cells/all/ellipse, imported, replicate), `threshold.boundaries`, `debris.cut`, `kde.extent`, `saturation.cut`, `convex.hull.polygon`, `ellipse.polygon`.
- `R/gate_application.R`: `inside.polygon` (even-odd), `edge.clamped`, `population.masks`, `shape.masks`, `resolve.replicates`, `boundary.drift`, `run.gating` (fixed / recalculate with drift), `train.boundaries`, `population.table`, `population.medians`.
- `tests/testthat.R`, `tests/testthat/test-gating-engine.R`.
- `scripts/engine_reference.py`, `scripts/engine_check_test.R`.

Edited
- `NAMESPACE`: exports anchor.population, calculate.gate, compute.threshold, population.masks, population.medians, population.table, recommend.algorithm, run.gating, threshold.kde.min, threshold.mixture, threshold.otsu, threshold.tail, train.boundaries; importFrom(AutoSpectral, match.fluorophores).

Conventions
- List keys use dots (gate.marker.x, gate.param$tail.fraction); string values are Honeychrome's ("kde_min", "x-y-|x+y+", "fixed"), so .agmodel I/O only renames keys.
- Missing thresholds are NA (Python None).
- Populations are stored as row indices of the data the gating ran on (`result$events`), not full-length masks.
- Gate parameter lists are read with `[[ ]]` (exact names): `$` partially matches, e.g. `threshold.factor` would return `threshold.factor.x`.

Differences from ag_core (deliberate)
- Mixture: own EM instead of scikit-learn; deterministic start, so last-digit differences in thresholds are expected.
- delta BIC subsample (> 5,000 events) uses R's RNG (seed 0, caller's stream restored), so values differ above 5,000 events; the recommendation rule is the same.
- No mixture-fit cache, no stop/pause (stop_check, keep_partial).
- Points exactly on a polygon edge may be classified differently from matplotlib.

Not in this phase
- Old engine (`do_gate_generic.R`, `calculate_gate.r`, `calculate_threshold*.R`, `apply_gate.R`, singlet wrappers) still present and used by define.gates / run.autogating; removed with the pipeline rewrite in Phase 3, together with the flowCore / flowStats / Biobase / modeest / fields / sp / mclust imports.

Check results (2026-11-11, Oliver's Mac; testthat suite passes)
- Algorithms vs snapshot ag_core: tail, kde_min, bimodal, otsu, kde_mode, extents <= 1.1e-16. Mixture (2 components) <= 2.7e-13 on two-population data; 0.029 (unimodal_tail) and 0.00084 (dbic_unimodal) on one-population data; 3 components 6e-6 to 2.1e-3. Debris cut follows the mixture (0.030 on unimodal_tail). delta BIC 0.063 on 3,000 events.
- Decision (Oliver): keep scikit-learn's stopping rule (tol 1e-3) in both; differences arise only when a mixture is fitted to one-population data.
- 10-gate hierarchy, trained on A (20,000 events) and applied to B: every gate runs in both; thresholds within 2e-5 (cd4cd8 mixture x), others <= 1e-6; counts within 2 events, except the convex hull ("all") at training, 5,029 (Python) vs 5,016 (R): the hull's own vertices are training events on the boundary, which matplotlib counts inside and the R even-odd test may not (1 event apart when applied to B). Sources, flags, stats.parent and anchor agree.
