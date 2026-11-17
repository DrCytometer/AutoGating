# Data layer (2026-11-11)

Plan: project doc `claude/AUTOGATING_R_REVIEW_PLAN.md`, Part E Phase 1.

New
- `R/transform.R`: `make.transform()` (Gating-ML flin / flog written out; flogicle is `AutoSpectral::biexp.transform(channelRange = 1)`, A = 0 only), `flowjo.biex.to.logicle()` (T = max.range, M = pos, A = neg, W = log10(-width)/2).
- `scripts/transform_reference.py` (FlowKit reference values), `scripts/data_layer_check_test.R`.
- `.Rbuildignore` (scripts/, dev/, .Rproj).

Rewritten
- `R/analyze_aux.r` `read.transformation()`: params per antigen for every panel channel (logicle fluorescence, linear scatter, identity time); CSV in Gating-ML or legacy FlowJo columns; returns list(forward, inverse, params). New trailing argument `fcs.panel.time.idx`.
- `R/setup_autogating.R`: `AutoSpectral::readFCS(columns = panel dyes)`, compensation X %*% solve(S), transforms per channel, one matrix without row names + `fcs.event.sample` factor. Removed from the return value: `fcs.event`, `fcs.event.regexp`. Added: `fcs.panel.time.idx`, `ap$transform.params`. No `$PnR` truncation; time not rescaled.

Edited
- `R/get_autogating_param*.r`: `transformation.defaults` (FlowJo biex) replaced by `transform.defaults` (Gating-ML, fluorescence + scatter).
- `R/create_panel_file.R`: keywords via `readFCS` instead of `flowCore::read.FCS`.
- `R/read_compensation_csv.R`: documents the spillover orientation.
- `DESCRIPTION`, `NAMESPACE`: + AutoSpectral (Remotes DrCytometer/AutoSpectral); - flowWorkspace; flowCore imports reduced to flowFrame, polygonGate (removed with the Phase 2 algorithms).

Behaviour changes
- Display scale 0-1 for all transformed channels (was 0-256 FlowJo channels for fluorescence, raw for scatter).
- Values above `$PnR` no longer truncated (as Honeychrome).
- Time channel left in raw units (was divided by 100).

Not yet updated (later phases): `scripts/analyze.r` (legacy copy of the pipeline, calls the old `read.transformation`), `inst/shiny/app.R`, roxygen references to CLAUDE.md in define_gates.R / create_experiment_file.R.

Check results (data_layer_check_test.R, 2026-11-11)
- Transforms vs FlowKit: forward <= 2.2e-16, round trip <= 2e-12 (relative), all nine parameter sets.
- FlowJo biex settings vs flowjo_biexp/256: max 0.033 (width -100), 0.018 (neg 1, width -20), 0.0065 (width -1000). FlowJo biex is not the Gating-ML logicle; the legacy CSV is mapped to the closest logicle (documented). FlowJo workspace import (later phase) needs FlowJo's biex itself to redraw FlowJo shapes, as Honeychrome does with FlowKit WSPBiexTransform.
- 3 Aurora spleen files (1,325,736 events, 47 channels): identical to flowCore read + compensate + same transforms (max difference 0); 9.0 s, 475 MB. Old $PnR truncation changed 3,672 to 4,308 values per file (of ~17-28 million).

Revision (2026-11-11, Oliver): logicle uses AutoSpectral::biexp.transform instead of a second implementation; flowWorkspace comparison removed from the check script (flowjo_biexp interpolates a spline over FlowJo's exported table, not a reference). The hand-written logicle (Halley solve, A != 0 supported) matched FlowKit to 2e-16 before it was removed; biexp.transform needs re-checking against the same FlowKit values.

Re-check with biexp.transform (2026-11-11): max |forward difference| vs FlowKit 4.5e-11 (W 1, M 4.42), 6.0e-12 (W 1.5, M 5.62), <= 1.5e-15 for the others; round trip <= 1.9e-12 relative; A = 1 refused by biexp.transform as expected. Data layer unchanged (max difference 0), 7.5 s. Phase 1 complete.
