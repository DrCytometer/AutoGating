# FlowJo import (2026-11-11)

Plan: project doc `claude/AUTOGATING_R_REVIEW_PLAN.md`, Part E Phase 4. Assumes the Phase 1-3 changes (`CHANGES_DATA_LAYER.md`, `CHANGES_GATE_ENGINE.md`, `CHANGES_PIPELINE.md`) landed.

Ported from snapshots of Honeychrome `ag_flowjo.py` (sha1 bd6bb1e4, unchanged since 2026-10-17) and `ag_review.py` (working tree 2026-11-11, sha1 31d09703) saved in `dev/reference/`.

Rename (Phase 3 follow-up): `upgrade.gate.definitions` -> `convert.gate.definitions` (roxygen registered it as an S3 method of utils::upgrade); file `R/convert_gate_definitions.R`.

New
- `R/flowjo_workspace.R`: `read.flowjo.workspace()` (groups and samples, labels, FlowJo transforms, matrix affixes, populations with rectangle/polygon/ellipse gates in raw units), `fj.default.source`, `fj.gated.parameters`, `fj.parameter.stem`, `fj.experiment` (panel side), `align.flowjo.parameters` (name -> fluorophore -> label -> closest name), `fj.plan.import`, `fj.plan.messages`.
- `R/flowjo_import.R`: FlowJo display scales (biex = AutoSpectral::biexp.transform with FlowJo's parameters, logicle/log via make.transform, arcsinh formula), polygon densification on FlowJo's scale, ellipses -> 10-vertex equal-area polygons, pieces and signs, sibling grouping, quadrant layout (80-cell grid), histogram rule, gate definitions with template / template.original / proposed.algorithm / proposal.rule / name.group / source, unique names, singlet gate moved to the top.
- `R/flowjo_review.R`: `convert.imported.gates()` (conversion evidence on the events: below-zero share, drawn-positive share, fit of tail / mixture / 3-component mixtures to the shapes as drawn, dropped non-separating axis, replicate on the anchor, debris side of cell gates; proposals applied; overlap check) and `import.flowjo.workspace()` (read, align with optional overrides, convert to gate definitions; messages for uncertain matches and dropped branches).
- `R/gate_model.R`: `join.quadrants`.
- `tests/testthat/test-flowjo.R`; `scripts/flowjo_reference.py`, `scripts/flowjo_check_test.R`.

Deleted: `R/wsp_to_autogating.R` (wsp.to.autogating, write.gate.definitions, get.cytometer.zero.params).

`NAMESPACE` regenerated from the roxygen tags (+ read.flowjo.workspace, import.flowjo.workspace, convert.imported.gates, importFrom AutoSpectral match.markers; - the three removed functions and the xml2 importFrom lines, xml2 is called with ::).

Deliberate differences from Honeychrome
- Fluorophore/marker recognition with AutoSpectral match.fluorophores / match.markers (same databases); closest-name fallback uses 1 - edit distance / longer length >= 0.8 instead of difflib.
- FlowJo hyperlog axes treated as linear (Honeychrome: FlowKit hyperlog).
- biexp.transform vs FlowKit WSPBiexTransform (FlowJo biex, width -250, pos 5.62, max 4194304): within 2 of 4096 display channels everywhere except far below zero, where FlowKit floors at 0 and biexp.transform continues.
- Review subsample of > 10,000 parent events drawn with R's RNG.
- No interactive split of 2dsep gates, no conversion-check flags (Builder UI features).

Cloud check (2026-11-11; scripts/flowjo_check_test.R vs flowjo_reference.py)
- Workspace_test.wsp (6 gates) and IMMPROVE template .wspt (1,638 gates, 1,678 shapes): alignment identical (10 and 90 parameters); gate names, types, channels, sep axis, parents, regions, proposed algorithms, proposal rules, name groups and numbers identical; shapes within 4.5e-11 (same vertex counts); no dropped populations in either.
- Conversion on synthetic events: Workspace_test 6/6 identical proposals (below-zero and positive shares within 0.006 from the subsample); IMMPROVE 1,637/1,638 identical; the exception, "BAFFR+ (Atypical)" (91 parent events), fits 3-component mixtures differently (R/Python EM on 91 events), so tail is chosen instead of mixture.
- FlowJo import 3.8 s, conversion 5.0 s for the IMMPROVE template; testthat suites pass.
