# Shiny app rebuild

The old app (`inst/shiny/app.R`, flowCore, hand-sourced helper files, R gate lists) is replaced by an app on the new engine, opened with `autogating.app()`.

## How it works, in lay terms

The app keeps one experiment, one gating model and the pooled training events (the same events `define.gates()` trains on: `ap$events.per.sample` from each sample marked `gate.calculation`). Every change you make (a gate edited, added, deleted or recalculated, a transform changed) is applied at once: gates that already have a boundary keep it, gates without one are calculated on the training events, and every population is gated again from the top, so child gates always see their parent's current events. The previous 20 versions of the model are kept for Undo.

- **Data.** Load an experiment from its files (the same inputs as `get.autogating.param()` + `setup.autogating()`) or start the app with one already read in R. Load a model: empty, an `.agmodel` file, a FlowJo workspace (`import.flowjo.workspace`, then optionally `convert.imported.gates` to turn drawn gates into calculated ones), or the earlier R gate lists (`convert.gate.definitions`).
- **Transforms.** Pick a channel and change its logicle (T, W, M), log (T, M) or linear (T, A) transform; the histogram shows the training events or any population with the new transform before it is applied. Applying puts the channel's events back in raw units and transforms them again, and stored boundaries on that channel are converted (`remap.model`), so a threshold stays at the same raw value. Save the transforms as a CSV that `transformation.filename` reads.
- **Gates.** The hierarchy (indented by depth) with each population's % of parent, the algorithm ("manual" for a manual threshold) and the status (trained, failed with the reason). The selected gate's figure is the Phase 6 gate figure on the training events. The editor shows: name, type, parent population, channels, threshold axis (1dsep shown against a second channel), algorithm per axis, population names, threshold settings (manual threshold, search range, tail fraction, threshold factor, replicate source), singlet band width, cells debris algorithm and edge density, and whether the gate is fixed or recalculated on each sample. On a 2dsep gate, quadrants given the same name form one population (an L, or two opposite quadrants); an empty name leaves a quadrant out. Clicking the figure fills the manual threshold(s); Apply uses them; clearing the field returns to the algorithm. "Recommend algorithm" runs `recommend.algorithm` on each threshold channel of the parent population.
- **Renaming.** Renaming a gate or a population (matched to the old one by region) updates every reference: children's parents, replicate origins and sources, statistics parents, application modes and the stored boundary. A gate whose calculation settings changed loses its stored boundary and is calculated again; a gate only renamed keeps it.
- **Save and run.** Save the `.agmodel` (`write.agmodel`), write the training figures (`write.gate.figures`), run the model on every sample (`run.autogating`, with optional per-sample figures; flagged gates are listed), and Finish, which closes the app and returns `list(model, fcs.setup)` to R.

## Code changes

- `R/autogating_app.R`: `autogating.app(fcs.setup, ap, model, launch.browser)` (exported). Checks for shiny, bslib and DT, evaluates `inst/shiny/app.R` in an environment whose parent is the package namespace, and runs it.
- `inst/shiny/app.R`: rewritten (UI and server only; all gating, editing and file work is done by package functions).
- `R/gate_editing.R` (new, internal unless stated): `population.keys`, `form.regions`, `gate.form.values`, `gate.from.form`, `put.gate`, `drop.gate`, `gate.table`, `retransform.channel`, `write.transformation.csv` (exported).
- `R/define_gates.R`: the body of `define.gates` is split into `model.from.files`, `prepare.model`, `training.set` and `record.training`, which the app calls too, so the app and the pipeline train identically. Behaviour of `define.gates` unchanged.
- Root gates built in the editor have parent population "root", as Honeychrome and the converters write them.
- `README.md` rewritten: it described the old app; it now covers the package, the pipeline and the app.
- DESCRIPTION Suggests: shiny, bslib, DT, htmltools.

## Checks (cloud, shiny 1.8.0, bslib 0.6.1, DT 0.31)

- `tests/testthat/test-gate-editing.R`: form round trip (incl. unions and replicate sources), gate and population renames reaching children, replicate sources and stored boundaries, boundary dropped when settings change, new gate numbering and parent checks, dropping descendants, retransform keeps raw values, transformation CSV read back by `read.transformation`. All test files pass.
- `shiny::testServer` on the simulated Aurora experiment: start with an experiment and model; edit a 2dsep gate (rename a population, manual X threshold, Y algorithm tail, recalculate mode); add a child gate; rename its parent population (child follows); Undo; change the CD8 transform (threshold keeps its raw value, 859); delete a gate with descendants; save and re-read the `.agmodel`; start empty, build a singlets gate, run on every sample; missing files reported by name.
- Headless Chromium on the running app: all four tabs render; selecting a gate draws its figure; clicking the figure fills the manual thresholds and Apply moves the lines.

Not checked: loading an experiment from FCS files inside the app (the cloud has no AutoSpectral FCS reader; the app calls `setup.autogating` as the pipeline does) and FlowJo import inside the app (the functions themselves were checked in Phase 4).

Speed: every change re-applies the whole hierarchy to the training events. With many gates or many training samples, lower "Training events per sample" on the Data tab ("Resample training events") while building the model, then retrain at full size with `define.gates(fcs.setup, ap, model, retrain = TRUE)` (manual thresholds are part of the gate definitions and are kept).
