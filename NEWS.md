# AutoGating 0.3.0

A rewrite on the model of Honeychrome Automated Gating.

## New

* Gating models are read and written as Honeychrome `.agmodel` files
  (`read.agmodel()`, `write.agmodel()`); models move between R and
  Honeychrome.
* Gate engine as in Honeychrome: singlets band, cells rectangle, 1dsep and
  2dsep thresholds (populations are open-ended regions, quadrants can be
  joined), free polygons, replicate gates and replicated thresholds;
  threshold algorithms mixture (one-dimensional Gaussian mixture, no
  mclust), tail, kde_min, otsu; `recommend.algorithm()`.
* Training on an equal-sized subsample of each training sample, and
  per-sample application with fixed or recalculated boundaries and drift
  flags (`define.gates()`, `run.autogating()`, `run.gating()`), in parallel
  with `worker.n`.
* FlowJo workspace import (`read.flowjo.workspace()`,
  `import.flowjo.workspace()`) and data-driven conversion of drawn gates to
  calculated gates (`convert.imported.gates()`).
* Differential statistics on log2 odds within the parent and on marker
  medians, with limma (`differential.statistics()`), and their figures
  (`write.differential.figures()`).
* Gate figures from the gating result (`gate.figure()`,
  `write.gate.figures()`), with axes from the panel type and transform.
* The AutoGating app (`autogating.app()`): load an experiment and a model,
  change transforms, edit the gate hierarchy with live figures, save and run.
* `write.transformation.csv()`.
* Vignette on the ST1 example.

## Changed

* FCS files are read with `AutoSpectral::readFCS()`; flowCore, flowStats,
  flowWorkspace, Biobase, modeest, fields and sp are no longer used.
* All channels are on the Gating-ML 0 to 1 display scale (logicle via
  `AutoSpectral::biexp.transform()`); the time channel is taken from the
  panel type and kept in raw units.
* Earlier R gate lists are converted on load (`convert.gate.definitions()`).
* Population counts, medians and gate status are written as tidy CSV files;
  missing populations are `NA`, not 0.
* Parameters: `figure.dir`, `figure.training` (default `TRUE`) and
  `figure.samples` (default `FALSE`) replace `figure.dir.basename`; settings
  no longer read by the package (gate and plotting defaults, data ranges,
  sentinel strings) are removed.
* Random subsamples no longer change the session's random number stream.

## Fixed

* Activation medians reported the first sample's events for every sample.
* Threshold gates ended at the training range, losing events in brighter
  samples.
* Statistics parents `"*"` and `NULL`, three-quadrant joins, replicate
  gates whose origin has an ignored population, and activation markers on
  the parent's axis.
* FlowJo conversion: sides of positive and scatter gates, 2dsep gates, the
  separating axis, channel matching with matrix prefixes, nested activation
  gates; Boolean (NOT, AND, OR) gates are reported and left out instead of
  being given an algorithm that does not exist.
* Mixture, Otsu and density-minimum thresholds.
* Scatter axes on Aurora and ID7000, and axis labels on the Aurora.

## Removed

* `calculate_differential_abundance.R`, the old plotting functions, the old
  Shiny prototype and the optional AutoSpectralRcpp density estimate.
