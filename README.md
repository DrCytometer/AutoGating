# AutoGating

Automated, reproducible gating of flow cytometry data in R, compatible with the
Automated Gating plugins of Honeychrome.

A gate hierarchy (singlets, cells, thresholds on one or two channels, polygons,
gates copied from other gates) is trained on pooled events from a set of training
samples and applied to every sample. Thresholds are calculated from the data: a
Gaussian mixture crossover, a symmetric tail, a density valley or Otsu's method,
or copied from another gate or a reference population. Models are read and
written as Honeychrome `.agmodel` files, so a model built in R opens in
Honeychrome and the other way round. FlowJo workspaces can be imported and their
drawn gates converted to calculated ones.

## Installation

```r
install.packages("pak")
pak::pak("DrCytometer/AutoGating")

# for the app
pak::pak(c("shiny", "bslib", "DT", "shinyFiles"))
```

## Pipeline

The FCS files are compensated or unmixed data: unmix spectral data first
(AutoSpectral, Honeychrome or the instrument software). Compensation of
conventional data can be done here with `compensation.filename`.

```r
library(AutoGating)

ap <- get.autogating.param(
    cytometer           = "aurora",            # or "symphony", "fortessa", "id7000"
    data.dir            = "fcs",
    panel.filename      = "panel.csv",         # dye, antigen, type
    experiment.filename = "experiment.csv",    # filename, sample, gate.calculation
    model.filename      = "model.agmodel",     # or population.gates.filename for R gate lists
    statistics.dir      = "statistics",
    figure.dir          = "figures"
)

fcs.setup <- setup.autogating(ap)           # read, compensate, transform
model     <- define.gates(fcs.setup, ap)    # train gates without a stored boundary
write.agmodel(model, "model.agmodel", fcs.setup$fcs.panel)
gating    <- run.autogating(fcs.setup, model, ap)   # counts, medians, gate status

stats <- differential.statistics(gating, "samples.csv", group.column = "group",
                                 output.dir = "differential")
```

`create.panel.file()` and `create.experiment.file()` draft the panel and
experiment files from the FCS files; check the draft panel and remove the
channels that are not gated (for example an event ID or autofluorescence
channel). Figures of every gate on the training events are written by
`define.gates()` (`figure.training = FALSE` skips them); `figure.samples = TRUE`
draws every gate on every sample as well.

## From a FlowJo workspace

```r
imported  <- import.flowjo.workspace("template.wspt", fcs.setup, ap)
imported$mapping                                    # FlowJo parameter -> panel channel
converted <- convert.imported.gates(imported$gate.definitions,
                                    training.set(fcs.setup, ap)$data,   # the events define.gates trains on
                                    transforms = fcs.setup$ap$transform.params,
                                    panel      = fcs.setup$fcs.panel)
converted$proposals                                 # algorithm chosen for each gate, with its agreement
model <- define.gates(fcs.setup, ap, model = new.agmodel(converted$gate.definitions))
```

Drawn gates are first applied as drawn; `convert.imported.gates()` replaces each
with the calculated gate (type and threshold algorithm) that best reproduces it
on the data. Gates with a low `agreement` are worth checking in the figures or
the app.

## App

```r
result <- autogating.app()                          # load an experiment in the app
result <- autogating.app(fcs.setup, ap, model)      # or edit one read in R
model  <- result$model
```

- **Data**: load an experiment and a model (empty, `.agmodel`, FlowJo workspace,
  or the earlier R gate lists); every path has a Browse button.
- **Transforms**: change the logicle, log or linear transform of each channel,
  previewed on any population; save the transforms as a CSV.
- **Gates**: the gate hierarchy with each population's share of its parent, the
  figure of the selected gate and an editor for its parent, channels,
  populations, algorithm and settings. Click the figure to set a manual
  threshold. Every change is recalculated at once and can be undone.
- **Save and run**: save the model, write the training figures, run the model on
  every sample, and return to R.

## Funding

This work has received funding from the KU Leuven C1 program, the European Union’s
Horizon 2020 research and innovation programme under grant agreement No 874707
(EXIMIOUS), Wellcome Investigator Award, 222442/A/21/Z and UKRI Proactive
Vaccinology Award, MR/Y004450/1 (IMMPROVE).
