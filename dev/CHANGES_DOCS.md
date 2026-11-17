# Clean-up and documentation

- Roxygen: `@importFrom` tags in `R/gate_figures.R` and `R/differential_figures.R` are one line each (roxygen2 rejects continuation lines). Roxygen run on a copy of the package (roxygen2 7.3.1, code sourced) gives no errors and the same NAMESPACE.
- App: every path field has a Browse button (shinyFiles: file chooser for the CSV, `.agmodel`, FlowJo and R gate files, folder chooser for the FCS, statistics and figure folders, save dialog for the model and the transformation CSV); roots are the home folder, the working folder and the mounted volumes. The text fields stay editable. The training-events button is now "Resample training events". `autogating.app()` checks for shinyFiles too.
- Parameters: removed settings nothing reads (`leukocyte.scatter.parameter`, `default.scatter.param(eter)`, `time.param`, `scatter.data.min/max`, `expr.data.min/max`, `ignore`, `parent.string`, `stats.separator`, `activation.marker.selected`, `activation.population.label`, `parent.regexp`, `parent.regexp.sub`). Kept: `singlet.scatter.parameter` (FlowJo import), `gate.number.width` (figure names). Fortessa file header corrected (file name and `cytometer = "fortessa"`).
- Non-ASCII characters in comments and documentation replaced (em dashes, arrows, a multiplication sign).
- DESCRIPTION: version 0.3.0; Suggests drop dplyr, parallel and AutoSpectralRcpp (none used), add shinyFiles.
- `NEWS.md` (new), `vignettes/autogating.Rmd` (new, ST1 example, code not run: the FCS files are not in the package; knits with rmarkdown 2.25).
- Removed `scripts/review_diagnostics_test.R` (the 2026-11-11 review diagnostics, written for the old functions; their results are in the review plan) and `scripts/__pycache__`.

