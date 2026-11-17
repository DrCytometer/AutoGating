# Differential statistics (2026-11-11)

Plan: project doc `claude/AUTOGATING_R_REVIEW_PLAN.md`, Part E Phase 5. Assumes Phases 1-4 landed (`CHANGES_*.md`).

Reference: snapshots of Honeychrome `ag_results.py` (sha1 78126188) and `controller_components/differential_stats.py`, `moderated_stats.py`, `anova_stats.py` (HEAD, unmodified) in `dev/reference/` (package layout so `from honeychrome.controller_components import` works).

Rd fix (Phase 4 follow-up): the internal `%|NA|%` operator removed (its Rd name was invalid); `man/grapes-or-NA-or-grapes.Rd` deleted.

New
- `R/differential_statistics.R`: `differential.statistics(gating, sample.table, ...)` and internals `share.prior`, `population.log2.odds`, `population.frequency.tables`, `ancestor.channels`, `marker.median.table`, `resolve.pairing`, `design.matrix`, `moderated.fit`, `moderated.test`, `stagewise.adjust`, `apply.significance`. limma lmFit / eBayes / duplicateCorrelation; TREAT p-values; BH per comparison and global; population-first stage-wise FDR for marker medians. Defaults as Honeychrome: min 20 parent events, 3 measured and 3 non-zero samples, 3 samples per group, thresholds 0.5 (log2 odds) and 0.03 (medians), FDR 0.05 global, reference mode with the first group as baseline.
- `tests/testthat/test-statistics.R`; `scripts/stats_reference.py`, `scripts/stats_check_test.R`.

Deleted: `R/calculate_differential_abundance.R` (limma on raw fractions, wide tables, heatmap/volcano). DESCRIPTION: - ggrepel, pheatmap. NAMESPACE regenerated (+ differential.statistics, - calculate.differential.abundance; limma via ::).

Fixes A2 #19 (missing/empty populations are NA and tested on the samples that have them, not 0) and #20 (log2 odds with the share prior; comparisons named "other vs baseline", any group names).

Not ported (v1): the counts test (moderated vs anchor, QL GLM), ANOVA variance tables, capped weights; samples with a missing covariate value are left out instead of fitting each comparison on its own samples. Figures of the results: see CHANGES_PLOTTING.md.

Cloud check (2026-11-11; limma 3.58.1; stats_check_test.R vs run_statistics on 15 synthetic samples, 3 groups, 11 populations, 6 channels)
- reference, pairwise, fixed pairing, numeric covariate: frequencies and marker medians identical to 1e-14 (logFC, SE, P, adj, adj.global, stage-wise); significance identical.
- random individual effect (nested mice): logFC within 7e-5, P within 7e-5, adj within 1.4e-4 (limma duplicateCorrelation vs Honeychrome's REML port); significance identical.
