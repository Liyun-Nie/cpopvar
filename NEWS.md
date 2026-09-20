# cpopvar 0.0.0.9000

## Development candidate

- Starts the independent public-package version lineage from the internal
  backend-only V4 baseline.
- Excludes the retired Shiny workbench and its runtime dependencies.
- Retains the M03/M04 legacy sharing, hypergeometric, driver, and two-speed
  paths for workflow compatibility; their statistical interpretation remains
  scheduled for a later method review.
- Adds corrected author, ORCID, license, repository, and citation metadata.
- Adds synthetic inputs, initial unit tests, a legacy contract smoke test, and
  a single cross-platform package-check workflow.
- Removes runtime package installation from backend helpers.
- Arranges M02 composite plots with patchwork first so ggplot2 4.x can keep
  ggpubr pairwise comparison panels instead of collapsing to a single plot.
- Draws M02 single-factor overall p-values with annotate() because ggplot2 4.x
  leaves `stat_compare_means()` overall tests blank.
- GitHub Actions checks use `--no-manual`; the PDF manual is validated on the
  WSL XeLaTeX host. Plot-layout PNG size gates accept smaller Windows files.

This four-component development version is not a stable release and must not be
tagged or submitted to CRAN.
