# cpopvar 0.9.0

## Backend beta pre-release

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
- Treats DESCRIPTION `Version:` as the single maintained version source and
  checks NEWS, CITATION, NAMESPACE imports, and tarball naming against it.
- Adds the public configuration helpers `cpopvar_config()` and
  `validate_cpopvar_config()`. `run_analysis()` accepts a YAML path or an
  in-memory config and sends both through the same validator, parser, and
  dispatcher.
- Uses the Cairo PDF device when available so Windows R does not silently drop
  plot titles, axes, legends, annotations, or composite panels. If Cairo is
  advertised but fails to open (for example `failed to load cairo DLL` on some
  macOS CI hosts), falls back to the native PDF device instead of aborting.
- Builds grid composites without drawing to an ambient graphics device, then
  renders them through the same cross-platform PDF path.
- Checks `magick` and `pdftools` before requested raster conversion and reports
  a structured PDF-only fallback instead of silently producing zero PNG files.
- Records the actual PDF filenames in M01/M02 logs, summaries, and manifests.
- Adds bilingual English and Simplified Chinese user guides covering scientific
  questions, usage, pre-release installation, inputs, and a synthetic quick
  start.
- Documents clean local installation from a source archive with
  `remotes::install_local(..., dependencies = NA, upgrade = "never")` and a
  base-R fallback when hard dependencies are already installed.
- Opens the intended PDF device before UpSetR construction so M01, M03, and
  M04 plots write the named PDF once and do not create a default `Rplots.pdf`.

This backend beta is a GitHub Pre-release. It is not a CRAN release.
