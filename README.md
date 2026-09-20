# cpopvar

`cpopvar` is a configuration-driven R package for chloroplast genome
population-variation analysis. It provides preprocessing, sample-frequency
filtering, per-kilobase normalisation, exploratory analysis (M01), comparative
analysis (M02), CDS hotspot analysis (M03), intergenic hotspot analysis (M04),
and supported plastome ideogram output (M05).

This package contains the programmatic backend only. The retired Shiny
workbench is not included.

## Development status

This directory is the unique live package source. The current private candidate
version is `0.0.0.9000`. It is undergoing package checks and API hardening and
is not a stable release. The first public beta is planned as `0.9.0`; the first
stable CRAN release is planned as `1.0.0`. Versioned artifacts and runtime
validation live in the sibling test/PM repository.

## Installation

After a source archive has passed the package gates, install it with:

```r
install.packages("cpopvar_VERSION.tar.gz", repos = NULL, type = "source")
library(cpopvar)
```

R 4.1 or later is required. Some plotting and annotation features have
additional system requirements.

## Synthetic example

The candidate includes a small synthetic bundle for installation and workflow
checks:

```r
example_dir <- system.file("extdata", "synthetic", package = "cpopvar")
config <- yaml::read_yaml(file.path(example_dir, "config.yml"))
config$input_files <- lapply(
  config$input_files,
  function(path) file.path(example_dir, path)
)

resolved_config <- tempfile(fileext = ".yml")
yaml::write_yaml(config, resolved_config)

result_dir <- cpopvar::run_analysis(
  config_file = resolved_config,
  output_dir = file.path(tempdir(), "cpopvar-example")
)
```

The example is intended for installation smoke tests. A successful check must
also confirm that enabled tasks completed without errors in the session log.

## Configuration and outputs

External column names are declared in the YAML configuration. The package maps
core inputs to stable internal names, while auxiliary grouping and ordering
files retain their configured external names.

`run_analysis()` creates a session directory under `output_dir` containing the
effective configuration, processed tables, module outputs, logs, and a run
summary. Exact files depend on the enabled tasks.

## Legacy statistical compatibility

The M03/M04 sharing, hypergeometric, driver, and two-speed paths are retained
from an earlier internal workflow so existing analyses remain runnable.
Compatibility tests cover invocation and output contracts only. Scientific
validity is not re-established by this release and will be reviewed and
replaced, where appropriate, in the next statistical-method iteration.

## Optional features

Raster conversion of PDF figures requires `magick`. Bioconductor-backed
annotation operations may require `GenomicRanges` and `rtracklayer`. Optional
features report missing dependencies instead of installing packages at runtime.

## Development checks

```bash
R CMD build .
R CMD check --as-cran cpopvar_VERSION.tar.gz
```

Package-check results must be reported from the corresponding archived logs.

## Citation

Use `citation("cpopvar")` after installation. Article metadata will be added
when a DOI is available.

## License

MIT © 2026 Liyun Nie. See `LICENSE` and `LICENSE.md`.
