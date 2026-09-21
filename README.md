# cpopvar

[简体中文](README.zh-CN.md)

`cpopvar` is an R package for analysing population variation in chloroplast
genomes (plastomes). After the input tables are prepared, it can preprocess
variants, filter them by sample-occurrence frequency, normalise counts per
kilobase, and then analyse variant distributions, group comparisons, CDS and
intergenic (IGS) hotspots, and plastome ideograms.

Analyses are run from R functions or an equivalent YAML configuration. There
is no graphical interface.

cpopvar is currently a pre-release and is not available from CRAN.

## Biological questions

Shared preprocessing and normalisation come first; later tasks are enabled by
question. You do not need to turn every module on in a single run. Results
from one run can narrow the inputs, groups, or candidate loci of the next.
Enabling any of M01–M05 also runs P01–P03.

- **Shared processing.** P01 annotates variants and compiles region
  information. P02 filters variants by a sample-occurrence frequency
  threshold. P03 normalises counts per kilobase under unified region types,
  variant types, and lengths.
- **Where does variation occur, and is it uneven?** Enable
  `M01_distribution`.
- **Which biological groups differ?** Enable `M02_comparative`.
- **Which CDS genes or intergenic (IGS) intervals are high-variation
  candidates?** Enable `M03_hotspot` and/or `M04_igs_hotspot`.
- **Where do those candidates fall on the plastome map?** Enable
  `M05_hotspot_ideogram` in the same run as `M03_hotspot` and
  `M04_igs_hotspot`.

## Usage

1. Prepare the CSV/TSV input tables (see Input data below).
2. Choose the M01–M05 tasks that match the questions you want to ask.
3. Build a configuration in memory with `cpopvar_config()`, or write an
   equivalent YAML file.
4. Check the configuration with `validate_cpopvar_config()` and, if it
   passes, call `run_analysis()`.
5. Inspect the session directory returned by the run. Adjust the inputs or
   tasks and rerun if needed.

## Pre-release installation

cpopvar is currently a pre-release and is not available from CRAN. Install it
from a local source archive. In the commands below, replace `VERSION` with
the version string in the archive file name.

**Recommended for a clean R library.** Install `remotes`, then install the
archive so that hard dependencies declared in DESCRIPTION (`Depends`,
`Imports`, and `LinkingTo`) are resolved from your configured repositories:

```r
install.packages("remotes")
remotes::install_local("cpopvar_VERSION.tar.gz", dependencies = NA, upgrade = "never")
library(cpopvar)
```

`dependencies = NA` installs those hard dependencies. It does not install
packages listed under `Suggests`. `upgrade = "never"` leaves already-installed
packages in place.

**Base-R fallback.** Use this path only when the required dependencies are
already installed:

```r
install.packages("cpopvar_VERSION.tar.gz", repos = NULL, type = "source")
library(cpopvar)
```

`repos = NULL` does not resolve or install missing dependencies from CRAN.
On a clean machine this command can fail if those packages are not already
present.

## Environment and dependencies

- R 4.1.0 or later is required (`Depends: R (>= 4.1.0)`).
- Core R packages are listed under `Imports` in DESCRIPTION. The recommended
  `remotes::install_local(..., dependencies = NA)` call installs or resolves
  those hard dependencies. cpopvar itself does not install packages at
  runtime.
- Snippy is an external upstream tool for preparing `main_data`. cpopvar
  neither installs nor invokes it.
- On Windows, Rtools is needed only if R must compile a dependency from
  source. It is not required for every installation.
- Packages in `Suggests` are feature-specific and are not installed by
  `dependencies = NA`. Install a missing optional package only for the
  feature you want to use.

## Quick start

The package includes a small synthetic example that walks through one
analysis. The code below reads the bundled files, builds a configuration, and
runs CDS hotspot analysis (`M03_hotspot`) after the configuration validates:

```r
example_dir <- system.file("extdata", "synthetic", package = "cpopvar")
bundled <- yaml::read_yaml(file.path(example_dir, "config.yml"))
bundled$input_files <- lapply(
  bundled$input_files,
  function(path) file.path(example_dir, path)
)

cfg <- cpopvar::cpopvar_config(
  main_data = bundled$input_files$main_data,
  annotations = bundled$input_files$annotations,
  genome_regions = bundled$input_files$genome_regions,
  group_info = bundled$input_files$group_info,
  species_order = bundled$input_files$species_order,
  gene_function_mapping = bundled$input_files$gene_function_mapping,
  tasks = "M03_hotspot",
  task_options = list(M03_hotspot = bundled$analysis_tasks$M03_hotspot$parameters),
  visualization_settings = bundled$visualization_settings,
  output_settings = bundled$output_settings
)
stopifnot(cpopvar::validate_cpopvar_config(cfg)$valid)

result_dir <- cpopvar::run_analysis(
  cfg,
  output_dir = file.path(tempdir(), "cpopvar-example")
)
```

`result_dir` is the session directory written for this run, under
`output_dir` at `app_data/sessions/<session_id>/`. Open that directory to
inspect input copies, processed tables, figures, and logs.

You can also pass a YAML configuration path to `run_analysis()` as `config`
or `config_file`. Relative paths in that YAML must be converted to absolute
paths, or the call must be made from the directory that contains the YAML
file.

## Input data

cpopvar currently starts from prepared CSV/TSV tables. It does not read
assembled FASTA or GenBank files directly.

External column names for the core tables can be mapped to internal names in
YAML `column_mappings` (or the corresponding arguments of
`cpopvar_config()`). Auxiliary grouping and ordering files keep the column
names given in the configuration; later tasks use those columns as they are.

Bundled examples live under
[`inst/extdata/synthetic/`](inst/extdata/synthetic/config.yml).

### Tables to prepare

- **`main_data` (required).** One row per sample-level variant. A single
  sample usually contributes many rows. Recommended column names:

  `sample_id`, `species`, `position`, `var_type`, `gene`, `region_type`,
  `ref_allele`, `alt_allele`

  Other headers can be mapped in the configuration. Example:
  [`inst/extdata/synthetic/main_data.csv`](inst/extdata/synthetic/main_data.csv).

- **`annotations` (needed for the full P01–P03 workflow).** The current
  input is CSV/TSV, not a `.gb` file. Expected fields are `Species`, `Type`,
  `Gene`, `Minimum`, `Maximum`, `Length`, and `Intervals`. If a Geneious
  export uses `Sequence Name` / `# Intervals`, rename those headers first.
  Example:
  [`inst/extdata/synthetic/annotations.csv`](inst/extdata/synthetic/annotations.csv).

- **`genome_regions` (needed for the full P01–P03 workflow).**
  One row per species, with start and end coordinates for the LSC, inverted
  repeat (IR), and SSC:

  `species`, `lsc_start`, `lsc_end`, `ir_start`, `ir_end`, `ssc_start`,
  `ssc_end`

  Include `total_length` when the genome length is available. Genomes that
  lack inverted repeats can be marked as `IR_lacking_genome` in
  `special_handling`. Example:
  [`inst/extdata/synthetic/genome_regions.csv`](inst/extdata/synthetic/genome_regions.csv).

- **Optional auxiliary tables.** Match the bundled examples rather than
  assuming a universal schema:
  - [`group_info.csv`](inst/extdata/synthetic/group_info.csv) — species-level
    grouping labels used by `M02_comparative`.
  - [`species_order.csv`](inst/extdata/synthetic/species_order.csv) — display
    order and optional plot labels.
  - [`gene_function_map.csv`](inst/extdata/synthetic/gene_function_map.csv) —
    gene-to-function lookup used by hotspot enrichment views.

YAML example:
[`inst/extdata/synthetic/config.yml`](inst/extdata/synthetic/config.yml).

### Preparing `main_data`

These steps are done outside cpopvar. The package does not read FASTA or
GenBank files and does not call Snippy.

For each species, prepare assembled plastome FASTA files for all samples and
one detailed annotated GenBank reference (usually from a selected sample).
Remove one inverted-repeat (IR) copy from both the assemblies and that
reference so the same variant is not counted twice and coordinates stay
aligned.

Then run Snippy for each assembly against that species reference, for
example:

```bash
snippy --outdir SAMPLE --ctgs SAMPLE.single_IR.fasta --ref SPECIES_REFERENCE.gb --cpus 1
```

Prepend species and sample identifiers to each sample’s `snps.tab`,
concatenate tables within species, merge the species tables, and map columns
to the `main_data` schema above. Concatenated Snippy intermediates are not
themselves a valid input; supply the mapped table.

## Outputs

`run_analysis()` writes results to a session directory:

`<output_dir>/app_data/sessions/<session_id>/`

Exact files depend on the tasks enabled for that run. A typical layout is
shown below.

- `raw/` — copies of the configured input tables (`annotations` under
  `raw/annotations/`).
- `results/processed_data/P01_preprocessed/` — annotated variant tables
  (whole-genome and single-IR-copy tables) and `region_info_complete.csv`.
- `results/processed_data/P02_filtered/` — `variants_filtered.csv` (variants
  after the sample-occurrence frequency filter).
- `results/processed_data/P03_normalized/` — `normalized_frequencies.csv` and
  `normalized_frequencies_summary.csv`.
- `results/plots/` — tables and figures from enabled M01–M05 tasks. These may
  include M01 distribution summaries and plots, M02 comparative summaries and
  figures, M03 CDS hotspot candidates, M04 IGS hotspot candidates, and the
  M05 ideogram plus marker/heatmap tables.
- `results/run_summary.rds` — summary of executed tasks. When reporting is
  enabled, `results/Final_Analysis_Report.md` is also written. The
  configuration actually used for the run is saved as
  `results/config_used.yml`.
- `logs/` — session logs.

## Citation

Use `citation("cpopvar")` after installation. Article metadata will be added
when a DOI is available.

## License

MIT © 2026 Liyun Nie. See `LICENSE` and `LICENSE.md`.
