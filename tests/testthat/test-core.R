test_that("task parameters use explicit precedence", {
  config <- list(visualization = list(frequency_percentile = 0.70))
  direct <- list(frequency_percentile = 0.90)
  nested <- list(parameters = list(frequency_percentile = 0.80))

  expect_equal(get_task_parameter(direct, config, "frequency_percentile"), 0.90)
  expect_equal(get_task_parameter(nested, config, "frequency_percentile"), 0.80)
  expect_equal(get_task_parameter(list(), config, "frequency_percentile"), 0.70)
})

test_that("YAML validation accepts the documented required sections", {
  path <- tempfile(fileext = ".yml")
  yaml::write_yaml(
    list(
      session_info = list(session_id = "test"),
      input_files = list(
        main_data = "variants.csv",
        cds_lengths = "cds.csv",
        region_info = "regions.csv",
        genome_regions = "genome.csv"
      ),
      column_mappings = list(
        main_data = list(
          species = "species",
          variant_type = "var_type",
          position = "position",
          gene = "gene"
        ),
        cds_lengths = list(gene = "gene"),
        region_info = list(region_name = "region_name")
      ),
      analysis_tasks = list(
        preprocessing = list(enabled = TRUE)
      )
    ),
    path
  )

  expect_true(validate_yaml_config(path))
})

test_that("hotspot thresholds use positive frequencies and retain ties", {
  normalized <- data.frame(
    species = rep("Demo", 6),
    gene = letters[1:6],
    frequency_per_kb = c(0, 1, 2, 3, 4, 4),
    region_type = "CDS",
    var_type = "snp"
  )

  ranked <- identify_candidate_hotspots(
    normalized_data = normalized,
    config = list(),
    task_params = list(frequency_percentile = 0.75)
  )

  expect_equal(unique(ranked$threshold), 4)
  expect_equal(sum(ranked$is_potential_hotspot), 2)
  expect_false(ranked$is_potential_hotspot[ranked$frequency_per_kb == 0])
})

test_that("sample-frequency filtering keeps qualifying sites once", {
  variants <- data.frame(
    sample_id = c("s1", "s2", "s1"),
    species = "Demo",
    position = c(10, 10, 20),
    var_type = "snp",
    gene = c("matK", "matK", "rbcL"),
    region_type = "CDS",
    genome_region = "LSC"
  )

  filtered <- cpopvar:::filter_and_dedup_species(
    variants,
    list(thresholds = list(snp = 2L))
  )

  expect_equal(nrow(filtered), 1)
  expect_equal(filtered$position, 10)
})

test_that("run_analysis reports a missing configuration", {
  expect_error(
    run_analysis(file.path(tempdir(), "not-present.yml")),
    "Configuration file not found"
  )
})

test_that("session paths expose the reporting stage contract", {
  output_dir <- tempfile("cpopvar-path-contract-")
  paths <- get_session_paths(
    session_id = "reporting_contract",
    output_dir = output_dir
  )

  expect_identical(
    paths$P01_preprocessed_out,
    file.path(paths$processed_data, "P01_preprocessed")
  )
  expect_identical(
    paths$P02_filtered_out,
    file.path(paths$processed_data, "P02_filtered")
  )
  expect_identical(
    paths$P03_normalized_out,
    file.path(paths$processed_data, "P03_normalized")
  )
  expect_true(all(dir.exists(unlist(paths))))

  unlink(output_dir, recursive = TRUE)
})

test_that("reporting receives prior task manifests", {
  manifest <- list(module = "M03", entries = list())
  manifests <- cpopvar:::collect_task_manifests(
    list(
      task_results = list(
        M03_hotspot = list(
          success = TRUE,
          result_data = list(manifest = manifest)
        )
      )
    )
  )

  expect_named(manifests, "M03_hotspot")
  expect_identical(manifests$M03_hotspot, manifest)
})

test_that("session statistics use current processing-stage paths", {
  output_dir <- tempfile("cpopvar-reporting-stats-")
  withr::local_options(cpopvar.output_dir = output_dir)

  stats <- cpopvar:::collect_session_statistics(
    session_id = "reporting_stats",
    config_data = list(analysis_tasks = list()),
    execution_context = NULL
  )

  expect_false(isTRUE(stats$error_occurred))
  expect_match(stats$data_quality_summary, "Data not found")

  unlink(output_dir, recursive = TRUE)
})
