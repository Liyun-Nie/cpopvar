synthetic_config_dir <- function() {
  system.file("extdata", "synthetic", package = "cpopvar")
}

load_synthetic_yaml_config <- function(absolute_paths = TRUE) {
  example_dir <- synthetic_config_dir()
  config <- yaml::read_yaml(file.path(example_dir, "config.yml"))
  if (isTRUE(absolute_paths)) {
    config$input_files <- lapply(config$input_files, function(path) {
      file.path(example_dir, path)
    })
  }
  config
}

synthetic_memory_config <- function(yaml_cfg) {
  enabled <- names(yaml_cfg$analysis_tasks)[vapply(
    yaml_cfg$analysis_tasks,
    function(task) isTRUE(task$enabled %||% TRUE),
    logical(1)
  )]
  task_options <- lapply(yaml_cfg$analysis_tasks, function(task) {
    task$parameters %||% list()
  })
  cpopvar_config(
    main_data = yaml_cfg$input_files$main_data,
    annotations = yaml_cfg$input_files$annotations,
    genome_regions = yaml_cfg$input_files$genome_regions,
    group_info = yaml_cfg$input_files$group_info,
    species_order = yaml_cfg$input_files$species_order,
    gene_function_mapping = yaml_cfg$input_files$gene_function_mapping,
    tasks = enabled,
    task_options = task_options,
    column_mappings = yaml_cfg$column_mappings,
    session_id = yaml_cfg$session_info$session_id,
    visualization_settings = yaml_cfg$visualization_settings,
    output_settings = yaml_cfg$output_settings
  )
}

test_that("cpopvar_config requires main_data and rejects unknown arguments", {
  expect_error(cpopvar_config(), "main_data is required")
  expect_error(
    cpopvar_config(main_data = "variants.csv", frequency_percentile = 0.8),
    "Unknown arguments"
  )
  expect_error(
    cpopvar_config(main_data = "variants.csv", tasks = "not_a_real_task"),
    "Unknown task"
  )
})

test_that("cpopvar_config builds YAML-isomorphic defaults and enables core tasks", {
  cfg <- cpopvar_config(
    main_data = "variants.csv",
    annotations = "annotations.csv",
    genome_regions = "genome_regions.csv",
    tasks = c("M01_distribution", "M03_hotspot")
  )

  expect_s3_class(cfg, "cpopvar_config")
  expect_named(
    cfg,
    c("session_info", "input_files", "column_mappings", "analysis_tasks"),
    ignore.order = TRUE
  )
  expect_identical(cfg$input_files$main_data, "variants.csv")
  expect_true(cfg$analysis_tasks$preprocessing$enabled)
  expect_true(cfg$analysis_tasks$filtering$enabled)
  expect_true(cfg$analysis_tasks$normalization$enabled)
  expect_true(cfg$analysis_tasks$M01_distribution$enabled)
  expect_true(cfg$analysis_tasks$M03_hotspot$enabled)
  expect_false(cfg$analysis_tasks$M02_comparative$enabled)
})

test_that("auxiliary files cannot have column mappings", {
  cfg <- cpopvar_config(main_data = "variants.csv", tasks = "preprocessing")
  raw <- unclass(cfg)
  raw$column_mappings$group_info <- list(species = "species")
  result <- validate_cpopvar_config(raw)
  expect_false(result$valid)
  expect_true(any(vapply(result$errors, function(err) err$code, character(1)) == "invalid_mapping"))
})

test_that("validate_cpopvar_config returns structured errors without a YAML file", {
  result <- validate_cpopvar_config(list(foo = 1))
  expect_false(result$valid)
  expect_null(result$config)
  expect_true(length(result$errors) >= 1)
  expect_true(all(c("code", "path", "message") %in% names(result$errors[[1]])))
  expect_error(
    validate_cpopvar_config(list(foo = 1), error = TRUE),
    "Configuration validation failed"
  )
})

test_that("YAML path and memory config share the same validator and execution plan", {
  yaml_cfg <- load_synthetic_yaml_config(absolute_paths = TRUE)
  mem_cfg <- synthetic_memory_config(yaml_cfg)

  yaml_path <- tempfile(fileext = ".yml")
  yaml::write_yaml(yaml_cfg, yaml_path)

  yaml_check <- validate_cpopvar_config(yaml_path)
  mem_check <- validate_cpopvar_config(mem_cfg)
  expect_true(yaml_check$valid)
  expect_true(mem_check$valid)

  parsed_yaml <- parse_task_configuration(yaml_cfg, validate_files = FALSE)
  parsed_mem <- parse_task_configuration(unclass(mem_cfg), validate_files = FALSE)

  expect_setequal(
    cpopvar:::collect_enabled_task_names(yaml_cfg$analysis_tasks),
    cpopvar:::collect_enabled_task_names(mem_cfg$analysis_tasks)
  )
  expect_setequal(
    cpopvar:::execution_plan_task_names(parsed_yaml$execution_plan),
    cpopvar:::execution_plan_task_names(parsed_mem$execution_plan)
  )

  enabled <- collect_enabled_task_names(yaml_cfg$analysis_tasks)
  for (task_name in enabled) {
    expect_equal(
      parsed_yaml$tasks[[task_name]]$parameters,
      parsed_mem$tasks[[task_name]]$parameters
    )
  }
})

test_that("run_analysis accepts config_file and in-memory config equivalently", {
  yaml_cfg <- load_synthetic_yaml_config(absolute_paths = TRUE)
  mem_cfg <- synthetic_memory_config(yaml_cfg)
  yaml_path <- tempfile(fileext = ".yml")
  yaml::write_yaml(yaml_cfg, yaml_path)

  yaml_out <- tempfile("cpopvar-yaml-")
  mem_out <- tempfile("cpopvar-mem-")
  dir.create(yaml_out)
  dir.create(mem_out)

  yaml_session <- run_analysis(config_file = yaml_path, output_dir = yaml_out)
  mem_session <- run_analysis(mem_cfg, output_dir = mem_out)

  yaml_used <- yaml::read_yaml(file.path(yaml_session, "results", "config_used.yml"))
  mem_used <- yaml::read_yaml(file.path(mem_session, "results", "config_used.yml"))
  expect_setequal(
    cpopvar:::collect_enabled_task_names(yaml_used$analysis_tasks),
    cpopvar:::collect_enabled_task_names(mem_used$analysis_tasks)
  )

  yaml_norm <- file.path(yaml_session, "results", "processed_data", "P03_normalized")
  mem_norm <- file.path(mem_session, "results", "processed_data", "P03_normalized")
  expect_true(dir.exists(yaml_norm))
  expect_true(dir.exists(mem_norm))
  expect_setequal(list.files(yaml_norm), list.files(mem_norm))

  unlink(c(yaml_out, mem_out), recursive = TRUE)
})
