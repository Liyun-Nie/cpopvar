test_that("task parameters honor the documented precedence", {
  config <- list(visualization = list(cutoff = 0.5))

  expect_equal(
    get_task_parameter(list(cutoff = 0.8), config, "cutoff", 0.1),
    0.8
  )
  expect_equal(
    get_task_parameter(list(parameters = list(cutoff = 0.7)), config, "cutoff", 0.1),
    0.7
  )
  expect_equal(get_task_parameter(NULL, config, "cutoff", 0.1), 0.5)
  expect_equal(get_task_parameter(NULL, list(), "cutoff", 0.1), 0.1)
})

test_that("missing YAML files fail validation", {
  path <- tempfile(fileext = ".yml")
  expect_error(validate_yaml_config(path), "Configuration file not found")
  result <- validate_yaml_config(path, return_errors = TRUE)
  expect_false(result$valid)
})

test_that("synthetic M03 config provides every required color contract", {
  config_path <- system.file(
    "extdata",
    "synthetic",
    "config.yml",
    package = "cpopvar"
  )
  config <- yaml::read_yaml(config_path)

  expect_named(
    get_color_palette(config, "global_heatmap"),
    "gradient"
  )
  expect_named(
    get_color_palette(config, "global_pca"),
    c("scree_bar", "scree_line", "loadings_points")
  )
  expect_named(
    get_color_palette(config, "upset"),
    c("sets_bar_color", "main_bar_color", "matrix_dot_color")
  )
  expect_named(
    get_color_palette(config, "hotspot_correlation"),
    "gradient"
  )
  expect_named(
    get_color_palette(config, "gene_category", c("A", "B")),
    c("A", "B")
  )
  expect_named(
    get_color_palette(config, "M03_hotspot"),
    "primary"
  )
})

test_that("module-declared failures propagate to task and command status", {
  expect_error(
    cpopvar:::assert_module_result_success(
      list(success = FALSE, error = "palette contract failed"),
      "M03_hotspot"
    ),
    "M03_hotspot module failed: palette contract failed"
  )
  expect_no_error(
    cpopvar:::assert_module_result_success(
      list(success = TRUE),
      "M03_hotspot"
    )
  )

  expect_error(
    cpopvar:::stop_on_failed_tasks(
      list(
        preprocessing = list(success = TRUE),
        M03_hotspot = list(
          success = FALSE,
          error_message = "palette contract failed"
        )
      )
    ),
    "M03_hotspot: palette contract failed"
  )

  task_result <- cpopvar:::execute_single_task(
    task_name = "broken_task",
    task_config = list(
      task_id = "broken_task",
      module = "unsupported_module",
      parameters = list()
    ),
    execution_context = list(
      metadata = list(
        global_parameters = list(),
        config_data = list()
      ),
      data_cache = list()
    )
  )
  expect_false(task_result$success)
  expect_match(task_result$error_message, "Unknown module: unsupported_module")
})

test_that("hotspot thresholds use positive frequencies and retain ties", {
  data <- data.frame(
    species = rep("demo", 6),
    gene = letters[1:6],
    frequency_per_kb = c(0, 1, 2, 3, 3, 3),
    region_type = "CDS",
    var_type = "snp"
  )

  ranked <- identify_candidate_hotspots(
    normalized_data = data,
    config = list(),
    task_params = list(frequency_percentile = 0.75)
  )

  expect_equal(unique(ranked$threshold), 3)
  expect_equal(sum(ranked$is_potential_hotspot), 3)
  expect_false(any(ranked$is_potential_hotspot[ranked$frequency_per_kb == 0]))
})
