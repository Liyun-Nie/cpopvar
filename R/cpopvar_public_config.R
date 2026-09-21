cpopvar_public_task_names <- function() {
  c(
    "preprocessing",
    "filtering",
    "normalization",
    "M01_distribution",
    "M02_comparative",
    "M03_hotspot",
    "M04_igs_hotspot",
    "M05_hotspot_ideogram",
    "reporting"
  )
}

cpopvar_core_task_names <- function() {
  c("preprocessing", "filtering", "normalization")
}

cpopvar_task_aliases <- function() {
  c(
    core_preprocessing = "preprocessing",
    core_filtering = "filtering",
    core_normalization = "normalization",
    M01_data_distribution = "M01_distribution",
    M02_comparative_analysis = "M02_comparative",
    M03_hotspot_analysis = "M03_hotspot",
    M04_igs_analysis = "M04_igs_hotspot",
    M05_ideogram = "M05_hotspot_ideogram"
  )
}

canonical_cpopvar_task_name <- function(task_name) {
  aliases <- cpopvar_task_aliases()
  if (task_name %in% names(aliases)) {
    return(unname(aliases[[task_name]]))
  }
  task_name
}

default_cpopvar_column_mappings <- function(input_files) {
  mappings <- list(
    main_data = list(
      sample_id = "sample_id",
      species = "species",
      position = "position",
      var_type = "var_type",
      gene = "gene",
      region_type = "region_type"
    )
  )
  if (!is.null(input_files$genome_regions)) {
    mappings$genome_regions <- list(
      species = "species",
      lsc_start = "lsc_start",
      lsc_end = "lsc_end",
      ir_start = "ir_start",
      ir_end = "ir_end",
      ssc_start = "ssc_start",
      ssc_end = "ssc_end",
      total_length = "total_length",
      special_handling = "special_handling"
    )
  }
  if (!is.null(input_files$cds_lengths)) {
    mappings$cds_lengths <- list(
      gene = "gene",
      length = "length"
    )
  }
  if (!is.null(input_files$region_info)) {
    mappings$region_info <- list(
      species = "species",
      region_name = "region_name",
      region_type = "region_type",
      start_pos = "start_pos",
      end_pos = "end_pos"
    )
  }
  mappings
}

new_cpopvar_config_error <- function(code, path, message) {
  list(code = code, path = path, message = message)
}

format_cpopvar_config_errors <- function(errors) {
  if (length(errors) == 0) {
    return("Configuration validation failed.")
  }
  lines <- vapply(errors, function(err) {
    sprintf("  [%s] %s: %s", err$code, err$path, err$message)
  }, character(1))
  paste(c("Configuration validation failed:", lines), collapse = "\n")
}

as_cpopvar_config <- function(config_data) {
  structure(config_data, class = c("cpopvar_config", "list"))
}

collect_enabled_task_names <- function(tasks) {
  if (is.null(tasks) || !is.list(tasks) || length(tasks) == 0) {
    return(character(0))
  }
  names(tasks)[vapply(tasks, function(task) {
    isTRUE(task$enabled %||% TRUE)
  }, logical(1))]
}

execution_plan_task_names <- function(execution_plan) {
  if (is.null(execution_plan) || is.null(execution_plan$stages)) {
    return(character(0))
  }
  unique(unlist(lapply(execution_plan$stages, function(stage) stage$tasks), use.names = FALSE))
}

#' Build an in-memory cpopvar configuration
#'
#' Construct a YAML-isomorphic configuration list from stable workflow inputs.
#' The YAML file remains the serialized configuration contract; this helper
#' only builds the same nested structure in memory.
#'
#' Core preprocessing, filtering, and normalisation tasks are enabled whenever
#' any visualisation task is requested. Module-specific options belong in
#' `task_options`, not in `run_analysis()`.
#'
#' @param main_data Path to the main variant table.
#' @param annotations Optional annotation table path.
#' @param genome_regions Optional genome-region table path.
#' @param cds_lengths Optional CDS length table path.
#' @param region_info Optional region-info table path.
#' @param group_info Optional grouping table path.
#' @param species_order Optional species-order table path.
#' @param gene_function_mapping Optional gene-function mapping path.
#' @param tasks Character vector of task names to enable. Aliases such as
#'   `M01_data_distribution` are accepted and mapped to the canonical YAML
#'   names.
#' @param task_options Named list of per-task parameter lists. Names must be
#'   canonical task names or aliases; values become `analysis_tasks[[task]]$parameters`.
#' @param column_mappings Optional nested column-mapping list. Auxiliary files
#'   such as `group_info` and `species_order` must not be mapped.
#' @param session_id Optional session identifier.
#' @param visualization_settings Optional visualization settings list copied
#'   into the config unchanged.
#' @param output_settings Optional output settings list copied into the config
#'   unchanged.
#' @param ... Must remain empty. Unknown arguments are rejected so module
#'   parameters cannot be flattened onto this constructor.
#'
#' @return A `cpopvar_config` list that can be passed to
#'   [validate_cpopvar_config()] or [run_analysis()].
#' @export
#' @examples
#' \dontrun{
#' cfg <- cpopvar_config(
#'   main_data = "variants.csv",
#'   annotations = "annotations.csv",
#'   genome_regions = "genome_regions.csv",
#'   tasks = c("M01_distribution", "M03_hotspot")
#' )
#' result <- run_analysis(cfg, output_dir = tempdir())
#' }
cpopvar_config <- function(main_data,
                           annotations = NULL,
                           genome_regions = NULL,
                           cds_lengths = NULL,
                           region_info = NULL,
                           group_info = NULL,
                           species_order = NULL,
                           gene_function_mapping = NULL,
                           tasks = c("M01_distribution", "M03_hotspot"),
                           task_options = list(),
                           column_mappings = NULL,
                           session_id = NULL,
                           visualization_settings = NULL,
                           output_settings = NULL,
                           ...) {
  extra <- list(...)
  if (length(extra) > 0) {
    extra_names <- names(extra)
    if (is.null(extra_names) || any(!nzchar(extra_names))) {
      stop(
        "Unknown positional arguments in cpopvar_config(). Module options belong in task_options.",
        call. = FALSE
      )
    }
    stop(
      paste0(
        "Unknown arguments in cpopvar_config(): ",
        paste(extra_names, collapse = ", "),
        ". Module-specific options belong in task_options."
      ),
      call. = FALSE
    )
  }

  if (missing(main_data) || is.null(main_data) || !nzchar(as.character(main_data)[1])) {
    stop("main_data is required.", call. = FALSE)
  }
  if (!is.character(tasks) || length(tasks) < 1) {
    stop("`tasks` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.list(task_options)) {
    stop("`task_options` must be a list.", call. = FALSE)
  }

  input_files <- list(main_data = as.character(main_data)[1])
  optional_inputs <- list(
    annotations = annotations,
    genome_regions = genome_regions,
    cds_lengths = cds_lengths,
    region_info = region_info,
    group_info = group_info,
    species_order = species_order,
    gene_function_mapping = gene_function_mapping
  )
  for (name in names(optional_inputs)) {
    value <- optional_inputs[[name]]
    if (!is.null(value) && nzchar(as.character(value)[1])) {
      input_files[[name]] <- as.character(value)[1]
    }
  }

  requested <- unique(vapply(tasks, canonical_cpopvar_task_name, character(1)))
  unknown <- setdiff(requested, cpopvar_public_task_names())
  if (length(unknown) > 0) {
    stop(
      paste0("Unknown task(s): ", paste(unknown, collapse = ", ")),
      call. = FALSE
    )
  }
  if (any(!requested %in% cpopvar_core_task_names())) {
    requested <- unique(c(cpopvar_core_task_names(), requested))
  }

  option_names <- names(task_options)
  if (is.null(option_names)) {
    option_names <- character(0)
  }
  canonical_options <- list()
  for (name in option_names) {
    canonical <- canonical_cpopvar_task_name(name)
    if (!canonical %in% cpopvar_public_task_names()) {
      stop(
        paste0("Unknown task_options name: ", name),
        call. = FALSE
      )
    }
    canonical_options[[canonical]] <- task_options[[name]]
  }

  analysis_tasks <- lapply(cpopvar_public_task_names(), function(task_name) {
    params <- canonical_options[[task_name]] %||% list()
    list(
      enabled = task_name %in% requested,
      parameters = params
    )
  })
  names(analysis_tasks) <- cpopvar_public_task_names()

  if (is.null(column_mappings)) {
    column_mappings <- default_cpopvar_column_mappings(input_files)
  }

  config_data <- list(
    session_info = list(
      session_id = session_id %||% paste0("session_", format(Sys.time(), "%Y%m%d_%H%M%S")),
      run_name = "cpopvar_config",
      description = "In-memory configuration generated by cpopvar_config()"
    ),
    input_files = input_files,
    column_mappings = column_mappings,
    analysis_tasks = analysis_tasks
  )
  if (!is.null(visualization_settings)) {
    config_data$visualization_settings <- visualization_settings
  }
  if (!is.null(output_settings)) {
    config_data$output_settings <- output_settings
  }

  validated <- validate_cpopvar_config(config_data, error = TRUE)
  validated$config
}

#' Validate a cpopvar configuration
#'
#' Validate a YAML path or in-memory configuration without writing a temporary
#' file. File-path and memory inputs use the same structure contract as the
#' parser and dispatcher.
#'
#' @param config A YAML file path or an in-memory configuration list.
#' @param error If `TRUE`, invalid configurations raise an error instead of
#'   returning the structured result.
#'
#' @return A list with `valid`, `errors`, `warnings`, and `config`. Each error
#'   is a list with `code`, `path`, and `message`. `config` is a
#'   `cpopvar_config` object when validation succeeds, otherwise `NULL`.
#' @export
#' @examples
#' \dontrun{
#' result <- validate_cpopvar_config(cfg)
#' result$valid
#' }
validate_cpopvar_config <- function(config, error = FALSE) {
  errors <- list()
  warnings <- list()
  config_data <- NULL

  append_error <- function(code, path, message) {
    errors <<- c(errors, list(new_cpopvar_config_error(code, path, message)))
  }
  append_warning <- function(code, path, message) {
    warnings <<- c(warnings, list(new_cpopvar_config_error(code, path, message)))
  }

  if (is.character(config)) {
    if (length(config) != 1 || is.na(config) || !nzchar(config)) {
      append_error("invalid_type", "config", "Configuration file path must be a single non-empty string.")
    } else if (!file.exists(config)) {
      append_error("file_not_found", "config", paste("Configuration file not found:", config))
    } else {
      loaded <- tryCatch(
        yaml::read_yaml(config),
        error = function(e) e
      )
      if (inherits(loaded, "error")) {
        append_error("parse_error", "config", paste("Failed to parse YAML:", loaded$message))
      } else {
        config_data <- loaded
      }
    }
  } else if (is.list(config)) {
    config_data <- unclass(config)
  } else {
    append_error(
      "invalid_type",
      "config",
      "config must be a YAML file path or an in-memory configuration list."
    )
  }

  if (!is.null(config_data)) {
    structure_result <- validate_task_config_structure(config_data)
    for (msg in structure_result$errors) {
      append_error("invalid_structure", "config", msg)
    }
    for (msg in structure_result$warnings) {
      append_warning("structure_warning", "config", msg)
    }

    if (is.null(config_data$input_files) || is.null(config_data$input_files$main_data)) {
      append_error("missing_input", "input_files.main_data", "main_data is required.")
    }

    task_names <- names(config_data$analysis_tasks)
    if (!is.null(task_names)) {
      canonical <- vapply(task_names, canonical_cpopvar_task_name, character(1))
      unknown <- unique(setdiff(canonical, cpopvar_public_task_names()))
      if (length(unknown) > 0) {
        append_error(
          "unknown_task",
          "analysis_tasks",
          paste("Unknown task(s):", paste(unknown, collapse = ", "))
        )
      }
    }

    mapped_aux <- intersect(names(config_data$column_mappings), c("group_info", "species_order"))
    if (length(mapped_aux) > 0) {
      append_error(
        "invalid_mapping",
        "column_mappings",
        paste(
          "Auxiliary files must not use column mappings:",
          paste(mapped_aux, collapse = ", ")
        )
      )
    }
  }

  valid <- length(errors) == 0
  result <- list(
    valid = valid,
    errors = errors,
    warnings = warnings,
    config = if (valid) as_cpopvar_config(config_data) else NULL
  )
  if (isTRUE(error) && !valid) {
    stop(format_cpopvar_config_errors(errors), call. = FALSE)
  }
  result
}

load_cpopvar_config_input <- function(config) {
  if (is.character(config)) {
    if (length(config) != 1 || is.na(config) || !nzchar(config)) {
      stop("Configuration file path must be a single non-empty string.", call. = FALSE)
    }
    if (!file.exists(config)) {
      stop(paste("Configuration file not found:", config), call. = FALSE)
    }
    list(
      data = yaml::read_yaml(config),
      source_path = normalizePath(config, winslash = "/", mustWork = TRUE)
    )
  } else if (is.list(config)) {
    list(data = unclass(config), source_path = NULL)
  } else {
    stop(
      "config must be a YAML file path or an in-memory cpopvar config list.",
      call. = FALSE
    )
  }
}
