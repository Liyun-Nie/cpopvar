############################################################
#### task_parser.R - Task-Driven Configuration Parser ####
############################################################
#
# Advanced task parsing and validation system for the task-driven
# analysis pipeline. Handles YAML task configuration parsing,
# dependency resolution, parameter validation, and execution planning.
#
# Key Features:
# - YAML task configuration parsing
# - Task dependency graph resolution
# - Parameter validation with type checking
# - Task execution planning and scheduling
# - Error handling and user-friendly validation messages
# - Task naming convention validation
#
############################################################

# Load required libraries

# Source required utilities
# source("src/utils/common_utils.R")
# source("src/utils/config_validator.R")

# Note: Using null coalescing operator (%||%) exported from annotation_processing.R

#' Parse and validate task-driven configuration with execution mode support
#' 
#' @param config_data Pre-loaded configuration data object (V3.14: Single Source of Truth)
#' @param validate_files Whether to validate file existence
#' @param validate_dependencies Whether to validate task dependencies
#' @param override_mode Optional execution mode override ("full_analysis" or "visualization_only")
#' @return List containing parsed tasks, execution plan, and validation results
#' @export
parse_task_configuration <- function(config_data,
                                   validate_files = TRUE,
                                   validate_dependencies = TRUE,
                                   override_mode = NULL) {
  
  log_message("=== Starting Task-Driven Configuration Parsing ===")
  log_message("Using provided config_data object (V3.14: Single Source of Truth)")
  
  # V3.14: Use provided config_data directly instead of re-reading file
  # This ensures config_file_path injected in run_analysis is preserved
  
  # V3.16.14: Simplified - always use full_analysis mode
  execution_mode_info <- list(
    mode = "full_analysis",
    mode_changed = FALSE,
    config_data = config_data,
    adjustments = character(0)
  )
  
  log_message("Execution mode: full_analysis (V3.16 simplified config)")
  
  # Validate configuration structure
  config_validation <- validate_task_config_structure(config_data)
  if (!config_validation$valid) {
    stop(paste("Configuration structure validation failed:", 
               paste(config_validation$errors, collapse = "; ")), call. = FALSE)
  }
  
  # V3.16.14: No execution mode adjustments needed - always full_analysis
  
  # Parse individual tasks
  task_parsing_result <- parse_analysis_tasks(config_data$analysis_tasks)
  
  # Validate task dependencies
  if (validate_dependencies) {
    dependency_validation <- validate_task_dependencies(task_parsing_result$tasks)
    if (!dependency_validation$valid) {
      stop(paste("Task dependency validation failed:", 
                 paste(dependency_validation$errors, collapse = "; ")), call. = FALSE)
    }
  } else {
    dependency_validation <- list(valid = TRUE, message = "Dependency validation skipped for visualization-only mode")
  }
  
  # Create execution plan
  execution_plan <- create_task_execution_plan(task_parsing_result$tasks, config_data$task_execution)
  
  # V3.16.14: File validation simplified - no mode-aware logic needed
  file_validation <- NULL
  if (validate_files) {
    file_validation <- list(valid = TRUE, message = "File validation simplified in V3.16")
  }
  
  log_message(sprintf("Successfully parsed %d analysis tasks", length(task_parsing_result$tasks)))
  log_message(sprintf("Execution plan created with %d execution stages", length(execution_plan$stages)))
  
  return(list(
    config_data = config_data,
    tasks = task_parsing_result$tasks,
    task_metadata = task_parsing_result$metadata,
    execution_plan = execution_plan,
    execution_mode = execution_mode_info$mode,
    mode_adjustments = execution_mode_info$adjustments,
    validations = list(
      config_structure = config_validation,
      dependencies = dependency_validation,
      files = file_validation
    ),
    parsing_time = Sys.time()
  ))
}

#' Validate task configuration structure
#' 
#' @param config_data Parsed YAML configuration data
#' @return List with validation results
validate_task_config_structure <- function(config_data) {
  
  log_message("Validating task configuration structure")
  
  errors <- character(0)
  warnings <- character(0)
  
  # Check required top-level sections
  required_sections <- c("session_info", "input_files", "column_mappings", "analysis_tasks")
  missing_sections <- setdiff(required_sections, names(config_data))
  if (length(missing_sections) > 0) {
    errors <- c(errors, paste("Missing required sections:", paste(missing_sections, collapse = ", ")))
  }
  
  # Validate analysis_tasks section
  if ("analysis_tasks" %in% names(config_data)) {
    if (!is.list(config_data$analysis_tasks) || length(config_data$analysis_tasks) == 0) {
      errors <- c(errors, "analysis_tasks section must be a non-empty list")
    } else {
      # V3.16: Minimal validation - task names are module names
      # Only check for 'enabled' field
      for (task_name in names(config_data$analysis_tasks)) {
        task <- config_data$analysis_tasks[[task_name]]
        
        if (!"enabled" %in% names(task)) {
          warnings <- c(warnings, paste("Task", task_name, "missing 'enabled' field, assuming enabled=TRUE"))
        }
      }
    }
  }
  
  # Validate task_execution section
  if ("task_execution" %in% names(config_data)) {
    execution_config <- config_data$task_execution
    
    if ("execution_mode" %in% names(execution_config)) {
      if (!execution_config$execution_mode %in% c("sequential", "parallel")) {
        errors <- c(errors, "execution_mode must be 'sequential' or 'parallel'")
      }
    }
  }
  
  return(list(
    valid = length(errors) == 0,
    errors = errors,
    warnings = warnings,
    sections_found = names(config_data)
  ))
}


#' Parse individual analysis tasks
#' 
#' @param tasks_config Analysis tasks configuration
#' @return List containing parsed tasks and metadata
parse_analysis_tasks <- function(tasks_config) {
  
  log_message("Parsing individual analysis tasks")
  
  parsed_tasks <- list()
  task_metadata <- list(
    total_tasks = length(tasks_config),
    enabled_tasks = 0,
    disabled_tasks = 0,
    modules_used = character(0),
    dependency_count = 0
  )
  
  # Handle both named list and array structures
  if (is.null(names(tasks_config))) {
    # Array structure: use task_id as name
    for (i in seq_along(tasks_config)) {
      task_config <- tasks_config[[i]]
      task_name <- task_config$task_id
      
      # Parse task
      parsed_task <- parse_single_task(task_name, task_config)
      parsed_tasks[[task_name]] <- parsed_task
      
      # Update metadata
      if (parsed_task$enabled) {
        task_metadata$enabled_tasks <- task_metadata$enabled_tasks + 1
      } else {
        task_metadata$disabled_tasks <- task_metadata$disabled_tasks + 1
      }
      
      task_metadata$modules_used <- unique(c(task_metadata$modules_used, parsed_task$module))
      
      if (!is.null(parsed_task$depends_on)) {
        task_metadata$dependency_count <- task_metadata$dependency_count + length(parsed_task$depends_on)
      }
    }
  } else {
    # Named list structure: use names directly
    for (task_name in names(tasks_config)) {
      task_config <- tasks_config[[task_name]]
      
      # Parse task
      parsed_task <- parse_single_task(task_name, task_config)
      parsed_tasks[[task_name]] <- parsed_task
      
      # Update metadata
      if (parsed_task$enabled) {
        task_metadata$enabled_tasks <- task_metadata$enabled_tasks + 1
      } else {
        task_metadata$disabled_tasks <- task_metadata$disabled_tasks + 1
      }
      
      task_metadata$modules_used <- unique(c(task_metadata$modules_used, parsed_task$module))
      
      if (!is.null(parsed_task$depends_on)) {
        task_metadata$dependency_count <- task_metadata$dependency_count + length(parsed_task$depends_on)
      }
    }
  }
  
  log_message(sprintf("Parsed %d tasks (%d enabled, %d disabled)", 
                     task_metadata$total_tasks, task_metadata$enabled_tasks, task_metadata$disabled_tasks))
  log_message(sprintf("Modules used: %s", paste(task_metadata$modules_used, collapse = ", ")))
  
  return(list(
    tasks = parsed_tasks,
    metadata = task_metadata
  ))
}

#' Parse a single task configuration
#' 
#' @param task_name Name of the task
#' @param task_config Task configuration
#' @return Parsed task object
parse_single_task <- function(task_name, task_config) {
  
  # V3.16: Simplified - task name IS the module name
  # No inference needed, direct mapping
  
  # V3.16.14: Extract parameters first, then check for data_scoping
  task_parameters <- task_config$parameters %||% list()
  
  # Create standardized task object
  parsed_task <- list(
    task_name = task_name,
    task_id = task_name,  # task_id = task_name = module
    enabled = task_config$enabled %||% TRUE,
    description = sprintf("%s module", task_name),  # Simple default description
    module = task_name,  # Task name IS the module name
    data_source = task_config$data_source %||% "filtered",
    parameters = task_parameters,
    # V3.16.14: data_scoping is INSIDE parameters, not at task level
    # Keep this for backward compatibility, but parameters$data_scoping takes precedence
    data_scoping = task_config$data_scoping %||% list(),
    output = task_config$output %||% list(),
    depends_on = task_config$depends_on %||% NULL,
    validation = task_config$validation %||% list(),
    priority = task_config$priority %||% "medium"
  )
  
  # Validate parameters based on module
  param_validation <- validate_task_parameters(parsed_task)
  parsed_task$parameter_validation = param_validation
  
  return(parsed_task)
}

#' Validate task parameters based on module type
#' 
#' @param task Parsed task object
#' @return Parameter validation results
validate_task_parameters <- function(task) {
  
  module <- task$module
  parameters <- task$parameters
  
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0)
  )
  
  # Module-specific parameter validation
  if (module == "filtering") {
    # Check threshold parameters
    threshold_params <- c("snp_threshold", "indel_threshold", "complex_threshold", "mnp_threshold")
    for (param in threshold_params) {
      if (param %in% names(parameters)) {
        value <- parameters[[param]]
        if (!is.numeric(value) || value < 0) {
          validation_results$errors <- c(validation_results$errors, 
                                       paste("Parameter", param, "must be a non-negative number"))
        }
      }
    }
  } else if (module == "normalization") {
    # Check frequency unit
    if ("frequency_unit" %in% names(parameters)) {
      if (!is.numeric(parameters$frequency_unit) || parameters$frequency_unit <= 0) {
        validation_results$errors <- c(validation_results$errors, 
                                     "frequency_unit must be a positive number")
      }
    }
  } else if (startsWith(module, "M0")) {
    # Visualization module parameters
    if ("plot_types" %in% names(parameters)) {
      valid_plot_types <- c("histogram", "boxplot", "violin", "heatmap", "stacked_bar", "correlation")
      invalid_types <- setdiff(parameters$plot_types, valid_plot_types)
      if (length(invalid_types) > 0) {
        validation_results$warnings <- c(validation_results$warnings,
                                       paste("Unknown plot types:", paste(invalid_types, collapse = ", ")))
      }
    }
  }
  
  validation_results$valid <- length(validation_results$errors) == 0
  
  return(validation_results)
}

#' Validate task dependencies
#' 
#' @param tasks List of parsed tasks
#' @return Dependency validation results
validate_task_dependencies <- function(tasks) {
  
  log_message("Validating task dependencies")
  
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    dependency_graph = NULL
  )
  
  # Check for circular dependencies and missing dependencies
  enabled_tasks <- names(tasks)[sapply(tasks, function(t) t$enabled)]
  
  # Build dependency graph
  edges <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
  
  for (task_name in enabled_tasks) {
    task <- tasks[[task_name]]
    if (!is.null(task$depends_on)) {
      for (dependency in task$depends_on) {
        # Check if dependency exists and is enabled
        if (!dependency %in% enabled_tasks) {
          if (dependency %in% names(tasks)) {
            validation_results$warnings <- c(validation_results$warnings,
                                           paste("Task", task_name, "depends on disabled task", dependency))
          } else {
            validation_results$errors <- c(validation_results$errors,
                                         paste("Task", task_name, "depends on non-existent task", dependency))
          }
        } else {
          edges <- rbind(edges, data.frame(from = dependency, to = task_name, stringsAsFactors = FALSE))
        }
      }
    }
  }
  
  # Check for circular dependencies using igraph
  if (nrow(edges) > 0) {
    tryCatch({
      dependency_graph <- igraph::graph_from_data_frame(edges, directed = TRUE)
      
      # Check for cycles
      if (!igraph::is_dag(dependency_graph)) {
        validation_results$errors <- c(validation_results$errors, "Circular dependencies detected in task graph")
      } else {
        validation_results$dependency_graph <- dependency_graph
      }
    }, error = function(e) {
      validation_results$warnings <- c(validation_results$warnings, 
                                     paste("Could not build dependency graph:", e$message))
    })
  }
  
  validation_results$valid <- length(validation_results$errors) == 0
  
  return(validation_results)
}

#' Create task execution plan
#' 
#' @param tasks List of parsed tasks
#' @param execution_config Task execution configuration
#' @return Execution plan with ordered task stages
create_task_execution_plan <- function(tasks, execution_config) {
  
  log_message("Creating task execution plan")
  
  enabled_tasks <- tasks[sapply(tasks, function(t) t$enabled)]
  
  if (length(enabled_tasks) == 0) {
    return(list(stages = list(), total_tasks = 0, execution_mode = "sequential"))
  }
  
  execution_mode <- execution_config$execution_mode %||% "sequential"
  
  if (execution_mode == "sequential") {
    execution_plan <- create_sequential_execution_plan(enabled_tasks)
  } else {
    execution_plan <- create_parallel_execution_plan(enabled_tasks)
  }
  
  log_message(sprintf("Execution plan created: %d stages, %d total tasks", 
                     length(execution_plan$stages), execution_plan$total_tasks))
  
  return(execution_plan)
}

#' Create sequential execution plan with dependency resolution
#' 
#' @param tasks List of enabled tasks
#' @return Sequential execution plan
create_sequential_execution_plan <- function(tasks) {
  
  # Topological sort based on dependencies
  task_names <- names(tasks)
  dependency_levels <- list()
  processed_tasks <- character(0)
  current_level <- 0
  
  while (length(processed_tasks) < length(task_names)) {
    current_level <- current_level + 1
    level_tasks <- character(0)
    
    for (task_name in task_names) {
      if (task_name %in% processed_tasks) next
      
      task <- tasks[[task_name]]
      dependencies <- task$depends_on %||% character(0)
      
      # Check if all dependencies are already processed
      if (all(dependencies %in% processed_tasks)) {
        level_tasks <- c(level_tasks, task_name)
      }
    }
    
    if (length(level_tasks) == 0) {
      # Should not happen if dependency validation passed
      stop("Cannot resolve task dependencies - possible circular dependency")
    }
    
    dependency_levels[[current_level]] <- level_tasks
    processed_tasks <- c(processed_tasks, level_tasks)
  }
  
  # Create execution stages
  stages <- list()
  for (level in seq_along(dependency_levels)) {
    stage_tasks <- dependency_levels[[level]]
    stages[[level]] <- list(
      stage_id = level,
      tasks = stage_tasks,
      execution_mode = "sequential",
      estimated_duration = sum(sapply(stage_tasks, function(tn) estimate_task_duration(tasks[[tn]])))
    )
  }
  
  return(list(
    stages = stages,
    total_tasks = length(task_names),
    execution_mode = "sequential",
    estimated_total_duration = sum(sapply(stages, function(s) s$estimated_duration))
  ))
}

#' Create parallel execution plan (placeholder)
#' 
#' @param tasks List of enabled tasks
#' @return Parallel execution plan
create_parallel_execution_plan <- function(tasks) {
  
  # For now, convert to sequential plan
  # Future enhancement: implement true parallel execution with dependency constraints
  log_message("Parallel execution mode not fully implemented, falling back to sequential")
  
  sequential_plan <- create_sequential_execution_plan(tasks)
  sequential_plan$execution_mode <- "parallel_fallback"
  
  return(sequential_plan)
}

#' Estimate task execution duration
#' 
#' @param task Task object
#' @return Estimated duration in minutes
estimate_task_duration <- function(task) {
  
  # Simple duration estimation based on module type
  module <- task$module
  
  duration_estimates <- list(
    "preprocessing" = 5,
    "filtering" = 10,
    "normalization" = 8,
    "M01_distribution" = 15,
    "M02_comparative" = 20,
    "M03_hotspot" = 12
  )
  
  return(duration_estimates[[module]] %||% 10) # Default 10 minutes
}

#' Validate input files configuration
#' 
#' @param input_files_config Input files configuration
#' @return File validation results
validate_input_files <- function(input_files_config) {
  
  log_message("Validating input files")
  
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    file_status = list()
  )
  
  # V3.16: Simplified validation for flat file structure
  # input_files now contains direct file paths, not nested structures
  
  # Skip validation - files will be validated during copy_auxiliary_files
  # This function is kept for backward compatibility but is now a no-op
  
  log_message("File validation skipped - will be performed during file copying")
  
  return(validation_results)
}

#' Get execution plan summary
#' 
#' @param execution_plan Task execution plan
#' @return Human-readable execution summary
#' @export
get_execution_plan_summary <- function(execution_plan) {
  
  summary_lines <- character(0)
  
  summary_lines <- c(summary_lines, 
                    paste("Execution Plan Summary:"))
  summary_lines <- c(summary_lines,
                    paste("- Total tasks:", execution_plan$total_tasks))
  summary_lines <- c(summary_lines,
                    paste("- Execution stages:", length(execution_plan$stages)))
  summary_lines <- c(summary_lines,
                    paste("- Execution mode:", execution_plan$execution_mode))
  
  if ("estimated_total_duration" %in% names(execution_plan)) {
    summary_lines <- c(summary_lines,
                      paste("- Estimated duration:", execution_plan$estimated_total_duration, "minutes"))
  }
  
  summary_lines <- c(summary_lines, "")
  summary_lines <- c(summary_lines, "Stage Details:")
  
  for (i in seq_along(execution_plan$stages)) {
    stage <- execution_plan$stages[[i]]
    summary_lines <- c(summary_lines,
                      paste("  Stage", i, ":", paste(stage$tasks, collapse = ", ")))
  }
  
  return(paste(summary_lines, collapse = "\n"))
}

#' Detect execution mode from configuration
#' 
#' @param config_data Parsed configuration data
#' @param override_mode Optional mode override
#' @return List with mode information and adjusted config
detect_execution_mode <- function(config_data, override_mode = NULL) {
  
  # V3.16: Simplified configuration - always use full_analysis mode
  # The simplified config removes execution_mode and visualization_only support
  # to reduce complexity and ensure users run the complete pipeline
  
  log_message("Detecting execution mode")
  
  # Initialize result structure
  mode_info <- list(
    mode = "full_analysis",  # V3.16: Always full_analysis
    mode_changed = FALSE,
    config_data = config_data,
    adjustments = character(0)
  )
  
  # Check for override mode (for backward compatibility with tests)
  if (!is.null(override_mode)) {
    if (override_mode %in% c("full_analysis", "visualization_only")) {
      mode_info$mode <- override_mode
      mode_info$mode_changed <- TRUE
      mode_info$adjustments <- c(mode_info$adjustments, paste("Mode overridden to", override_mode))
      log_message(sprintf("Execution mode overridden to: %s", override_mode))
    } else {
      log_message(sprintf("Invalid override mode '%s', using default 'full_analysis'", override_mode))
    }
    } else {
    log_message("Execution mode: full_analysis (V3.16 simplified config)")
  }
  
  # Apply mode-specific configuration adjustments
  mode_info$config_data <- apply_mode_to_config(config_data, mode_info$mode)
  
  return(mode_info)
}

#' Auto-detect execution mode from task configuration
#' 
#' @param tasks_config Analysis tasks configuration
#' @return List with detection results
auto_detect_mode_from_tasks <- function(tasks_config) {
  
  if (is.null(tasks_config) || length(tasks_config) == 0) {
    return(list(detected = FALSE, mode = "full_analysis", reason = "No tasks configured"))
  }
  
  # Count enabled core processing tasks
  core_tasks <- c("core_preprocessing", "core_filtering", "core_normalization")
  enabled_core_count <- 0
  total_core_count <- 0
  
  # Count visualization tasks with use_existing_data parameter
  viz_with_existing_data <- 0
  total_viz_tasks <- 0
  
  for (task_name in names(tasks_config)) {
    task <- tasks_config[[task_name]]
    
    # Check core processing tasks
    if (task_name %in% core_tasks || 
        (is.list(task) && "module" %in% names(task) && 
         task$module %in% c("preprocessing", "filtering", "normalization"))) {
      total_core_count <- total_core_count + 1
      if (is.list(task) && "enabled" %in% names(task) && task$enabled == TRUE) {
        enabled_core_count <- enabled_core_count + 1
      }
    }
    
    # Check visualization tasks
    if (is.list(task) && "module" %in% names(task) && 
        startsWith(task$module, "M0")) {
      total_viz_tasks <- total_viz_tasks + 1
      if ("parameters" %in% names(task) && 
          "use_existing_data" %in% names(task$parameters) &&
          task$parameters$use_existing_data == TRUE) {
        viz_with_existing_data <- viz_with_existing_data + 1
      }
    }
  }
  
  # Detection logic
  if (total_core_count > 0 && enabled_core_count == 0 && viz_with_existing_data > 0) {
    return(list(
      detected = TRUE, 
      mode = "visualization_only",
      reason = "Core processing tasks disabled, visualization tasks have use_existing_data"
    ))
  }
  
  if (viz_with_existing_data > 0 && viz_with_existing_data == total_viz_tasks) {
    return(list(
      detected = TRUE,
      mode = "visualization_only", 
      reason = "All visualization tasks configured with use_existing_data"
    ))
  }
  
  return(list(detected = FALSE, mode = "full_analysis", reason = "Default mode"))
}

#' Apply execution mode to configuration
#' 
#' @param config_data Configuration data
#' @param mode Execution mode
#' @return Adjusted configuration data
apply_mode_to_config <- function(config_data, mode) {
  
  # Ensure workflow section exists
  if (!"workflow" %in% names(config_data)) {
    config_data$workflow <- list()
  }
  
  # Set workflow flags based on mode
  if (mode == "visualization_only") {
    config_data$workflow$visualization_only <- TRUE
    
    # Adjust output settings for visualization-only mode
    if ("output" %in% names(config_data)) {
      config_data$output$save_csv <- FALSE  # Don't save CSV in viz-only mode
    }
    
    # Adjust task execution settings
    if (!"task_execution" %in% names(config_data)) {
      config_data$task_execution <- list()
    }
    if (!"visualization_mode" %in% names(config_data$task_execution)) {
      config_data$task_execution$visualization_mode <- list()
    }
    config_data$task_execution$visualization_mode$enabled <- TRUE
    
  } else {
    config_data$workflow$visualization_only <- FALSE
    
    # Ensure task execution settings for full analysis
    if (!"task_execution" %in% names(config_data)) {
      config_data$task_execution <- list()
    }
    if ("visualization_mode" %in% names(config_data$task_execution)) {
      config_data$task_execution$visualization_mode$enabled <- FALSE
    }
  }
  
  return(config_data)
}

#' Apply execution mode adjustments to tasks configuration
#' 
#' @param tasks_config Analysis tasks configuration
#' @param mode Execution mode
#' @return Adjusted tasks configuration
apply_execution_mode_to_tasks <- function(tasks_config, mode) {
  
  # V3.16: Simplified configuration - mode is always full_analysis
  # No need for mode-specific adjustments
  
  log_message(sprintf("Applying execution mode adjustments: %s", mode))
  
  if (is.null(tasks_config) || length(tasks_config) == 0) {
    return(tasks_config)
  }
  
  adjusted_tasks <- tasks_config
  
  # V3.16: In simplified config, all tasks are treated equally
  # No execution_modes field to check, no mode-specific adjustments needed
  for (task_name in names(tasks_config)) {
    task <- tasks_config[[task_name]]
    
    if (!is.list(task)) next
    
    # Apply mode-specific adjustments (for backward compatibility)
    if (mode == "visualization_only") {
      adjusted_tasks[[task_name]] <- apply_visualization_only_adjustments(task, task_name)
    } else if (mode == "full_analysis") {
      adjusted_tasks[[task_name]] <- apply_full_analysis_adjustments(task, task_name)
    }
  }
  
  return(adjusted_tasks)
}

#' Apply visualization-only mode adjustments to a task
#' 
#' @param task Task configuration
#' @param task_name Task name
#' @return Adjusted task configuration
apply_visualization_only_adjustments <- function(task, task_name) {
  
  # Get module type
  module <- task$module %||% "unknown"
  
  # Disable core processing tasks
  if (module %in% c("preprocessing", "filtering", "normalization") ||
      task_name %in% c("core_preprocessing", "core_filtering", "core_normalization")) {
    task$enabled <- FALSE
    task$description <- paste("SKIPPED -", task$description %||% "Using existing processed data")
    return(task)
  }
  
  # Adjust visualization tasks
  if (startsWith(module, "M0") || module %in% c("M01_distribution", "M02_comparative", "M03_hotspot")) {
    # Add use_existing_data parameter
    if (!"parameters" %in% names(task)) {
      task$parameters <- list()
    }
    task$parameters$use_existing_data <- TRUE
    
    # Remove dependencies for visualization-only mode
    task$depends_on <- NULL
    
    # Ensure task is enabled if it supports visualization-only mode
    if ("execution_modes" %in% names(task) && "visualization_only" %in% task$execution_modes) {
      if (!"force_enabled" %in% names(task) || task$force_enabled != FALSE) {
        task$enabled <- TRUE
      }
    }
  }
  
  return(task)
}

#' Apply full analysis mode adjustments to a task
#' 
#' @param task Task configuration
#' @param task_name Task name
#' @return Adjusted task configuration
apply_full_analysis_adjustments <- function(task, task_name) {
  
  # Remove visualization-only parameters
  if ("parameters" %in% names(task) && "use_existing_data" %in% names(task$parameters)) {
    task$parameters$use_existing_data <- NULL
  }
  
  # Ensure core processing tasks are enabled (unless explicitly disabled)
  module <- task$module %||% "unknown"
  if (module %in% c("preprocessing", "filtering", "normalization") ||
      task_name %in% c("core_preprocessing", "core_filtering", "core_normalization")) {
    if (!"force_enabled" %in% names(task) || task$force_enabled != FALSE) {
      task$enabled <- TRUE
    }
  }
  
  return(task)
}

#' Mode-aware input files validation
#' 
#' @param input_files_config Input files configuration
#' @param mode Execution mode
#' @return Validation results
validate_input_files_mode_aware <- function(input_files_config, mode) {
  
  log_message(sprintf("Validating input files for mode: %s", mode))
  
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    file_status = list(),
    mode = mode
  )
  
  if (mode == "visualization_only") {
    # In visualization-only mode, we mainly need auxiliary files
    # Main data and region data files are optional since we use existing processed data
    
    # Check auxiliary files
    if ("auxiliary" %in% names(input_files_config)) {
      auxiliary_files <- input_files_config$auxiliary
      for (file_type in names(auxiliary_files)) {
        file_config <- auxiliary_files[[file_type]]
        if ("path" %in% names(file_config)) {
          file_path <- file_config$path
          if (!file.exists(file_path)) {
            if (file_config$required %||% FALSE) {
              validation_results$errors <- c(validation_results$errors,
                                           paste("Required auxiliary file not found:", file_path))
            } else {
              validation_results$warnings <- c(validation_results$warnings,
                                             paste("Optional auxiliary file not found:", file_path))
            }
          } else {
            validation_results$file_status[[file_type]] <- "found"
          }
        }
      }
    }
    
    # Check for existing processed data (this would be implemented elsewhere)
    validation_results$warnings <- c(validation_results$warnings,
                                   "Note: Visualization-only mode requires existing processed data in session directory")
    
  } else {
    # Full analysis mode - use standard validation
    validation_results <- validate_input_files(input_files_config)
    validation_results$mode <- mode
  }
  
  validation_results$valid <- length(validation_results$errors) == 0
  
  return(validation_results)
}

# log_message("Task-driven configuration parser with execution mode support loaded successfully")