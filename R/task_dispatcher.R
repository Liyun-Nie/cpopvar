############################################################
#### task_dispatcher.R - Task Execution Dispatcher ####
############################################################
#
# Task execution engine for the task-driven analysis pipeline.
# Handles task routing, parameter application, result aggregation,
# and error management across multiple analysis tasks.
#
# Key Features:
# - Task-specific parameter routing
# - Module function mapping and execution
# - Result aggregation and cross-task analysis
# - Error handling with task isolation
# - Progress tracking and reporting
# - Task result validation and quality control
#
############################################################

# Load required libraries

# Source required utilities
# source("src/utils/common_utils.R")
# source("src/utils/task_parser.R")
# source("src/functions/task_specific_visualizations.R")
# Dependencies automatically loaded in R package context

#' Dispatch and execute analysis tasks
#' 
#' @param parsed_config Parsed task configuration from task_parser
#' @param session_id Session ID for file paths and logging
#' @param output_dir Optional base directory for session outputs
#' @param progress_callback Optional progress callback function
#' @return List containing task execution results and metadata
#' @importFrom magrittr %>%
#' @export
dispatch_analysis_tasks <- function(parsed_config,
                                   session_id,
                                   output_dir = NULL,
                                   progress_callback = NULL) {
  
  log_message("=== Starting Task-Driven Analysis Execution ===")
  log_message(sprintf("Session ID: %s", session_id))
  
  # Initialize execution context
  execution_context <- initialize_execution_context(parsed_config, session_id, output_dir)
  
  # Execute tasks according to execution plan
  # V3.16.3: Pass parsed_tasks to ensure standardized task structure
  execution_results <- execute_task_stages(
    execution_context,
    parsed_config$execution_plan,
    progress_callback,
    parsed_config$tasks
  )
  
  # Aggregate and validate results
  aggregated_results <- aggregate_task_results(execution_results, parsed_config)
  
  # Generate execution summary
  execution_summary <- generate_execution_summary(execution_results, parsed_config)

  stop_on_failed_tasks(execution_results$task_results)
  
  log_message("=== Task-Driven Analysis Execution Completed ===")
  log_message(sprintf("Executed %d tasks across %d stages", 
                     execution_summary$total_tasks_executed,
                     execution_summary$total_stages))
  
  return(list(
    task_results = execution_results$task_results,
    aggregated_results = aggregated_results,
    execution_summary = execution_summary,
    execution_context = execution_context,
    completion_time = Sys.time()
  ))
}

# Convert module-declared failures into regular R errors so that task and CLI
# status cannot report success when a module returned success = FALSE.
assert_module_result_success <- function(module_result, module_name) {
  if (!is.list(module_result) || !identical(module_result$success, FALSE)) {
    return(invisible(module_result))
  }

  failure_message <- module_result$error_message %||%
    module_result$error %||%
    module_result$summary$error %||%
    "Module returned success = FALSE without an error message"

  stop(
    sprintf("%s module failed: %s", module_name, failure_message),
    call. = FALSE
  )
}

# Preserve task isolation during execution, then fail the overall command after
# all enabled tasks have had a chance to record their diagnostics.
stop_on_failed_tasks <- function(task_results) {
  failed_names <- names(task_results)[vapply(
    task_results,
    function(result) !isTRUE(result$success),
    logical(1)
  )]

  if (length(failed_names) == 0L) {
    return(invisible(task_results))
  }

  failed_details <- vapply(
    failed_names,
    function(task_name) {
      message <- task_results[[task_name]]$error_message %||% "unknown error"
      sprintf("%s: %s", task_name, message)
    },
    character(1)
  )

  stop(
    paste(
      "Analysis failed for one or more tasks:",
      paste(failed_details, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

#' Initialize execution context
#' 
#' @param parsed_config Parsed task configuration
#' @param session_id Session ID
#' @param output_dir Optional base directory for session outputs
#' @return Execution context
initialize_execution_context <- function(parsed_config, session_id, output_dir = NULL) {
  
  log_message("Initializing task execution context")
  
  # Setup session-specific logging
  setup_logging(session_id = session_id, output_dir = output_dir)
  log_message(sprintf("Session-specific logging initialized for session: %s", session_id))
  
  # Get session paths
  session_paths <- get_session_paths(session_id, create_dirs = TRUE, output_dir = output_dir)
  
  # Copy auxiliary files to session directory before any tasks are run
  copy_auxiliary_files(parsed_config$config_data, session_paths)
  
  # SHARED DEPENDENCY INJECTION: Load species genome data once for all tasks
  log_message("Loading shared species genome region data...")
  config_data <- parsed_config$config_data
  
  # V3.16: Simplified configuration - direct file path
  # Get genome_regions file from simplified input_files structure
  species_genome_regions_path <- config_data$input_files$genome_regions
  
  # V3.16.1: Safer check - handle character(0), NULL, NA cases
  has_valid_genome_regions <- !is.null(species_genome_regions_path) && 
                               length(species_genome_regions_path) > 0 && 
                               !is.na(species_genome_regions_path) && 
                               nzchar(species_genome_regions_path)
  
  # If path is NULL or doesn't exist, try to find it in session directory
  if (!has_valid_genome_regions || !file.exists(species_genome_regions_path)) {
    # File should have been copied to session's raw directory
    if (has_valid_genome_regions) {
      session_raw_path <- file.path(session_paths$raw, basename(species_genome_regions_path))
        } else {
      # Try default filename
      session_raw_path <- file.path(session_paths$raw, "species_genome_regions.csv")
    }
    
    # V3.16.1: Safer file.exists check
    if (length(session_raw_path) > 0 && file.exists(session_raw_path)) {
      species_genome_regions_path <- session_raw_path
    } else {
      stop(sprintf("Genome regions file not found. Expected at: %s", session_raw_path))
    }
  }
  
  species_genome_data <- readr::read_csv(species_genome_regions_path, show_col_types = FALSE)
  log_message(sprintf("Successfully loaded %d species configurations from %s", 
                      nrow(species_genome_data), species_genome_regions_path))
  
  # Initialize data cache for cross-task data sharing
  data_cache <- list(
    preprocessing_results = NULL,
    filtering_results = NULL,
    normalization_results = NULL,
    species_genome_data = species_genome_data, # Shared dependency
    session_paths = session_paths
  )
  
  # Task execution metadata
  execution_metadata <- list(
    session_id = session_id,
    start_time = Sys.time(),
    config_data = parsed_config$config_data,
    global_parameters = parsed_config$config_data$global_parameters %||% list(),
    task_count = length(parsed_config$tasks),
    enabled_task_count = sum(sapply(parsed_config$tasks, function(t) t$enabled))
  )
  
  return(list(
    data_cache = data_cache,
    metadata = execution_metadata,
    session_paths = session_paths
  ))
}

#' Execute task stages according to execution plan
#' 
#' @param execution_context Execution context
#' @param execution_plan Task execution plan
#' @param progress_callback Progress callback function
#' @param parsed_tasks Parsed task configurations
#' @return Stage execution results
execute_task_stages <- function(execution_context,
                               execution_plan,
                               progress_callback = NULL,
                               parsed_tasks = NULL) {
  
  log_message("Executing task stages")
  
  stage_results <- list()
  task_results <- list()
  total_stages <- length(execution_plan$stages)
  
  for (stage_idx in seq_along(execution_plan$stages)) {
    stage <- execution_plan$stages[[stage_idx]]
    
    log_message(sprintf("=== Executing Stage %d/%d ===", stage_idx, total_stages))
    log_message(sprintf("Tasks in stage: %s", paste(stage$tasks, collapse = ", ")))
    
    # Update progress
    if (!is.null(progress_callback)) {
      progress <- (stage_idx - 1) / total_stages
      progress_callback(sprintf("Executing Stage %d: %s", stage_idx, paste(stage$tasks, collapse = ", ")), progress)
    }
    
    # Execute tasks in current stage
    # V3.16.3: Pass parsed_tasks to execute_stage_tasks
    stage_task_results <- execute_stage_tasks(
      stage,
      execution_context,
      progress_callback,
      parsed_tasks
    )
    
    # Update execution context with results
    execution_context <- update_execution_context(execution_context, stage_task_results)
    
    stage_results[[stage_idx]] <- stage_task_results
    task_results <- c(task_results, stage_task_results$individual_results)
    
    log_message(sprintf("Stage %d completed successfully", stage_idx))
  }
  
  return(list(
    stage_results = stage_results,
    task_results = task_results,
    execution_context = execution_context
  ))
}

#' Execute tasks within a single stage
#' 
#' @param stage Stage configuration
#' @param execution_context Execution context
#' @param progress_callback Progress callback function
#' @param parsed_tasks Parsed task configurations
#' @return Stage task results
execute_stage_tasks <- function(stage,
                               execution_context,
                               progress_callback = NULL,
                               parsed_tasks = NULL) {
  
  stage_task_results <- list()
  individual_results <- list()
  
  for (task_idx in seq_along(stage$tasks)) {
    task_name <- stage$tasks[[task_idx]]
    
    log_message(sprintf("Executing task: %s", task_name))
    
    # V3.16.3: Get task configuration from parsed_tasks (not raw config_data)
    # Parsed tasks have standardized structure with module field
    task_config <- parsed_tasks[[task_name]]
    
    if (is.null(task_config)) {
      log_message(sprintf("ERROR: Task configuration not found for: %s", task_name), level = "error")
      stop(sprintf("Task configuration not found for: %s", task_name))
    }
    
    # Execute individual task
    task_result <- execute_single_task(
      task_name,
      task_config,
      execution_context,
      progress_callback
    )
    
    stage_task_results[[task_name]] <- task_result
    individual_results[[task_name]] <- task_result

    # Make completed task results available to later tasks in the same stage.
    # Reporting depends on these manifests before the stage-level merge occurs.
    if (is.null(execution_context$task_results)) {
      execution_context$task_results <- list()
    }
    execution_context$task_results[[task_name]] <- task_result
    
    if (task_result$success) {
      log_message(sprintf("Task %s completed successfully", task_name))
      
      # V3.16.5: CRITICAL FIX - Update execution_context immediately after each successful task
      # This ensures subsequent tasks in the same stage can access previous task results
      if (task_result$module == "preprocessing") {
        execution_context$data_cache$preprocessing_results <- task_result$result_data
        log_message("[DEBUG] Updated data_cache with preprocessing results")
      } else if (task_result$module == "filtering") {
        execution_context$data_cache$filtering_results <- task_result$result_data
        log_message("[DEBUG] Updated data_cache with filtering results")
      } else if (task_result$module == "normalization") {
        execution_context$data_cache$normalization_results <- task_result$result_data
        log_message("[DEBUG] Updated data_cache with normalization results")
      }
    } else {
      log_message(sprintf("Task %s failed: %s", task_name, task_result$error_message), level = "error")
    }
  }
  
  return(list(
    stage_id = stage$stage_id,
    stage_results = stage_task_results,
    individual_results = individual_results,
    stage_success = all(sapply(stage_task_results, function(r) r$success))
  ))
}

#' Execute a single analysis task
#' 
#' @param task_name Name of the task
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param progress_callback Progress callback function
#' @return Task execution result
execute_single_task <- function(task_name,
                               task_config,
                               execution_context,
                               progress_callback = NULL) {
  
  start_time <- Sys.time()
  
  # Create task-specific result structure
  task_result <- list(
    task_name = task_name,
    task_id = task_config$task_id,
    module = task_config$module,
    success = FALSE,
    result_data = NULL,
    error_message = NULL,
    execution_time = NULL,
    warnings = character(0)
  )
  
  module_error <- NULL
  module_result <- tryCatch({
    
    # Route task to appropriate module function
    module_result <- route_task_to_module(
      task_config,
      execution_context,
      progress_callback
    )

    assert_module_result_success(module_result, task_config$module)
    module_result
  }, error = function(e) {
    module_error <<- e
    
    # V3.16.2: Enhanced error logging with full traceback
    log_message(sprintf("Task %s error: %s", task_name, e$message), level = "error")
    log_message("[DEBUG] Full error traceback:", level = "error")
    
    # Capture and log the call stack
    call_stack <- sys.calls()
    for (i in seq_along(call_stack)) {
      call_str <- deparse(call_stack[[i]])[1]  # Get first line of each call
      if (nchar(call_str) > 100) {
        call_str <- paste0(substr(call_str, 1, 97), "...")
      }
      log_message(sprintf("  [%d] %s", i, call_str), level = "error")
    }
    NULL
  })

  if (is.null(module_error)) {
    task_result$success <- TRUE
    task_result$result_data <- module_result
  } else {
    task_result$success <- FALSE
    task_result$error_message <- conditionMessage(module_error)
  }
  
  task_result$execution_time <- difftime(Sys.time(), start_time, units = "secs")
  
  return(task_result)
}

#' Route task to appropriate module function
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param progress_callback Progress callback function
#' @return Module execution result
route_task_to_module <- function(task_config,
                                execution_context,
                                progress_callback = NULL) {
  
  module <- task_config$module
  parameters <- task_config$parameters
  data_source <- task_config$data_source %||% "filtered"
  
  # V3.16.15: Debug logging for data_scoping parameter passing
  log_message(sprintf("[DEBUG] route_task_to_module - module: %s", module))
  log_message(sprintf("[DEBUG] route_task_to_module - parameters keys: %s", 
                     paste(names(parameters), collapse = ", ")))
  if (!is.null(parameters$data_scoping)) {
    log_message(sprintf("[DEBUG] route_task_to_module - data_scoping keys: %s", 
                       paste(names(parameters$data_scoping), collapse = ", ")))
    if (!is.null(parameters$data_scoping$target_var_types)) {
      log_message(sprintf("[DEBUG] route_task_to_module - target_var_types: %s", 
                         paste(parameters$data_scoping$target_var_types, collapse = ", ")))
    }
  } else {
    log_message("[DEBUG] route_task_to_module - data_scoping is NULL in parameters")
  }
  
  # CRITICAL FIX: Include data_scoping and other task-level fields in parameters
  # NOTE: data_scoping is already in parameters from YAML, but plot_types etc are at task level
  # Only override if task-level fields exist AND are not empty (for backward compatibility)
  # V3.16.15: Check length > 0 to avoid overriding with empty list()
  if (!is.null(task_config$data_scoping) && length(task_config$data_scoping) > 0) {
    log_message("[DEBUG] route_task_to_module - Overriding data_scoping from task_config level")
    parameters$data_scoping <- task_config$data_scoping
  }
  if (!is.null(task_config$plot_types)) {
    parameters$plot_types <- task_config$plot_types
  }
  if (!is.null(task_config$output_settings)) {
    parameters$output_settings <- task_config$output_settings
  }
  if (!is.null(task_config$factors)) {
    parameters$factors <- task_config$factors
  }
  
  # Merge global parameters with task-specific parameters
  global_params <- execution_context$metadata$global_parameters %||% list()
  merged_params <- merge_parameters(global_params, parameters)
  
  # Check for visualization-only mode
  config_data <- execution_context$metadata$config_data
  # V3.16: Simplified configuration - workflow section may not exist
  is_visualization_mode <- FALSE  # Always FALSE in V3.16 (full_analysis only)
  use_existing_data <- merged_params$use_existing_data %||% FALSE
  
  # In visualization-only mode, load existing data for visualization tasks
  if (is_visualization_mode || use_existing_data) {
    
    # For data processing tasks in visualization mode, either skip or load existing data
    if (module %in% c("preprocessing", "filtering", "normalization")) {
      if (task_config$enabled %||% TRUE) {
        log_message(sprintf("Skipping %s task in visualization-only mode", module))
      }
      # Return empty result to indicate task was skipped
      return(list(
        skipped = TRUE,
        module = module,
        reason = "visualization_only_mode",
        message = sprintf("%s task skipped in visualization-only mode", module)
      ))
    }
    
    # For visualization tasks, ensure we have loaded data
    if (module %in% c("M01_distribution", "M02_comparative", "M03_hotspot")) {
      # Check if we already have normalization data in context
      if (is.null(execution_context$data_cache$normalization_results)) {
        log_message("Loading existing data for visualization tasks")
        existing_data <- load_existing_data_for_visualization(task_config, execution_context, merged_params)
        
        # Update execution context with loaded data
        execution_context$data_cache$normalization_results <- existing_data
      }
    }
  }
  
  # ARCHITECTURE V6.2: Module-specific data loading - each execute function handles its own data needs
  # The execution context provides raw normalization data, modules apply scoping as needed
  log_message("=== MODULE-SPECIFIC DATA LOADING ARCHITECTURE ===")
  log_message("Raw normalization data available - modules will handle scoping and type conversions")

  # Route to appropriate module
  if (module == "preprocessing") {
    return(execute_preprocessing_task(task_config, execution_context, merged_params))
    
  } else if (module == "filtering") {
    return(execute_filtering_task(task_config, execution_context, merged_params))
    
  } else if (module == "normalization") {
    return(execute_normalization_task(task_config, execution_context, merged_params))
    
  } else if (module == "M01_distribution") {
    return(execute_M01_task(task_config, execution_context, merged_params))
    
  } else if (module == "M02_comparative") {
    return(execute_M02_task(task_config, execution_context, merged_params))
    
  } else if (module == "M03_hotspot") {
    return(execute_M03_task(task_config, execution_context, merged_params))
    
  } else if (module == "M04_igs_hotspot") {
    return(execute_M04_task(task_config, execution_context, merged_params))
    
  } else if (module == "M05_hotspot_ideogram") {
    # M05 chromosome hotspot visualization
    # Integrates M03 CDS hotspots and M04 poiGS hotspots
    return(execute_M05_task(task_config, execution_context, merged_params))
    
  } else if (module == "M05_integration") {
    # Deprecated in the internal V4 workflow: M05 module has been migrated to M01.
    # All M05 functionality is now available through M01 sub-modules
    stop("M05_integration module is deprecated. Please use M01_standard_distribution with enabled_submodules instead. See configuration file for migration guide.")
    
  } else if (module == "generate_report" || module == "reporting") {
    return(execute_generate_report_task(task_config, execution_context, merged_params))
    
  } else {
    stop(paste("Unknown module:", module))
  }
}

#' Execute preprocessing task
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param parameters Merged parameters
#' @return Preprocessing result
execute_preprocessing_task <- function(task_config, execution_context, parameters) {

  # V3.16.2: Debug logging
  log_message("[DEBUG] execute_preprocessing_task started")
  log_message(sprintf("[DEBUG] execution_context$metadata exists: %s", !is.null(execution_context$metadata)))
  log_message(sprintf("[DEBUG] execution_context$session_paths exists: %s", !is.null(execution_context$session_paths)))

  # Step 1: Use helper function to load and map data
  config_data <- execution_context$metadata$config_data
  log_message("[DEBUG] config_data extracted from execution_context")
  log_message(sprintf("[DEBUG] config_data$input_files exists: %s", !is.null(config_data$input_files)))
  
  variant_df <- load_main_variant_data(config_data, execution_context$session_paths)
  log_message("[DEBUG] load_main_variant_data completed")

  # Step 2 (key fix): Inject loaded data frame into parameter list (parameters) passed to module
  parameters$data <- variant_df
  
  # Step 2.5 (shared dependency injection): Also inject shared species_genome_data into parameters
  parameters$species_genome_data <- execution_context$data_cache$species_genome_data

  # Step 3: Prepare runtime configuration
  runtime_config <- apply_task_parameters_to_config(parameters, config_data)
  runtime_config$session_info <- list(session_id = execution_context$metadata$session_id)

  # Step 4 (key fix): Use R package's new standard function signature to call processor
  result <- preprocess_variants(
    config = runtime_config,
    task_params = parameters,
    save_intermediate = TRUE
  )

  return(result)
}

#' Execute filtering task
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context  
#' @param parameters Merged parameters
#' @return Filtering result
execute_filtering_task <- function(task_config, execution_context, parameters) {
  
  # GOLDEN RULE: IRA-only data is the absolute processing mainline from P02 onwards
  # All core processing (filtering, normalization) MUST use IRA-only data exclusively
  preprocessing_results <- execution_context$data_cache$preprocessing_results
  if (is.null(preprocessing_results)) {
    stop("Filtering task requires preprocessing results, but none found in data cache")
  }
  
  preprocessing_ira_data <- preprocessing_results$data_ira_only  # IRA-only data for ALL processing
  
  if (is.null(preprocessing_ira_data)) {
    stop("Filtering task requires IRA-only preprocessing data, but none found in preprocessing results")
  }
  
  log_message(sprintf("GOLDEN RULE ENFORCEMENT: Using %d rows from IRA-only data as the sole processing mainline", nrow(preprocessing_ira_data)))
  log_message("Full Genome data is excluded from core processing per architectural golden rules")
  
  # Apply parameters to runtime config and include session info
  config_data <- execution_context$metadata$config_data
  # V25 FINAL FIX: Use the config_data loaded from YAML as the base config,
  # instead of the obsolete SHINY_CONFIG object.
  runtime_config <- apply_task_parameters_to_config(parameters, config_data)
  runtime_config$session_info <- list(session_id = execution_context$metadata$session_id)
  
  # Get shared species_genome_data dependency  
  species_genome_data <- execution_context$data_cache$species_genome_data
  
  # Execute filtering with IRA-only data as the sole processing mainline
  result <- filter_variants(
    input_data = preprocessing_ira_data, # Pass IRA-only data frame exclusively
    config = runtime_config,
    species_genome_data = species_genome_data, # Inject shared dependency
    save_results = TRUE
  )
  
  return(result)
}

#' Execute normalization task
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param parameters Merged parameters
#' @return Normalization result
execute_normalization_task <- function(task_config, execution_context, parameters) {
  
  # Get filtering results
  filtering_results <- execution_context$data_cache$filtering_results
  if (is.null(filtering_results)) {
    stop("Normalization task requires filtering results, but none found in data cache")
  }
  
  # The data source is now linear; always use the output from the filtering stage
  variant_data <- filtering_results$data
  
  log_message(sprintf("Normalizing frequencies for %d variants from the filtering stage", nrow(variant_data)))
  
  # Apply parameters to runtime config and include session info  
  config_data <- execution_context$metadata$config_data
  # V25 FINAL FIX: Use the config_data loaded from YAML as the base config,
  # instead of the obsolete SHINY_CONFIG object.
  runtime_config <- apply_task_parameters_to_config(parameters, config_data)
  runtime_config$session_info <- list(session_id = execution_context$metadata$session_id)
  
  # Get configuration data
  session_paths <- execution_context$session_paths
  
  # Load auxiliary data using configuration-driven data loading helpers
  log_message("Loading auxiliary data for normalization using data loading helpers")
  
  # Load species genome data for genome_region assignment
  species_genome_data <- execution_context$data_cache$species_genome_data
  
  # Load CDS lengths data (pass species_genome_data for consistency)
  gene_lengths_df <- load_cds_lengths_data(config_data, session_paths, species_genome_data)
  if (is.null(gene_lengths_df)) {
    stop("Failed to load gene lengths data for normalization")
  }
  
  # Load region information, passing species_genome_data for genome_region assignment
  region_info_df <- load_region_info_data(config_data, session_paths, species_genome_data)
  if (is.null(region_info_df)) {
    log_message("Region info data not available, some region types may not be processed", level = "warning")
  }
  
  log_message(sprintf("Auxiliary data loaded: %d gene entries, %d region info entries", 
                     nrow(gene_lengths_df), ifelse(is.null(region_info_df), 0, nrow(region_info_df))))
  
  # Execute normalization with the simplified, linear data flow
  result <- normalize_frequencies(
    variant_data = variant_data,
    gene_lengths_df = gene_lengths_df,
    region_info_df = region_info_df,
    config = runtime_config
  )
  
  # Saving is now handled within normalize_frequencies, so we just log it here.
  log_message("Normalization results have been saved by the module.")
  
  # Region coverage validation is currently disabled
  # This validation is a diagnostic tool and does not affect core analysis
  # Disabled due to persistent R package loading issues that require further investigation
  # P02/P03 discrepancies can be analyzed manually by comparing CSV files
  log_message("=== Region coverage validation DISABLED ===")
  log_message("Coverage validation temporarily disabled. Core analysis pipeline unaffected.")
  log_message("P02/P03 discrepancies can be analyzed manually if needed.")
  
  return(result)
}

#' Load existing data for visualization-only mode
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param parameters Merged parameters
#' @return Existing normalized data
load_existing_data_for_visualization <- function(task_config, execution_context, parameters) {
  
  log_message("Loading existing data for visualization-only mode")
  
  # Get session information
  session_id <- execution_context$metadata$session_id
  data_source <- task_config$data_source %||% "filtered"
  
  # Try to load existing normalized data using helper function
  existing_data <- load_existing_normalized_data(session_id, data_source)
  
  if (is.null(existing_data)) {
    # Try alternative data sources or locations
    log_message("Primary data source not found, trying alternative locations", level = "warning")
    
    # Check if there's cached data in execution context
    if (!is.null(execution_context$data_cache$normalization_results)) {
      log_message("Using normalization results from execution context cache")
      existing_data <- execution_context$data_cache$normalization_results$frequencies
    } else {
      stop(sprintf("No existing normalized data found for session '%s' and data source '%s'. Please run data processing first or check session ID.", session_id, data_source))
    }
  }
  
  # Validate data structure
  required_cols <- c("species", "frequency_per_kb", "variant_count", "region_type", "gene", "region_length")
  missing_cols <- setdiff(required_cols, names(existing_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Existing data is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  log_message(sprintf("Successfully loaded existing data: %d rows, %d columns", nrow(existing_data), ncol(existing_data)))
  
  # Return in same format as normalization results
  return(list(
    frequencies = existing_data,
    frequencies_summary = existing_data %>%
      dplyr::group_by(species, gene, region_type) %>%
      dplyr::summarise(
        variant_count = sum(variant_count, na.rm = TRUE),
        region_length = dplyr::first(region_length),
        frequency_per_kb = sum(frequency_per_kb, na.rm = TRUE), 
        .groups = "drop"
      ) %>%
      dplyr::mutate(var_type = "all_types") %>%
      dplyr::select(species, var_type, region_type, gene, variant_count, region_length, frequency_per_kb),
    metadata = list(
      data_source = data_source,
      session_id = session_id,
      loaded_from_existing = TRUE,
      load_timestamp = Sys.time()
    )
  ))
}

#' Execute M01 distribution analysis task
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param parameters Merged parameters
#' @return M01 result
execute_M01_task <- function(task_config, execution_context, parameters) {
  
  # ARCHITECTURE V6.2: Module-specific data loading for M01
  # M01 is the only module requiring dual data sources and will handle its own scoping
  normalization_data <- execution_context$data_cache$normalization_results
  if (is.null(normalization_data)) {
    stop("M01 task requires normalization results, but none found in data cache")
  }
  
  log_message("=== M01 MODULE-SPECIFIC DATA LOADING ===")
  log_message("Passing raw normalization cache to M01 module for internal data handling")
  log_message(sprintf("Raw frequencies data: %d rows", nrow(normalization_data$frequencies)))
  log_message(sprintf("Raw frequencies_summary data: %d rows", nrow(normalization_data$frequencies_summary)))
  
  # V4.0 ENHANCEMENT: Pass additional data for M01 sub-modules
  # M01-2 (chromosome_distribution) requires filtered_data
  # M01-3 (cds_snp_annotation, igs_variant_features) requires both filtered and preprocessed data
  filtering_results <- execution_context$data_cache$filtering_results
  preprocessing_results <- execution_context$data_cache$preprocessing_results
  
  if (!is.null(filtering_results)) {
    log_message(sprintf("Filtered data available for M01 sub-modules: %d rows", 
                       nrow(filtering_results$data)))
  }
  if (!is.null(preprocessing_results)) {
    log_message(sprintf("Preprocessed data available for M01 sub-modules: IRA=%d rows, Full=%d rows", 
                       nrow(preprocessing_results$data_ira_only),
                       nrow(preprocessing_results$data_full_genome)))
  }
  
  # Apply task-specific parameters and merge YAML config data
  config_data <- execution_context$metadata$config_data
  # V25 FINAL FIX: Use the config_data loaded from YAML as the base config,
  # instead of the obsolete SHINY_CONFIG object.
  runtime_config <- apply_task_parameters_to_config(parameters, config_data)
  
  # ARCHITECTURE PRINCIPLE #10: Inject session_id into config object
  runtime_config$session_info <- list(session_id = execution_context$metadata$session_id)

  # Generate M01 visualizations using Phase 2 orchestrator with raw normalization data
  # The M01 module will handle its own dual data loading and scoping internally
  log_message("M01 module will handle dual data loading internally per architect's elegant solution")
  result <- run_m01_distribution(
    pie_data = normalization_data$frequencies,      # Raw site-level data
    freq_data = normalization_data$frequencies_summary, # Raw gene-level data
    normalized_data = normalization_data$frequencies,   # V4.0: For poiGS correlation
    filtered_data = if (!is.null(filtering_results)) filtering_results$data else NULL,  # V4.0: For chromosome & sequence analyses
    preprocessing_results = preprocessing_results,  # V4.0: For sample-level annotations
    config = runtime_config,
    output_dir = execution_context$session_paths$plots,
    task_name = task_config$task_id,
    task_params = parameters
  )
  
  return(result)
}

#' Execute M02 comparative analysis task
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param parameters Merged parameters
#' @return M02 result
execute_M02_task <- function(task_config, execution_context, parameters) {
  
  # ARCHITECTURE V6.2: Single data mode for M02 - module handles its own data cleanup
  normalization_data <- execution_context$data_cache$normalization_results
  if (is.null(normalization_data)) {
    stop("M02 task requires normalization results, but none found in data cache")
  }
  
  # M02 always starts from detailed (per-variant-type) data
  # Architecture principle: M02 module handles its own aggregation internally
  # This ensures data_scoping can filter by var_type before aggregation
  normalized_data_to_use <- normalization_data$frequencies
  log_message("=== M02 DATA SOURCE: frequencies (detailed data) ===")
  log_message(sprintf("Detailed data: %d rows", nrow(normalized_data_to_use)))
  log_message("M02 module will handle aggregation internally based on analysis_mode")
  
  # Apply task-specific parameters and merge YAML config data
  config_data <- execution_context$metadata$config_data
  # V25 FINAL FIX: Use the config_data loaded from YAML as the base config,
  # instead of the obsolete SHINY_CONFIG object.
  runtime_config <- apply_task_parameters_to_config(parameters, config_data)
  
  # ARCHITECTURE PRINCIPLE #10: Inject session_id into config object
  runtime_config$session_info <- list(session_id = execution_context$metadata$session_id)
  
  # Generate M02 visualizations using Phase 2 orchestrator with appropriate data source
  # M02 module will handle data scoping and cleanup internally per "who uses, who cleans" principle
  result <- run_m02_comparative(
    normalized_data_to_use,  # Data source selected based on analysis_mode
    config = runtime_config,
    output_dir = execution_context$session_paths$plots,
    task_name = task_config$task_id,
    task_params = parameters
  )
  
  return(result)
}

#' Execute M03 hotspot analysis task
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param parameters Merged parameters
#' @return M03 result
execute_M03_task <- function(task_config, execution_context, parameters) {
  
  # ARCHITECTURE V6.2: Single data mode for M03 - module handles its own data cleanup
  normalization_data <- execution_context$data_cache$normalization_results
  if (is.null(normalization_data)) {
    stop("M03 task requires normalization results, but none found in data cache")
  }
  
  log_message("=== M03 SINGLE DATA MODE ===")
  log_message("Passing raw normalization data to M03 module for internal cleanup")
  log_message(sprintf("Raw normalized frequencies: %d rows", nrow(normalization_data$frequencies)))
  
  # Apply task-specific parameters and merge YAML config data
  config_data <- execution_context$metadata$config_data
  # V25 FINAL FIX: Use the config_data loaded from YAML as the base config,
  # instead of the obsolete SHINY_CONFIG object.
  runtime_config <- apply_task_parameters_to_config(parameters, config_data)
  
  # ARCHITECTURE PRINCIPLE #10: Inject session_id into config object
  runtime_config$session_info <- list(session_id = execution_context$metadata$session_id)
  
  # Generate M03 visualizations using Phase 2 orchestrator with raw data
  # M03 module will handle data scoping and cleanup internally per "who uses, who cleans" principle
  result <- run_m03_hotspot(
    normalization_data$frequencies,  # Raw frequencies data
    config = runtime_config,
    output_dir = execution_context$session_paths$plots,
    task_name = task_config$task_id,
    task_params = parameters
  )
  
  return(result)
}

#' Execute M04 IGS hotspot analysis task
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param parameters Merged parameters
#' @return M04 result
execute_M04_task <- function(task_config, execution_context, parameters) {
  
  # M04 uses normalized frequency data for IGS-specific hotspot analysis
  # This ensures analysis is based on standardized frequency values across regions
  normalization_results <- execution_context$data_cache$normalization_results
  if (is.null(normalization_results)) {
    stop("M04 task requires normalization results, but none found in data cache")
  }
  
  log_message("=== M04 IGS HOTSPOT MODE ===")
  log_message("Using normalized frequency data for IGS-specific hotspot analysis")
  log_message(sprintf("Normalized frequency data: %d rows", nrow(normalization_results$frequencies)))
  
  # Apply task-specific parameters and merge YAML config data
  config_data <- execution_context$metadata$config_data
  runtime_config <- apply_task_parameters_to_config(parameters, config_data)
  
  # Inject session_id into config object
  runtime_config$session_info <- list(session_id = execution_context$metadata$session_id)
  
  # Generate M04 IGS hotspot visualization
  # M04 module handles IGS-specific data scoping internally
  result <- run_m04_igs_analysis(
    normalized_data = normalization_results$frequencies,  # Normalized frequency data
    config = runtime_config,
    output_dir = execution_context$session_paths$plots,
    task_name = task_config$task_id,
    task_params = parameters
  )
  
  return(result)
}

#' Execute M05 Hotspot Ideogram task
#' 
#' @description
#' Chromosome-level hotspot visualization (M05).
#' Integrates CDS hotspots from M03 and poiGS hotspots from M04 into
#' unified chromosome ideogram plots.
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param parameters Merged parameters
#' @return M05 hotspot ideogram result
#' 
#' @keywords internal
execute_M05_task <- function(task_config, execution_context, parameters) {
  
  log_message("=== M05 HOTSPOT IDEOGRAM MODE ===")
  log_message("Integrating M03 and M04 hotspot results for chromosome visualization")
  
  # M05 depends on M03 and M04 outputs, not on normalized data
  # The module will auto-detect and load the required CSV files
  
  # Apply task-specific parameters and merge YAML config data
  config_data <- execution_context$metadata$config_data
  runtime_config <- apply_task_parameters_to_config(parameters, config_data)
  
  # Inject session_id into config object
  runtime_config$session_info <- list(session_id = execution_context$metadata$session_id)
  
  # Generate M05 hotspot ideogram visualization
  # M05 module handles dependency validation and file loading internally
  result <- run_M05_hotspot_ideogram(
    normalized_data = NULL,  # Not used, M05 loads M03/M04 outputs directly
    config = runtime_config,
    output_dir = execution_context$session_paths$plots,
    task_name = task_config$task_id,
    task_params = parameters
  )
  
  return(result)
}

#' Execute M05 integration analysis task (DEPRECATED)
#' 
#' @description
#' Deprecated in the internal V4 workflow; retained for backward compatibility only.
#' All M05 functionality has been migrated to M01 EDA Hub.
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param parameters Merged parameters
#' @return M05 integration result
#' 
#' @keywords internal
execute_M05_integration_task <- function(task_config, execution_context, parameters) {
  
  # Deprecated: this function should not be called by the current backend.
  # The dispatcher now returns an error message before reaching this point
  
  stop(paste(
    "M05_integration module is deprecated in the current backend.",
    "All M05 functionality has been migrated to M01 EDA Hub.",
    "Please update your configuration file:",
    "  1. Set M05_integration_analysis.enabled = false",
    "  2. Enable M01 sub-modules: enabled_submodules = ['foundational', 'structural_relational', 'sequence_variants']",
    "See full_refactor_test.yml for the complete migration guide.",
    sep = "\n"
  ))
}

#' Execute report generation task
#' 
#' @param task_config Task configuration
#' @param execution_context Execution context
#' @param parameters Merged parameters
#' @return Report generation result
execute_generate_report_task <- function(task_config, execution_context, parameters) {
  
  log_message("=== GENERATING FINAL REPORT ===")
  
  # Get session information
  session_id <- execution_context$metadata$session_id
  config_data <- execution_context$metadata$config_data
  session_paths <- execution_context$session_paths
  
  log_message(sprintf("Generating final report for session: %s", session_id))
  
  # Convert plots to additional formats before reporting
  log_message("=== PRE-REPORTING: FORMAT CONVERSION ===")
  plots_dir <- session_paths$plots
  
  if (dir.exists(plots_dir)) {
    conversion_result <- convert_plots_for_reporting(
      plots_dir = plots_dir,
      config = config_data,
      session_id = session_id
    )
    
    if (conversion_result$success) {
      log_message(sprintf("Plot conversion completed: %d files converted", 
                         conversion_result$total_converted))
    } else {
      conversion_error <- conversion_result$error %||%
        "Requested plot conversion failed."
      log_message(
        sprintf("[ERROR] Plot conversion failed: %s", conversion_error),
        level = "error"
      )
      return(list(
        success = FALSE,
        error_message = conversion_error,
        conversion = conversion_result
      ))
    }
  } else {
    log_message("Plots directory not found, skipping format conversion", level = "warning")
  }
  
  # Collect manifests from all task results
  manifests <- collect_task_manifests(execution_context)
  
  # Generate the final report
  result <- generate_final_report(
    session_id = session_id,
    config_data = config_data,
    execution_context = execution_context,
    manifests = manifests,
    output_path = parameters$output_path %||% NULL
  )
  
  if (result$success) {
    log_message(sprintf("[SUCCESS] Final report generated successfully: %s", result$output_path))
  } else {
    log_message(sprintf("[ERROR] Failed to generate final report: %s", result$error_message), level = "error")
  }
  
  return(result)
}

#' Update execution context with task results
#' 
#' @param execution_context Current execution context
#' @param stage_results Stage execution results
#' @return Updated execution context
update_execution_context <- function(execution_context, stage_results) {
  
  # V6.5 FIX: Accumulate task results for manifest collection
  if (is.null(execution_context$task_results)) {
    execution_context$task_results <- list()
  }
  
  # Add individual task results to execution context
  for (task_name in names(stage_results$individual_results)) {
    execution_context$task_results[[task_name]] <- stage_results$individual_results[[task_name]]
  }
  
  # Update data cache with successful task results
  for (task_name in names(stage_results$individual_results)) {
    task_result <- stage_results$individual_results[[task_name]]
    
    if (task_result$success) {
      # Cache results based on module type, using the new simplified structure
      if (task_result$module == "preprocessing") {
        # The result_data itself is the list containing the 'data' dataframe
        execution_context$data_cache$preprocessing_results <- task_result$result_data
      } else if (task_result$module == "filtering") {
        # The result_data is the list containing the 'data' dataframe
        execution_context$data_cache$filtering_results <- task_result$result_data
      } else if (task_result$module == "normalization") {
        # The result_data is the full list from normalization (frequencies, summary, etc.)
        execution_context$data_cache$normalization_results <- task_result$result_data
      }
    }
  }
  
  return(execution_context)
}

#' Merge global and task-specific parameters
#' 
#' @param global_params Global parameters
#' @param task_params Task-specific parameters
#' @return Merged parameters
merge_parameters <- function(global_params, task_params) {
  
  # Deep merge of parameter lists, with task params taking precedence
  merged <- global_params
  
  for (param_name in names(task_params)) {
    merged[[param_name]] <- task_params[[param_name]]
  }
  
  return(merged)
}

# V29 Simplify the function to only merge task-specific
# parameters into the base configuration loaded from YAML.

apply_task_parameters_to_config <- function(parameters, base_config) {

  runtime_config <- base_config

  # Merge visualization_settings from base_config for V4 color system support
  if (!is.null(base_config$visualization_settings)) {
    runtime_config$visualization_settings <- base_config$visualization_settings
    log_message("Successfully merged visualization_settings from YAML config")
  } else {
    log_message("WARNING: No visualization_settings found in YAML config", level = "warning")
  }

  # Apply task-specific filtering parameters if they exist
  if (!is.null(parameters)) {
      # Ensure filtering object exists
      if (is.null(runtime_config$filtering)) {
        runtime_config$filtering <- list()
      }
      
      # V3.6 DYNAMIC CONFIGURATION: Support new thresholds structure
      if (!is.null(parameters$thresholds) && is.list(parameters$thresholds)) {
        # New dynamic thresholds structure - this is the modern approach
        runtime_config$filtering$thresholds <- parameters$thresholds
        log_message("Applied dynamic threshold configuration:")
        for (vtype in names(parameters$thresholds)) {
          log_message(sprintf("  %s: >= %d samples", vtype, parameters$thresholds[[vtype]]))
        }
      } else {
        # Legacy individual threshold parameters support (backward compatibility)
        has_filtering_params <- any(c("snp_threshold", "indel_threshold", "complex_threshold", "mnp_threshold") %in% names(parameters))
        
        if ("snp_threshold" %in% names(parameters)) {
          runtime_config$filtering$snp_threshold <- parameters$snp_threshold
        }
        if ("indel_threshold" %in% names(parameters)) {
          runtime_config$filtering$indel_threshold <- parameters$indel_threshold
        }
        if ("complex_threshold" %in% names(parameters)) {
          runtime_config$filtering$complex_threshold <- parameters$complex_threshold
        }
        if ("mnp_threshold" %in% names(parameters)) {
          runtime_config$filtering$mnp_threshold <- parameters$mnp_threshold
        }
        
        # Ensure all required thresholds have default values (using most lenient default value 1)
        if (is.null(runtime_config$filtering$complex_threshold)) {
          runtime_config$filtering$complex_threshold <- 1
        }
        if (is.null(runtime_config$filtering$mnp_threshold)) {
          runtime_config$filtering$mnp_threshold <- 1
        }
        
        # Only log filtering configuration for tasks that actually have filtering parameters
        if (has_filtering_params) {
          log_message(sprintf("Applied legacy filtering configuration: SNP=%s, INDEL=%s, complex=%s, mnp=%s", 
                             runtime_config$filtering$snp_threshold %||% "default",
                             runtime_config$filtering$indel_threshold %||% "default", 
                             runtime_config$filtering$complex_threshold,
                             runtime_config$filtering$mnp_threshold))
        }
      }
  }

  # Ensure normalization config exists
  if (is.null(runtime_config$normalization)) {
    runtime_config$normalization <- list()
  }
  if (is.null(runtime_config$normalization$frequency_unit)) {
    runtime_config$normalization$frequency_unit <- 1000  # Default: per kilobase
  }

  return(runtime_config)
}

#' Aggregate task results across all executed tasks
#' 
#' @param execution_results Task execution results
#' @param parsed_config Parsed configuration
#' @return Aggregated results
aggregate_task_results <- function(execution_results, parsed_config) {
  
  log_message("Aggregating task results")
  
  task_results <- execution_results$task_results
  
  # Organize results by module
  results_by_module <- list(
    preprocessing = list(),
    filtering = list(),
    normalization = list(),
    M01_distribution = list(),
    M02_comparative = list(),
    M03_hotspot = list()
  )
  
  # Group results
  for (task_name in names(task_results)) {
    task_result <- task_results[[task_name]]
    if (task_result$success) {
      module <- task_result$module
      results_by_module[[module]][[task_name]] <- task_result$result_data
    }
  }
  
  # Create cross-task analysis summary
  cross_task_summary <- create_cross_task_summary(task_results)
  
  return(list(
    by_module = results_by_module,
    cross_task_summary = cross_task_summary,
    successful_tasks = sum(sapply(task_results, function(r) r$success)),
    failed_tasks = sum(sapply(task_results, function(r) !r$success))
  ))
}

#' Create cross-task analysis summary
#' 
#' @param task_results Task execution results
#' @return Cross-task summary
create_cross_task_summary <- function(task_results) {
  
  summary <- list(
    total_tasks = length(task_results),
    successful_tasks = sum(sapply(task_results, function(r) r$success)),
    failed_tasks = sum(sapply(task_results, function(r) !r$success)),
    modules_executed = unique(sapply(task_results, function(r) r$module)),
    total_execution_time = sum(sapply(task_results, function(r) as.numeric(r$execution_time)), na.rm = TRUE)
  )
  
  return(summary)
}

#' Generate execution summary
#' 
#' @param execution_results Execution results
#' @param parsed_config Parsed configuration
#' @return Execution summary
generate_execution_summary <- function(execution_results, parsed_config) {
  
  task_results <- execution_results$task_results
  
  summary <- list(
    total_tasks_executed = length(task_results),
    successful_tasks = sum(sapply(task_results, function(r) r$success)),
    failed_tasks = sum(sapply(task_results, function(r) !r$success)),
    total_stages = length(execution_results$stage_results),
    total_execution_time = sum(sapply(task_results, function(r) as.numeric(r$execution_time)), na.rm = TRUE),
    modules_used = unique(sapply(task_results, function(r) r$module))
  )
  
  return(summary)
}

#' Collect manifests from all task results
#' 
#' @param execution_context Execution context containing task results
#' @return List of manifests from all successful tasks
collect_task_manifests <- function(execution_context) {
  
  log_message("Collecting task manifests from execution context")
  
  manifests <- list()
  
  # Check if task_results exists in execution_context
  if (is.null(execution_context$task_results)) {
    log_message("No task_results found in execution_context", level = "warning")
    return(manifests)
  }
  
  # Iterate through all task results
  for (task_name in names(execution_context$task_results)) {
    task_result <- execution_context$task_results[[task_name]]
    
    # Check if task was successful and has manifest data
    if (task_result$success && !is.null(task_result$result_data$manifest)) {
      manifests[[task_name]] <- task_result$result_data$manifest
      log_message(sprintf("[SUCCESS] Collected manifest from task: %s", task_name))
    }
  }
  
  log_message(sprintf("Collected %d manifests from %d tasks", length(manifests), length(execution_context$task_results)))
  
  return(manifests)
}

# log_message("Task execution dispatcher loaded successfully")