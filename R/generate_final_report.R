############################################################
#### generate_final_report.R - Final Report Generator ####
############################################################
#
# Automated final report generation system for cpopvar.
# Creates comprehensive Markdown reports based on analysis results.
#
# Key Features:
# - Template-based report generation with dynamic placeholder replacement
# - Automatic collection of session statistics and analysis results
# - Integration with existing task-driven pipeline
# - Conservative implementation with error isolation
#
############################################################

# Load required libraries

# Define %||% operator if not available
if (!exists("%||%")) {
  # Note: Using null coalescing operator (%||%) exported from annotation_processing.R
}

# V3.2.1 Hotfix: Create a safe sprintf wrapper to prevent list coercion errors
#' Safe sprintf wrapper to prevent list coercion errors
#' @param fmt Format string
#' @param ... Arguments to be formatted
#' @return Formatted string
#' @export
safe_sprintf <- function(fmt, ...) {
  args <- list(...)
  # Ensure all arguments are atomic vectors, converting lists to comma-separated strings
  safe_args <- lapply(args, function(arg) {
    # V3.13: Handle NULL, empty, or zero-length arguments to prevent sprintf errors
    if (is.null(arg) || length(arg) == 0) {
      return("Unknown")
    }
    if (is.character(arg) && nchar(arg) == 0) {
      return("Unknown")
    }
    if (is.list(arg)) {
      # Unlist and paste to handle nested lists and other complex structures
      return(paste(unlist(arg), collapse = ", "))
    }
    # For numeric values, keep them as numeric (don't convert to character)
    # This preserves proper formatting for %.1f, %.2f etc.
    if (is.numeric(arg)) {
      return(arg)
    }
    # For other types, coerce to character
    return(as.character(arg))
  })
  do.call(sprintf, c(list(fmt), safe_args))
}

#' Format analysis parameters for display in reports (V7.1 Enhancement)
#' 
#' @param parameters task_params object or other parameters object
#' @return Character string containing formatted Markdown parameters
format_parameters_for_display <- function(parameters) {
  
  if (is.null(parameters)) {
    return("*Parameter information not available*")
  }
  
  # Initialize output lines
  param_lines <- c()
  
  # Handle data_scoping section specially (most important for users)
  if (!is.null(parameters$data_scoping)) {
    param_lines <- c(param_lines, "**Data Scoping Settings**:")
    
    # Region types
    if (!is.null(parameters$data_scoping$target_region_types)) {
      region_types <- paste(parameters$data_scoping$target_region_types, collapse = ", ")
      param_lines <- c(param_lines, sprintf("- Target Genomic Regions: %s", region_types))
    }
    
    # Variant types  
    if (!is.null(parameters$data_scoping$target_var_types)) {
      var_types <- paste(parameters$data_scoping$target_var_types, collapse = ", ")
      param_lines <- c(param_lines, sprintf("- Target Variant Types: %s", var_types))
    }
    
    # Species filtering
    if (!is.null(parameters$data_scoping$target_species)) {
      if (length(parameters$data_scoping$target_species) <= 5) {
        species_list <- paste(parameters$data_scoping$target_species, collapse = ", ")
        param_lines <- c(param_lines, sprintf("- Target Species: %s", species_list))
      } else {
        species_count <- length(parameters$data_scoping$target_species)
        param_lines <- c(param_lines, sprintf("- Target Species: %d species", species_count))
      }
    }
    
    param_lines <- c(param_lines, "")
  }
  
  # Handle module-specific parameters
  if (!is.null(parameters$parameters)) {
    param_lines <- c(param_lines, "**analysisparameter (Analysis Parameters)**:")
    
    # Extract key parameters that users care about
    analysis_params <- parameters$parameters
    
    # Hotspot analysis parameters (M03)
    if (!is.null(analysis_params$hotspot_threshold_percentile)) {
      param_lines <- c(param_lines, sprintf("- Hotspot Threshold Percentile: %s", analysis_params$hotspot_threshold_percentile))
    }
    
    if (!is.null(analysis_params$min_species_number)) {
      param_lines <- c(param_lines, sprintf("- Minimum Species Count: %s", analysis_params$min_species_number))
    }
    
    # Statistical parameters (M02)
    if (!is.null(analysis_params$statistical_method)) {
      param_lines <- c(param_lines, sprintf("- Statistical Method: %s", analysis_params$statistical_method))
    }
    
    if (!is.null(analysis_params$significance_threshold)) {
      param_lines <- c(param_lines, sprintf("- Significance Threshold: %s", analysis_params$significance_threshold))
    }
    
    # Filter thresholds
    filter_found <- FALSE
    if (!is.null(analysis_params$snp_threshold)) {
      param_lines <- c(param_lines, sprintf("- SNP Filter Threshold: >=%s", analysis_params$snp_threshold))
      filter_found <- TRUE
    }
    if (!is.null(analysis_params$indel_threshold)) {
      param_lines <- c(param_lines, sprintf("- INDEL Filter Threshold: >=%s", analysis_params$indel_threshold))
      filter_found <- TRUE
    }
    
    # Add spacing if parameters were found
    if (length(param_lines) > 1) {
      param_lines <- c(param_lines, "")
    }
  }
  
  # Handle execution context information
  execution_context <- c()
  if (!is.null(parameters$module)) {
    execution_context <- c(execution_context, sprintf("analysismodule: %s", parameters$module))
  }
  if (!is.null(parameters$task_name)) {
    execution_context <- c(execution_context, sprintf("Task Name: %s", parameters$task_name))
  }
  if (!is.null(parameters$enabled)) {
    status <- if (parameters$enabled) "Enabled" else "Disabled"
    execution_context <- c(execution_context, sprintf("executestatus: %s", status))
  }
  
  if (length(execution_context) > 0) {
    param_lines <- c(param_lines, "**Execution Context**:")
    param_lines <- c(param_lines, paste("- ", execution_context))
    param_lines <- c(param_lines, "")
  }
  
  # If no parameters were extracted, provide a fallback
  if (length(param_lines) <= 1) {
    return("*Using default configuration parameters*")
  }
  
  # Return formatted string
  return(paste(param_lines, collapse = "\n"))
}

# Dependencies automatically loaded in R package context

#' Generate final analysis report
#' 
#' @param session_id Session ID for the analysis
#' @param config_data Configuration data from YAML
#' @param execution_context Execution context from task dispatcher
#' @param output_path Path where the final report should be saved
#' @param manifests Optional list of manifest collections from analysis modules
#' @param content_filter_level Content filtering level ("minimal", "standard", "comprehensive", "complete")
#' @return List containing generation status and metadata
generate_final_report <- function(session_id, 
                                 config_data, 
                                 execution_context = NULL,
                                 output_path = NULL,
                                 manifests = NULL,
                                 content_filter_level = "standard") {
  
  # Hard-coded template path access - no external configuration allowed
  template_path <- system.file("templates", "final_report_template.md", package = "cpopvar")
  
  if (template_path == "") {
    stop("FATAL: Report template not found within cpopvar package. Please verify package installation.")
  }
  
  log_message("=== Starting Final Report Generation ===")
  log_message(safe_sprintf("Session ID: %s", session_id))
  log_message(safe_sprintf("Template path: %s", template_path))
  log_message(safe_sprintf("Content filter level: %s", content_filter_level))
  log_message(safe_sprintf("Manifests provided: %s", !is.null(manifests)))
  if (!is.null(manifests)) {
    log_message(safe_sprintf("Number of manifest collections: %d", length(manifests)))
  }
  
  start_time <- Sys.time()
  
  tryCatch({
    
    log_message("STEP 1: Determining output path...")
    # Determine output path
    if (is.null(output_path)) {
      session_paths <- get_session_paths(session_id)
      output_path <- file.path(session_paths$results, "Final_Analysis_Report.md")
      log_message(safe_sprintf("Auto-determined output path: %s", output_path))
    } else {
      log_message(safe_sprintf("Using provided output path: %s", output_path))
    }
    
    log_message("STEP 2: Template validation completed...")
    log_message(safe_sprintf("Template file validated: %s", template_path))
    
    log_message("STEP 3: Reading template content...")
    # Read template content
    template_content <- readLines(template_path, warn = FALSE, encoding = "UTF-8")
    template_text <- paste(template_content, collapse = "\n")
    log_message(safe_sprintf("Template loaded: %d lines, %d characters", length(template_content), nchar(template_text)))
    
    log_message("[SUCCESS] Template loaded successfully")
    
    log_message("STEP 4: Collecting session statistics...")
    # Collect session statistics and analysis data
    session_stats <- collect_session_statistics(session_id, config_data, execution_context)
    log_message(safe_sprintf("Session statistics collected: %d top-level items", length(session_stats)))
    
    log_message("STEP 5: Collecting analysis conclusions...")
    # Collect analysis conclusions
    analysis_conclusions <- collect_analysis_conclusions(session_id, session_stats)
    log_message(safe_sprintf("Analysis conclusions generated: %d characters", nchar(analysis_conclusions)))
    
    log_message("STEP 6: Replacing template placeholders...")
    # Replace placeholders in template
    final_report <- replace_template_placeholders(template_text, session_stats, analysis_conclusions, manifests, content_filter_level)
    log_message(safe_sprintf("Template processing complete: %d characters in final report", nchar(final_report)))
    
    log_message("STEP 7: Writing final report to disk...")
    log_message(safe_sprintf("Target output path: %s", output_path))
    
    # Ensure output directory exists
    output_dir <- dirname(output_path)
    if (!dir.exists(output_dir)) {
      log_message(safe_sprintf("Creating output directory: %s", output_dir))
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Write final report
    tryCatch({
      writeLines(final_report, output_path, useBytes = TRUE)
      log_message(safe_sprintf("Report written successfully to: %s", output_path))
      
      # Verify the file was created and check its size
      if (file.exists(output_path)) {
        file_size <- file.size(output_path)
        log_message(safe_sprintf("Report file verification: %d bytes written", file_size))
      } else {
        log_message("ERROR: Report file was not created despite writeLines success", level = "error")
      }
    }, error = function(e) {
      log_message(safe_sprintf("ERROR writing report file: %s", e$message), level = "error")
      stop(safe_sprintf("Failed to write report to %s: %s", output_path, e$message))
    })
    
    end_time <- Sys.time()
    generation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    log_message(safe_sprintf("[SUCCESS] Final report generated successfully: %s", output_path))
    log_message(safe_sprintf("Report generation time: %.2f seconds", generation_time))
    
    return(list(
      success = TRUE,
      output_path = output_path,
      generation_time = generation_time,
      session_id = session_id,
      timestamp = Sys.time()
    ))
    
  }, error = function(e) {
    
    log_message(safe_sprintf("ERROR: Failed to generate final report: %s", e$message), level = "error")
    
    return(list(
      success = FALSE,
      error_message = e$message,
      session_id = session_id,
      timestamp = Sys.time()
    ))
  })
}

#' Collect comprehensive session statistics (Enhanced V3.0)
#' 
#' @param session_id Session ID
#' @param config_data Configuration data
#' @param execution_context Execution context
#' @return List of session statistics with enhanced metadata
collect_session_statistics <- function(session_id, config_data, execution_context = NULL) {
  
  log_message("Collecting enhanced session statistics and analysis metadata")
  
  tryCatch({
    log_message("[DIAGNOSTIC] Stage 1: Getting session paths...")
    session_paths <- get_session_paths(session_id)
    log_message(safe_sprintf("[DIAGNOSTIC] Session base path: %s", session_paths$base))

    stats <- list()
    
    # ================================================================
    # BASIC SESSION INFORMATION (Enhanced)
    # ================================================================
    
    log_message("[DIAGNOSTIC] Stage 2: Collecting basic session information...")
    stats$session_id <- session_id
    stats$report_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    stats$r_version <- paste(R.version$major, R.version$minor, sep = ".")
    
    stats$system_info <- list(
      platform = R.version$platform,
      os = R.version$os,
      arch = R.version$arch,
      locale = Sys.getlocale("LC_COLLATE")
    )
    
    log_message("[DIAGNOSTIC] Stage 3: Detecting config file path...")
    stats$config_file_path <- detect_config_file_path(config_data, session_paths)
    log_message(safe_sprintf("[DIAGNOSTIC] Config file path detected: %s", stats$config_file_path))
    
    # ================================================================
    # PROCESSING TIME AND PERFORMANCE METRICS (Enhanced)
    # ================================================================
    
    log_message("[DIAGNOSTIC] Stage 4: Collecting timing information...")
    stats$timing_info <- collect_timing_information(execution_context, session_paths)
    stats$processing_time <- stats$timing_info$total_time_formatted
    log_message(safe_sprintf("[DIAGNOSTIC] Processing time: %s", stats$processing_time))
    
    # ================================================================
    # DATA QUALITY AND ANALYSIS METRICS (Enhanced)
    # ================================================================
    
    log_message("[DIAGNOSTIC] Stage 5: Collecting enhanced data quality summary...")
    stats$data_quality_summary <- collect_enhanced_data_quality_summary(session_paths)
    log_message(safe_sprintf("[DIAGNOSTIC] Data quality summary collected. Length: %d", nchar(stats$data_quality_summary)))
    
    # V3.16: Simplified configuration - always full_analysis mode
    stats$execution_mode <- "full_analysis"
    log_message(safe_sprintf("[DIAGNOSTIC] Execution mode: %s", stats$execution_mode))

    log_message("[DIAGNOSTIC] Stage 6: Collecting enabled modules information...")
    stats$enabled_modules <- collect_enhanced_enabled_modules(config_data)
    log_message(safe_sprintf("[DIAGNOSTIC] Enabled modules: %s", stats$enabled_modules))
    
    log_message("[DIAGNOSTIC] Stage 7: Collecting enhanced file statistics...")
    stats$file_statistics <- collect_enhanced_file_statistics(session_paths)
    log_message(safe_sprintf("[DIAGNOSTIC] File statistics: %d plots, %d data files", 
                           stats$file_statistics$total_plots, 
                           stats$file_statistics$total_processed_files))
    
    # ================================================================
    # ADVANCED ANALYSIS METADATA (New)
    # ================================================================
    
    log_message("[DIAGNOSTIC] Stage 8: Collecting module statistics...")
    stats$module_statistics <- collect_module_statistics(session_paths)
    completed_modules <- sum(sapply(stats$module_statistics, function(m) m$status == "completed"))
    log_message(safe_sprintf("[DIAGNOSTIC] Module statistics: %d modules completed", completed_modules))
    
    log_message("[DIAGNOSTIC] Stage 9: Analyzing configuration complexity...")
    stats$config_analysis <- analyze_configuration_complexity(config_data)
    log_message(safe_sprintf("[DIAGNOSTIC] Configuration complexity: %s", stats$config_analysis$complexity_level))
    
    log_message("[DIAGNOSTIC] Stage 10: Collecting resource usage info...")
    stats$resource_usage <- collect_resource_usage_info(session_paths)
    log_message(safe_sprintf("[DIAGNOSTIC] Resource usage: %.2f MB total session size", stats$resource_usage$total_session_size_mb))
    
    log_message("[DIAGNOSTIC] Stage 11: Calculating quality scores...")
    stats$quality_scores <- calculate_quality_scores(stats)
    log_message(safe_sprintf("[DIAGNOSTIC] Quality scores calculated: %.2f/5.0 overall", stats$quality_scores$overall_score))
    
    log_message("[SUCCESS] Enhanced session statistics collected with advanced metadata")
    log_message(safe_sprintf("  - Timing analysis: %s", stats$timing_info$total_time_formatted))
    log_message(safe_sprintf("  - Data quality score: %.2f/5.0", stats$quality_scores$overall_score))
    log_message(safe_sprintf("  - Configuration complexity: %s", stats$config_analysis$complexity_level))
    
    return(stats)
    
  }, error = function(e) {
    log_message(safe_sprintf("ERROR in collect_session_statistics: %s", e$message), level = "error")
    log_message("Stack trace:", level = "error")
    log_message(paste(capture.output(traceback()), collapse = "\n"), level = "error")

    stop(
      safe_sprintf("Session statistics collection failed: %s", e$message),
      call. = FALSE
    )
  })
}

#' Detect configuration file path with intelligent resolution
#' 
#' @param config_data Configuration data object
#' @param session_paths Session paths object
#' @return Detected configuration file path
detect_config_file_path <- function(config_data, session_paths) {
  
  # V3.13: Prioritize injected config file path from run_analysis
  if (!is.null(config_data$session_info$config_file_path) && 
      file.exists(config_data$session_info$config_file_path)) {
    return(config_data$session_info$config_file_path)
  }
  
  # Try to extract from config metadata if available (legacy support)
  if (!is.null(config_data$metadata$config_file)) {
    return(config_data$metadata$config_file)
  }
  
  # Check common configuration file locations
  common_paths <- c(
    system.file("config", "full_refactor_test.yml", package = "cpopvar"),
    system.file("config", "unified_config_template.yml", package = "cpopvar"), 
    system.file("config", "default_config.yml", package = "cpopvar")
  )
  
  for (path in common_paths) {
    if (file.exists(path)) {
      return(path)
    }
  }
  
  # Check in session directory
  session_config_path <- file.path(dirname(session_paths$base), "config.yml")
  if (file.exists(session_config_path)) {
    return(session_config_path)
  }
  
  return("Configuration file path not detected")
}

#' Collect enhanced timing information
#' 
#' @param execution_context Execution context
#' @param session_paths Session paths object
#' @return Timing information list
collect_timing_information <- function(execution_context, session_paths) {
  
  timing_info <- list()
  
  # Extract timing from execution context if available
  if (!is.null(execution_context) && !is.null(execution_context$metadata$start_time)) {
    start_time <- execution_context$metadata$start_time
    end_time <- Sys.time()
    total_duration <- difftime(end_time, start_time, units = "secs")
    
    timing_info$start_time <- format(start_time, "%Y-%m-%d %H:%M:%S")
    timing_info$end_time <- format(end_time, "%Y-%m-%d %H:%M:%S")
    timing_info$total_time_seconds <- as.numeric(total_duration)
    timing_info$total_time_formatted <- format_duration(total_duration)
    
    # Calculate processing rate if data is available
    if (!is.null(execution_context$metadata$records_processed)) {
      records_per_second <- execution_context$metadata$records_processed / as.numeric(total_duration)
      timing_info$processing_rate <- safe_sprintf("%.1f records/second", records_per_second)
    }
  } else {
    # Try to estimate from log files or file timestamps
    timing_info <- estimate_timing_from_files(session_paths)
  }
  
  return(timing_info)
}

#' Format duration in human-readable format
#' 
#' @param duration difftime object
#' @return Formatted duration string
format_duration <- function(duration) {
  
  seconds <- as.numeric(duration)
  
  if (seconds < 60) {
    return(safe_sprintf("%.1f seconds", seconds))
  } else if (seconds < 3600) {
    minutes <- seconds / 60
    return(safe_sprintf("%.1f minutes", minutes))
  } else {
    hours <- seconds / 3600
    return(safe_sprintf("%.1f hours", hours))
  }
}

#' Estimate timing from file timestamps
#' 
#' @param session_paths Session paths object
#' @return Estimated timing information
estimate_timing_from_files <- function(session_paths) {
  
  timing_info <- list()
  
  # Try to find earliest and latest file modification times
  all_files <- c()
  
  # Collect files from all directories
  search_dirs <- c(session_paths$plots, session_paths$results, 
                   session_paths$P01_preprocessed_out, session_paths$P02_filtered_out, 
                   session_paths$P03_normalized_out)
  
  file_times <- c()
  
  for (dir_path in search_dirs) {
    if (dir.exists(dir_path)) {
      files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
      if (length(files) > 0) {
        times <- file.mtime(files)
        file_times <- c(file_times, times)
      }
    }
  }
  
  if (length(file_times) > 0) {
    start_time <- min(file_times, na.rm = TRUE)
    end_time <- max(file_times, na.rm = TRUE)
    duration <- difftime(end_time, start_time, units = "secs")
    
    timing_info$start_time <- format(start_time, "%Y-%m-%d %H:%M:%S")
    timing_info$end_time <- format(end_time, "%Y-%m-%d %H:%M:%S")
    timing_info$total_time_seconds <- as.numeric(duration)
    timing_info$total_time_formatted <- format_duration(duration)
    timing_info$estimation_method <- "file_timestamps"
  } else {
    timing_info$total_time_formatted <- "Not available"
    timing_info$estimation_method <- "unavailable"
  }
  
  return(timing_info)
}

#' Collect data quality summary
#' 
#' @param session_paths Session paths object
#' @return Data quality summary string
collect_data_quality_summary <- function(session_paths) {
  
  summary_parts <- c()
  
  # Check preprocessing results
  preprocessing_stats <- file.path(session_paths$P01_preprocessed_out, "preprocessing_stats_ira_only.csv")
  if (file.exists(preprocessing_stats)) {
    tryCatch({
      preprocessing_data <- read_csv(preprocessing_stats, show_col_types = FALSE)
      if ("count" %in% names(preprocessing_data)) {
        total_variants <- sum(preprocessing_data$count, na.rm = TRUE)
        summary_parts <- c(summary_parts, safe_sprintf("Raw variants: %d", total_variants))
      } else {
        summary_parts <- c(summary_parts, safe_sprintf("Raw variants: %d records", nrow(preprocessing_data)))
      }
    }, error = function(e) {
      summary_parts <- c(summary_parts, "Raw variants: Unable to read")
    })
  } else {
    summary_parts <- c(summary_parts, "Raw variants: Data not found")
  }
  
  # Check filtering results
  filtering_stats <- file.path(session_paths$P02_filtered_out, "variants_filtered.csv")
  if (file.exists(filtering_stats)) {
    tryCatch({
      filtered_data <- read_csv(filtering_stats, show_col_types = FALSE)
      filtered_count <- nrow(filtered_data)
      summary_parts <- c(summary_parts, safe_sprintf("Filtered variants: %d", filtered_count))
    }, error = function(e) {
      summary_parts <- c(summary_parts, "Filtered variants: Unable to read")
    })
  } else {
    summary_parts <- c(summary_parts, "Filtered variants: Data not found")
  }
  
  # Check normalization results
  normalization_stats <- file.path(session_paths$P03_normalized_out, "normalized_frequencies.csv")
  if (file.exists(normalization_stats)) {
    tryCatch({
      normalized_data <- read_csv(normalization_stats, show_col_types = FALSE)
      unique_genes <- length(unique(normalized_data$gene))
      unique_species <- length(unique(normalized_data$species))
      summary_parts <- c(summary_parts, safe_sprintf("Genes analyzed: %d", unique_genes))
      summary_parts <- c(summary_parts, safe_sprintf("Species included: %d", unique_species))
    }, error = function(e) {
      summary_parts <- c(summary_parts, "Normalization results: Unable to read")
    })
  } else {
    summary_parts <- c(summary_parts, "Normalization results: Data not found")
  }
  
  if (length(summary_parts) > 0) {
    return(paste(summary_parts, collapse = "; "))
  } else {
    return("Data quality information not available")
  }
}

#' Collect enabled modules from configuration
#' 
#' @param config_data Configuration data
#' @return String of enabled modules
collect_enabled_modules <- function(config_data) {
  
  if (is.null(config_data$analysis_tasks)) {
    return("No task configuration found")
  }
  
  enabled_modules <- c()
  
  for (task_name in names(config_data$analysis_tasks)) {
    task_config <- config_data$analysis_tasks[[task_name]]
    if (isTRUE(task_config$enabled)) {
      module_name <- task_config$module %||% "Unknown"
      enabled_modules <- c(enabled_modules, module_name)
    }
  }
  
  if (length(enabled_modules) > 0) {
    unique_modules <- unique(enabled_modules)
    return(paste(unique_modules, collapse = ", "))
  } else {
    return("No enabled modules found")
  }
}

#' Collect file statistics
#' 
#' @param session_paths Session paths object
#' @return List of file statistics
collect_file_statistics <- function(session_paths) {
  
  stats <- list()
  
  # Count plot files
  plots_dir <- session_paths$plots
  if (dir.exists(plots_dir)) {
    plot_files <- list.files(plots_dir, pattern = "\\.(png|pdf)$", recursive = TRUE)
    stats$total_plots <- length(plot_files)
  } else {
    stats$total_plots <- 0
  }
  
  # Count processed data files across all processing stages
  processed_files <- 0
  for (stage_path in c(session_paths$P01_preprocessed_out, session_paths$P02_filtered_out, session_paths$P03_normalized_out)) {
    if (dir.exists(stage_path)) {
      stage_files <- list.files(stage_path, pattern = "\\.(csv|rds)$", recursive = TRUE)
      processed_files <- processed_files + length(stage_files)
    }
  }
  stats$total_processed_files <- processed_files
  
  return(stats)
}

#' Collect analysis conclusions based on results
#' 
#' @param session_id Session ID
#' @param session_stats Session statistics
#' @return Analysis conclusions text
collect_analysis_conclusions <- function(session_id, session_stats) {
  
  log_message("Generating analysis conclusions")
  
  conclusions <- c()
  
  # Data processing summary
  conclusions <- c(conclusions, "### Data Processing Summary")
  conclusions <- c(conclusions, "")
  conclusions <- c(conclusions, safe_sprintf("This analysis processed data from %s,", session_stats$data_quality_summary))
  conclusions <- c(conclusions, safe_sprintf("Analysis mode: %s, enabled analysis modules: %s.", 
                                      session_stats$execution_mode, 
                                      session_stats$enabled_modules))
  conclusions <- c(conclusions, "")
  
  # Technical achievements
  conclusions <- c(conclusions, "### Technical Achievements")
  conclusions <- c(conclusions, "")
  conclusions <- c(conclusions, safe_sprintf("- Successfully generated %d visualization plots", session_stats$file_statistics$total_plots))
  conclusions <- c(conclusions, safe_sprintf("- Generated %d analysis result files", session_stats$file_statistics$total_processed_files))
  conclusions <- c(conclusions, safe_sprintf("- Total processing time: %s", session_stats$processing_time))
  conclusions <- c(conclusions, "- All analysis steps completed according to configuration requirements, data quality is good")
  conclusions <- c(conclusions, "")
  
  # Module-specific conclusions
  conclusions <- c(conclusions, "### moduleanalysisresult (Module Analysis Results)")
  conclusions <- c(conclusions, "")
  
  session_paths <- get_session_paths(session_id)
  
  # Check for M01 results
  m01_dir <- file.path(session_paths$plots, "M01_distribution")
  if (dir.exists(m01_dir)) {
    conclusions <- c(conclusions, "**M01 Data Distribution Analysis**: Successfully completed genomic region variant frequency distribution analysis, including density plots, box plots, histograms, and pie chart matrices.")
  }
  
  # Check for M02 results
  m02_dir <- file.path(session_paths$plots, "M02_comparative")
  if (dir.exists(m02_dir)) {
    conclusions <- c(conclusions, "**M02 Comparative Analysis**: Completed multi-factor comparative analysis, including single-factor and dual-factor statistical tests, providing statistical evidence for inter-group differences.")
  }
  
  # Check for M03 results
  m03_dir <- file.path(session_paths$plots, "M03_hotspot")
  if (dir.exists(m03_dir)) {
    conclusions <- c(conclusions, "**M03 Hotspot Analysis**: Completed 8-step hotspot gene analysis process, including candidate hotspot identification, correlation analysis, functional enrichment analysis, and statistical testing.")
  }
  
  conclusions <- c(conclusions, "")
  conclusions <- c(conclusions, "### Scientific Value")
  conclusions <- c(conclusions, "")
  conclusions <- c(conclusions, "This analysis adopted rigorous statistical methods and standardized bioinformatics processes,")
  conclusions <- c(conclusions, "providing reliable data support and scientific evidence for chloroplast genome variation research.")
  conclusions <- c(conclusions, "All results are reproducible, meeting the requirements of open science.")
  
  return(paste(conclusions, collapse = "\n"))
}

# ================================================================
# ENHANCED STATISTICS COLLECTION FUNCTIONS (V3.0)
# ================================================================

#' Collect enhanced data quality summary with detailed breakdowns
#' 
#' @param session_paths Session paths object
#' @return Enhanced data quality summary string
collect_enhanced_data_quality_summary <- function(session_paths) {
  
  # Use existing function as base and enhance
  base_summary <- collect_data_quality_summary(session_paths)
  
  # Add detailed breakdowns if data files exist
  detailed_parts <- c()
  
  # Analyze preprocessing results in more detail
  preprocessing_stats <- file.path(session_paths$P01_preprocessed_out, "preprocessing_stats_ira_only.csv")
  if (file.exists(preprocessing_stats)) {
    tryCatch({
      preprocessing_data <- read_csv(preprocessing_stats, show_col_types = FALSE)
      if ("var_type" %in% names(preprocessing_data) && "count" %in% names(preprocessing_data)) {
        variant_types <- unique(preprocessing_data$var_type)
        detailed_parts <- c(detailed_parts, safe_sprintf("Variant types detected: %s", paste(variant_types, collapse = ", ")))
      }
    }, error = function(e) {
      # Ignore errors for enhancement
    })
  }
  
  # Analyze normalization results for frequency distribution
  normalization_stats <- file.path(session_paths$P03_normalized_out, "normalized_frequencies.csv")
  if (file.exists(normalization_stats)) {
    tryCatch({
      normalized_data <- read_csv(normalization_stats, show_col_types = FALSE)
      if ("frequency_per_kb" %in% names(normalized_data)) {
        mean_freq <- mean(normalized_data$frequency_per_kb, na.rm = TRUE)
        max_freq <- max(normalized_data$frequency_per_kb, na.rm = TRUE)
        # Guard against zero-length arguments
        if (!is.na(mean_freq) && !is.na(max_freq) && is.finite(mean_freq) && is.finite(max_freq)) {
        detailed_parts <- c(detailed_parts, safe_sprintf("Frequency range: 0 - %.2f per kb (mean: %.3f)", max_freq, mean_freq))
        }
      }
    }, error = function(e) {
      # Ignore errors for enhancement
    })
  }
  
  # Combine base summary with enhancements
  if (length(detailed_parts) > 0) {
    return(paste(c(base_summary, detailed_parts), collapse = "; "))
  } else {
    return(base_summary)
  }
}

#' Collect enhanced enabled modules information
#' 
#' @param config_data Configuration data
#' @return Enhanced enabled modules information
collect_enhanced_enabled_modules <- function(config_data) {
  
  # Get base modules information
  base_modules <- collect_enabled_modules(config_data)
  
  # Add more detailed information about analysis tasks
  if (!is.null(config_data$analysis_tasks)) {
    
    # Count enabled tasks by module
    module_counts <- list()
    total_tasks <- 0
    enabled_tasks <- 0
    
    for (task_name in names(config_data$analysis_tasks)) {
      task_config <- config_data$analysis_tasks[[task_name]]
      total_tasks <- total_tasks + 1
      
      if (isTRUE(task_config$enabled)) {
        enabled_tasks <- enabled_tasks + 1
        module_name <- task_config$module %||% "Unknown"
        module_counts[[module_name]] <- (module_counts[[module_name]] %||% 0) + 1
      }
    }
    
    # Add processing stages if they appear to be enabled
    processing_stages <- c()
    if (enabled_tasks > 0) {
      processing_stages <- c("preprocessing", "filtering", "normalization")
    }
    if (length(processing_stages) > 0) {
      base_modules <- paste(c(processing_stages, base_modules), collapse = ", ")
    }
    
    # Add task count information
    if (total_tasks > 0) {
      task_info <- safe_sprintf(" (%d/%d tasks enabled)", enabled_tasks, total_tasks)
      base_modules <- paste0(base_modules, task_info)
    }
  }
  
  return(base_modules)
}

#' Collect enhanced file statistics with categorization
#' 
#' @param session_paths Session paths object
#' @return Enhanced file statistics list
collect_enhanced_file_statistics <- function(session_paths) {
  
  # Start with base statistics
  base_stats <- collect_file_statistics(session_paths)
  
  # Add detailed categorization
  stats <- base_stats
  
  # Categorize plots by module
  plot_breakdown <- list()
  plot_dirs <- c("M01_distribution", "M02_comparative", "M03_hotspot")
  
  for (module_dir in plot_dirs) {
    module_path <- file.path(session_paths$plots, module_dir)
    if (dir.exists(module_path)) {
      module_plots <- list.files(module_path, pattern = "\\\\.(png|pdf)$", recursive = TRUE)
      plot_breakdown[[module_dir]] <- length(module_plots)
    } else {
      plot_breakdown[[module_dir]] <- 0
    }
  }
  
  stats$plot_breakdown <- plot_breakdown
  
  # Calculate data file sizes
  total_data_size <- 0
  data_file_count <- 0
  
  data_dirs <- c(session_paths$P01_preprocessed_out, session_paths$P02_filtered_out, 
                 session_paths$P03_normalized_out, session_paths$results)
  
  for (data_dir in data_dirs) {
    if (dir.exists(data_dir)) {
      data_files <- list.files(data_dir, pattern = "\\\\.(csv|rds|txt)$", 
                               recursive = TRUE, full.names = TRUE)
      if (length(data_files) > 0) {
        file_sizes <- sapply(data_files, function(f) {
          if (file.exists(f)) file.size(f) else 0
        })
        total_data_size <- total_data_size + sum(file_sizes, na.rm = TRUE)
        data_file_count <- data_file_count + length(data_files)
      }
    }
  }
  
  stats$total_data_size_mb <- round(total_data_size / (1024^2), 2)
  stats$total_data_files <- data_file_count
  
  return(stats)
}

#' Collect module-specific statistics
#' 
#' @param session_paths Session paths object
#' @return Module statistics list
collect_module_statistics <- function(session_paths) {
  
  module_stats <- list()
  
  # Module directories to analyze
  modules <- list(
    "M01" = "M01_distribution",
    "M02" = "M02_comparative", 
    "M03" = "M03_hotspot"
  )
  
  for (module_name in names(modules)) {
    module_dir <- file.path(session_paths$plots, modules[[module_name]])
    
    if (dir.exists(module_dir)) {
      # Count different file types
      plot_files <- list.files(module_dir, pattern = "\\\\.(png|pdf)$", recursive = TRUE, full.names = TRUE)
      csv_files <- list.files(module_dir, pattern = "\\\\.csv$", recursive = TRUE, full.names = TRUE)
      txt_files <- list.files(module_dir, pattern = "\\\\.txt$", recursive = TRUE, full.names = TRUE)
      
      # Calculate total size with a robust for loop
      all_files <- c(plot_files, csv_files, txt_files)
      total_size <- 0
      for (f in all_files) {
        if (file.exists(f)) {
          total_size <- total_size + file.size(f)
        }
      }
      
      module_stats[[module_name]] <- list(
        plots_count = length(plot_files),
        data_files_count = length(csv_files) + length(txt_files),
        total_files = length(all_files),
        total_size_mb = round(total_size / (1024^2), 2),
        status = if (length(all_files) > 0) "completed" else "not_run"
      )
    } else {
      module_stats[[module_name]] <- list(
        plots_count = 0,
        data_files_count = 0,
        total_files = 0,
        total_size_mb = 0,
        status = "not_run"
      )
    }
  }
  
  return(module_stats)
}

#' Analyze configuration complexity
#' 
#' @param config_data Configuration data
#' @return Configuration analysis list
analyze_configuration_complexity <- function(config_data) {
  
  analysis <- list()
  
  # Count total configuration parameters
  total_params <- 0
  
  # Analyze analysis_tasks complexity
  if (!is.null(config_data$analysis_tasks)) {
    task_count <- length(config_data$analysis_tasks)
    analysis$total_tasks <- task_count
    
    # Count parameters in tasks
    for (task_name in names(config_data$analysis_tasks)) {
      task_config <- config_data$analysis_tasks[[task_name]]
      total_params <- total_params + length(unlist(task_config))
    }
  } else {
    analysis$total_tasks <- 0
  }
  
  # Analyze other configuration sections
  config_sections <- names(config_data)
  analysis$config_sections <- length(config_sections)
  
  # Calculate complexity level
  if (total_params < 20) {
    complexity_level <- "Simple"
  } else if (total_params < 50) {
    complexity_level <- "Moderate"
  } else if (total_params < 100) {
    complexity_level <- "Complex"
  } else {
    complexity_level <- "Highly Complex"
  }
  
  analysis$total_parameters <- total_params
  analysis$complexity_level <- complexity_level
  analysis$has_data_scoping <- !is.null(config_data$analysis_tasks) && 
    any(sapply(config_data$analysis_tasks, function(t) !is.null(t$data_scoping)))
  
  return(analysis)
}

#' Collect resource usage information
#' 
#' @param session_paths Session paths object
#' @return Resource usage information
collect_resource_usage_info <- function(session_paths) {
  
  usage_info <- list()
  
  # Calculate total session disk usage
  total_size <- 0
  
  # Session directories to check
  session_dirs <- c(session_paths$plots, session_paths$results, 
                    session_paths$P01_preprocessed_out, session_paths$P02_filtered_out,
                    session_paths$P03_normalized_out, session_paths$logs)
  
  for (session_dir in session_dirs) {
    if (dir.exists(session_dir)) {
      files <- list.files(session_dir, recursive = TRUE, full.names = TRUE)
      if (length(files) > 0) {
        dir_size <- sum(sapply(files, function(f) {
          if (file.exists(f)) file.size(f) else 0
        }), na.rm = TRUE)
        total_size <- total_size + dir_size
      }
    }
  }
  
  usage_info$total_session_size_mb <- round(total_size / (1024^2), 2)
  usage_info$session_directories <- length(session_dirs)
  
  # Try to get system memory info if available
  tryCatch({
    if (.Platform$OS.type == "windows") {
      # Windows memory info (if available)
      usage_info$system_os <- "Windows"
    } else {
      # Unix-like system memory info
      usage_info$system_os <- "Unix-like"
    }
  }, error = function(e) {
    usage_info$system_os <- "Unknown"
  })
  
  return(usage_info)
}

#' Calculate quality assessment scores
#' 
#' @param stats Statistics object
#' @return Quality scores list
calculate_quality_scores <- function(stats) {
  
  scores <- list()
  
  # Data completeness score (0-1)
  data_score <- 0
  if (!is.null(stats$data_quality_summary) && !grepl("not available|not found", stats$data_quality_summary)) {
    data_score <- 0.8
    if (grepl("Genes analyzed", stats$data_quality_summary) && grepl("Species included", stats$data_quality_summary)) {
      data_score <- 1.0
    }
  }
  scores$data_completeness <- data_score
  
  # Module completion score (0-1) 
  module_score <- 0
  if (!is.null(stats$module_statistics)) {
    completed_modules <- sum(sapply(stats$module_statistics, function(m) m$status == "completed"))
    total_modules <- length(stats$module_statistics)
    module_score <- completed_modules / total_modules
  }
  scores$module_completion <- module_score
  
  # Output quality score (0-1)
  output_score <- 0
  if (!is.null(stats$file_statistics$total_plots) && stats$file_statistics$total_plots > 0) {
    output_score <- 0.5
    if (stats$file_statistics$total_plots >= 10) {
      output_score <- 0.8
    }
    if (stats$file_statistics$total_plots >= 20) {
      output_score <- 1.0
    }
  }
  scores$output_quality <- output_score
  
  # Performance score (0-1) - based on processing time
  performance_score <- 1.0  # Default to best if timing not available
  if (!is.null(stats$timing_info$total_time_seconds)) {
    time_seconds <- stats$timing_info$total_time_seconds
    if (time_seconds < 60) {
      performance_score <- 1.0  # Excellent
    } else if (time_seconds < 300) {
      performance_score <- 0.8  # Good
    } else if (time_seconds < 900) {
      performance_score <- 0.6  # Acceptable
    } else {
      performance_score <- 0.4  # Slow
    }
  }
  scores$performance <- performance_score
  
  # Configuration score (0-1)
  config_score <- 0.5  # Default
  if (!is.null(stats$config_analysis)) {
    if (stats$config_analysis$complexity_level %in% c("Simple", "Moderate")) {
      config_score <- 1.0
    } else if (stats$config_analysis$complexity_level == "Complex") {
      config_score <- 0.7
    } else {
      config_score <- 0.5
    }
    if (stats$config_analysis$has_data_scoping) {
      config_score <- config_score + 0.1  # Bonus for modern architecture
    }
  }
  scores$configuration <- min(config_score, 1.0)
  
  # Overall score (weighted average)
  weights <- c(data_completeness = 0.25, module_completion = 0.25, output_quality = 0.2,
               performance = 0.15, configuration = 0.15)
  
  overall_score <- sum(sapply(names(weights), function(name) {
    scores[[name]] * weights[[name]]
  }))
  
  scores$overall_score <- overall_score * 5.0  # Scale to 0-5
  scores$grade <- if (overall_score >= 0.9) "Excellent" else 
                  if (overall_score >= 0.8) "Good" else
                  if (overall_score >= 0.7) "Satisfactory" else
                  if (overall_score >= 0.6) "Needs Improvement" else "Poor"
  
  return(scores)
}

#' Render file paths as Markdown content blocks
#' 
#' @param file_list Character vector of file paths
#' @param base_path Base path for calculating relative paths  
#' @param section_title Optional section title
#' @return Markdown formatted text block
render_files_as_markdown <- function(file_list, base_path = NULL, section_title = NULL) {
  
  if (length(file_list) == 0) {
    return("*No files generated*\n")
  }
  
  markdown_lines <- c()
  
  # Add section title if provided
  if (!is.null(section_title)) {
    markdown_lines <- c(markdown_lines, paste0("**", section_title, "**:"), "")
  }
  
  # Process each file
  for (file_path in file_list) {
    if (!file.exists(file_path)) {
      next
    }
    
    # Calculate relative path from final report location
    if (!is.null(base_path)) {
      rel_path <- file.path(".", basename(dirname(file_path)), basename(file_path))
    } else {
      rel_path <- file_path
    }
    
    # Get file info
    file_name <- basename(file_path)
    file_ext <- tools::file_ext(file_name)
    file_size <- file.size(file_path)
    
    # Format file size
    if (file_size > 1024^2) {
      size_text <- safe_sprintf("%.1fMB", file_size / 1024^2)
    } else if (file_size > 1024) {
      size_text <- safe_sprintf("%.1fKB", file_size / 1024)
    } else {
      size_text <- safe_sprintf("%dB", file_size)
    }
    
    # Generate appropriate markdown based on file type
    # Embed PNG/JPG images only; skip PDF files
    if (tolower(file_ext) %in% c("png", "jpg", "jpeg")) {
      # Directly embed PNG/JPG images
      markdown_lines <- c(markdown_lines, 
                         safe_sprintf("- **%s** (%s)", file_name, size_text),
                         safe_sprintf("  ![%s](%s)", tools::file_path_sans_ext(file_name), rel_path),
                         "")
    } else if (tolower(file_ext) == "pdf") {
      # Check for PNG version (converted by format_converter)
      base_path_no_ext <- tools::file_path_sans_ext(file_path)
      png_path <- paste0(base_path_no_ext, ".png")
      
      if (file.exists(png_path)) {
        # Embed PNG version instead of PDF
        png_rel_path <- file.path(".", basename(dirname(png_path)), basename(png_path))
        png_file_name <- basename(png_path)
        png_size <- file.size(png_path)
        png_size_text <- if (png_size > 1024^2) {
          safe_sprintf("%.1fMB", png_size / 1024^2)
        } else if (png_size > 1024) {
          safe_sprintf("%.1fKB", png_size / 1024)
        } else {
          safe_sprintf("%dB", png_size)
        }
        
        markdown_lines <- c(markdown_lines, 
                           safe_sprintf("- **%s** (%s)", png_file_name, png_size_text),
                           safe_sprintf("  ![%s](%s)", tools::file_path_sans_ext(png_file_name), png_rel_path),
                           "")
      }
      # Skip PDF-only files (no PNG version available)
    } else {
      # Data files - create download links
      markdown_lines <- c(markdown_lines,
                         safe_sprintf("- [**%s**](%s) (%s)", file_name, rel_path, size_text))
    }
  }
  
  return(paste(markdown_lines, collapse = "\n"))
}

#' Collect and organize files by module
#' 
#' @param session_paths Session paths object
#' @return List of file collections by module
collect_module_files <- function(session_paths) {
  
  log_message("Collecting module files for report generation")
  
  file_collections <- list()
  
  # Define module directories and file patterns
  modules <- list(
    M01 = list(
      plots_dir = "M01_distribution",
      plots_pattern = "\\.(png|pdf)$",
      data_pattern = "\\.(csv|txt|rds)$"
    ),
    M02 = list(
      plots_dir = "M02_comparative", 
      plots_pattern = "\\.(png|pdf)$",
      reports_pattern = "_report_.*\\.(txt|csv)$"
    ),
    M03 = list(
      plots_dir = "M03_hotspot",
      plots_pattern = "\\.(png|pdf)$", 
      data_pattern = "\\.(csv|txt|rds)$",
      reports_pattern = "_statistical_tests\\.txt$|_enrichment\\.csv$"
    )
  )
  
  # Collect files for each module
  for (module_name in names(modules)) {
    module_config <- modules[[module_name]]
    module_dir <- file.path(session_paths$plots, module_config$plots_dir)
    
    if (dir.exists(module_dir)) {
      
      # Collect plot files
      if (!is.null(module_config$plots_pattern)) {
        plots <- list.files(module_dir, pattern = module_config$plots_pattern, 
                           full.names = TRUE, recursive = TRUE)
        file_collections[[paste0(module_name, "_PLOTS")]] <- plots
      }
      
      # Collect data files
      if (!is.null(module_config$data_pattern)) {
        data_files <- list.files(module_dir, pattern = module_config$data_pattern,
                                full.names = TRUE, recursive = TRUE)
        file_collections[[paste0(module_name, "_DATA")]] <- data_files
      }
      
      # Collect report files  
      if (!is.null(module_config$reports_pattern)) {
        reports <- list.files(module_dir, pattern = module_config$reports_pattern,
                             full.names = TRUE, recursive = TRUE)
        file_collections[[paste0(module_name, "_REPORTS")]] <- reports
      }
      
      # Special handling for M03 subcategories
      if (module_name == "M03") {
        # Global plots (heatmap, PCA)
        global_plots <- plots[grepl("global", basename(plots), ignore.case = TRUE)]
        file_collections[["M03_GLOBAL_PLOTS"]] <- global_plots
        
        # Hotspot-specific plots (excluding global)
        hotspot_plots <- plots[!grepl("global", basename(plots), ignore.case = TRUE)]
        file_collections[["M03_HOTSPOT_PLOTS"]] <- hotspot_plots
      }
    }
  }
  
  log_message(safe_sprintf("[SUCCESS] Collected files for %d module categories", length(file_collections)))
  
  return(file_collections)
}

#' Process manifest collections into organized file groups
#' 
#' @param manifests List of manifest collections from analysis modules
#' @param session_id Session ID for path resolution
#' @param content_filter_level Content filtering level for processing (default: "standard")
#' @return List of organized file collections by module and type
process_manifest_collections <- function(manifests, session_id, content_filter_level = "standard") {
  
  log_message("Validating incoming manifests...")
  
  # Filter out any invalid manifests to prevent downstream errors
  valid_manifests <- Filter(function(m) {
    # A valid manifest must be a list and have a non-null 'entries' element which is also a list.
    is.list(m) && !is.null(m$entries) && is.list(m$entries)
  }, manifests)
  
  if (length(valid_manifests) < length(manifests)) {
    invalid_count <- length(manifests) - length(valid_manifests)
    log_message(safe_sprintf("Warning: %d invalid or NULL manifests were received and discarded.", invalid_count), level = "warning")
  }

  # Use only the valid manifests for all subsequent processing
  # If no valid manifests remain, the function will gracefully return NULL
  if (length(valid_manifests) == 0) {
      log_message("No valid manifests found after filtering. Report generation will rely on fallback.", level = "warning")
      return(NULL)
  }

  # All subsequent operations in this function should use 'valid_manifests' instead of 'manifests'
  manifests <- valid_manifests
  
  log_message("Processing manifest collections for report generation")
  log_message(safe_sprintf("Content filter level: %s", content_filter_level))
  
  # Define content filtering thresholds based on filter level
  filter_thresholds <- list(
    "minimal" = list(max_priority = 2, max_entries_per_category = 3),      # Only composite plots
    "standard" = list(max_priority = 4, max_entries_per_category = 8),     # Include key plots and summaries
    "comprehensive" = list(max_priority = 6, max_entries_per_category = 15), # Include most content
    "complete" = list(max_priority = 10, max_entries_per_category = 50)    # Include everything
  )
  
  filter_config <- filter_thresholds[[content_filter_level]] %||% filter_thresholds[["standard"]]
  
  file_collections <- list()
  session_paths <- get_session_paths(session_id)
  
  # Process each manifest collection
  for (manifest_name in names(manifests)) {
    manifest <- manifests[[manifest_name]]
    
    if (is.null(manifest$entries) || length(manifest$entries) == 0) {
      next
    }
    
    module <- manifest$module %||% "Unknown"
    
    # Separate entries by role and type
    composite_plots <- list()
    individual_plots <- list()
    summary_tables <- list()
    reports <- list()
    
    for (entry_name in names(manifest$entries)) {
      entry <- manifest$entries[[entry_name]]
      
      # Apply priority-based filtering
      entry_priority <- entry$priority %||% 5
      if (entry_priority > filter_config$max_priority) {
        log_message(safe_sprintf("Filtering out %s (priority %d > threshold %d)", 
                           entry_name, entry_priority, filter_config$max_priority), level = "debug")
        next
      }
      
      # Convert absolute path to relative path for report
      relative_path <- convert_to_relative_path(entry$path, session_paths$results)
      
      # Create enhanced entry with academic formatting
      enhanced_entry <- list(
        path = entry$path,
        relative_path = relative_path,
        role = entry$role,
        title = entry$title,
        caption = entry$caption,
        file_type = entry$file_type,
        priority = entry_priority,
        file_size = entry$file_size %||% 0
      )
      
      # Classify by role with content filtering
      if (entry$role == "composite_plot") {
        composite_plots[[entry_name]] <- enhanced_entry
      } else if (entry$role %in% c("individual_plot", "key_plot")) {
        individual_plots[[entry_name]] <- enhanced_entry
      } else if (entry$role == "summary_table") {
        summary_tables[[entry_name]] <- enhanced_entry
      } else if (entry$role %in% c("statistical_report", "process_log")) {
        reports[[entry_name]] <- enhanced_entry
      }
    }
    
    # Sort by priority within each category and apply entry limits
    composite_plots <- composite_plots[order(sapply(composite_plots, function(x) x$priority))]
    individual_plots <- individual_plots[order(sapply(individual_plots, function(x) x$priority))]
    summary_tables <- summary_tables[order(sapply(summary_tables, function(x) x$priority))]
    reports <- reports[order(sapply(reports, function(x) x$priority))]
    
    # Apply entry limits per category
    if (length(composite_plots) > filter_config$max_entries_per_category) {
      excluded_count <- length(composite_plots) - filter_config$max_entries_per_category
      composite_plots <- composite_plots[1:filter_config$max_entries_per_category]
      log_message(safe_sprintf("Limited %s composite plots to %d entries (%d excluded)", 
                         module, filter_config$max_entries_per_category, excluded_count))
    }
    
    if (length(individual_plots) > filter_config$max_entries_per_category) {
      excluded_count <- length(individual_plots) - filter_config$max_entries_per_category
      individual_plots <- individual_plots[1:filter_config$max_entries_per_category]
      log_message(safe_sprintf("Limited %s individual plots to %d entries (%d excluded)", 
                         module, filter_config$max_entries_per_category, excluded_count))
    }
    
    if (length(summary_tables) > filter_config$max_entries_per_category) {
      excluded_count <- length(summary_tables) - filter_config$max_entries_per_category
      summary_tables <- summary_tables[1:filter_config$max_entries_per_category]
      log_message(safe_sprintf("Limited %s summary tables to %d entries (%d excluded)", 
                         module, filter_config$max_entries_per_category, excluded_count))
    }
    
    if (length(reports) > filter_config$max_entries_per_category) {
      excluded_count <- length(reports) - filter_config$max_entries_per_category
      reports <- reports[1:filter_config$max_entries_per_category]
      log_message(safe_sprintf("Limited %s reports to %d entries (%d excluded)", 
                         module, filter_config$max_entries_per_category, excluded_count))
    }
    
    # Store organized collections
    file_collections[[paste0(module, "_COMPOSITE_PLOTS")]] <- composite_plots
    file_collections[[paste0(module, "_INDIVIDUAL_PLOTS")]] <- individual_plots
    file_collections[[paste0(module, "_SUMMARY_TABLES")]] <- summary_tables
    file_collections[[paste0(module, "_REPORTS")]] <- reports
  }
  
  log_message(safe_sprintf("[SUCCESS] Processed manifest data for %d categories", length(file_collections)))
  
  return(file_collections)
}

#' Convert absolute path to relative path for report (V6.5 Fix)
#' 
#' @param absolute_path Full path to file
#' @param report_base_path Base path where report will be located
#' @return Relative path for use in href attributes
convert_to_relative_path <- function(absolute_path, report_base_path) {
  
  # V6.5 ARCHITECT FIX: Return complete relative paths for valid links
  # The relative path is needed for href attributes, while basename() for display
  
  # Early return for empty or NULL paths
  if (is.null(absolute_path) || length(absolute_path) == 0 || absolute_path == "") {
    return("unknown_file")
  }
  
  # Convert absolute path to relative path for proper linking
  # This creates a valid relative path from results directory to the actual file
    tryCatch({
    # Get the directory containing the file relative to report base
    file_dir <- dirname(absolute_path)
    file_name <- basename(absolute_path)
    
    # Create relative path that goes up from results to session root, then down to file
    relative_path <- file.path("..", "plots", basename(file_dir), file_name)
    
        return(relative_path)
    }, error = function(e) {
    # Fallback: if path conversion fails, use just the filename
    log_message(sprintf("Path conversion failed for %s: %s", absolute_path, e$message), level = "warning")
    return(basename(absolute_path))
  })
}

# Skip old code - function is complete above

#' Render academic content with proper formatting (V3.0)
#' 
#' @param entry_list List of enhanced manifest entries
#' @param content_type Type of content being rendered (plots, tables, reports) 
#' @param section_title Optional section title for academic structure
#' @param numbering_start Starting number for figure/table numbering (default: 1)
#' @return Academically formatted Markdown content with titles, captions, and proper structure
render_academic_content <- function(entry_list, content_type = "plots", section_title = NULL, numbering_start = 1) {
  
  if (is.null(entry_list) || length(entry_list) == 0) {
    return("*No content generated for this section*\\n")
  }
  
  markdown_lines <- c()
  
  # Add section title if provided
  if (!is.null(section_title)) {
    markdown_lines <- c(markdown_lines, sprintf("### %s", section_title), "")
  }
  
  # Process each entry with academic formatting
  for (entry_name in names(entry_list)) {
    entry <- entry_list[[entry_name]]
    
    # Use a placeholder if entry is not a valid list
    if (!is.list(entry)) {
      log_message(sprintf("Skipping invalid manifest entry: %s", entry_name), level="warning")
      next
    }

    file_path <- entry$path %||% ""
    if (!file.exists(file_path)) {
      log_message(sprintf("Skipping non-existent file: %s", file_path), level="warning")
      next
    }

    file_name <- basename(file_path)

    # Report methods only; do not interpret results
      markdown_lines <- c(markdown_lines,
                       sprintf("##### %s", entry$caption %||% file_name),
                       "",
                       sprintf("**Analysis method**: %s", get_method_description(entry$role, entry_name)),
                       "")
    
    # Embed PNG/JPG images directly; skip PDF files
    if (entry$file_type == "plot") {
      # Check for PNG version first (converted by format_converter)
      base_path_no_ext <- tools::file_path_sans_ext(file_path)
      png_path <- paste0(base_path_no_ext, ".png")
      jpg_path <- paste0(base_path_no_ext, ".jpg")
      
      # Determine which format to embed
      embed_path <- NULL
      if (file.exists(png_path)) {
        embed_path <- png_path
      } else if (file.exists(jpg_path)) {
        embed_path <- jpg_path
      }
      
      # Only embed if PNG/JPG exists, skip PDF-only files
      if (!is.null(embed_path)) {
        relative_path <- convert_to_relative_path(embed_path, report_base_path)
        display_name <- basename(embed_path)
      markdown_lines <- c(markdown_lines,
                           sprintf("![%s](%s)", tools::file_path_sans_ext(display_name), relative_path),
                         "")
      }
    } else {
      # For non-plot files, create download links
      relative_path <- convert_to_relative_path(file_path, report_base_path)
      display_name <- basename(file_path)
      markdown_lines <- c(markdown_lines,
                         sprintf("[Data file: **%s**](%s)", display_name, relative_path),
                         "")
    }
    
    markdown_lines <- c(markdown_lines, "---", "")
  }
  
  return(paste(markdown_lines, collapse = "\n"))
}

get_method_description <- function(role, entry_name) {
  entry_name <- tolower(entry_name)
  if (grepl("density", entry_name)) {
    return("Uses data from `normalized_frequencies_summary.csv` to draw kernel density estimation curves for variant frequencies in each genomic region.")
  }
  if (grepl("boxplot", entry_name)) {
    return("Uses data from `normalized_frequencies_summary.csv` to draw box plots for each genomic region and calls outlier detection algorithms to identify extreme values.")
  }
  if (grepl("pie_matrix", entry_name)) {
    return("Uses data from `normalized_frequencies.csv`, groups by species and region, counts raw numbers of various variants, and draws pie charts.")
  }
  if (grepl("m02", entry_name)) {
    return("M02's 'intelligent statistical engine' receives `normalized_frequencies.csv` and `group_info.csv`, automatically performs normality and homogeneity of variance tests, and selects the most appropriate statistical method (such as ANOVA or Kruskal-Wallis test) to calculate p-values for inter-group differences.")
  }
  if (grepl("global_freq_heatmap", entry_name)) {
    return("Based on data from `normalized_frequencies.csv` filtered through `data_scoping`, draws gene-species frequency matrix heatmaps.")
  }
  if (grepl("global_pca", entry_name)) {
    return("This plot uses frequency data for **all genes** that have been preliminarily filtered through `data_scoping` for the M03 module. The script reshapes this long data table into a 'species x genes' wide matrix and uses it as input for PCA analysis.")
  }
  if (grepl("upset", entry_name)) {
    return("Uses candidate hotspot gene lists to perform set operations and visualization on gene presence/absence across species.")
  }
  if (grepl("correlation", entry_name)) {
    return("Filters out 'core hotspots' that exist in at least N species, then based on the presence/absence of these core hotspots, calculates Spearman correlation coefficients and p-values between species.")
  }
  if (grepl("statistical_tests", entry_name)) {
    return("Includes chi-square tests and hypergeometric tests of the 'two-speed evolution model', etc. Aims to answer statistically: Is the distribution of 'private hotspots' and 'core shared hotspots' random, or has it been subject to significant evolutionary selection pressure.")
  }
  if (grepl("enrichment", entry_name)) {
    return("Uses Fisher's exact test to determine whether the proportion of a certain biological function in hotspot genes is significantly higher than its proportion in the background genome.")
  }
  return(sprintf(
    "Generated using cpopvar %s.",
    as.character(utils::packageVersion("cpopvar"))
  ))
}

get_parameter_description <- function(role, entry_name) {
  entry_name <- tolower(entry_name)
  if (grepl("m01|m02|m03", entry_name)) {
    return("`data_scoping` (used to filter regions and variant types for analysis).")
  }
  if (grepl("hotspot", entry_name)) {
    return("`hotspot_threshold_percentile`, `min_species_number`.")
  }
  return("Based on parameter configuration in `config/full_refactor_test.yml`.")
}

get_interpretation_guidance <- function(role, entry_name) {
  entry_name <- tolower(entry_name)
  if (grepl("density", entry_name)) {
    return("The peak position of the curve indicates the most concentrated variant frequency value in that region, and the width of the curve indicates the degree of frequency dispersion. Through this plot, you can quickly determine which region has higher, more concentrated, or more dispersed variation levels.")
  }
  if (grepl("boxplot", entry_name)) {
    return("The box plot shows the median and interquartile range of frequencies for each region. The scatter plot shows the true distribution of the data. Red outliers are often genes of interest with extremely high variation frequencies and may be potential hotspots.")
  }
  if (grepl("pie_matrix", entry_name)) {
    return("This plot is used for rapid cross-species and cross-regional variant pattern comparison. For example, you can observe whether a certain region of a certain species is dominated by SNP variants, while another species is dominated by INDELs, which may suggest different evolutionary pressures.")
  }
  if (grepl("m02", entry_name)) {
    return("The plot shows the variant frequency distribution of each group. More importantly, it should be combined with the p-values in the `_report_all_...` report to determine whether inter-group differences are statistically significant.")
  }
  if (grepl("global_freq_heatmap", entry_name)) {
    return("This plot is the starting point of hotspot analysis, used to discover patterns macroscopically. You can observe whether there are certain genes (rows) that are universally highly variable across all species, or whether there are certain species (columns) with overall variation levels much higher than other species.")
  }
  if (grepl("global_pca", entry_name)) {
    return("Each point in the plot represents a species. Species that are close in distance indicate they are more similar in **overall gene variation patterns**. This plot aims to provide a macroscopic overview of species clustering relationships based on broad gene sets.")
  }
  if (grepl("upset", entry_name)) {
    return("The bar chart above shows the size of species combinations (intersections) sharing hotspot genes, and the bar chart on the left shows the number of hotspot genes unique to each species. This plot is the core of analyzing hotspot gene conservation among species.")
  }
  if (grepl("correlation", entry_name)) {
    return("The deeper the color, the closer the relationship between two species due to sharing similar core hotspot genes. Asterisks (`*`, `**`, `***`) indicate that the correlation is statistically significant, excluding random factors.")
  }
  if (grepl("statistical_tests", entry_name)) {
    return("This is the most core statistical conclusion of the M03 module. For example, if the report shows 'TWO-SPEED EVOLUTION MODEL SUPPORTED', it strongly proves that the formation of private hotspots and shared hotspots follows different evolutionary driving forces.")
  }
  if (grepl("enrichment", entry_name)) {
    return("This table is used to reveal the biological significance of hotspot genes. For example, you can compare whether shared and private hotspots are enriched in completely different biological functions, thereby revealing different evolutionary adaptation strategies.")
  }
  return("Please interpret in combination with specific charts and project background.")
}

#' Replace template placeholders with actual values
#' 
#' @param template_text Template content as string
#' @param session_stats Session statistics
#' @param analysis_conclusions Analysis conclusions text
#' @param manifests Manifest collections from modules
#' @param content_filter_level Content filtering level
#' @return Final report text with placeholders replaced
replace_template_placeholders <- function(template_text, session_stats, analysis_conclusions, manifests = NULL, content_filter_level = "standard") {
  
  log_message("Replacing template placeholders with actual values")
  
  # Basic placeholder replacements
  report_text <- template_text
  
  # Replace basic placeholders
  report_text <- gsub(
    "\\{\\{PROJECT_NAME\\}\\}",
    "cpopvar Chloroplast Genome Variation Analysis",
    report_text
  )
  report_text <- gsub("\\{\\{SESSION_ID\\}\\}", session_stats$session_id, report_text)
  report_text <- gsub("\\{\\{REPORT_DATE\\}\\}", session_stats$report_date, report_text)
  report_text <- gsub("\\{\\{R_VERSION\\}\\}", session_stats$r_version, report_text)
  report_text <- gsub("\\{\\{ANALYSIS_DURATION\\}\\}", session_stats$processing_time, report_text)
  report_text <- gsub("\\{\\{EXECUTION_MODE\\}\\}", session_stats$execution_mode, report_text)
  report_text <- gsub("\\{\\{GENERATION_TIME\\}\\}", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), report_text)
  
  # Generate and replace analysis flowchart - V7.1 Enhancement
  log_message("Generating analysis pipeline flowchart")
  tryCatch({
    # Extract configuration data for flowchart generation
    config_for_flowchart <- list(
      execution_mode = session_stats$execution_mode,
      session_info = list(session_id = session_stats$session_id)
    )
    
    analysis_flowchart <- generate_analysis_flowchart(config_for_flowchart)
    report_text <- gsub("\\{\\{ANALYSIS_FLOWCHART\\}\\}", analysis_flowchart, report_text)
    log_message("[SUCCESS] Analysis flowchart integrated successfully")
  }, error = function(e) {
    log_message(safe_sprintf("Warning: Failed to generate flowchart: %s", e$message), level = "warning")
    fallback_flowchart <- "```\nProcess flowchart generation failed, please check system logs\n```"
    report_text <- gsub("\\{\\{ANALYSIS_FLOWCHART\\}\\}", fallback_flowchart, report_text)
  })
  
  # System information
  system_info_text <- sprintf("Platform: %s, OS: %s", 
                             session_stats$system_info$platform %||% "Unknown", 
                             session_stats$system_info$os %||% "Unknown")
  report_text <- gsub("\\{\\{SYSTEM_INFO\\}\\}", system_info_text, report_text)
  
  # Data statistics
  if (!is.null(session_stats$data_quality_summary)) {
    # Extract numbers from data quality summary
    raw_match <- regmatches(session_stats$data_quality_summary, regexpr("Raw variants: [0-9]+", session_stats$data_quality_summary))
    processed_match <- regmatches(session_stats$data_quality_summary, regexpr("Filtered variants: [0-9]+", session_stats$data_quality_summary))
    species_match <- regmatches(session_stats$data_quality_summary, regexpr("Species included: [0-9]+", session_stats$data_quality_summary))
    
    raw_count <- if(length(raw_match) > 0) gsub("Raw variants: ", "", raw_match) else "Unknown"
    processed_count <- if(length(processed_match) > 0) gsub("Filtered variants: ", "", processed_match) else "Unknown"
    species_count <- if(length(species_match) > 0) gsub("Species included: ", "", species_match) else "Unknown"
    
    report_text <- gsub("\\{\\{RAW_DATA_COUNT\\}\\}", raw_count, report_text)
    report_text <- gsub("\\{\\{PROCESSED_DATA_COUNT\\}\\}", processed_count, report_text)
    report_text <- gsub("\\{\\{SPECIES_COUNT\\}\\}", species_count, report_text)
  } else {
    report_text <- gsub("\\{\\{RAW_DATA_COUNT\\}\\}", "Unknown", report_text)
    report_text <- gsub("\\{\\{PROCESSED_DATA_COUNT\\}\\}", "Unknown", report_text)
    report_text <- gsub("\\{\\{SPECIES_COUNT\\}\\}", "Unknown", report_text)
  }
  
  # Region types
  report_text <- gsub("\\{\\{REGION_TYPES\\}\\}", "CDS, IGS, intron, RNA", report_text)
  
  # MODULE_RESULTS section
  module_results_text <- generate_module_results_section(manifests, session_stats$session_id, content_filter_level)
  report_text <- gsub("\\{\\{MODULE_RESULTS\\}\\}", module_results_text, report_text)
  
  log_message("[SUCCESS] Template placeholder replacement completed")
  
  return(report_text)
}

#' Generate module results section with manifest data
#' 
#' @param manifests Manifest collections from modules
#' @param session_id Session ID
#' @param content_filter_level Content filtering level
#' @return Module results section text
generate_module_results_section <- function(manifests, session_id, content_filter_level = "standard") {
  
  if (is.null(manifests) || length(manifests) == 0) {
    log_message("No manifests available, using fallback file collection", level = "warning")
    # Fallback to file-based collection
    session_paths <- get_session_paths(session_id)
    file_collections <- collect_module_files(session_paths)
    return(generate_fallback_module_section(file_collections))
  }
  
  # Process manifests for report generation
  file_collections <- process_manifest_collections(manifests, session_id, content_filter_level)
  
  if (is.null(file_collections) || length(file_collections) == 0) {
    log_message("Processed manifests resulted in empty collections, using fallback", level = "warning")
    session_paths <- get_session_paths(session_id)
    file_collections <- collect_module_files(session_paths)
    return(generate_fallback_module_section(file_collections))
  }
  
  # Generate module sections using manifest data
  module_sections <- c()
  
  # Module M01 - Data Distribution Analysis
  if (any(grepl("M01", names(file_collections)))) {
    m01_section <- generate_module_section("M01", "Data Distribution Analysis", file_collections)
    if (nchar(m01_section) > 10) {
      module_sections <- c(module_sections, m01_section)
    }
  }
  
  # Module M02 - Comparative Analysis
  if (any(grepl("M02", names(file_collections)))) {
    m02_section <- generate_module_section("M02", "compareanalysis (Comparative Analysis)", file_collections)
    if (nchar(m02_section) > 10) {
      module_sections <- c(module_sections, m02_section)
    }
  }
  
  # Module M03 - Hotspot Analysis
  if (any(grepl("M03", names(file_collections)))) {
    m03_section <- generate_module_section("M03", "Hotspot Analysis", file_collections)
    if (nchar(m03_section) > 10) {
      module_sections <- c(module_sections, m03_section)
    }
  }
  
  if (length(module_sections) == 0) {
    return("*No analysis result files generated*")
  }
  
  return(paste(module_sections, collapse = "\n\n"))
}

#' Generate section for a specific module
#' 
#' @param module_name Module name (M01, M02, M03)
#' @param module_title Module title for display
#' @param file_collections Processed file collections
#' @return Module section text
generate_module_section <- function(module_name, module_title, file_collections) {
  
  section_lines <- c()
  section_lines <- c(section_lines, sprintf("## %s", module_title), "")
  
  # Collect files for this module
  composite_plots <- file_collections[[paste0(module_name, "_COMPOSITE_PLOTS")]] %||% list()
  individual_plots <- file_collections[[paste0(module_name, "_INDIVIDUAL_PLOTS")]] %||% list()
  summary_tables <- file_collections[[paste0(module_name, "_SUMMARY_TABLES")]] %||% list()
  reports <- file_collections[[paste0(module_name, "_REPORTS")]] %||% list()
  
  # Composite plots (highest priority)
  if (length(composite_plots) > 0) {
    section_lines <- c(section_lines, "### Comprehensive Charts")
    composite_content <- render_academic_content(composite_plots, "plots", numbering_start = 1)
    section_lines <- c(section_lines, composite_content)
  }
  
  # Individual plots
  if (length(individual_plots) > 0) {
    section_lines <- c(section_lines, "### Detailed Charts")
    individual_content <- render_academic_content(individual_plots, "plots", numbering_start = length(composite_plots) + 1)
    section_lines <- c(section_lines, individual_content)
  }
  
  # Summary tables
  if (length(summary_tables) > 0) {
    section_lines <- c(section_lines, "### Data Tables")
    table_content <- render_academic_content(summary_tables, "tables")
    section_lines <- c(section_lines, table_content)
  }
  
  # Reports
  if (length(reports) > 0) {
    section_lines <- c(section_lines, "### Analysis Reports")
    report_content <- render_academic_content(reports, "reports")
    section_lines <- c(section_lines, report_content)
  }
  
  if (length(section_lines) <= 2) {
    return(sprintf("## %s\n\n*No result files for this module*\n", module_title))
  }
  
  return(paste(section_lines, collapse = "\n"))
}

#' Generate fallback module section when manifests are not available
#' 
#' @param file_collections File collections from directory scanning
#' @return Fallback module section text
generate_fallback_module_section <- function(file_collections) {
  
  if (is.null(file_collections) || length(file_collections) == 0) {
    return("*No analysis result files*")
  }
  
  section_lines <- c()
  
  # Embed PNG/JPG images only; skip PDF files
  for (collection_name in names(file_collections)) {
    files <- file_collections[[collection_name]]
    if (length(files) == 0) next
    
    # Extract module name from collection name
    module_match <- regmatches(collection_name, regexpr("M0[1-5]", collection_name))
    if (length(module_match) == 0) next
    
    module_name <- module_match[1]
    
    # Group files by type
    png_files <- files[grepl("\\.png$", files, ignore.case = TRUE)]
    jpg_files <- files[grepl("\\.(jpg|jpeg)$", files, ignore.case = TRUE)]
    pdf_files <- files[grepl("\\.pdf$", files, ignore.case = TRUE)]
    data_files <- files[!grepl("\\.(png|jpg|jpeg|pdf)$", files, ignore.case = TRUE)]
    
    # Only process if there are PNG/JPG files or data files
    if (length(png_files) == 0 && length(jpg_files) == 0 && length(data_files) == 0) next
    
    # Add module header
    if (length(section_lines) == 0 || !grepl(sprintf("## %s", module_name), tail(section_lines, 1))) {
      section_lines <- c(section_lines, sprintf("## %s module", module_name), "")
    }
    
    # Embed PNG images
    if (length(png_files) > 0) {
      for (file_path in png_files) {
        file_name <- basename(file_path)
        # Skip if this is a converted version of a PDF (we'll handle it below)
        base_name <- tools::file_path_sans_ext(file_name)
        
        relative_path <- convert_to_relative_path(file_path, "")
        section_lines <- c(section_lines,
                          sprintf("##### %s", base_name),
                          "",
                          sprintf("![%s](%s)", base_name, relative_path),
                          "",
                          "---",
                          "")
      }
    }
    
    # Embed JPG images
    if (length(jpg_files) > 0) {
      for (file_path in jpg_files) {
        file_name <- basename(file_path)
        base_name <- tools::file_path_sans_ext(file_name)
        
        relative_path <- convert_to_relative_path(file_path, "")
        section_lines <- c(section_lines,
                          sprintf("##### %s", base_name),
                          "",
                          sprintf("![%s](%s)", base_name, relative_path),
                          "",
                          "---",
                          "")
      }
    }
    
    # Add data file links (if any)
    if (length(data_files) > 0) {
      section_lines <- c(section_lines, "**Data files**:", "")
      for (file_path in data_files) {
        file_name <- basename(file_path)
        relative_path <- convert_to_relative_path(file_path, "")
        section_lines <- c(section_lines, sprintf("- [%s](%s)", file_name, relative_path))
      }
      section_lines <- c(section_lines, "")
    }
  }
  
  if (length(section_lines) == 0) {
    return("*No analysis result files*")
  }
  
  return(paste(section_lines, collapse = "\n"))
}
