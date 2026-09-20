############################################################
#### path_manager.R - Centralized Path Management System ####
############################################################
#
# Unified path configuration and management for the entire application
# Replaces scattered path definitions with centralized system
# Provides both session-aware and global path management
#
############################################################

# Initialize global PATH_CONFIG when this module is loaded
PATH_CONFIG <- NULL

# Load required libraries (here package no longer needed in R package context)

#' Get session-specific path structure
#' 
#' @param session_id Unique session identifier
#' @param create_dirs Whether to create directories
#' @param output_dir Base directory for output. Defaults to current working directory.
#' @return List of session-specific paths
#' @export
get_session_paths <- function(session_id, create_dirs = TRUE, output_dir = NULL) {
  
  if (is.null(output_dir)) {
    output_dir <- getOption("cpopvar.output_dir", NULL)
  }

  # Determine base directory
  base_dir <- if (!is.null(output_dir)) {
    # If user specified output directory, use it
    normalizePath(output_dir, mustWork = FALSE)
  } else {
    # Otherwise, use current working directory
    getwd()
  }
  
  # Validate directory existence and create if necessary
  if (!dir.exists(base_dir)) {
    tryCatch({
      dir.create(base_dir, recursive = TRUE)
    }, warning = function(w) {
      stop(paste("Output directory does not exist and cannot be created:", base_dir, "Error:", w$message))
    }, error = function(e) {
      stop(paste("Output directory does not exist and cannot be created:", base_dir, "Error:", e$message))
    })
  }
  
  # Use app_data structure for session storage under the base directory
  base_session_dir <- file.path(base_dir, "app_data", "sessions", session_id)
  
  session_paths <- list(
    # Base session directory
    base = base_session_dir,
    
    # Input data storage (all user uploads)
    raw = file.path(base_session_dir, "raw"),
    
    # Final results and outputs
    results = file.path(base_session_dir, "results"),
    
    # Processed data (intermediate results) - positioned under results as per golden standard
    processed_data = file.path(base_session_dir, "results", "processed_data"),

    # Named processing-stage paths retained as part of the reporting contract
    P01_preprocessed_out = file.path(base_session_dir, "results", "processed_data", "P01_preprocessed"),
    P02_filtered_out = file.path(base_session_dir, "results", "processed_data", "P02_filtered"),
    P03_normalized_out = file.path(base_session_dir, "results", "processed_data", "P03_normalized"),
    
    # Generated plots and visualizations
    plots = file.path(base_session_dir, "results", "plots"),
    
    # Session-specific logs
    logs = file.path(base_session_dir, "logs"),
    
    # Annotation files (subfolder of raw)
    annotations = file.path(base_session_dir, "raw", "annotations"),
    
    # Region statistics files (subfolder of raw)
    region_stats = file.path(base_session_dir, "raw", "region_stats"),
    
    # Exported results
    exports = file.path(base_session_dir, "exports")
  )
  
  # Create directories if requested
  if (create_dirs) {
    for (path in session_paths) {
      if (!dir.exists(path)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
      }
    }
  }
  
  return(session_paths)
}

#' Get path for a specific processing stage
#'
#' This function now supports a testing override via an environment variable.
#'
#' @param session_id Session ID
#' @param stage_name Name of the processing stage (e.g., "P01_preprocessed_out")
#' @param create_dir Whether to create the directory if it doesn't exist
#' @param output_dir Base directory for output. Defaults to current working directory.
#' @return Path to the processing stage directory
#' @export
get_processing_stage_path <- function(session_id, stage_name, create_dir = FALSE, output_dir = NULL) {

  # --- NEW: Testing Override Logic ---
  # Check if a special testing environment variable is set.
  test_data_path_override <- Sys.getenv("CPOP_TEST_DATA_PATH")

  if (nzchar(test_data_path_override)) {
    # If the variable is set, we are in "test mode".
    # Ignore the session_id and use the provided path directly.
    base_path <- file.path(test_data_path_override, stage_name)
  } else {
    # --- Fixed Logic: Use unified path management system ---
    # Call get_session_paths() to get standardized paths
    session_paths <- get_session_paths(session_id, create_dirs = FALSE, output_dir = output_dir)
    # Append stage_name to the correct processed_data base path
    base_path <- file.path(session_paths$processed_data, stage_name)
  }

  if (create_dir && !dir.exists(base_path)) {
    dir.create(base_path, recursive = TRUE, showWarnings = FALSE)
  }

  return(base_path)
}
