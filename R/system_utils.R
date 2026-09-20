############################################################
#### system_utils.R - System and Logging Utilities ####
############################################################
#
# Provides utilities for logging, session management, and
# environment control. Migrated from src/utils and made compliant.
#
############################################################

# Note: Using null coalescing operator (%||%) exported from annotation_processing.R

#' Setup logging for the analysis
#' @param log_file Path to log file (if NULL, auto-generated)
#' @param level Logging level (default: "INFO")
#' @param session_id Session identifier for session-specific logging
#' @param output_dir Base directory for output. Defaults to current working directory.
setup_logging <- function(log_file = NULL, level = "INFO", session_id = NULL, output_dir = NULL) {
  if (is.null(log_file)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    
    # Determine base directory
    base_dir <- if (!is.null(output_dir)) {
      normalizePath(output_dir, mustWork = FALSE)
    } else {
      getwd()
    }
    
    if (!is.null(session_id)) {
      # Session-specific logging
      session_log_dir <- file.path(base_dir, "app_data", "sessions", session_id, "logs")
      log_file <- file.path(session_log_dir, paste0("session_", timestamp, ".log"))
    } else {
      # Global logging
      log_file <- file.path(base_dir, "app_data", "global_logs", paste0("shiny_workflow_", timestamp, ".log"))
    }
  }
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  options(shiny_workflow_log_file = log_file)
}

#' Log messages with timestamp and session support
#' @param message The message to log
#' @param level Log level (default: "INFO")
#' @param to_console Whether to print to console (default: TRUE)
#' @param session_id Session identifier for session-specific logging
#' @export
log_message <- function(message, level = "INFO", to_console = TRUE, session_id = NULL) {
  log_file <- getOption("shiny_workflow_log_file")
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] %s: %s", timestamp, level, message)
  if (!is.null(session_id)) {
      log_entry <- sprintf("[%s] [%s] %s: %s", timestamp, session_id, level, message)
  }
  if (to_console) cat(log_entry, "\n")
  if (!is.null(log_file) && file.exists(dirname(log_file))) {
      write(log_entry, file = log_file, append = TRUE)
  }
}

#' Progress callback function for Shiny
#' @param message Progress message to display
#' @param value Progress value (0-1 or NULL)
#' @param detail Additional progress detail text
progress_callback <- function(message, value = NULL, detail = NULL) {
  log_message(sprintf("Progress: %s", message))
}

#' Memory cleanup function
#' 
#' @title Memory cleanup function
#' @description Performs garbage collection to free up memory
#' @param verbose Logical, whether to print verbose output (default: FALSE)
#' @return Invisible garbage collection results
#' @export
cleanup_memory <- function(verbose = FALSE) {
  gc(verbose = verbose)
}

#' Check required packages
#' 
#' @title Check required packages
#' @description Checks for missing packages without modifying the user's library.
#' @param packages Character vector of package names to check
#' @param quiet Logical, whether to suppress the success message (default: TRUE)
#' @return Invisibly returns TRUE when all packages are available.
#' @export
check_and_install_packages <- function(packages, quiet = TRUE) {
  missing_packages <- packages[!packages %in% installed.packages()[, "Package"]]
  if (length(missing_packages) > 0) {
    stop(
      sprintf(
        "Missing required packages: %s. Install them before running this feature.",
        paste(missing_packages, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!quiet) {
    message("All required packages are available.")
  }
  invisible(TRUE)
}

#' Initialize session-specific logging with unified path management
#' 
#' @param session_id Session identifier
#' @param output_dir Base directory for output. Defaults to current working directory.
#' @return Path to session log file
#' @export
initialize_session_logging <- function(session_id, output_dir = NULL) {
  if (is.null(session_id) || session_id == "") {
    warning("Session ID is NULL or empty, using global logging")
    return(NULL)
  }
  
  # Create session log directory using unified path management
  if (exists("get_session_paths") && is.function(get_session_paths)) {
    tryCatch({
      session_paths <- get_session_paths(session_id, create_dirs = TRUE, output_dir = output_dir)
      session_log_dir <- session_paths$logs
    }, error = function(e) {
      # Fallback to path structure with dynamic output directory
      base_dir <- if (!is.null(output_dir)) {
        normalizePath(output_dir, mustWork = FALSE)
      } else {
        getwd()
      }
      session_log_dir <<- file.path(base_dir, "app_data", "sessions", session_id, "logs")
      dir.create(session_log_dir, recursive = TRUE, showWarnings = FALSE)
    })
  } else {
    # Fallback to path structure with dynamic output directory
    base_dir <- if (!is.null(output_dir)) {
      normalizePath(output_dir, mustWork = FALSE)
    } else {
      getwd()
    }
    session_log_dir <- file.path(base_dir, "app_data", "sessions", session_id, "logs")
    dir.create(session_log_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Generate session log file
  timestamp_file <- format(Sys.time(), "%Y-%m-%d")
  session_log_file <- file.path(session_log_dir, paste0("session_", session_id, "_", timestamp_file, ".log"))
  
  # Log session initialization
  log_message("=== SESSION LOGGING INITIALIZED ===", session_id = session_id)
  log_message(paste("Session ID:", session_id), session_id = session_id)
  log_message(paste("Log file:", session_log_file), session_id = session_id)
  log_message(paste("Timestamp:", Sys.time()), session_id = session_id)
  
  return(session_log_file)
}
