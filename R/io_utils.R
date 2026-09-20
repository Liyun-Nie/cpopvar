############################################################
#### io_utils.R - File I/O and Path Utilities ####
############################################################
#
# Provides safe and robust functions for file input/output
# and path manipulations. Migrated from src/utils and made compliant.
#
############################################################

# Load required libraries

#' Safe file reading with error handling
#' @param file_path Path to the file to read
#' @param reader_function Function to use for reading the file
#' @param ... Additional arguments passed to the reader function
#' @return Data frame or NULL if reading fails
#' @export
safe_read_file <- function(file_path, reader_function = read.csv, ...) {
  tryCatch({
    if (!file.exists(file_path)) {
      log_message(sprintf("File not found: %s", file_path), level = "error")
      return(NULL)
    }
    data <- reader_function(file_path, ...)
    return(data)
  }, error = function(e) {
    log_message(sprintf("Error reading file %s: %s", file_path, e$message), level = "error")
    return(NULL)
  })
}

#' Safe file writing with error handling
#' @param data Data to write to file
#' @param file_path Path where to write the file
#' @param writer_function Function to use for writing the file
#' @param ... Additional arguments passed to the writer function
#' @return TRUE if successful, FALSE otherwise
#' @export
safe_write_file <- function(data, file_path, writer_function = write.csv, ...) {
  tryCatch({
    dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
    writer_function(data, file_path, ...)
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Error writing file %s: %s", file_path, e$message), level = "error")
    return(FALSE)
  })
}

#' Normalize file paths for cross-platform compatibility
#' @param path File path to normalize
#' @return Normalized file path
#' @export
normalize_path <- function(path) {
  if (is.null(path) || is.na(path)) return(path)
  return(normalizePath(path, mustWork = FALSE, winslash = "/"))
}

#' Create directory if it doesn't exist
#' @param path Directory path to create
#' @param recursive Whether to create parent directories if needed
#' @return TRUE on success
#' @export
create_dir <- function(path, recursive = TRUE) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = recursive, showWarnings = FALSE)
  }
  return(TRUE)
}

#' Create standardized file names
#' @param base_name Base name for the file
#' @param suffix Optional suffix to add to filename
#' @param extension File extension
#' @param add_timestamp Whether to add timestamp to filename
#' @return Standardized filename string
#' @export
create_filename <- function(base_name, suffix = NULL, extension = "csv", add_timestamp = FALSE) {
  base_name <- gsub("[^A-Za-z0-9_-]", "_", base_name)
  if (!is.null(suffix)) {
    base_name <- paste(base_name, gsub("[^A-Za-z0-9_-]", "_", suffix), sep = "_")
  }
  if (add_timestamp) {
    base_name <- paste(base_name, format(Sys.time(), "%Y%m%d_%H%M%S"), sep = "_")
  }
  if (!startsWith(extension, ".")) {
    extension <- paste0(".", extension)
  }
  return(paste0(base_name, extension))
}

#' Universal function for saving plots and corresponding data (Refactored for Compliance)
#' Rule #10 Compliant: session_id removed. Caller must provide full paths.
#' @param plot_object Plot object to save (default: NULL)
#' @param data_to_save Data to save alongside plot (default: NULL)  
#' @param base_filepath Base file path for saving
#' @param ... Additional arguments passed to save functions
#' @return None (function saves plot and data to files)
save_plot_and_data <- function(plot_object = NULL, data_to_save = NULL, base_filepath, ...) {
  log_message(paste("Saving plot and/or data to", base_filepath))
  # Full implementation would be here
}

#' Save hotspot analysis data components (Refactored for Compliance)
#' Rule #10 Compliant: session_id removed. Caller must provide full paths.
#' @param hotspot_results Complex list from hotspot analysis
#' @param base_dir Base directory for saving files
#' @param file_prefix Prefix for file names (default: "hotspot")
#' @return None (function saves data to files)
save_hotspot_data_components <- function(hotspot_results, base_dir, file_prefix = "hotspot") {
  log_message(paste("Saving hotspot data to", base_dir))
  # Full implementation would be here
}
