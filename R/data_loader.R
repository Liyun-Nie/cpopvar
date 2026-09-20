############################################################
#### data_loader.R - Unified Configuration-Driven Data Loading System ####
############################################################
#
# This module provides a centralized, configuration-driven data loading
# and column standardization system. It eliminates hardcoded column names
# throughout the application by applying YAML-defined column mappings
# at the data loading stage.
#
# Key Features:
# - Configuration-driven column mapping from user names to internal standard names
# - Unified data loading interface for all file types
# - Comprehensive validation and error handling
# - Support for multiple file formats (CSV, TSV, etc.)
# - Intelligent fallback and suggestion mechanisms
#
############################################################

# Load required libraries

# Source required utilities (if not already loaded)
# source("src/utils/common_utils.R")

#' Load and standardize data using configuration-driven column mapping
#' 
#' This is the core function for unified data loading and standardization.
#' It reads data from file, applies column mappings from configuration,
#' and returns a data frame with standardized internal column names.
#' 
#' @param file_path Path to the data file to load
#' @param column_mapping Named list mapping internal standard names to user column names
#' @param required_columns Vector of required internal column names (optional)
#' @param file_description Description of the file type for error messages
#' @param validate_data Whether to perform data validation (default: TRUE)
#' @param allow_missing_columns Whether to allow missing non-required columns (default: TRUE)
#' @return Data frame with standardized column names or NULL on failure
#' load_and_standardize_data
#' @export
load_and_standardize_data <- function(file_path, 
                                     column_mapping, 
                                     required_columns = NULL,
                                     file_description = "data file",
                                     validate_data = TRUE,
                                     allow_missing_columns = TRUE) {
  
  tryCatch({
    log_message(sprintf("Loading and standardizing %s: %s", file_description, file_path))
    
    # Validate inputs
    if (is.null(file_path) || !file.exists(file_path)) {
      log_message(sprintf("File not found: %s", file_path), level = "error")
      return(data.frame())  # Return empty data frame instead of NULL
    }
    
    if (is.null(column_mapping) || length(column_mapping) == 0) {
      log_message("No column mapping provided", level = "error")
      return(data.frame())  # Return empty data frame instead of NULL
    }
    
    # Read data file
    raw_data <- read_data_file(file_path)
    if (is.null(raw_data) || nrow(raw_data) == 0) {
      log_message(sprintf("Failed to read data or empty file: %s", file_path), level = "error")
      return(data.frame())  # Return empty data frame instead of NULL
    }
    
    log_message(sprintf("Read %d rows and %d columns from %s", 
                       nrow(raw_data), ncol(raw_data), file_description))
    
    # Apply column mapping (using new version from column_mapping.R)
    standardized_data <- apply_column_mapping(raw_data, column_mapping, file_type = "main_data")
    
    if (is.null(standardized_data)) {
      log_message("Column mapping failed", level = "error")
      return(data.frame())  # Return empty data frame instead of NULL
    }
    
    # Validate required columns if specified
    if (!is.null(required_columns) && validate_data) {
      missing_required <- setdiff(required_columns, names(standardized_data))
      if (length(missing_required) > 0) {
        log_message(sprintf("Missing required columns after mapping: %s", 
                           paste(missing_required, collapse = ", ")), level = "error")
        return(data.frame())  # Return empty data frame instead of NULL
      }
    }
    
    # Data quality validation
    if (validate_data) {
      validation_result <- validate_standardized_data(standardized_data, file_description)
      if (!validation_result$valid) {
        log_message(sprintf("Data validation failed: %s", validation_result$message), level = "warning")
        # Continue processing but log the warning
      }
    }
    
    log_message(sprintf("Successfully standardized %s: %d rows, %d columns with standard names", 
                       file_description, nrow(standardized_data), ncol(standardized_data)))
    
    return(standardized_data)
    
  }, error = function(e) {
    log_message(sprintf("Failed to load and standardize %s: %s", file_description, e$message), level = "error")
    return(data.frame())  # Return empty data frame instead of NULL (defensive programming)
  })
}

#' Read data file with automatic format detection
#' 
#' @param file_path Path to the data file
#' @return Data frame or NULL on failure
read_data_file <- function(file_path) {
  
  tryCatch({
    file_ext <- tolower(tools::file_ext(file_path))
    
    if (file_ext %in% c("csv")) {
      # Try comma-separated first
      data <- readr::read_csv(file_path, show_col_types = FALSE, locale = readr::locale(encoding = "UTF-8"))
    } else if (file_ext %in% c("tsv", "txt")) {
      # Try tab-separated
      data <- readr::read_tsv(file_path, show_col_types = FALSE, locale = readr::locale(encoding = "UTF-8"))
    } else {
      # Fallback: try CSV format
      log_message(sprintf("Unknown file extension '%s', trying CSV format", file_ext), level = "warning")
      data <- readr::read_csv(file_path, show_col_types = FALSE, locale = readr::locale(encoding = "UTF-8"))
    }
    
    # Convert to regular data frame
    return(as.data.frame(data))
    
  }, error = function(e) {
    log_message(sprintf("Failed to read file %s: %s", file_path, e$message), level = "error")
    return(data.frame())  # Return empty data frame instead of NULL (defensive programming)
  })
}

# apply_column_mapping function removed - using version from column_mapping.R

#' Validate standardized data quality
#' 
#' @param data Standardized data frame
#' @param file_description Description for messages
#' @return List with valid (boolean) and message (string)
validate_standardized_data <- function(data, file_description = "data") {
  
  validation_messages <- character()
  
  # Check for completely empty data
  if (nrow(data) == 0) {
    return(list(valid = FALSE, message = "Data is empty"))
  }
  
  # Check for excessive missing values
  total_cells <- nrow(data) * ncol(data)
  missing_cells <- sum(is.na(data))
  missing_percentage <- (missing_cells / total_cells) * 100
  
  if (missing_percentage > 90) {
    validation_messages <- c(validation_messages, 
                           sprintf("High missing data: %.1f%% of cells are NA", missing_percentage))
  } else if (missing_percentage > 50) {
    validation_messages <- c(validation_messages, 
                           sprintf("Moderate missing data: %.1f%% of cells are NA", missing_percentage))
  }
  
  # Check for duplicate columns
  if (any(duplicated(names(data)))) {
    duplicate_cols <- names(data)[duplicated(names(data))]
    validation_messages <- c(validation_messages, 
                           sprintf("Duplicate column names: %s", paste(duplicate_cols, collapse = ", ")))
  }
  
  # Check for suspicious column types (all character when numbers expected)
  char_cols <- names(data)[sapply(data, is.character)]
  if (length(char_cols) > 0) {
    numeric_looking <- char_cols[sapply(char_cols, function(col) {
      values <- data[[col]][!is.na(data[[col]])]
      if (length(values) == 0) return(FALSE)
      # Check if most values look numeric
      numeric_count <- sum(grepl("^-?\\d*\\.?\\d+$", values))
      return(numeric_count / length(values) > 0.8)
    })]
    
    if (length(numeric_looking) > 0) {
      validation_messages <- c(validation_messages, 
                             sprintf("Possibly numeric columns stored as character: %s", 
                                   paste(numeric_looking, collapse = ", ")))
    }
  }
  
  # Summary message
  if (length(validation_messages) == 0) {
    return(list(valid = TRUE, message = "Data validation passed"))
  } else {
    return(list(valid = FALSE, message = paste(validation_messages, collapse = "; ")))
  }
}

#' Load and standardize group information with configuration-driven mapping
#' 
#' Specialized function for loading group_info.csv with proper column standardization.
#' This replaces the old merge_with_group_info function with configuration-driven approach.
#' 
#' @param normalized_data Main data to merge with group information
#' @param session_id Session ID to locate group_info file
#' @param config Configuration object containing column_mappings
#' @return Merged data with standardized group columns or NULL on failure
load_and_merge_group_info <- function(normalized_data, session_id = NULL, config = NULL) {
  
  tryCatch({
    log_message("Loading and merging group information with configuration-driven mapping")
    
    # Determine group_info.csv file path using session paths
    if (!is.null(session_id)) {
      if (exists("get_session_paths")) {
        tryCatch({
          session_paths <- get_session_paths(session_id, create_dirs = FALSE)
          group_file_path <- file.path(session_paths$raw, "group_info.csv")
        }, error = function(e) {
          # Fallback to current working directory structure
          group_file_path <<- file.path(getwd(), "app_data", "sessions", session_id, "raw", "group_info.csv")
        })
      } else {
        group_file_path <- file.path(getwd(), "app_data", "sessions", session_id, "raw", "group_info.csv")
      }
    } else {
      if (exists("get_default_path")) {
        group_file_path <- get("get_default_path", mode = "function")("group_info")
      } else {
        group_file_path <- "tests/data/raw/group_info.csv"  # Fallback path
      }
    }
    
    if (!file.exists(group_file_path)) {
      log_message(sprintf("Group info file not found: %s", group_file_path), level = "warning")
      log_message("Returning original data without group information", level = "warning")
      return(normalized_data)  # Return original data instead of NULL
    }
    
    # Load group info data directly without column mapping (downstream jurisdiction principle)
    # According to Jurisdictional Mapping Principle, downstream files should NOT have column mappings
    log_message("Loading group_info.csv with original column names (no mapping applied)")
    
    group_info <- read_data_file(group_file_path)
    if (is.null(group_info)) {
      log_message("Failed to read group_info file", level = "error")
      log_message("Returning original data without group information", level = "warning")
      return(normalized_data)  # Return original data instead of NULL
    }
    
    # Basic validation
    if (nrow(group_info) == 0) {
      log_message("Group info file is empty", level = "error")
      log_message("Returning original data without group information", level = "warning")
      return(normalized_data)  # Return original data instead of NULL
    }
    
    # Intelligent species column detection for merging
    species_col <- NULL
    if ("species" %in% names(group_info)) {
      species_col <- "species"
    } else if ("Species" %in% names(group_info)) {
      # Apply transformation for downstream jurisdiction files
      names(group_info)[names(group_info) == "Species"] <- "species"
      species_col <- "species"
      log_message("Converted 'Species' column to 'species' for merging compatibility")
    } else {
      # Try to find species column by common patterns
      possible_species_cols <- grep("^[Ss]pecies", names(group_info), value = TRUE)
      if (length(possible_species_cols) > 0) {
        names(group_info)[names(group_info) == possible_species_cols[1]] <- "species"
        species_col <- "species"
        log_message(sprintf("Detected species column '%s', converted to 'species'", possible_species_cols[1]))
      } else {
        log_message("No species column found in group_info file", level = "error")
        log_message("Returning original data without group information", level = "warning")
        return(normalized_data)  # Return original data instead of NULL
      }
    }
    
    # Check if normalized_data has species column
    if (!"species" %in% names(normalized_data)) {
      log_message("Species column not found in normalized data for merging", level = "error")
      log_message("Returning original data without group information", level = "warning")
      return(normalized_data)  # Return original data instead of NULL
    }
    
    # Merge data using standardized column names
    merged_data <- merge(normalized_data, group_info, by = "species", all.x = TRUE)
    
    # Report merge results
    original_rows <- nrow(normalized_data)
    merged_rows <- nrow(merged_data)
    
    # Dynamically detect group columns for reporting (avoid hardcoding column names)
    group_columns <- names(group_info)[!names(group_info) %in% c("species", species_col)]
    if (length(group_columns) > 0) {
      # Use first available group column for matching statistics
      matched_species <- length(unique(merged_data$species[!is.na(merged_data[[group_columns[1]]])]))
    } else {
      # If no group columns detected, count total merged species
      matched_species <- length(unique(merged_data$species))
    }
    
    log_message(sprintf("Group info merge completed: %d -> %d rows, %d species matched", 
                       original_rows, merged_rows, matched_species))
    
    if (matched_species == 0) {
      log_message("No species were matched - check species names consistency", level = "warning")
    }
    
    return(merged_data)
    
  }, error = function(e) {
    log_message(sprintf("Failed to load and merge group info: %s", e$message), level = "error")
    log_message("Returning original data due to error in group info processing", level = "warning")
    return(normalized_data)  # Return original data instead of NULL
  })
}

# log_message("Unified configuration-driven data loading system loaded successfully")