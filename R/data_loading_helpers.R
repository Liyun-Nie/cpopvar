############################################################
#### data_loading_helpers.R - Data Loading Auxiliary Functions ####
############################################################
#
# Helper functions to support the task dispatcher with data loading
# operations. These functions bridge the gap between file-based and
# data-frame-based function signatures.
#
############################################################

#' @importFrom magrittr %>%
NULL

# Load required libraries

# Dependencies automatically loaded in R package context

#' Load main variant data using configuration-driven column mapping
#' 
#' @title Load main variant data using configuration-driven column mapping
#' @description Loads and processes main variant data with automatic column mapping
#' @param config_data Configuration data containing input file paths and column mappings
#' @param session_paths Session paths object
#' @return Data frame with standardized column names
#' @export
load_main_variant_data <- function(config_data, session_paths) {
  
  # V3.16: Simplified configuration - direct file path
  # Get the main data file path from the session directory (it should have been copied)
  original_path <- config_data$input_files$main_data
  
  # V3.16.2: Debug logging
  log_message(sprintf("[DEBUG] original_path type: %s, length: %d", class(original_path), length(original_path)))
  if (length(original_path) > 0) {
    log_message(sprintf("[DEBUG] original_path value: %s", original_path))
  }
  log_message(sprintf("[DEBUG] session_paths$raw: %s", session_paths$raw))
  
  # V3.16.1: Safer check - handle character(0), NULL, NA cases
  if (is.null(original_path) || length(original_path) == 0 || is.na(original_path) || !nzchar(original_path)) {
    stop("Configuration error: input_files$main_data is empty/null. Please check your configuration file.")
  }
  
  main_data_path <- file.path(session_paths$raw, basename(original_path))
  log_message(sprintf("[DEBUG] main_data_path: %s", main_data_path))
  
  # V3.16.1: Safer file.exists check
  if (length(main_data_path) == 0 || !file.exists(main_data_path)) {
      stop(sprintf("Main data file not found in session directory: %s. It should have been copied at startup.", main_data_path))
  }
  
  # Get column mappings
  column_mapping <- config_data$column_mappings$main_data
  if (is.null(column_mapping)) {
    stop("Column mappings for main_data not found in configuration")
  }
  
  # Required columns for variant data (using internal standard names)
  required_columns <- c("species", "var_type", "position", "sample_id")
  
  # Load and standardize the data
  log_message("Loading main variant data with configuration-driven column mapping")
  variant_data <- load_and_standardize_data(
    file_path = main_data_path,
    column_mapping = column_mapping,
    required_columns = required_columns,
    file_description = "main variant data",
    validate_data = TRUE,
    allow_missing_columns = FALSE
  )
  
  if (is.null(variant_data)) {
    stop("Failed to load main variant data")
  }
  
  log_message(sprintf("Successfully loaded %d rows of variant data", nrow(variant_data)))
  
  return(variant_data)
}

#' Load CDS lengths data from unified region_info_complete.csv
#' 
#' @description
#' Extracts CDS length information from the unified region_info_complete.csv file.
#' This function replaces the old dual-data-source approach (cds_lengths + region_info).
#' 
#' @param config_data Configuration data
#' @param session_paths Session paths object
#' @param species_genome_data Species genome configuration data (optional, passed through to load_region_info_data)
#' @return Data frame with CDS length data (species, gene_name, gene_length)
load_cds_lengths_data <- function(config_data, session_paths, species_genome_data = NULL) {
  
  log_message("Loading CDS lengths from unified region_info_complete.csv")
  
  # Load unified region info, passing species_genome_data
  region_info_complete <- load_region_info_data(config_data, session_paths, species_genome_data)

  if (is.null(region_info_complete)) {
    stop("Failed to load region_info_complete.csv")
  }
  
  # Extract CDS records and aggregate by gene
  # For multi-exon genes, sum all exon lengths
  cds_lengths <- region_info_complete %>%
    dplyr::filter(region_type == "CDS") %>%
    dplyr::group_by(species, region_name) %>%
    dplyr::summarise(
      gene_length = sum(region_length, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(gene_name = region_name)
  
  log_message(sprintf("Extracted CDS lengths: %d genes from %d species", 
                     nrow(cds_lengths), 
                     length(unique(cds_lengths$species))))
  
  return(cds_lengths)
  }
  
#' Load unified region info data (region_info_complete.csv)
#' 
#' @description
#' Loads the unified region_info_complete.csv file which contains:
#' - CDS records (with positions)
#' - IGS records (intergenic spacers)
#' - Intron records
#' 
#' Auto-detects mode based on provided input files:
#' - If "annotations" file is provided: Generate from source
#' - Otherwise: Error (annotations file is required)
#' 
#' @param config_data Configuration data
#' @param session_paths Session paths object
#' @param species_genome_data Species genome configuration data (optional, for genome_region assignment)
#' @return Data frame with complete region information (species, region_name, region_type, gene, region_start, region_end, region_length, genome_region)
load_region_info_data <- function(config_data, session_paths, species_genome_data = NULL) {
  
  # Auto-detect mode based on what files are provided
  annotations_path <- config_data$input_files$annotations
  
  # V3.16: Safer check - handle character(0) case
  has_valid_annotations <- !is.null(annotations_path) && 
                           length(annotations_path) > 0 && 
                           !is.na(annotations_path) && 
                           nzchar(annotations_path)
  
  if (has_valid_annotations) {
    # Mode: Generate from annotations
    log_message("Auto-detected mode: annotations (generating from source)")
    
    # Source file should be in raw/annotations directory
    source_file <- file.path(session_paths$raw, "annotations", basename(annotations_path))
    
    # V3.16.1: Safer file.exists check
    if (length(source_file) == 0 || !file.exists(source_file)) {
      stop(sprintf("Source annotations file not found: %s", source_file))
    }
    
    # V3.16.6: Use standardized path manager for P01 output
    # Extract session_id from session_paths (it should be available in the calling context)
    session_id <- basename(session_paths$base)
    output_dir <- get_processing_stage_path(session_id, "P01_preprocessed", create_dir = TRUE)
    log_message(sprintf("Using standardized P01 path: %s", output_dir))
    
    output_file <- file.path(output_dir, "region_info_complete.csv")
    
    # V3.16.8: Check if file already exists to avoid redundant processing
    if (file.exists(output_file)) {
      log_message(sprintf("Loading existing region_info_complete.csv: %s", output_file))
      region_data <- readr::read_csv(output_file, show_col_types = FALSE)
    } else {
      # Run preprocessing module
      log_message(sprintf("Preprocessing annotations: %s -> %s", source_file, output_file))
      
      # Load preprocessing module
      if (!exists("preprocess_annotations_to_region_info")) {
        stop("Preprocessing module not loaded. Please ensure annotation_preprocessing_helpers.R is sourced.")
      }
      
      result <- preprocess_annotations_to_region_info(
        annotations_path = source_file,
        species_genome_data = species_genome_data,  # used for genome_region assignment
        output_path = output_file,
        min_igs_length = 3,
        gap_threshold = 5000
      )
      
      region_data <- result$region_info
    }
    
  } else {
    stop("No valid annotation source found. Please provide 'annotations' file path in configuration (input_files$annotations).")
  }
  
  # Validate required columns
  required_cols <- c("species", "region_name", "region_type", "region_start", "region_end", "region_length")
  missing_cols <- setdiff(required_cols, colnames(region_data))
  
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns in region_info_complete.csv: %s", 
                paste(missing_cols, collapse = ", ")))
  }
  
  # V3.16.7: Standardize region_type (preserve case for CDS and intron)
      # V3.16.12: Simplified standardization - preprocessing already handles this
      # Only ensure consistent case for standard terms
      region_data <- region_data %>%
        dplyr::mutate(
          region_type = dplyr::case_when(
            region_type == "IGS" ~ "IGS",
            region_type == "CDS" ~ "CDS",
            region_type == "intron" ~ "intron",
            TRUE ~ region_type
          )
        )
  
  log_message(sprintf("Loaded %d region records from %d species", 
                     nrow(region_data), 
                     length(unique(region_data$species))))
  log_message(sprintf("  - CDS:     %d", sum(region_data$region_type == "CDS")))
  log_message(sprintf("  - IGS:     %d", sum(region_data$region_type == "IGS")))
  log_message(sprintf("  - Intron:  %d", sum(region_data$region_type == "intron")))
  
  return(region_data)
}

#' Load existing normalized data for visualization-only mode
#' 
#' @param session_id Session identifier
#' @param data_source Data source type ("filtered", "unfiltered", "normalized")
#' @return Existing normalized data or NULL if not found
load_existing_normalized_data <- function(session_id, data_source = "filtered") {
  
  log_message(sprintf("Attempting to load existing normalized data for session: %s", session_id))
  
  # Build expected file paths using the new P01/P02/P03 structure
  session_processed_dir <- get_processing_stage_path(session_id, "P03_normalized", create_dir = FALSE)
  
  if (!dir.exists(session_processed_dir)) {
    log_message(sprintf("Normalized data directory not found: %s", session_processed_dir), level="warning")
    return(NULL)
  }
  
  # Search for the standard file name
  file_path <- file.path(session_processed_dir, "normalized_frequencies.csv")
  if (file.exists(file_path)) {
    log_message(sprintf("Found existing normalized data: %s", file_path))
    tryCatch({
      data <- readr::read_csv(file_path, show_col_types = FALSE)
      log_message(sprintf("Successfully loaded %d rows of existing data", nrow(data)))
      return(data)
    }, error = function(e) {
      log_message(sprintf("Error loading file %s: %s", file_path, e$message), level = "warning")
    })
  }
  
  log_message("No existing normalized data found", level = "warning")
  return(NULL)
}

#' Copy all configured input files to the session directory
#'
#' This function recursively scans the entire `input_files` configuration
#' and copies any item that has a `path` and `target_location` to the
#' appropriate session directory.
#'
#' @param config_data Configuration data containing input_files section
#' @param session_paths Session paths object from get_session_paths()
#' @return Invisibly returns a list of copied files
copy_auxiliary_files <- function(config_data, session_paths) {
  log_message("Scanning configuration and copying all required input files...")
  
  input_files_config <- config_data$input_files
  
  if (is.null(input_files_config) || length(input_files_config) == 0) {
    log_message("No input files configured.")
    return(invisible(NULL))
  }
  
  # Hardcoded internal mapping of file types to target locations
  # This ensures consistent directory structure and prevents user misconfiguration
  TARGET_LOCATION_MAP <- list(
    main_data = "raw",
    annotations = "raw/annotations",
    group_info = "raw",
    genome_regions = "raw",
    species_order = "raw",
    gene_function_mapping = "raw"
  )
  
  # Required files list
  REQUIRED_FILES <- c("main_data", "annotations", "group_info")
  
  copied_files <- list()
  
  # Process each file type in input_files (simplified format)
  for (file_key in names(input_files_config)) {
    file_path <- input_files_config[[file_key]]
    
    # Skip if not a simple file path (i.e., nested structures for legacy support)
    if (!is.character(file_path) || length(file_path) != 1) {
      next
    }
    
    # Get target location from internal map
    target_location <- TARGET_LOCATION_MAP[[file_key]]
    
    if (is.null(target_location)) {
      log_message(sprintf("Unknown file type '%s', skipping", file_key), level = "warning")
      next
    }
    
    # Construct target directory path
    target_dir <- file.path(session_paths$base, target_location)
    
    # Ensure target directory exists
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    
    target_path <- file.path(target_dir, basename(file_path))
    
    if (!file.exists(target_path)) {
      if (file.exists(file_path)) {
        log_message(sprintf("Copying '%s' to: %s", file_key, target_location))
        file.copy(file_path, target_path, overwrite = FALSE)
        copied_files[[file_key]] <- target_path
      } else {
        # Check if this is a required file
        if (file_key %in% REQUIRED_FILES) {
          stop(sprintf("Required file '%s' not found at: %s", file_key, file_path))
        } else {
          log_message(sprintf("Optional file '%s' not found, skipping: %s", file_key, file_path), level = "warning")
        }
      }
    } else {
      log_message(sprintf("File '%s' already exists in session directory", file_key))
      copied_files[[file_key]] <- target_path
    }
  }
  
  # Legacy support: also handle old nested structure if present
  traverse_and_copy_legacy <- function(config_node, path_prefix = "") {
    if (is.list(config_node)) {
      # Check if this node represents a file to be copied (old format)
      if (!is.null(config_node$path) && !is.null(config_node$target_location)) {
        file_type_name <- paste(path_prefix, collapse = ".")
        
        original_path <- config_node$path
        target_dir_key <- config_node$target_location
        target_dir <- session_paths[[target_dir_key]]
        
        if (is.null(target_dir)) {
          log_message(sprintf("Target location '%s' for file '%s' is not a valid session path. Skipping.", 
                          target_dir_key, file_type_name), level = "warning")
          return()
        }
        
        # Ensure target directory exists
        dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
        
        target_path <- file.path(target_dir, basename(original_path))
        
        if (!file.exists(target_path)) {
          if (file.exists(original_path)) {
            log_message(sprintf("Copying file '%s' (legacy format) to session directory: %s", file_type_name, target_path))
            file.copy(original_path, target_path, overwrite = FALSE)
            copied_files[[file_type_name]] <<- target_path
          } else {
            log_message(sprintf("File not found (legacy format), skipping: %s", original_path), level = "warning")
          }
        }
      } else {
        # If not a file node, recurse through its children
        for (name in names(config_node)) {
          traverse_and_copy_legacy(config_node[[name]], c(path_prefix, name))
        }
      }
    }
  }
  
  # Also process legacy nested structures for backward compatibility
  traverse_and_copy_legacy(input_files_config, "input_files")
  
  log_message("File copying complete.")
  return(invisible(copied_files))
}

# log_message("Data loading helper functions loaded successfully")
