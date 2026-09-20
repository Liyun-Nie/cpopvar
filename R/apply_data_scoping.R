############################################################
#### apply_data_scoping.R - Data Scoping Utility Function ####
############################################################
#
# Pure utility function for applying data_scoping filters to datasets
# Implements modular data scoping logic for visualization modules
#
# Phase 5 Task J: Core utility for modular data_scoping implementation
# This function enables visualization modules to dynamically filter
# complete normalized datasets based on their task-specific requirements
#
############################################################

#' Apply data scoping filters to a dataset
#'
#' This is a pure utility function that applies filtering based on task parameters.
#' It follows the dependency injection principle and does not access any global
#' configuration or state. All filtering parameters are passed explicitly.
#'
#' @param data Data frame to be filtered (typically normalized frequency data)
#' @param task_params Task parameters list containing data_scoping configuration
#' @return Filtered data frame based on data_scoping criteria
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # Example task_params structure:
#' task_params <- list(
#'   data_scoping = list(
#'     target_region_types = c("CDS", "IGS", "intron"),
#'     target_var_types = c("snp", "INDEL"),
#'     target_species = c("Species_A", "Species_B")
#'   )
#' )
#' # Example usage (requires data frame with appropriate columns):
#' # filtered_data <- apply_data_scoping(your_data, task_params)
#' @importFrom magrittr %>%
#' @export
apply_data_scoping <- function(data, task_params) {
  
  # Input validation
  if (is.null(data) || !is.data.frame(data)) {
    log_message("apply_data_scoping: Invalid data input (NULL or not a data.frame)", level = "error")
    return(data.frame())
  }
  
  if (nrow(data) == 0) {
    log_message("apply_data_scoping: Input data is empty, returning empty data.frame", level = "warning")
    return(data)
  }
  
  # Extract data_scoping configuration
  data_scoping <- NULL
  if (!is.null(task_params) && is.list(task_params)) {
    data_scoping <- task_params$data_scoping
  }
  
  # V3.16.14: Debug logging for data_scoping configuration
  if (is.null(data_scoping)) {
    log_message("[DEBUG] data_scoping is NULL in task_params", level = "warning")
    log_message(sprintf("[DEBUG] task_params keys: %s", paste(names(task_params), collapse = ", ")))
  } else if (!is.list(data_scoping)) {
    log_message(sprintf("[DEBUG] data_scoping is not a list, type: %s", class(data_scoping)), level = "warning")
  } else {
    log_message(sprintf("[DEBUG] data_scoping found with keys: %s", paste(names(data_scoping), collapse = ", ")))
  }
  
  # If no data_scoping configuration is provided, return original data
  if (is.null(data_scoping) || !is.list(data_scoping)) {
    log_message("apply_data_scoping: No data_scoping configuration found, returning original data")
    return(data)
  }
  
  log_message("=== APPLYING DATA_SCOPING FILTERS ===")
  
  # Track original data dimensions
  original_rows <- nrow(data)
  original_variants <- if("variant_count" %in% names(data)) sum(data$variant_count, na.rm = TRUE) else NA
  
  log_message(sprintf("ORIGINAL DATA: %d rows, %s variants", 
                     original_rows, 
                     if(is.na(original_variants)) "unknown" else as.character(original_variants)))
  
  # Make a copy of the data for filtering
  filtered_data <- data
  
  # Apply region_type filtering
  if (!is.null(data_scoping$target_region_types) && 
      is.character(data_scoping$target_region_types) && 
      length(data_scoping$target_region_types) > 0) {
    
    # CRITICAL: Check for "all_types" special value
    if (length(data_scoping$target_region_types) == 1 && 
        data_scoping$target_region_types == "all_types") {
      log_message("REGION_TYPE FILTER: all_types -> 0 rows (0 filtered out)")
      # Skip filtering - accept all region types
    } else if ("region_type" %in% names(filtered_data)) {
      before_rows <- nrow(filtered_data)
      
      filtered_data <- filtered_data %>%
        dplyr::filter(region_type %in% data_scoping$target_region_types)
      
      after_rows <- nrow(filtered_data)
      log_message(sprintf("REGION_TYPE FILTER: %s -> %d rows (%d filtered out)", 
                         paste(data_scoping$target_region_types, collapse=", "), 
                         after_rows, before_rows - after_rows))
    } else {
      log_message("REGION_TYPE FILTER: 'region_type' column not found, skipping filter", level = "warning")
    }
  }
  
  # Apply var_type filtering
  if (!is.null(data_scoping$target_var_types) && 
      is.character(data_scoping$target_var_types) && 
      length(data_scoping$target_var_types) > 0) {
    
    # CRITICAL: Check for "all_types" special value
    if (length(data_scoping$target_var_types) == 1 && 
        data_scoping$target_var_types == "all_types") {
      log_message("VAR_TYPE FILTER: all_types -> 0 rows (0 filtered out)")
      # Skip filtering - accept all var types
    } else if ("var_type" %in% names(filtered_data)) {
      before_rows <- nrow(filtered_data)
      
      filtered_data <- filtered_data %>%
        dplyr::filter(var_type %in% data_scoping$target_var_types)
      
      after_rows <- nrow(filtered_data)
      log_message(sprintf("VAR_TYPE FILTER: %s -> %d rows (%d filtered out)", 
                         paste(data_scoping$target_var_types, collapse=", "), 
                         after_rows, before_rows - after_rows))
    } else {
      log_message("VAR_TYPE FILTER: 'var_type' column not found, skipping filter", level = "warning")
    }
  }
  
  # Apply genome_region filtering
  if (!is.null(data_scoping$target_genome_regions) && 
      is.character(data_scoping$target_genome_regions) && 
      length(data_scoping$target_genome_regions) > 0) {
    
    if ("genome_region" %in% names(filtered_data)) {
      before_rows <- nrow(filtered_data)
      
      filtered_data <- filtered_data %>%
        dplyr::filter(genome_region %in% data_scoping$target_genome_regions)
      
      after_rows <- nrow(filtered_data)
      log_message(sprintf("GENOME_REGION FILTER: %s -> %d rows (%d filtered out)", 
                         paste(data_scoping$target_genome_regions, collapse=", "), 
                         after_rows, before_rows - after_rows))
    } else {
      log_message("GENOME_REGION FILTER: 'genome_region' column not found, skipping filter", level = "warning")
    }
  }
  
  # Apply species filtering (optional, less commonly used)
  if (!is.null(data_scoping$target_species) && 
      is.character(data_scoping$target_species) && 
      length(data_scoping$target_species) > 0) {
    
    if ("species" %in% names(filtered_data)) {
      before_rows <- nrow(filtered_data)
      
      filtered_data <- filtered_data %>%
        dplyr::filter(species %in% data_scoping$target_species)
      
      after_rows <- nrow(filtered_data)
      log_message(sprintf("SPECIES FILTER: %d species -> %d rows (%d filtered out)", 
                         length(data_scoping$target_species), 
                         after_rows, before_rows - after_rows))
    } else {
      log_message("SPECIES FILTER: 'species' column not found, skipping filter", level = "warning")
    }
  }
  
  # Apply gene filtering (optional, for gene-specific analyses)
  if (!is.null(data_scoping$target_genes) && 
      is.character(data_scoping$target_genes) && 
      length(data_scoping$target_genes) > 0) {
    
    if ("gene" %in% names(filtered_data)) {
      before_rows <- nrow(filtered_data)
      
      filtered_data <- filtered_data %>%
        dplyr::filter(gene %in% data_scoping$target_genes)
      
      after_rows <- nrow(filtered_data)
      log_message(sprintf("GENE FILTER: %d genes -> %d rows (%d filtered out)", 
                         length(data_scoping$target_genes), 
                         after_rows, before_rows - after_rows))
    } else {
      log_message("GENE FILTER: 'gene' column not found, skipping filter", level = "warning")
    }
  }
  
  # Calculate final statistics
  final_rows <- nrow(filtered_data)
  final_variants <- if("variant_count" %in% names(filtered_data)) sum(filtered_data$variant_count, na.rm = TRUE) else NA
  
  # Summary of filtering results
  retention_rate <- if(original_rows > 0) round((final_rows / original_rows) * 100, 2) else 0
  
  log_message("DATA_SCOPING RESULTS:")
  log_message(sprintf("  Rows: %d -> %d (%.2f%% retained)", 
                     original_rows, final_rows, retention_rate))
  
  if (!is.na(original_variants) && !is.na(final_variants)) {
    variant_retention_rate <- if(original_variants > 0) round((final_variants / original_variants) * 100, 2) else 0
    log_message(sprintf("  Variants: %d -> %d (%.2f%% retained)", 
                       original_variants, final_variants, variant_retention_rate))
  }
  
  # Validation checks
  if (final_rows == 0) {
    log_message("WARNING: Data_scoping resulted in empty dataset", level = "warning")
  }
  
  # Log unique values in key columns for debugging
  if (final_rows > 0) {
    if ("region_type" %in% names(filtered_data)) {
      unique_regions <- unique(filtered_data$region_type)
      log_message(sprintf("  Final region_types: %s", paste(unique_regions, collapse=", ")))
    }
    
    if ("var_type" %in% names(filtered_data)) {
      unique_var_types <- unique(filtered_data$var_type)
      log_message(sprintf("  Final var_types: %s", paste(unique_var_types, collapse=", ")))
    }
    
    # Log genome_region distribution
    if ("genome_region" %in% names(filtered_data)) {
      unique_genome_regions <- unique(filtered_data$genome_region)
      log_message(sprintf("  Final genome_regions: %s", paste(unique_genome_regions, collapse=", ")))
    }
    
    if ("species" %in% names(filtered_data)) {
      unique_species <- length(unique(filtered_data$species))
      log_message(sprintf("  Final species count: %d", unique_species))
    }
  }
  
  log_message("=== DATA_SCOPING COMPLETED ===")
  
  return(filtered_data)
}

# Load required dependencies

# Data scoping utility function loaded successfully