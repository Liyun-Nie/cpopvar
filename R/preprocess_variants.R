############################################################
#### preprocess_variants.R - Variant Data Preprocessing Functions ####
############################################################
#
# Functions for preprocessing variant data from raw input
# Includes region annotation, duplicate removal, and IR coordinate conversion
# Designed for Shiny compatibility with progress callbacks
#
############################################################

#' @importFrom magrittr %>%
#' @importFrom dplyr rowwise
NULL

# Load required libraries

# Source required configurations
# source("config/default_config.R")
# source("config/species_genome_regions.R")
# source("src/utils/common_utils.R")

# Dependencies automatically loaded in R package context

#' Preprocess variant data from task parameters
#' 
#' @title Preprocess variant data from task parameters
#' @description Processes variant data including region annotation and IR coordinate conversion
#' @param config Configuration list (required, no default - follows Task Controller pattern)
#' @param task_params Task parameters containing pre-loaded data in task_params$data
#' @param progress_callback Optional function for progress updates
#' @param remove_duplicates Whether to remove duplicate positions per species
#' @param save_intermediate Whether to save intermediate results
#' @return List containing processed data (both full genome and IRA-only) and metadata
#' @export
preprocess_variants <- function(config,
                               task_params,
                               progress_callback = NULL,
                               remove_duplicates = TRUE,
                               save_intermediate = TRUE) {
  
  if (!is.null(progress_callback)) {
    progress_callback("Starting variant preprocessing", 0.1)
  }
  
  log_message("Starting variant preprocessing")
  
  # Key fix: Get pre-loaded data from task_params, no file reading needed
  if (is.null(task_params$data)) {
    stop("Preprocessing task requires a 'data' object in task_params, but it was not found.")
  }
  raw_data <- task_params$data
  log_message(sprintf("Processing variant data frame provided by the task dispatcher (%d rows)", nrow(raw_data)))

  # Get species genome data from shared dependency injection
  if (is.null(task_params$species_genome_data)) {
    stop("Preprocessing task requires 'species_genome_data' in task_params, but it was not found.")
  }
  species_genome_data <- task_params$species_genome_data
  log_message(sprintf("Using shared species genome configuration data (%d species)", nrow(species_genome_data)))

  # Column mapping already completed by load_main_variant_data upstream
  # standardize_column_names and column_mapped parameter deprecated
  
  # Note: Duplicate removal is now performed at the end of the pipeline
  # to ensure proper handling of variant types and IRB conversions
  
  # Add genome region annotations
  if (!is.null(progress_callback)) {
    progress_callback("Adding genome region annotations", 0.4)
  }
  
  # Inject the loaded data as a dependency into the annotation function
  annotated_data <- add_genome_region_annotations(raw_data, species_genome_data, progress_callback)
  
  # Handle IR coordinate conversion
  if (!is.null(progress_callback)) {
    progress_callback("Converting IR coordinates", 0.8)
  }
  
  # Inject the loaded data as a dependency into the coordinate conversion function
  final_data <- convert_ir_coordinates(annotated_data, species_genome_data, progress_callback)
  
  # V3.16.12: MANDATORY STANDARDIZATION - Enforce standard terminology
  # This ensures global consistency regardless of user input format
  log_message("Applying mandatory standardization for region_type and var_type")
  
  # Standardize region_type: Convert all variations to standard terms
  if ("region_type" %in% names(final_data)) {
    final_data <- final_data %>%
      dplyr::mutate(
        region_type = dplyr::case_when(
          tolower(trimws(region_type)) %in% c("intergenic", "igs") ~ "IGS",
          tolower(trimws(region_type)) == "cds" ~ "CDS",
          tolower(trimws(region_type)) == "intron" ~ "intron",
          tolower(trimws(region_type)) %in% c("rrna", "trna-cds", "trna-intron", "rna") ~ "RNA",
          TRUE ~ "others"  # Other types categorized as "others"
        )
      )
    log_message(sprintf("Standardized region_type: %s", paste(unique(final_data$region_type), collapse = ", ")))
  }
  
  if ("region_type" %in% names(annotated_data)) {
    annotated_data <- annotated_data %>%
      dplyr::mutate(
        region_type = dplyr::case_when(
          tolower(trimws(region_type)) %in% c("intergenic", "igs") ~ "IGS",
          tolower(trimws(region_type)) == "cds" ~ "CDS",
          tolower(trimws(region_type)) == "intron" ~ "intron",
          tolower(trimws(region_type)) %in% c("rrna", "trna-cds", "trna-intron", "rna") ~ "RNA",
          TRUE ~ "others"
        )
      )
  }
  
  # Standardize var_type: Convert all variations to standard terms
  if ("var_type" %in% names(final_data)) {
    final_data <- final_data %>%
      dplyr::mutate(
        var_type = dplyr::case_when(
          tolower(trimws(var_type)) == "snp" ~ "snp",
          tolower(trimws(var_type)) %in% c("indel", "del", "ins") ~ "INDEL",
          tolower(trimws(var_type)) == "complex" ~ "complex",
          tolower(trimws(var_type)) == "mnp" ~ "mnp",
          TRUE ~ var_type  # Keep unknown types as-is
        )
      )
    log_message(sprintf("Standardized var_type: %s", paste(unique(final_data$var_type), collapse = ", ")))
  }
  
  if ("var_type" %in% names(annotated_data)) {
    annotated_data <- annotated_data %>%
      dplyr::mutate(
        var_type = dplyr::case_when(
          tolower(trimws(var_type)) == "snp" ~ "snp",
          tolower(trimws(var_type)) %in% c("indel", "del", "ins") ~ "INDEL",
          tolower(trimws(var_type)) == "complex" ~ "complex",
          tolower(trimws(var_type)) == "mnp" ~ "mnp",
          TRUE ~ var_type
        )
      )
  }
  
  # Universal value mapping: apply user-defined preprocessing rules
  # Replace hardcoded gsub logic with flexible configuration-driven mapping
  log_message("Applying user-defined value mapping rules for data preprocessing")
  
  # Check if config contains preprocessing rules
  if (!is.null(config) && !is.null(config$preprocessing) && !is.null(config$preprocessing$value_mapping_rules)) {
    log_message("Using configuration-defined preprocessing rules")
    preprocessing_rules <- config$preprocessing$value_mapping_rules
    
    # Apply variant type mapping rules
    variant_mapping <- convert_preprocessing_rules_to_mapping(preprocessing_rules, "variant_types")
    if (!is.null(variant_mapping)) {
      final_data <- apply_value_mapping(final_data, "var_type", variant_mapping)
      annotated_data <- apply_value_mapping(annotated_data, "var_type", variant_mapping)
    }
    
    # Apply region type mapping rules (if var_location column exists)
    region_mapping <- convert_preprocessing_rules_to_mapping(preprocessing_rules, "region_types")
    if (!is.null(region_mapping) && "var_location" %in% names(final_data)) {
      final_data <- apply_value_mapping(final_data, "var_location", region_mapping)
      annotated_data <- apply_value_mapping(annotated_data, "var_location", region_mapping)
    }
    
  } else {
    # Fallback to legacy hardcoded logic for backward compatibility
    log_message("No preprocessing rules found in config - using legacy ins/del -> INDEL conversion")
    
    # Ensure var_type is character type to support gsub operations
    if (is.factor(final_data$var_type)) {
      final_data$var_type <- as.character(final_data$var_type)
    }
    if (is.factor(annotated_data$var_type)) {
      annotated_data$var_type <- as.character(annotated_data$var_type)
    }
    
    # Apply legacy hardcoded conversion
    final_data$var_type <- gsub("ins", "INDEL", final_data$var_type)
    final_data$var_type <- gsub("del", "INDEL", final_data$var_type)
    annotated_data$var_type <- gsub("ins", "INDEL", annotated_data$var_type)
    annotated_data$var_type <- gsub("del", "INDEL", annotated_data$var_type)
  }

  # Note: Duplicate removal is now handled by the filtering module
  # This ensures that frequency filtering is based on original sample counts
  if (remove_duplicates) {
    log_message("WARNING: remove_duplicates=TRUE ignored in preprocessing. Deduplication now handled by filtering module.")
  }
  
  # Return unfiltered data with IRB variants included for downstream processing
  # Filtering module will handle deduplication based on sample frequencies
  data_full_genome <- final_data
  data_ira_only <- annotated_data
  
  # Create summary statistics using the IRA-only data to match original script's scope
  summary_stats <- create_preprocessing_summary(data_ira_only, raw_data)
  
  # Log the species-specific summary table as requested
  if (!is.null(summary_stats$by_species)) {
    log_message("--- Preprocessing Summary (IRA-only data) ---")
    # Use capture.output to format the data frame for logging
    summary_text <- paste(utils::capture.output(print(summary_stats$by_species)), collapse = "\n")
    log_message(summary_text)
  }
  
  # Save intermediate results if requested
  if (isTRUE(save_intermediate)) {
    log_message("DEBUG: Entering save_intermediate block", level = "debug")
    
    # Get session_id from config according to Rule #10 - with safer access
    session_id <- NULL
    if (!is.null(config) && !is.null(config$session_info)) {
      session_id <- config$session_info$session_id
    }
    log_message(sprintf("DEBUG: session_id = %s", as.character(session_id)), level = "debug")

    # V3.16: Safer check for session_id - handle character(0) case
    has_valid_session_id <- !is.null(session_id) && length(session_id) > 0 && !is.na(session_id) && nchar(as.character(session_id)) > 0
    
    if (has_valid_session_id) {
      # Use new structured P01 preprocessing output path
      log_message("DEBUG: Calling get_processing_stage_path", level = "debug")
      output_dir <- get_processing_stage_path(session_id, "P01_preprocessed", create_dir = TRUE)
      log_message(sprintf("DEBUG: output_dir = %s", output_dir), level = "debug")
    } else {
      # Fallback to legacy path for backward compatibility
      log_message("DEBUG: Using fallback path", level = "debug")
      if (exists("get_default_path")) {
        output_dir <- get("get_default_path", mode = "function")("processed_data")
      } else {
        output_dir <- "app_data/processed"
      }
    }
    
    # Directory is already created by get_processing_stage_path
    
    # =============================================
    # Compatibility output: write processed_genome_regions.csv for M04
    # =============================================
    
    log_message("Generating processed_genome_regions.csv for M04 module compatibility")
    
    # Create processed data output directory if needed  
    if (has_valid_session_id) {
      processed_output_dir <- get_processing_stage_path(session_id, "P01_preprocessed", create_dir = TRUE)
    } else {
      processed_output_dir <- output_dir
    }
    
    # Copy genome regions file to processed directory for M04 module access
    tryCatch({
      # V3.16: Simplified configuration - direct file path
      # Determine source genome regions file path from config
      genome_regions_source <- config$input_files$genome_regions
      
      if (!is.null(genome_regions_source) && file.exists(genome_regions_source)) {
        processed_genome_regions_path <- file.path(processed_output_dir, "processed_genome_regions.csv")
        file.copy(genome_regions_source, processed_genome_regions_path, overwrite = TRUE)
        log_message(sprintf("Successfully generated processed_genome_regions.csv: %s", processed_genome_regions_path))
      } else {
        log_message("Warning: Could not locate genome_regions source file for M04 processing", level = "warning")
      }
    }, error = function(e) {
      log_message(sprintf("Error generating processed_genome_regions.csv: %s", e$message), level = "warning")
    })
    
    # Save the full genome data
    output_file_full <- file.path(output_dir, "annotated_variants_full_genome.csv")
    if (exists("safe_write_file")) {
      safe_write_file(data_full_genome, output_file_full)
    } else {
      utils::write.csv(data_full_genome, output_file_full, row.names = FALSE)
    }
    log_message(sprintf("Full genome data with IRB variants saved to: %s", output_file_full), level = "INFO")
    
    # Save the IRA-only data for downstream analysis
    output_file_ira <- file.path(output_dir, "annotated_variants_ira_only.csv")
    if (exists("safe_write_file")) {
      safe_write_file(data_ira_only, output_file_ira)
    } else {
      utils::write.csv(data_ira_only, output_file_ira, row.names = FALSE)
    }
    log_message(sprintf("IRA-only variant data saved to: %s", output_file_ira), level = "INFO")
    
    # Save enhanced preprocessing statistics tables
    # Full genome statistics
    log_message("DEBUG: About to call create_preprocessing_summary for full genome", level = "debug")
    full_genome_summary <- create_preprocessing_summary(data_full_genome, raw_data)
    log_message("DEBUG: create_preprocessing_summary for full genome completed", level = "debug")
    output_stats_full <- file.path(output_dir, "preprocessing_stats_full_genome.csv")
    if (exists("safe_write_file")) {
      safe_write_file(full_genome_summary$by_species, output_stats_full)
    } else {
      utils::write.csv(full_genome_summary$by_species, output_stats_full, row.names = FALSE)
    }
    log_message(sprintf("Full genome preprocessing statistics saved to: %s", output_stats_full), level = "INFO")
    
    # IRA-only statistics
    log_message("DEBUG: About to call create_preprocessing_summary for IRA-only", level = "debug")
    ira_only_summary <- create_preprocessing_summary(data_ira_only, raw_data)
    log_message("DEBUG: create_preprocessing_summary for IRA-only completed", level = "debug")
    output_stats_ira <- file.path(output_dir, "preprocessing_stats_ira_only.csv")
    if (exists("safe_write_file")) {
      safe_write_file(ira_only_summary$by_species, output_stats_ira)
    } else {
      utils::write.csv(ira_only_summary$by_species, output_stats_ira, row.names = FALSE)
    }
    log_message(sprintf("IRA-only preprocessing statistics saved to: %s", output_stats_ira), level = "INFO")
  }
  
  if (!is.null(progress_callback)) {
    progress_callback("Preprocessing completed", 1.0)
  }
  
  log_message("Variant preprocessing completed successfully")

  # V3.5.1 Hotfix: Add defensive check to ensure result is a valid list
  result <- list(
    data = data_full_genome,
    data_ira_only = data_ira_only,
    summary = summary_stats,
    metadata = list(
      task_params = task_params,
      processing_time = Sys.time(),
      original_rows = nrow(raw_data),
      final_rows_with_irb = nrow(data_full_genome),
      final_rows_ira_only = nrow(data_ira_only),
      deduplication_deferred = TRUE,
      config_version = if(!is.null(config) && !is.null(config$version)) config$version else "unknown"
    )
  )

  # V3.6 Enhanced: More robust result validation
  if (is.null(result) || !is.list(result) || is.null(result$data) || !is.data.frame(result$data) || nrow(result$data) == 0) {
    log_message("FATAL: Preprocessing returned an invalid result. Returning NULL.", level = "error")
    return(NULL)
  }

  return(result)
}

#' Standardize column names in raw data
#' 
#' @param data Raw variant data (either with V columns or already mapped columns)
#' @return Data with standardized column names
standardize_column_names <- function(data) {
  
  # Define mapping between standard names and possible column names
  # This supports both legacy V column structure and configuration-driven mapped columns
  column_mapping <- list(
    sample_id = c("sample_id", "V1", "sample", "id"),
    species = c("species", "V2", "organism", "species_name"),
    position = c("position", "V3", "pos", "location"),
    var_type = c("var_type", "variant_type", "V4", "mutation_type"),  # Updated order
    ref_allele = c("ref_allele", "V5", "ref", "reference"),
    alt_allele = c("alt_allele", "V6", "alt", "alternative")
  )
  
  # Check if data already has standard column names (from column mapping)
  standard_names <- names(column_mapping)
  has_standard_names <- all(standard_names %in% names(data))
  
  if (has_standard_names) {
    log_message("Data already has standard column names (from column mapping)")
    return(data)
  }
  
  # Check for V column structure (backward compatibility)
  expected_v_cols <- c("V1", "V2", "V3", "V4", "V5", "V6")
  has_v_columns <- all(expected_v_cols %in% names(data))
  
  if (has_v_columns) {
    log_message("Standardizing V column names to descriptive names")
    data <- data %>%
      dplyr::rename(
        sample_id = V1,
        species = V2,
        position = V3,
        var_type = V4,
        ref_allele = V5,
        alt_allele = V6
      )
  } else {
    # Try to map columns based on available names
    log_message("Attempting to map columns based on available names")
    
    for (standard_name in standard_names) {
      possible_names <- column_mapping[[standard_name]]
      found_column <- NULL
      
      for (possible_name in possible_names) {
        if (possible_name %in% names(data)) {
          found_column <- possible_name
          break
        }
      }
      
      if (is.null(found_column)) {
        stop(paste("Could not find column for", standard_name, ". Expected one of:", paste(possible_names, collapse = ", ")))
      }
      
      # Rename the column if it's not already the standard name
      if (found_column != standard_name) {
        names(data)[names(data) == found_column] <- standard_name
      }
    }
  }
  
  # The gsub call for sample_id is no longer needed as we now handle this during file read.
  # data$sample_id <- gsub("\"", "", data$sample_id)
  
  # Standardize species names
  if (exists("clean_species_names")) {
    data$species <- clean_species_names(data$species)
  } else {
    # Basic species name cleaning
    data$species <- gsub("\"", "", data$species)
    data$species <- trimws(data$species)
  }
  
  # The standardize_variant_types function introduces logic that conflicts with
  # the original script. The correct handling (ins/del -> INDEL) is already
  # performed later in the pipeline, so we remove this redundant call.
  # data$var_type <- standardize_variant_types(data$var_type)
  
  log_message("Column names standardized")
  
  return(data)
}

#' Remove duplicate variants based on species, position, and variant type
#' 
#' @param data Variant data
#' @return Data with duplicates removed
#' @details This function removes duplicates by considering species, position, and variant type
#'          This ensures that different variant types at the same position are preserved
remove_duplicates <- function(data) {
  
  original_rows <- nrow(data)
  
  # Remove duplicates by species, position, and variant type, keeping the first occurrence
  # This matches the logic from the original 01_filter_raw_data.R script
  deduplicated_data <- data %>%
    dplyr::group_by(species, position, var_type) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
  
  duplicates_removed <- original_rows - nrow(deduplicated_data)
  
  log_message(sprintf("Removed %d duplicate variants (species+position+type)", duplicates_removed))
  
  return(deduplicated_data)
}

#' Add genome region annotations to variant data (Refactored for Dependency Injection)
#' 
#' @param data Variant data
#' @param species_genome_data The data frame loaded from user-provided species_genome_regions.csv
#' @param progress_callback Optional progress callback
#' @return Data with genome region annotations
#' @details This function is optimized to use vectorized operations from dplyr,
#'          avoiding slow row-by-row loops for much faster performance.
add_genome_region_annotations <- function(data, species_genome_data, progress_callback = NULL) {
  
  log_message("Starting vectorized genome region annotation...")
  
  # Ensure 'position' is numeric for comparisons, as it was read as character.
  data$position <- as.numeric(data$position)
  
  # Use dplyr's group_by and mutate for a fast, vectorized operation.
  annotated_data <- data %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      genome_region = {
        # Get the specific configuration for the current species being processed
        species_name <- dplyr::cur_group()$species
        # Use injected data, not a global variable
        species_config <- species_genome_data[species_genome_data$species == species_name, ]
        
        if(nrow(species_config) == 0) return("unknown")
        
        # Check for special handling cases like IR-lacking genomes
        is_ir_lacking <- !is.null(species_config$special_handling) && 
                         !is.na(species_config$special_handling) && 
                         species_config$special_handling == "IR_lacking_genome"

        # Use case_when for clear, vectorized conditional logic
        # IMPORTANT: Use the vectorized '&' operator, not '&&'
        dplyr::case_when(
          is_ir_lacking ~ "whole_genome",
          !is.na(position) & position >= species_config$lsc_start & position <= species_config$lsc_end ~ "LSC",
          !is.na(position) & position >= species_config$ir_start & position <= species_config$ir_end ~ "IRA",
          !is.na(position) & position >= species_config$ssc_start & position <= species_config$ssc_end ~ "SSC",
          TRUE ~ "unknown"
        )
      }
    ) %>%
    dplyr::ungroup()

  missing_regions <- sum(annotated_data$genome_region == "unknown", na.rm = TRUE)
  if (missing_regions > 0) {
    log_message(sprintf("Warning: %d variants were classified as 'unknown' region", missing_regions), level = "warning")
  }
  
  log_message("Vectorized genome region annotation completed successfully.")
  
  return(annotated_data)
}

#' Convert IR coordinates to include IRB positions (Refactored for Dependency Injection)
#' 
#' @param data Variant data with genome regions
#' @param species_genome_data The data frame loaded from user-provided species_genome_regions.csv
#' @param progress_callback Optional progress callback
#' @return Data with IR coordinates converted
convert_ir_coordinates <- function(data, species_genome_data, progress_callback = NULL) {
  
  # Fix: Intelligently identify variants that need conversion (located in IR regions defined in config file)
  # Use more stable method for IRA region filtering (already standardized as IRA)
  ir_candidates <- data %>%
    dplyr::filter(genome_region == "IRA" & !is.na(genome_region))
  
  if (nrow(ir_candidates) == 0) {
    log_message("No variants found in IRA regions for coordinate conversion.")
    return(data)
  }
  
  # Check each variant row-by-row to see if it's in the corresponding species' primary IR region
  variants_in_primary_ir <- ir_candidates %>%
    dplyr::rowwise() %>%
    dplyr::filter({
      # Use injected data
      current_species_config <- species_genome_data[species_genome_data$species == species, ]
      # Ensure only variants located in primary IR region are selected for conversion
      # Fix: Use correct CSV column names to access IR region boundaries
      nrow(current_species_config) > 0 &&
      !is.na(current_species_config$ir_start) && !is.na(current_species_config$ir_end) &&
      position >= current_species_config$ir_start && 
      position <= current_species_config$ir_end
    }) %>%
    dplyr::ungroup()

  if (nrow(variants_in_primary_ir) == 0) {
    log_message("No variants found in primary IRA regions for coordinate conversion.")
    return(data)
  }
  
  log_message(sprintf("Converting IRA coordinates for %d variants from primary IRA region", nrow(variants_in_primary_ir)))
  
  # Mark original IR variants as IRA specifically (now handled at annotation stage)
  # data$genome_region[data$genome_region == "IR"] <- "IRA"  # No longer needed - IRA assigned during annotation
  
  # Create IRB copies
  irb_variants <- variants_in_primary_ir
  
  # Prepare genome regions data for joining - vectorized approach for optimal performance
  # Use injected data for join
  species_info <- species_genome_data %>%
    dplyr::select(species, total_length, ir_start)

  # Vectorized conversion - fast and efficient physical coordinate transformation
  irb_variants <- irb_variants %>%
    dplyr::left_join(species_info, by = "species") %>%
    dplyr::mutate(
      position = total_length - (position - ir_start)
    ) %>%
    dplyr::select(-total_length, -ir_start) # Clean up helper columns
    
  log_message(sprintf("Successfully applied vectorized physical coordinate conversion to %d IRB variants.", nrow(irb_variants)))
  
  # Explicitly mark all copied variants as IRB
  irb_variants$genome_region <- "IRB"
  
  # Add IRB variants to original data
  combined_data <- rbind(data, irb_variants)
  
  log_message(sprintf("Added %d IRB variants, total variants: %d", nrow(irb_variants), nrow(combined_data)))
  
  return(combined_data)
}

#' Create preprocessing summary statistics
#' 
#' @param processed_data Final processed data
#' @param original_data Original raw data
#' @return Summary statistics data frame
create_preprocessing_summary <- function(processed_data, original_data) {
  # Overall statistics
  overall_stats <- data.frame(
    metric = c("original_variants", "processed_variants", "species_count", "variant_types"),
    value = c(
      nrow(original_data),
      nrow(processed_data),
      length(unique(processed_data$species)),
      length(unique(processed_data$var_type))
    ),
    stringsAsFactors = FALSE
  )
  
  # Enhanced species-wise statistics with detailed breakdowns
  species_stats <- processed_data %>% 
    dplyr::group_by(species) %>% 
    dplyr::summarise(
      total_variants = dplyr::n(),
      # Variant type counts
      snp_count = sum(var_type == "snp", na.rm = TRUE),
      indel_count = sum(var_type %in% c("del", "ins", "INDEL"), na.rm = TRUE),
      complex_count = sum(var_type == "complex", na.rm = TRUE),
      mnp_count = sum(var_type == "mnp", na.rm = TRUE),
      # Genome region counts (LSC, IRA, SSC, IRB, whole_genome)
      lsc_variants = sum(genome_region == "LSC", na.rm = TRUE),
      ira_variants = sum(genome_region == "IRA", na.rm = TRUE),
      ssc_variants = sum(genome_region == "SSC", na.rm = TRUE),
      irb_variants = sum(genome_region == "IRB", na.rm = TRUE),
      whole_genome_variants = sum(genome_region == "whole_genome", na.rm = TRUE),
      # Region type counts (CDS, IGS, intron, RNA, others)
      cds_variants = sum(region_type == "CDS", na.rm = TRUE),
      igs_variants = sum(region_type == "IGS", na.rm = TRUE),
      intron_variants = sum(region_type == "intron", na.rm = TRUE),
      rna_variants = sum(region_type == "RNA", na.rm = TRUE),
      other_variants = sum(region_type == "others", na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    # V3.5.4 Final Hotfix: More robust NA replacement across all count columns
    # V3.6 Enhancement: Use safer column selection to avoid potential NA issues
    dplyr::mutate(dplyr::across(
      .cols = dplyr::where(is.numeric),
      .fns = ~ as.numeric(tidyr::replace_na(., 0))
    ))
  
  # Variant type statistics
  var_type_stats <- processed_data %>%
    dplyr::group_by(var_type) %>%
    dplyr::summarise(
      count = dplyr::n(),
      species_count = length(unique(species)),
      .groups = "drop"
    )
  
  # Genome region statistics - using safer NA handling
  region_stats <- processed_data %>%
    # Replace NA and NULL values with "unknown" before filtering
    dplyr::mutate(genome_region = ifelse(is.na(genome_region) | is.null(genome_region), "unknown", genome_region)) %>%
    dplyr::filter(genome_region != "unknown") %>%
    dplyr::group_by(genome_region) %>%
    dplyr::summarise(
      count = dplyr::n(),
      species_count = length(unique(species)),
      .groups = "drop"
    )
  
  # Return a list of all summary tables for flexibility.
  # The calling function can then access whichever part it needs.
  return(list(
    overall = overall_stats,
    by_species = species_stats,
    by_var_type = var_type_stats,
    by_region = region_stats
  ))
}

#' Validate preprocessing results
#' 
#' @param results Results from preprocess_variants
#' @param config Configuration list
#' @return Logical, TRUE if validation passes
validate_preprocessing_results <- function(results, config) {
  
  data <- results$data
  
  # Check required columns
  required_cols <- c("sample_id", "species", "position", "var_type", "genome_region")
  if (!all(required_cols %in% names(data))) {
    log_message("Missing required columns in processed data", level = "error")
    return(FALSE)
  }
  
  # Check for minimum data requirements
  if (nrow(data) < config$validation$min_rows_filtered_data) {
    log_message(sprintf("Too few rows in processed data: %d", nrow(data)), level = "error")
    return(FALSE)
  }
  
  # Check species coverage
  if (exists("get_all_species")) {
    expected_species <- get("get_all_species", mode = "function")()
    missing_species <- setdiff(expected_species, unique(data$species))
    if (length(missing_species) > 0) {
      log_message(sprintf("Missing species in processed data: %s", 
                         paste(missing_species, collapse = ", ")), level = "warning")
    }
  } else {
    log_message("Species coverage validation skipped - get_all_species() not available")
  }
  
  # Check variant types
  valid_var_types <- c("snp", "del", "ins", "indel", "complex", "mnp")
  invalid_var_types <- setdiff(unique(data$var_type), valid_var_types)
  if (length(invalid_var_types) > 0) {
    log_message(sprintf("Invalid variant types found: %s", 
                       paste(invalid_var_types, collapse = ", ")), level = "warning")
  }
  
  log_message("Preprocessing validation completed")
  return(TRUE)
}

#' Preprocess mapped variant data (convenience wrapper for UI)
#' 
#' @param mapped_data Data frame with already mapped columns
#' @param config Configuration list
#' @param progress_callback Optional progress callback
#' @param remove_duplicates Whether to remove duplicates
#' @return Preprocessed data with standardized structure
preprocess_mapped_variants <- function(mapped_data, 
                                      config,
                                      progress_callback = NULL,
                                      remove_duplicates = TRUE) {
  
  log_message("Starting preprocessing of mapped variant data")
  
  # Validate that required columns exist (handle both var_type and variant_type)
  required_cols_primary <- c("sample_id", "species", "position", "ref_allele", "alt_allele")
  variant_type_cols <- c("var_type", "variant_type")
  
  missing_primary <- setdiff(required_cols_primary, names(mapped_data))
  has_variant_type <- any(variant_type_cols %in% names(mapped_data))
  
  if (length(missing_primary) > 0) {
    stop(paste("Missing required columns in mapped data:", paste(missing_primary, collapse = ", ")))
  }
  
  if (!has_variant_type) {
    stop("Missing variant type column in mapped data. Expected one of: var_type, variant_type")
  }
  
  # Standardize variant_type to var_type if needed
  if ("variant_type" %in% names(mapped_data) && !"var_type" %in% names(mapped_data)) {
    names(mapped_data)[names(mapped_data) == "variant_type"] <- "var_type"
    log_message("Renamed 'variant_type' column to 'var_type' for consistency")
  }
  
  # Ensure position is numeric
  if (!is.numeric(mapped_data$position)) {
    mapped_data$position <- as.numeric(mapped_data$position)
  }
  
  # Call main preprocessing function with modernized signature
  task_params <- list(data = mapped_data)
  result <- preprocess_variants(
    config = config,
    task_params = task_params,
    progress_callback = progress_callback,
    remove_duplicates = remove_duplicates,
    save_intermediate = FALSE  # Don't save intermediate results for UI processing
  )
  
  log_message("Mapped variant data preprocessing completed successfully")
  return(result)
}

# log_message("Variant preprocessing functions loaded successfully")