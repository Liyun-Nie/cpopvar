############################################################
#### filter_variants.R - Variant Filtering Functions (Redesigned) ####
############################################################

#' @importFrom magrittr %>%
NULL
# Dependencies automatically loaded in R package context
#
# Functions for filtering variant data based on frequency thresholds
# Creates filtered and unfiltered data branches for analysis
# Implements dual-branch processing with proper sample-based frequency filtering
#
# Key Design Changes:
# 1. Accepts unfiltered data from preprocessing (with all sample records)
# 2. Implements frequency filtering based on original sample counts
# 3. Applies deduplication only after filtering (species+var_type+position)
# 4. Returns both filtered and unfiltered branches
#
############################################################

# Load required libraries

# Source required configurations
# source("config/default_config.R")
# source("src/utils/common_utils.R")
# source("src/utils/visualization_themes.R")

#' Filter variants using configuration-driven thresholds (Refactored for Dependency Injection)
#' 
#' @title Filter variants using configuration-driven thresholds
#' @description Applies frequency-based filtering to variant data using configurable thresholds
#' @param input_data Preprocessed variant data from P01 stage
#' @param config Configuration list containing filtering parameters
#' @param species_genome_data A data frame containing species genome region info, injected by the dispatcher
#' @param progress_callback Optional function for progress updates
#' @param save_results Whether to save filtering results to P02 directory
#' @return Filtered variant data with statistics
#' @export
filter_variants <- function(input_data, 
                           config,
                           species_genome_data, # Injected Dependency
                           progress_callback = NULL,
                           save_results = TRUE) {
  
  # Initialize
  if (!is.null(progress_callback)) {
    progress_callback("Starting configuration-driven variant filtering", 0.1)
  }
  
  log_message("=== Starting Configuration-Driven Variant Filtering ===")
  log_message("Design: Unified filtering based on config parameters")
  
  # Load data if file path provided
  if (is.character(input_data)) {
    input_data <- safe_read_file(input_data)
    if (is.null(input_data)) {
      stop("Failed to read input data for filtering")
    }
  }
  
  # Validate input data
  required_columns <- c("species", "position", "var_type", "genome_region")
  if (!validate_data_frame(input_data, required_columns, min_rows = 1, data_name = "input data for filtering")) {
    stop("Invalid input data structure for filtering")
  }
  
  if (!is.null(progress_callback)) {
    progress_callback("Input data validated", 0.1)
  }
  
  log_message(sprintf("Input data: %d variants across %d species", 
                     nrow(input_data), length(unique(input_data$species))))
  
  # === Apply configuration-driven filtering ===
  if (!is.null(progress_callback)) {
    progress_callback("Applying frequency-based filtering", 0.3)
  }
  
  log_message("Applying configuration-driven frequency filtering")
  filtered_data <- apply_frequency_filters_and_dedup(input_data, config, progress_callback)
  
  # === Generate filtering statistics ===
  if (!is.null(progress_callback)) {
    progress_callback("Generating filtering statistics", 0.7)
  }
  
  # Generate statistics, passing down the injected dependency
  filtering_stats <- generate_filtering_statistics(input_data, filtered_data, config, species_genome_data)
  
  # === Create filtering visualizations ===
  plots <- NULL
  # Skip visualization generation for now - not critical for data processing
  if (!is.null(progress_callback)) {
    progress_callback("Skipping filtering visualizations (not critical)", 0.9)
  }
  
  # === Save results if requested ===
  if (save_results) {
    save_filtering_results(filtered_data, filtering_stats, plots, config)
  }
  
  if (!is.null(progress_callback)) {
    progress_callback("Configuration-driven filtering completed", 1.0)
  }
  
  log_message("=== Configuration-Driven Variant Filtering Completed ===")
  log_message(sprintf("Input variants: %d", nrow(input_data)))
  log_message(sprintf("Filtered variants: %d (%.1f%% retention)", 
                     nrow(filtered_data), 
                     (nrow(filtered_data) / nrow(input_data)) * 100))
  
  # Return simplified filtering results
  return(list(
    data = filtered_data,
    statistics = filtering_stats,
    plots = plots,
    metadata = list(
      filtering_time = Sys.time(),
      original_variants = nrow(input_data),
      filtered_variants = nrow(filtered_data),
      retention_rate = nrow(filtered_data) / nrow(input_data),
      config_version = config$version,
      filtering_thresholds = config$filtering$thresholds
    )
  ))
}

#' --- NEW: Simplified Helper Functions for Unified Filtering ---

#' Apply configuration-driven frequency filters and deduplication
#'
#' This is the core filtering function. It iterates through species, applies
#' frequency filters based on config thresholds, and then deduplicates the results.
#'
#' @param data Unfiltered variant data from P01 stage
#' @param config Configuration list
#' @param progress_callback Optional progress callback
#' @return A single data frame of filtered and deduplicated variants.
apply_frequency_filters_and_dedup <- function(data, config, progress_callback = NULL) {
  
  thresholds <- config$filtering
  log_message("Applying frequency filters with thresholds:")
  
  # Dynamic threshold configuration - exclusive approach
  log_message("Dynamic threshold configuration:")
  if (!is.null(thresholds$thresholds) && is.list(thresholds$thresholds)) {
    for (vtype in names(thresholds$thresholds)) {
      log_message(sprintf("  %s: >= %d samples", vtype, thresholds$thresholds[[vtype]]))
    }
  } else {
    log_message("Warning: No dynamic thresholds found in config$filtering$thresholds")
    log_message("System will use default threshold (1) for all variant types")
  }
  
  species_list <- unique(data$species)
  total_species <- length(species_list)
  
  # Process each species, filter, and collect results
  filtered_species_list <- lapply(seq_along(species_list), function(i) {
    species <- species_list[i]
    if (!is.null(progress_callback)) {
      progress <- 0.3 + (i / total_species) * 0.4 # Progress from 0.3 to 0.7
      progress_callback(sprintf("Filtering species: %s", species), progress)
    }
    
    species_data <- data[data$species == species, ]
    if (nrow(species_data) == 0) return(NULL)
    
    # This function now contains the complete logic for one species
    filter_and_dedup_species(species_data, thresholds)
  })
  
  # Combine results and return
  filtered_data <- do.call(rbind, filtered_species_list)
  if (is.null(filtered_data)) {
    return(data.frame())
  }
  
  log_message(sprintf("Total variants after filtering and deduplication: %d", nrow(filtered_data)))
  return(filtered_data)
}

#' Filter and deduplicate variants for a single species
#'
#' @param species_data Variant data for one species
#' @param thresholds Filtering thresholds from config
#' @return A filtered and deduplicated data frame for the species.
filter_and_dedup_species <- function(species_data, thresholds) {
  
  # Dynamic threshold configuration - exclusive approach
  final_thresholds <- list()
  
  # Extract thresholds from the dynamic configuration structure
  if (!is.null(thresholds$thresholds) && is.list(thresholds$thresholds)) {
    log_message("Using dynamic threshold configuration")
    final_thresholds <- thresholds$thresholds
    
    # Set default threshold for any missing variant types (most permissive = 1)
    unique_variant_types <- unique(species_data$var_type)
    for (vtype in unique_variant_types) {
      if (is.null(final_thresholds[[vtype]])) {
        final_thresholds[[vtype]] <- 1
        log_message(sprintf("Using default threshold (1) for variant type: %s", vtype), level = "WARN")
      }
    }
    
  } else {
    # Fallback: no dynamic configuration found, use default threshold for all types
    log_message("No dynamic threshold configuration found, using default threshold (1) for all variant types")
    unique_variant_types <- unique(species_data$var_type)
    for (vtype in unique_variant_types) {
      final_thresholds[[vtype]] <- 1
    }
  }
  
  # 1. Count samples per position/variant type
  position_sample_counts <- species_data %>%
    dplyr::group_by(position, var_type) %>%
    dplyr::summarise(sample_count = dplyr::n(), .groups = "drop")
  
  # 2. Identify positions that pass the frequency threshold
  # Dynamic threshold application: use configuration-driven thresholds
  passing_positions <- position_sample_counts %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      # Dynamic threshold lookup: get threshold for each variant type
      threshold = {
        if (!is.null(final_thresholds[[var_type]])) {
          final_thresholds[[var_type]]
        } else {
          # Default threshold for any missing variant types
          log_message(sprintf("No threshold found for variant type '%s', using default (1)", var_type), level = "WARN")
          1
        }
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(sample_count >= threshold) %>%
    dplyr::select(position, var_type)
  
  if (nrow(passing_positions) == 0) {
    return(NULL)
  }
  
  # 3. Select all original records for the passing positions
  filtered_variants <- species_data %>%
    dplyr::inner_join(passing_positions, by = c("position", "var_type"))
  
  # 4. Apply deduplication (one record per species/position/var_type)
  deduplicated_data <- filtered_variants %>%
    dplyr::group_by(species, position, var_type) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()
  
  return(deduplicated_data)
}

#' Generate simplified filtering statistics (Refactored for Dependency Injection)
#'
#' @param original_data The data before filtering.
#' @param filtered_data The data after filtering.
#' @param config Configuration list.
#' @param species_genome_data The data frame with species genome region info.
#' @return A list of statistics data frames.
generate_filtering_statistics <- function(original_data, filtered_data, config, species_genome_data) {
  log_message("Generating enhanced filtering statistics with P01 baseline data")
  
  session_id <- config$session_info$session_id
  
  # Try to load P01 preprocessing statistics as baseline
  p01_stats_full <- NULL
  p01_stats_ira <- NULL
  
  if (!is.null(session_id)) {
    p01_dir <- get_processing_stage_path(session_id, "P01_preprocessed", create_dir = FALSE)
    
    # Load full genome stats
    full_stats_file <- file.path(p01_dir, "preprocessing_stats_full_genome.csv")
    if (file.exists(full_stats_file)) {
      p01_stats_full <- read.csv(full_stats_file)
      log_message("Loaded P01 full genome statistics for baseline comparison")
    }
    
    # Load IRA-only stats
    ira_stats_file <- file.path(p01_dir, "preprocessing_stats_ira_only.csv")
    if (file.exists(ira_stats_file)) {
      p01_stats_ira <- read.csv(ira_stats_file)
      log_message("Loaded P01 IRA-only statistics for baseline comparison")
    }
  }
  
  # Calculate basic current statistics
  species_stats_orig <- original_data %>% 
    dplyr::group_by(species) %>% 
    dplyr::summarise(unfiltered_variants = dplyr::n(), .groups = "drop")
  species_stats_filt <- filtered_data %>% 
    dplyr::group_by(species) %>% 
    dplyr::summarise(filtered_variants = dplyr::n(), .groups = "drop")
  
  # Calculate sample counts per species
  species_sample_counts <- original_data %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(sample_count = dplyr::n_distinct(sample_id), .groups = "drop")
  
  # Calculate detailed filtered statistics by region type and variant type
  filtered_region_stats <- filtered_data %>%
    dplyr::group_by(species, region_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = region_type, values_from = count, values_fill = 0, names_prefix = "filtered_")
  
  # V3.16.13: After preprocessing standardization, RNA types should already be "RNA"
  # No need to aggregate RNA subtypes - they should not exist at this stage
  # This code block is kept for validation but should be a no-op in normal operation
  if (nrow(filtered_region_stats) > 0) {
    # Check if any legacy RNA subtypes exist (they shouldn't after preprocessing)
    legacy_rna_cols <- c("filtered_rRNA", "filtered_tRNA-CDS", "filtered_tRNA-intron")
    existing_legacy_cols <- intersect(legacy_rna_cols, names(filtered_region_stats))
    
    if (length(existing_legacy_cols) > 0) {
      log_message(sprintf("[WARNING] Found legacy RNA subtype columns in filtered data: %s", 
                         paste(existing_legacy_cols, collapse = ", ")), level = "warning")
      log_message("These should have been standardized to 'RNA' in preprocessing", level = "warning")
      
      # Aggregate for backward compatibility, but this indicates a preprocessing issue
    filtered_region_stats <- filtered_region_stats %>%
      dplyr::mutate(
          filtered_RNA = rowSums(dplyr::select(., dplyr::any_of(legacy_rna_cols)), na.rm = TRUE) +
                         dplyr::if_else("filtered_RNA" %in% names(.), filtered_RNA, 0)
      ) %>%
        dplyr::select(-dplyr::any_of(legacy_rna_cols))
    }
  }
  
  filtered_vartype_stats <- filtered_data %>%
    dplyr::group_by(species, var_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = var_type, values_from = count, values_fill = 0, names_prefix = "filtered_")

  # NEW: Calculate genome region statistics (LSC/IRA/SSC) for comprehensive reporting
  filtered_genome_stats <- filtered_data %>%
    dplyr::group_by(species, genome_region) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = genome_region, values_from = count, values_fill = 0, names_prefix = "filtered_")
  
  # Note: Keep the correct column names as created by pivot_wider (filtered_IRA with uppercase A)
  
  # Create comprehensive species statistics
  species_stats <- species_stats_orig %>%
    dplyr::left_join(species_stats_filt, by = "species") %>%
    dplyr::left_join(species_sample_counts, by = "species") %>%
    dplyr::mutate(
      unfiltered_variants = ifelse(is.na(unfiltered_variants), 0, unfiltered_variants),
      filtered_variants_IRA_only = ifelse(is.na(filtered_variants), 0, filtered_variants),
      sample_count = ifelse(is.na(sample_count), 0, sample_count)
    ) %>%
    dplyr::select(-dplyr::any_of("filtered_variants"))  # Remove the old column name
  
  # Add P01 baseline data if available
  if (!is.null(p01_stats_full)) {
    # Merge with full genome baseline
    p01_baseline <- p01_stats_full %>%
      dplyr::select(species, total_variants) %>%
      dplyr::rename(original_variants = total_variants)
    
    species_stats <- species_stats %>%
      dplyr::left_join(p01_baseline, by = "species") %>%
      dplyr::mutate(original_variants = ifelse(is.na(original_variants), unfiltered_variants, original_variants))
  } else {
    species_stats$original_variants <- species_stats$unfiltered_variants
  }
  
  # Add IRA-only baseline if available
  if (!is.null(p01_stats_ira)) {
    ira_baseline <- p01_stats_ira %>%
      dplyr::select(species, total_variants) %>%
      dplyr::rename(unfiltered_variants_ira_only = total_variants)
    
    species_stats <- species_stats %>%
      dplyr::left_join(ira_baseline, by = "species")
  } else {
    species_stats$unfiltered_variants_ira_only <- 0
  }
  
  # Add full genome baseline if available (this will be the deduplicated count before filtering)
  if (!is.null(p01_stats_full)) {
    full_baseline <- p01_stats_full %>%
      dplyr::select(species, total_variants) %>%
      dplyr::rename(deduplicated_variants = total_variants)
    
    species_stats <- species_stats %>%
      dplyr::left_join(full_baseline, by = "species") %>%
      dplyr::mutate(deduplicated_variants = ifelse(is.na(deduplicated_variants), unfiltered_variants, deduplicated_variants))
  } else {
    # If no P01 full genome stats, assume deduplicated equals unfiltered for now
    species_stats$deduplicated_variants <- species_stats$unfiltered_variants
  }
  
  # Add detailed region and variant type statistics
  if (nrow(filtered_region_stats) > 0) {
    species_stats <- species_stats %>%
      dplyr::left_join(filtered_region_stats, by = "species")
  }
  
  if (nrow(filtered_vartype_stats) > 0) {
    species_stats <- species_stats %>%
      dplyr::left_join(filtered_vartype_stats, by = "species")
  }
  
  # NEW: Add genome region statistics and calculate theoretical full genome count
  if (nrow(filtered_genome_stats) > 0) {
    species_stats <- species_stats %>%
      dplyr::left_join(filtered_genome_stats, by = "species")
  }
  
  # --- COLUMN UNIFICATION: Merge 'filtered_IR' into 'filtered_IRA' ---
  if ("filtered_IR" %in% names(species_stats)) {
    if (!"filtered_IRA" %in% names(species_stats)) {
      species_stats$filtered_IRA <- 0
    }
    species_stats$filtered_IRA <- species_stats$filtered_IRA + species_stats$filtered_IR
    species_stats <- species_stats %>% dplyr::select(-filtered_IR)
    log_message("Unified filtered_IR into filtered_IRA column")
  }
  
  # --- ROBUSTNESS: Ensure all genome region columns exist ---
  required_region_cols <- c("filtered_LSC", "filtered_IRA", "filtered_SSC")
  for (col in required_region_cols) {
    if (!col %in% names(species_stats)) {
      species_stats[[col]] <- 0
    }
  }

  # --- Handle potential NA values from join operations ---
  species_stats <- species_stats %>%
    dplyr::mutate(
      filtered_LSC = ifelse(is.na(filtered_LSC), 0, filtered_LSC),
      filtered_IRA = ifelse(is.na(filtered_IRA), 0, filtered_IRA),
      filtered_SSC = ifelse(is.na(filtered_SSC), 0, filtered_SSC)
    )
  
  # --- SPECIAL CASE HANDLING: Process IR-lacking species correctly ---
  if ("filtered_whole_genome" %in% names(species_stats)) {
    # Method A: Precise identification using injected species_genome_data
    species_handling_info <- species_genome_data %>%
      dplyr::select(species, special_handling)
    
    species_stats <- species_stats %>%
      dplyr::left_join(species_handling_info, by = "species") %>%
      dplyr::mutate(
        filtered_variants_full_genome = dplyr::case_when(
          special_handling == "IR_lacking_genome" & !is.na(filtered_whole_genome) ~ filtered_whole_genome,
          TRUE ~ filtered_LSC + (filtered_IRA * 2) + filtered_SSC
        ),
        # Calculate overall retention rate based on IRA-only baseline (GOLDEN RULE compliance)
        # The denominator MUST be unfiltered_variants_ira_only to ensure statistical consistency
        overall_retention_rate = ifelse(unfiltered_variants_ira_only > 0, filtered_variants_IRA_only / unfiltered_variants_ira_only, 0)
      ) %>%
      dplyr::select(-dplyr::any_of(c("filtered_whole_genome", "special_handling")))
      
    log_message("Applied specialized handling for IR-lacking species")
  } else {
    # Standard calculation: For cases without IR-lacking species
    species_stats <- species_stats %>%
      dplyr::mutate(
        filtered_variants_full_genome = filtered_LSC + (filtered_IRA * 2) + filtered_SSC,
        # Calculate overall retention rate based on IRA-only baseline (GOLDEN RULE compliance)  
        # The denominator MUST be unfiltered_variants_ira_only to ensure statistical consistency
        overall_retention_rate = ifelse(unfiltered_variants_ira_only > 0, filtered_variants_IRA_only / unfiltered_variants_ira_only, 0)
      )
  }
  
  # Ensure all required columns exist (expanded per user requirements)
  required_columns <- c("species", "sample_count", "original_variants", "unfiltered_variants_ira_only", 
                       "filtered_variants_full_genome", "filtered_variants_IRA_only", "overall_retention_rate")
  
  for (col in required_columns) {
    if (!col %in% names(species_stats)) {
      species_stats[[col]] <- 0
    }
  }
  
  # Select and order columns in precise sequence per user requirements
  species_stats <- species_stats %>%
    dplyr::select(
      # Core identification (with sample count)
      species,
      sample_count,  # NEW: Insert after species
      original_variants,
      unfiltered_variants_ira_only,
      # Enhanced filtering statistics
      filtered_variants_full_genome,
      filtered_variants_IRA_only, 
      overall_retention_rate,
      # Genome region breakdown
      dplyr::any_of(c("filtered_LSC", "filtered_IRA", "filtered_SSC")),
      # Region type breakdown (CDS, IGS, intron, etc.)
      dplyr::starts_with("filtered_")
    ) %>%
    dplyr::arrange(species)
  
  log_message(sprintf("Generated enhanced filtering statistics for %d species", nrow(species_stats)))
  
  return(list(
    by_species = species_stats
  ))
}

#' Save simplified filtering results to the P02 directory
#'
#' @param filtered_data The filtered data frame to save.
#' @param statistics The statistics list to save.
#' @param plots (Currently unused) placeholder for future plots.
#' @param config The application configuration.
save_filtering_results <- function(filtered_data, statistics, plots, config) {
  session_id <- config$session_info$session_id
  if (is.null(session_id) || session_id == "") {
    stop("FATAL: session_id is missing, cannot save filtering results.")
  }
  
  # Use the new path manager to get the correct P02 output directory
  output_dir <- get_processing_stage_path(session_id, "P02_filtered", create_dir = TRUE)
  log_message(sprintf("Saving filtering results to: %s", output_dir))
  
  # Save main filtered data
  output_file_filtered <- file.path(output_dir, "variants_filtered.csv")
  utils::write.csv(filtered_data, output_file_filtered, row.names = FALSE)
  
  # Save statistics (overall_stats removed - redundant information)
  safe_write_file(statistics$by_species, file.path(output_dir, "filtering_species_stats.csv"))
  
  # Placeholder for future plot saving logic
  # if (!is.null(plots)) { ... }
  
  log_message("Filtering results saved successfully.")
}

# log_message("Dual-branch variant filtering functions loaded successfully")