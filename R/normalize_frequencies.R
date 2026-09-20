############################################################
#### normalize_frequencies.R - Frequency Normalization Functions ####
############################################################

#' @importFrom magrittr %>%
NULL
#
# Functions for calculating and normalizing variant frequencies
# Handles CDS, IGS, and intron regions with proper length normalization
# Designed for Shiny compatibility with flexible parameters
#
############################################################

# Load required libraries

# Source required configurations
# source("config/default_config.R")
# source("src/utils/common_utils.R")

#' Normalize variant frequencies across genomic regions
#' 
#' @title Normalize variant frequencies across genomic regions
#' @description Calculates per-kilobase variant frequencies for genomic regions
#' @param variant_data Filtered or unfiltered variant data
#' @param gene_lengths_df Gene lengths data frame (preloaded)
#' @param region_info_df Region information data frame (preloaded)
#' @param config Configuration list
#' @param progress_callback Optional function for progress updates
#' @return List containing normalized frequency data and metadata
#' @export
normalize_frequencies <- function(variant_data, 
                                 gene_lengths_df,
                                 region_info_df,
                                 config,
                                 progress_callback = NULL) {
  
  # Initialize
  if (!is.null(progress_callback)) {
    progress_callback("Starting frequency normalization", 0.1)
  }
  
  log_message("Starting frequency normalization")
  
  # Ensure region_type is character to prevent factor level errors
  if (is.factor(variant_data$region_type)) {
    variant_data$region_type <- as.character(variant_data$region_type)
  }
  
  # ================== DATA TRANSFORMATION TRACKING ==================
  # Track original data structure and counts
  original_variant_count <- nrow(variant_data)
  original_species_count <- length(unique(variant_data$species))
  original_var_types <- unique(variant_data$var_type)
  # Use standardized column names (data should already be mapped)
  var_location_col <- "region_type"
  original_region_types <- unique(variant_data[[var_location_col]][!is.na(variant_data[[var_location_col]])])
  
  log_message("=== NORMALIZATION DATA TRANSFORMATION TRACKING ===", level = "debug")
  log_message(sprintf("ORIGINAL DATA: %d variants, %d species, %d variant types (%s)", 
                     original_variant_count, original_species_count, 
                     length(original_var_types), paste(original_var_types, collapse=", ")), level = "debug")
  log_message(sprintf("REGION TYPES: %s", paste(original_region_types, collapse=", ")), level = "debug")
  
  # Count variants by region type
  region_counts <- table(variant_data[[var_location_col]], useNA = "ifany")
  for (region in names(region_counts)) {
    log_message(sprintf("  - %s variants: %d", ifelse(is.na(region), "NA/Unknown", region), region_counts[region]), level = "debug")
  }
  
  # Validate input data frames
  if (is.null(gene_lengths_df) || nrow(gene_lengths_df) == 0) {
    stop("Gene lengths data frame is required but not provided or empty")
  }
  
  if (is.null(region_info_df)) {
    log_message("Region info data frame is NULL, some region types may not be processed", level = "warning")
  }
  
  gene_lengths <- gene_lengths_df
  
  # CRITICAL FIX: Structure region_info_df into the expected format for IGS/intron functions
  region_info <- NULL
  if (!is.null(region_info_df) && nrow(region_info_df) > 0) {
    log_message("Structuring region info data for IGS/intron processing", level = "debug")
    
    # Ensure required columns exist (these should be mapped by the caller)
    required_region_cols <- c("species", "region_name", "region_type", "region_start", "region_end", "region_length")
    missing_cols <- setdiff(required_region_cols, names(region_info_df))
    
    if (length(missing_cols) > 0) {
      log_message(sprintf("Missing required columns in region_info_df: %s", paste(missing_cols, collapse=", ")), level = "warning")
      log_message("Available columns: %s", paste(names(region_info_df), collapse=", "), level = "debug")
    } else {
      # Clean and structure the data
      # V3.16.12: No longer apply tolower() - preserve standard case from preprocessing
      region_info_clean <- region_info_df %>%
        dplyr::mutate(
          region_start = as.numeric(region_start),
          region_end = as.numeric(region_end),
          region_length = as.numeric(region_length)
        ) %>%
        dplyr::filter(!is.na(region_start) & !is.na(region_end) & !is.na(region_length))
      
      # Create structured list with separate components
      # V3.16.12: Use only "IGS" - no backward compatibility needed
      region_info <- list(
        all_regions = region_info_clean,
        intergenic = region_info_clean %>% dplyr::filter(.data$region_type == "IGS"),
        intron = region_info_clean %>% dplyr::filter(.data$region_type == "intron")
      )
      
      log_message(sprintf("Region info structured: %d total regions (%d intergenic, %d intron)", 
                         nrow(region_info$all_regions), 
                         nrow(region_info$intergenic), 
                         nrow(region_info$intron)), level = "debug")
    }
  }
  
  # Validate input data
  required_columns <- c("species", "position", "var_type", "genome_region")
  if (!validate_data_frame(variant_data, required_columns, 1, "variant data")) {
    stop("Invalid variant data structure for normalization")
  }
  
  if (!is.null(progress_callback)) {
    progress_callback("Calculating CDS frequencies", 0.3)
  }
  
  # Calculate CDS frequencies with region_info for exon-level processing
  cds_frequencies <- calculate_cds_frequencies(variant_data, gene_lengths, config, region_info)
  
  if (!is.null(progress_callback)) {
    progress_callback("Calculating IGS frequencies", 0.6)
  }
  
  # Calculate IGS frequencies
  igs_frequencies <- calculate_igs_frequencies(variant_data, region_info, config)
  
  if (!is.null(progress_callback)) {
    progress_callback("Calculating intron frequencies", 0.8)
  }
  
  # Calculate intron frequencies
  intron_frequencies <- calculate_intron_frequencies(variant_data, region_info, config)
  
  # EMERGENCY HOTFIX PHASE 1: Process unhandled variants (RNA-related types)
  # This implements "pass-through" logic for variants not processed by calculate_* functions
  if (!is.null(progress_callback)) {
    progress_callback("Processing unhandled variants (RNA types)", 0.85)
  }
  
  log_message("=== EMERGENCY HOTFIX: PROCESSING UNHANDLED VARIANTS ===")
  
  # Identify variants not processed by the three calculate functions
  # V3.16.12: Use only standard terms after preprocessing standardization
  processed_region_types <- c("CDS", "IGS", "intron", "others")
  unprocessed_variants <- variant_data %>%
    dplyr::filter(!.data$region_type %in% processed_region_types)
  
  log_message(sprintf("Found %d unprocessed variants with region types: %s", 
                     nrow(unprocessed_variants),
                     paste(unique(unprocessed_variants$region_type), collapse = ", ")))
  
  # Create pass-through data frame for unprocessed variants
  other_variants_df <- data.frame()
  
  if (nrow(unprocessed_variants) > 0) {
    # Group by species, var_type, region_type, and gene to get variant counts
    # Add genome_region and region_name for consistency
    other_variants_df <- unprocessed_variants %>%
      dplyr::group_by(species, var_type, region_type, gene) %>%
      dplyr::summarise(
        variant_count = dplyr::n(),  # Count raw variants in each group
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        genome_region = NA_character_,  # RNA variants don't have genome_region
        region_name = gene,              # For RNA, region_name is the same as gene
        region_length = NA_real_,        # Set to NA as specified
        frequency_per_kb = NA_real_      # Set to NA as specified
      ) %>%
      # Ensure column order matches calculate_* function outputs
      dplyr::select(.data$species, .data$var_type, .data$region_type, .data$genome_region, .data$region_name, .data$gene, .data$variant_count, .data$region_length, .data$frequency_per_kb)
    
    log_message(sprintf("Created pass-through data frame: %d groups from %d raw variants", 
                       nrow(other_variants_df), nrow(unprocessed_variants)), level = "debug")
    
    # Log summary by region type
    region_summary <- other_variants_df %>%
      dplyr::group_by(region_type) %>%
      dplyr::summarise(
        groups = dplyr::n(),
        total_variants = sum(variant_count),
        .groups = "drop"
      )
    
    log_message("UNPROCESSED VARIANTS SUMMARY:", level = "debug")
    for (i in 1:nrow(region_summary)) {
      log_message(sprintf("  %s: %d groups, %d variants", 
                         region_summary$region_type[i],
                         region_summary$groups[i],
                         region_summary$total_variants[i]), level = "debug")
    }
  } else {
    log_message("No unprocessed variants found - creating empty pass-through data frame", level = "debug")
  }
  
  # Combine all frequencies including unprocessed variants
  if (!is.null(progress_callback)) {
    progress_callback("Combining frequency data", 0.9)
  }
  
  combined_frequencies <- combine_frequency_data(cds_frequencies, igs_frequencies, intron_frequencies, other_variants_df)
  
  # Generate summary frequencies table (all variant types combined per gene/region)
  # This matches the format of the old script output for comparison purposes
  if (!is.null(progress_callback)) {
    progress_callback("Generating summary frequencies table", 0.95)
  }
  
  log_message("Generating summary frequencies table (all variant types combined)", level = "debug")
  
  # Include region_name in grouping so IGS identity is preserved
  # Previous logic only grouped by gene, causing NA genes (IGS) to be merged incorrectly
  summary_frequencies <- combined_frequencies %>%
    dplyr::group_by(species, region_type, region_name, gene, region_length) %>%
    dplyr::summarise(
      total_variant_count = sum(variant_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      var_type = "all_types",
      variant_count = total_variant_count,
      frequency_per_kb = (total_variant_count / region_length) * config$normalization$frequency_unit
    ) %>%
    dplyr::select(.data$species, .data$var_type, .data$region_type, .data$region_name, .data$gene, .data$variant_count, .data$region_length, .data$frequency_per_kb) %>%
    dplyr::filter(.data$variant_count > 0)  # Only keep genes/regions with variants
  
  log_message(sprintf("Generated summary frequencies: %d entries (genes/regions with variants)", nrow(summary_frequencies)))
  
  # Generate normalization statistics
  normalization_stats <- generate_normalization_statistics(combined_frequencies, variant_data, config)
  
  # Perform region coverage validation
  # This validation is performed AFTER normalization to check for any missing variants
  # between P02 (filtered variants) and P03 (normalized frequencies)
  coverage_validation_results <- NULL
  
  if (!is.null(progress_callback)) {
    progress_callback("Performing coverage validation", 0.95)
  }
  
  log_message("=== Starting region coverage validation ===")
  
  # Coverage validation requires access to:
  # 1. region_info (from region_info_df parameter)
  # 2. processed_genome_regions.csv (need to load)
  # 3. variants_filtered.csv (original variant_data before normalization)
  
  # Note: Coverage validation will be performed in the task_dispatcher
  # where all required data sources are available
  # Here we just log a placeholder
  log_message("Coverage validation will be performed after normalization completes")
  
  if (!is.null(progress_callback)) {
    progress_callback("Normalization completed", 1.0)
  }
  
  log_message(sprintf("Frequency normalization completed"))
  log_message(sprintf("Generated frequencies: %d entries", nrow(combined_frequencies)))
  
  # Enhanced metadata tracking with comprehensive transformation details
  detailed_variants <- sum(combined_frequencies$variant_count)
  summary_variants <- sum(summary_frequencies$variant_count)
  zero_entries <- sum(combined_frequencies$variant_count == 0)
  entries_with_variants <- sum(combined_frequencies$variant_count > 0)
  
  # Calculate region-wise statistics
  region_breakdown <- combined_frequencies %>%
    dplyr::group_by(region_type) %>%
    dplyr::summarise(
      entries = dplyr::n(),
      variants = sum(variant_count),
      entries_with_variants = sum(variant_count > 0),
      .groups = "drop"
    )
  
  # Calculate species and variant type counts
  species_analyzed <- length(unique(combined_frequencies$species))
  variant_types_analyzed <- length(unique(combined_frequencies$var_type))
  genes_regions_analyzed <- length(unique(combined_frequencies$gene))
  
  # Prepare results object
  results <- list(
    frequencies = combined_frequencies,
    frequencies_summary = summary_frequencies,
    statistics = normalization_stats,
    metadata = list(
      # Basic information
      normalization_time = Sys.time(),
      config_version = config$version,
      frequency_unit = config$normalization$frequency_unit,
      
      # Input data tracking
      original_input_variants = nrow(variant_data),
      original_species_count = length(unique(variant_data$species)),
      original_variant_types = length(unique(variant_data$var_type)),
      
      # Matrix expansion tracking
      detailed_entries = nrow(combined_frequencies),
      summary_entries = nrow(summary_frequencies),
      matrix_expansion_entries = nrow(combined_frequencies) - nrow(variant_data),
      
      # Variant preservation tracking
      detailed_variants_preserved = detailed_variants,
      summary_variants_preserved = summary_variants,
      variants_lost = nrow(variant_data) - detailed_variants,
      
      # Analysis structure
      species_analyzed = species_analyzed,
      variant_types_analyzed = variant_types_analyzed,
      genes_regions_analyzed = genes_regions_analyzed,
      
      # Matrix composition
      entries_with_variants = entries_with_variants,
      entries_with_zero_variants = zero_entries,
      zero_entry_percentage = round((zero_entries / nrow(combined_frequencies)) * 100, 2),
      
      # Region breakdown
      region_breakdown = region_breakdown,
      
      # Transformation explanations
      transformation_summary = list(
        matrix_expansion_reason = "Complete analytical matrix includes all species x gene x variant_type combinations",
        variant_preservation_status = if(detailed_variants == summary_variants) "All variants preserved" else "Some variants lost",
        zero_entries_purpose = "Matrix placeholders for proper normalization and comparative analysis",
        data_integrity_check = if(detailed_variants == summary_variants) "PASSED" else "FAILED"
      )
    )
  )
  
  # Save normalization results to P03 directory
  save_normalization_results(results, config)
  
  return(results)
}

#' Load gene lengths data
#' 
#' @param gene_lengths_file Path to gene lengths file
#' @return Data frame with gene lengths
load_gene_lengths <- function(gene_lengths_file) {
  
  log_message(sprintf("Loading gene lengths from: %s", gene_lengths_file))
  
  gene_lengths <- safe_read_file(gene_lengths_file)
  
  if (is.null(gene_lengths)) {
    stop("Failed to load gene lengths data")
  }
  
  # Validate structure
  required_columns <- c("Species", "Gene", "CDS_length")
  if (!validate_data_frame(gene_lengths, required_columns, 1, "gene lengths")) {
    stop("Invalid gene lengths data structure")
  }
  
  # Standardize column names
  gene_lengths <- gene_lengths %>%
    dplyr::rename(
      species = Species,
      gene = Gene,
      cds_length = CDS_length
    )
  
  # Clean species names
  gene_lengths$species <- clean_species_names(gene_lengths$species)
  
  log_message(sprintf("Gene lengths loaded: %d entries for %d species", 
                     nrow(gene_lengths), length(unique(gene_lengths$species))))
  
  return(gene_lengths)
}

#' Load region information data (intergenic and intron regions)
#' 
#' @param region_info_file Path to region info file
#' @param config Configuration list
#' @return Data frame with region information including coordinates
load_region_info <- function(region_info_file, config) {
  
  log_message(sprintf("Loading region information from: %s", region_info_file))
  
  # Try primary file first
  region_info <- safe_read_file(region_info_file, read.table, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  
  # If primary file fails, try fallback
  if (is.null(region_info)) {
    fallback_file <- config$data_paths$region_info_fallback
    log_message(sprintf("Trying fallback region info file: %s", fallback_file))
    region_info <- safe_read_file(fallback_file, read.table, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  }
  
  if (is.null(region_info)) {
    log_message("No region information file available, using default lengths", level = "warning")
    return(NULL)
  }
  
  # Validate expected columns
  expected_columns <- c("Species", "Region_Name", "Type", "Start_Position", "End_Position", "Length")
  if (!all(expected_columns %in% names(region_info))) {
    missing_cols <- setdiff(expected_columns, names(region_info))
    log_message(sprintf("Missing columns in region info file: %s", paste(missing_cols, collapse = ", ")), level = "error")
    return(NULL)
  }
  
  # Standardize column names and clean data
  region_info <- region_info %>%
    dplyr::rename(
      species = Species,
      region_name = Region_Name,
      region_type = Type,
      start_position = Start_Position,
      end_position = End_Position,
      region_length = Length
    ) %>%
    dplyr::mutate(
      species = clean_species_names(species),
      region_type = tolower(region_type),
      start_position = as.numeric(start_position),
      end_position = as.numeric(end_position),
      region_length = as.numeric(region_length)
    ) %>%
    dplyr::filter(!is.na(start_position) & !is.na(end_position) & !is.na(region_length))
  
  # V3.16.12: Extract IGS regions using standard terminology
  intergenic_regions <- region_info %>%
    dplyr::filter(.data$region_type == "IGS")
  
  intron_regions <- region_info %>%
    dplyr::filter(.data$region_type == "intron")
  
  log_message(sprintf("Region information loaded: %d entries (%d intergenic, %d intron)", 
                     nrow(region_info), nrow(intergenic_regions), nrow(intron_regions)))
  
  # Return as list with separate components
  return(list(
    all_regions = region_info,
    intergenic = intergenic_regions,
    intron = intron_regions
  ))
}

#' Calculate CDS frequencies (coordinate-based exon-level normalization)
#' 
#' Coordinate matching is used so CDS exons are processed like IGS/intron intervals.
#' This aligns CDS processing with IGS/intron methodology for consistency
#' 
#' @param variant_data Variant data with annotation and genome_region
#' @param gene_lengths Gene lengths data (legacy parameter, kept for compatibility)
#' @param config Configuration list
#' @param region_info Region information containing CDS exon coordinates
#' @return Data frame with CDS frequencies by exon with genome_region
calculate_cds_frequencies <- function(variant_data, gene_lengths, config, region_info = NULL) {
  
  log_message("Calculating CDS frequencies using coordinate-based exon-level normalization")
  
  # Filter for CDS variants
  var_location_col <- "region_type"
  
  cds_variants <- variant_data %>%
    dplyr::filter(!is.na(!!rlang::sym(var_location_col)) & 
                  !!rlang::sym(var_location_col) == "CDS") %>%
    dplyr::select(species, var_type, position, genome_region) %>%
    dplyr::mutate(position = as.numeric(position)) %>%
    dplyr::filter(!is.na(position))
  
  if (nrow(cds_variants) == 0) {
    log_message("No CDS variants found in the data")
    return(data.frame())
  }
  
  log_message(sprintf("Found %d CDS variants for exon-level analysis", nrow(cds_variants)))
  
  # Get CDS regions from region_info (exon-level)
  if (is.null(region_info) || is.null(region_info$all_regions)) {
    log_message("WARNING: No region_info available, falling back to legacy gene-level calculation", level = "warning")
    # Fallback to legacy logic (not recommended)
    return(calculate_cds_frequencies_legacy(variant_data, gene_lengths, config))
  }
  
  cds_regions <- region_info$all_regions %>%
    dplyr::filter(region_type == "CDS")
  
  if (nrow(cds_regions) == 0) {
    log_message("WARNING: No CDS regions found in region_info", level = "warning")
    return(data.frame())
  }
  
  log_message(sprintf("Processing %d CDS regions (exon-level)", nrow(cds_regions)))
  
  # Validate CDS regions have genome_region field
  if (!"genome_region" %in% names(cds_regions)) {
    log_message("WARNING: genome_region field missing in CDS regions, adding 'unknown'", level = "warning")
    cds_regions$genome_region <- "unknown"
  }
  
  # Initialize result data frame
  cds_frequencies <- data.frame()
  
  # For each species, match variants to CDS regions (exons)
  for (species_name in unique(cds_variants$species)) {
    species_variants <- cds_variants %>%
      dplyr::filter(species == species_name)
    
    species_regions <- cds_regions %>%
      dplyr::filter(species == species_name)
    
    if (nrow(species_regions) == 0) {
      log_message(sprintf("No CDS regions found for species: %s", species_name), level = "warning")
      next
    }
    
    # For each variant type
    for (var_type_name in unique(species_variants$var_type)) {
      var_type_variants <- species_variants %>%
        dplyr::filter(var_type == var_type_name)
      
      # For each CDS region (exon)
      for (i in 1:nrow(species_regions)) {
        region_info_row <- species_regions[i, ]
        
        # Find variants within this exon's coordinates
        variants_in_region <- var_type_variants %>%
          dplyr::filter(position >= region_info_row$region_start & 
                       position <= region_info_row$region_end)
        
        variant_count <- nrow(variants_in_region)
        
        # Calculate normalized frequency
        normalized_freq <- (variant_count / region_info_row$region_length) * config$normalization$frequency_unit
        
        # Add to results
        # Carry forward region_name for downstream P03 tables
        # Provide both region_name (exon-level) and gene (gene-level) for downstream flexibility
        cds_frequencies <- rbind(cds_frequencies, data.frame(
          species = species_name,
          var_type = var_type_name,
          region_type = "CDS",
          genome_region = region_info_row$genome_region,
          region_name = region_info_row$region_name,  # Exon-level identifier (e.g., atpF_exon1)
          gene = region_info_row$gene,                # Gene-level identifier (e.g., atpF)
          variant_count = variant_count,
          region_length = region_info_row$region_length,
          frequency_per_kb = normalized_freq,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  log_message(sprintf("Generated CDS frequencies for %d exon-species-variant_type combinations", nrow(cds_frequencies)))
  
  # Log summary statistics
  if (nrow(cds_frequencies) > 0) {
    exon_summary <- cds_frequencies %>%
      dplyr::group_by(species, gene) %>%
      dplyr::summarise(total_variants = sum(variant_count), .groups = "drop") %>%
      dplyr::filter(total_variants > 0)
    
    if (nrow(exon_summary) > 0) {
      log_message(sprintf("CDS exons with variants: %d (max variants per exon: %d)", 
                         nrow(exon_summary), max(exon_summary$total_variants)))
    }
    
    # Log genome_region distribution
    if ("genome_region" %in% names(cds_frequencies)) {
      genome_region_counts <- table(cds_frequencies$genome_region)
      log_message(sprintf("Genome region distribution: %s", 
                         paste(names(genome_region_counts), "=", genome_region_counts, collapse = ", ")))
    }
  }
  
  return(cds_frequencies)
}

#' Legacy CDS frequency calculation (gene-level)
#' 
#' @description Fallback function for backward compatibility
#' @keywords internal
calculate_cds_frequencies_legacy <- function(variant_data, gene_lengths, config) {
  log_message("Using legacy gene-level CDS frequency calculation")
  
  # Original gene-level logic preserved for fallback
  var_location_col <- "region_type"
  gene_name_col <- "gene"
  
  cds_variants <- variant_data %>%
    dplyr::filter(!is.na(!!rlang::sym(var_location_col)) & 
                  !!rlang::sym(var_location_col) == "CDS" & 
                  !is.na(!!rlang::sym(gene_name_col))) %>%
    dplyr::select(.data$species, .data$var_type, .data$position, !!rlang::sym(gene_name_col)) %>%
    dplyr::rename(gene = !!rlang::sym(gene_name_col))
  
  if (nrow(cds_variants) == 0) {
    return(data.frame())
  }
  
  cds_variant_counts <- cds_variants %>%
    dplyr::group_by(species, gene, var_type) %>%
    dplyr::summarise(variant_count = dplyr::n(), .groups = "drop")
  
  gene_lengths <- gene_lengths %>%
    dplyr::rename(gene = gene_name, cds_length = gene_length)
  
  all_species <- unique(gene_lengths$species)
  all_genes <- unique(gene_lengths$gene)
  all_var_types <- unique(cds_variants$var_type)
  
  complete_matrix <- expand.grid(
    species = all_species,
    gene = all_genes,
    var_type = all_var_types,
    stringsAsFactors = FALSE
  )
  
  cds_frequencies <- complete_matrix %>%
    dplyr::left_join(cds_variant_counts, by = c("species", "gene", "var_type")) %>%
    dplyr::mutate(variant_count = ifelse(is.na(variant_count), 0, variant_count)) %>%
    dplyr::left_join(gene_lengths, by = c("species", "gene")) %>%
    dplyr::filter(!is.na(cds_length) & cds_length > 0) %>%
    dplyr::mutate(
      region_type = "CDS",
      genome_region = "unknown",  # Legacy mode doesn't have genome_region
      region_length = cds_length,
      frequency_per_kb = (variant_count / cds_length) * config$normalization$frequency_unit
    ) %>%
    dplyr::select(species, var_type, region_type, genome_region, gene, variant_count, region_length, frequency_per_kb)
  
  return(cds_frequencies)
}

#' Calculate IGS frequencies (region-level normalization with coordinate matching)
#' 
#' @param variant_data Variant data with position information
#' @param region_info Region information list containing intergenic regions
#' @param config Configuration list
#' @return Data frame with IGS frequencies by region
calculate_igs_frequencies <- function(variant_data, region_info, config) {
  
  log_message("Calculating IGS frequencies using coordinate-based region matching")
  
  # Check if region_info is available
  if (is.null(region_info) || is.null(region_info$intergenic)) {
    log_message("No intergenic region information available", level = "warning")
    return(data.frame())
  }
  
  # Filter for IGS variants using standardized column name
  # V3.16.12: Use only "IGS" standard term
  var_location_col <- "region_type"
  
  igs_variants <- variant_data %>%
    dplyr::filter(!is.na(!!rlang::sym(var_location_col)) & 
                  !!rlang::sym(var_location_col) == "IGS") %>%
    dplyr::select(species, var_type, position) %>%
    dplyr::mutate(position = as.numeric(position)) %>%
    dplyr::filter(!is.na(position))
  
  if (nrow(igs_variants) == 0) {
    log_message("No intergenic variants found in the data")
    return(data.frame())
  }
  
  log_message(sprintf("Found %d intergenic variants for region-level analysis", nrow(igs_variants)))
  
  # Get intergenic regions - STRICT DATA SOURCE PURITY ENFORCEMENT
  log_message("ENFORCING IGS DATA SOURCE PURITY: Only using intergenic region data for intergenic variants")
  intergenic_regions <- region_info$intergenic
  
  # Validate intergenic region data structure
  if (!all(c("species", "region_name", "region_length") %in% names(intergenic_regions))) {
    stop("Invalid intergenic regions data structure - missing required intergenic-specific columns")
  }
  
  # Validate that we're working with intergenic data only
  if ("region_type" %in% names(intergenic_regions)) {
    # V3.16.12: Use only "IGS" standard term
    non_intergenic_entries <- sum(intergenic_regions$region_type != "IGS", na.rm = TRUE)
    if (non_intergenic_entries > 0) {
      log_message(sprintf("WARNING: Found %d non-IGS entries in intergenic regions data", non_intergenic_entries), level = "warning")
      intergenic_regions <- intergenic_regions %>% dplyr::filter(region_type == "IGS")
    }
  }
  
  log_message("IGS DATA SOURCE VALIDATION: Confirmed intergenic regions data integrity")
  
  # V3.16.10: Debug logging for IGS region filtering
  log_message(sprintf("[DEBUG] intergenic_regions structure: %d rows, columns: %s", 
                     nrow(intergenic_regions), 
                     paste(names(intergenic_regions), collapse = ", ")), level = "debug")
  log_message(sprintf("[DEBUG] Unique species in intergenic_regions: %s", 
                     paste(unique(intergenic_regions$species), collapse = ", ")), level = "debug")
  log_message(sprintf("[DEBUG] Unique species in igs_variants: %s", 
                     paste(unique(igs_variants$species), collapse = ", ")), level = "debug")
  
  # Initialize result data frame
  igs_frequencies <- data.frame()
  
  # For each species, match variants to intergenic regions
  for (species_name in unique(igs_variants$species)) {
    species_variants <- igs_variants %>%
      dplyr::filter(species == species_name)
    
    species_regions <- intergenic_regions %>%
      dplyr::filter(species == species_name)
    
    if (nrow(species_regions) == 0) {
      log_message(sprintf("No intergenic regions found for species: %s", species_name), level = "warning")
      next
    }
    
    # For each variant type
    for (var_type_name in unique(species_variants$var_type)) {
      var_type_variants <- species_variants %>%
        dplyr::filter(var_type == var_type_name)
      
      # For each intergenic region
      for (i in 1:nrow(species_regions)) {
        region_info_row <- species_regions[i, ]
        
        # Find variants within this region's coordinates
        variants_in_region <- var_type_variants %>%
          dplyr::filter(position >= region_info_row$region_start & 
                       position <= region_info_row$region_end)
        
        variant_count <- nrow(variants_in_region)
        
        # Calculate normalized frequency
        normalized_freq <- (variant_count / region_info_row$region_length) * config$normalization$frequency_unit
        
        # Get genome_region from region_info
        genome_region <- if ("genome_region" %in% names(region_info_row)) {
          region_info_row$genome_region
        } else {
          "unknown"
        }
        
        # Add to results
        # Carry forward region_name for downstream P03 tables
        # IGS regions do not belong to any gene, so gene = NA
        igs_frequencies <- rbind(igs_frequencies, data.frame(
          species = species_name,
          var_type = var_type_name,
          region_type = "IGS",
          genome_region = genome_region,
          region_name = region_info_row$region_name,  # IGS identifier (e.g., IGS_psbA-matK)
          gene = NA,                                   # IGS does not belong to any gene
          variant_count = variant_count,
          region_length = region_info_row$region_length,
          frequency_per_kb = normalized_freq,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Remove rows with zero variants (optional - keep for complete matrix)
  # igs_frequencies <- igs_frequencies %>%dplyr::filter(variant_count > 0)
  
  log_message(sprintf("Generated IGS frequencies for %d region-species-variant_type combinations", nrow(igs_frequencies)))
  
  # Log summary statistics
  region_summary <- igs_frequencies %>%
    dplyr::group_by(species, gene) %>%
    dplyr::summarise(total_variants = sum(variant_count), .groups = "drop") %>%
    dplyr::filter(total_variants > 0)
  
  if (nrow(region_summary) > 0) {
    log_message(sprintf("IGS regions with variants: %d (max variants per region: %d)", 
                       nrow(region_summary), max(region_summary$total_variants)))
  }
  
  return(igs_frequencies)
}

#' Calculate intron frequencies (region-level normalization with coordinate matching)
#' 
#' @param variant_data Variant data with position information
#' @param region_info Region information list containing intron regions
#' @param config Configuration list
#' @return Data frame with intron frequencies by region
calculate_intron_frequencies <- function(variant_data, region_info, config) {
  
  log_message("Calculating intron frequencies using coordinate-based region matching")
  
  # Check if region_info is available
  if (is.null(region_info) || is.null(region_info$intron)) {
    log_message("No intron region information available", level = "warning")
    return(data.frame())
  }
  
  # Filter for intron variants using standardized column name
  var_location_col <- "region_type"
  
  intron_variants <- variant_data %>%
dplyr::filter(!is.na(!!rlang::sym(var_location_col)) & !!rlang::sym(var_location_col) == "intron") %>%
dplyr::select(species, var_type, position) %>%
dplyr::mutate(position = as.numeric(position)) %>%
dplyr::filter(!is.na(position))
  
  if (nrow(intron_variants) == 0) {
    log_message("No intron variants found in the data")
    return(data.frame())
  }
  
  log_message(sprintf("Found %d intron variants for region-level analysis", nrow(intron_variants)))
  
  # Get intron regions - STRICT DATA SOURCE PURITY ENFORCEMENT
  log_message("ENFORCING INTRON DATA SOURCE PURITY: Only using intron region data for intron variants")
  intron_regions <- region_info$intron
  
  # Validate intron region data structure
  if (!all(c("species", "region_name", "region_length") %in% names(intron_regions))) {
    stop("Invalid intron regions data structure - missing required intron-specific columns")
  }
  
  # Validate that we're working with intron data only
  if ("region_type" %in% names(intron_regions)) {
    non_intron_entries <- sum(intron_regions$region_type != "intron", na.rm = TRUE)
    if (non_intron_entries > 0) {
      log_message(sprintf("WARNING: Found %d non-intron entries in intron regions data", non_intron_entries), level = "warning")
      intron_regions <- intron_regions %>%dplyr::filter(region_type == "intron")
    }
  }
  
  log_message("INTRON DATA SOURCE VALIDATION: Confirmed intron regions data integrity")
  
  # Initialize result data frame
  intron_frequencies <- data.frame()
  
  # For each species, match variants to intron regions
  for (species_name in unique(intron_variants$species)) {
    species_variants <- intron_variants %>%
      dplyr::filter(species == species_name)
    
    species_regions <- intron_regions %>%
      dplyr::filter(species == species_name)
    
    if (nrow(species_regions) == 0) {
      log_message(sprintf("No intron regions found for species: %s", species_name), level = "warning")
      next
    }
    
    # For each variant type
    for (var_type_name in unique(species_variants$var_type)) {
      var_type_variants <- species_variants %>%
        dplyr::filter(var_type == var_type_name)
      
      # For each intron region
      for (i in 1:nrow(species_regions)) {
        region_info_row <- species_regions[i, ]
        
        # Find variants within this region's coordinates
        variants_in_region <- var_type_variants %>%
          dplyr::filter(position >= region_info_row$region_start & 
                       position <= region_info_row$region_end)
        
        variant_count <- nrow(variants_in_region)
        
        # Calculate normalized frequency
        normalized_freq <- (variant_count / region_info_row$region_length) * config$normalization$frequency_unit
        
        # Get genome_region from region_info
        genome_region <- if ("genome_region" %in% names(region_info_row)) {
          region_info_row$genome_region
        } else {
          "unknown"
        }
        
        # Add to results
        # Carry forward region_name for downstream P03 tables
        # Provide both region_name (intron-level) and gene (gene-level) for downstream flexibility
        intron_frequencies <- rbind(intron_frequencies, data.frame(
          species = species_name,
          var_type = var_type_name,
          region_type = "intron",
          genome_region = genome_region,
          region_name = region_info_row$region_name,  # Intron identifier (e.g., atpF_intron1)
          gene = region_info_row$gene,                # Gene-level identifier (e.g., atpF)
          variant_count = variant_count,
          region_length = region_info_row$region_length,
          frequency_per_kb = normalized_freq,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Remove rows with zero variants (optional - keep for complete matrix)
  # intron_frequencies <- intron_frequencies %>%dplyr::filter(variant_count > 0)
  
  log_message(sprintf("Generated intron frequencies for %d region-species-variant_type combinations", nrow(intron_frequencies)))
  
  # Log summary statistics
  region_summary <- intron_frequencies %>%
    dplyr::group_by(species, gene) %>%
    dplyr::summarise(total_variants = sum(variant_count), .groups = "drop") %>%
    dplyr::filter(total_variants > 0)
  
  if (nrow(region_summary) > 0) {
    log_message(sprintf("Intron regions with variants: %d (max variants per region: %d)", 
                       nrow(region_summary), max(region_summary$total_variants)))
  }
  
  return(intron_frequencies)
}

#' Combine frequency data from different regions
#' 
#' @param cds_frequencies CDS frequency data
#' @param igs_frequencies IGS frequency data
#' @param intron_frequencies Intron frequency data
#' @param other_frequencies Additional frequency data from other regions
#' @return Combined frequency data with standardized structure
combine_frequency_data <- function(cds_frequencies, igs_frequencies, intron_frequencies, other_frequencies = data.frame()) {
  
  log_message("=== COMBINING FREQUENCY DATA FROM ALL REGIONS ===")
  log_message(sprintf("INPUT DATA: CDS=%d, IGS=%d, Intron=%d, Other=%d entries", 
                     nrow(cds_frequencies), nrow(igs_frequencies), nrow(intron_frequencies), nrow(other_frequencies)))
  
  # Log variant counts for each region type
  cds_variants <- if(nrow(cds_frequencies) > 0) sum(cds_frequencies$variant_count) else 0
  igs_variants <- if(nrow(igs_frequencies) > 0) sum(igs_frequencies$variant_count) else 0
  intron_variants <- if(nrow(intron_frequencies) > 0) sum(intron_frequencies$variant_count) else 0
  other_variants <- if(nrow(other_frequencies) > 0) sum(other_frequencies$variant_count) else 0
  
  log_message(sprintf("VARIANT COUNTS: CDS=%d, IGS=%d, Intron=%d, Other=%d (Total=%d)", 
                     cds_variants, igs_variants, intron_variants, other_variants,
                     cds_variants + igs_variants + intron_variants + other_variants))
  
  # Standardize column names and add data type
  standardize_freq_data <- function(data, region_name) {
    if (nrow(data) == 0) {
      log_message(sprintf("No data for %s region", region_name), level = "warning")
      return(data.frame())
    }
    
    # Log details about this region's data
    region_variants <- sum(data$variant_count)
    region_entries_with_variants <- sum(data$variant_count > 0)
    log_message(sprintf("%s REGION: %d entries, %d with variants, %d total variants", 
                       region_name, nrow(data), region_entries_with_variants, region_variants))
    
    # Include genome_region and region_name in the required columns
    required_cols <- c("species", "var_type", "region_type", "genome_region", "region_name", "gene", "variant_count", "region_length", "frequency_per_kb")
    
    for (col in required_cols) {
      if (!col %in% names(data)) {
        if (col == "genome_region") {
          data[[col]] <- "unknown"
        } else if (col == "region_name") {
          data[[col]] <- data$gene  # Fallback: use gene as region_name if missing
        } else {
          data[[col]] <- NA
        }
        log_message(sprintf("Missing column '%s' in %s data, filled with fallback value", col, region_name), level = "warning")
      }
    }
    
    # Keep a consistent column order, including genome_region and region_name
    standardized_data <- data[, c("species", "var_type", "region_type", "genome_region", "region_name", "gene", "variant_count", "region_length", "frequency_per_kb")]
    
    log_message(sprintf("Standardized %s data: %d rows", region_name, nrow(standardized_data)))
    
    return(standardized_data)
  }
  
  # Standardize all datasets
  cds_std <- standardize_freq_data(cds_frequencies, "CDS")
  igs_std <- standardize_freq_data(igs_frequencies, "IGS")
  intron_std <- standardize_freq_data(intron_frequencies, "intron")
  other_std <- standardize_freq_data(other_frequencies, "OTHER")
  
  # Combine all data including unprocessed variants (HOTFIX: RNA pass-through)
  combined_data <- rbind(cds_std, igs_std, intron_std, other_std)
  
  if (nrow(combined_data) == 0) {
    log_message("No frequency data to combine", level = "warning")
    return(data.frame())
  }
  
  # Clean up data types
  combined_data <- combined_data %>%
dplyr::mutate(
      species = as.character(species),
      var_type = as.character(var_type),
      region_type = as.character(region_type),
      gene = as.character(gene),
      variant_count = as.numeric(variant_count),
      region_length = as.numeric(region_length),
      frequency_per_kb = as.numeric(frequency_per_kb)
    )
  
  # V3.16.12: RNA Data Validation (not aggregation - preprocessing already handled this)
  # After preprocessing standardization, all RNA types should already be "RNA"
  # This is a validation check to ensure data integrity
  log_message("=== Validating RNA Data Integrity ===")
  
  # Check for any legacy RNA subtypes that shouldn't exist after preprocessing
  legacy_rna_types <- c("rRNA", "tRNA-intron", "tRNA-CDS")
  legacy_rna_found <- combined_data %>%
    dplyr::filter(region_type %in% legacy_rna_types)
  
  if (nrow(legacy_rna_found) > 0) {
    log_message(sprintf("[WARNING] Found %d entries with legacy RNA subtypes (should have been standardized in preprocessing)", 
                       nrow(legacy_rna_found)), level = "warning")
    log_message("Legacy RNA types found:", level = "warning")
    legacy_summary <- legacy_rna_found %>%
      dplyr::group_by(region_type) %>%
      dplyr::summarise(entries = dplyr::n(), .groups = "drop")
    for (i in 1:nrow(legacy_summary)) {
      log_message(sprintf("  %s: %d entries", legacy_summary$region_type[i], legacy_summary$entries[i]), level = "warning")
    }
  }
  
  # Count standardized RNA data
  rna_data <- combined_data %>%
    dplyr::filter(region_type == "RNA") %>%
    dplyr::summarise(
      entries = dplyr::n(),
      variants = sum(variant_count, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(rna_data) > 0 && rna_data$entries > 0) {
    log_message(sprintf("RNA data validated: %d entries, %d variants", 
                       rna_data$entries, rna_data$variants))
  } else {
    log_message("No RNA data found in normalized results")
  }
  
  # Log final region type summary
  final_region_summary <- combined_data %>%
    dplyr::group_by(region_type) %>%
    dplyr::summarise(
      entries =dplyr::n(),
      variants = sum(variant_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(variants))
  
  log_message("FINAL REGION TYPE SUMMARY AFTER RNA AGGREGATION:")
  for (i in 1:nrow(final_region_summary)) {
    log_message(sprintf("  %s: %d entries, %d variants", 
                       final_region_summary$region_type[i], 
                       final_region_summary$entries[i], 
                       final_region_summary$variants[i]))
  }
  
  # Validation: Ensure only expected region types remain
  # Include common region type variations to prevent false warnings
  # V3.16.12: Use only standard terms after preprocessing
  expected_types <- c("CDS", "IGS", "intron", "others", "RNA")
  actual_types <- unique(combined_data$region_type)
  unexpected_types <- setdiff(actual_types, expected_types)
  
  if (length(unexpected_types) > 0) {
    log_message(sprintf("  [WARNING] VALIDATION WARNING: Unexpected region types found: %s", 
                       paste(unexpected_types, collapse=", ")), level = "warning")
  } else {
    log_message("  [SUCCESS] VALIDATION SUCCESS: Only expected region types remain (CDS, IGS/intergenic, intron, RNA)")
  }
  
  # Final combination statistics
  final_total_variants <- sum(combined_data$variant_count)
  final_entries_with_variants <- sum(combined_data$variant_count > 0)
  final_zero_entries <- sum(combined_data$variant_count == 0)
  
  log_message(sprintf("FINAL COMBINED DATA: %d total entries (%d CDS, %d IGS, %d intron)", 
                     nrow(combined_data), nrow(cds_std), nrow(igs_std), nrow(intron_std)))
  log_message(sprintf("FINAL STATISTICS: %d entries with variants, %d zero entries, %d total variants preserved", 
                     final_entries_with_variants, final_zero_entries, final_total_variants))
  
  return(combined_data)
}

#' Generate comprehensive normalization statistics
#' 
#' @param frequency_data Combined frequency data
#' @param variant_data Original variant data
#' @param config Configuration list
#' @return List of detailed normalization statistics
generate_normalization_statistics <- function(frequency_data, variant_data, config) {
  
  log_message("Generating comprehensive normalization statistics")
  
  if (nrow(frequency_data) == 0) {
    log_message("No frequency data for statistics generation", level = "warning")
    return(list(
      overall = data.frame(),
      by_species = data.frame(),
      by_region_type = data.frame(),
      by_var_type = data.frame(),
      frequency_distribution = list()
    ))
  }
  
  # Calculate key statistics for clear labeling
  observed_variants <- sum(frequency_data$variant_count, na.rm = TRUE)
  analysis_entries <- nrow(frequency_data)
  entries_with_variants <- sum(frequency_data$variant_count > 0, na.rm = TRUE)
  zero_entries <- sum(frequency_data$variant_count == 0, na.rm = TRUE)
  
  # Overall statistics with improved labeling
  overall_stats <- data.frame(
    metric = c("analysis_entries_total", "entries_with_observed_variants", "entries_with_zero_variants", 
               "observed_variants_total", "original_input_variants", "species_count", 
               "region_types_analyzed", "variant_types_analyzed", "genes_regions_analyzed", 
               "matrix_expansion_entries", "frequency_unit"),
    value = c(
      analysis_entries,
      entries_with_variants,
      zero_entries,
      observed_variants,
      nrow(variant_data),
      length(unique(frequency_data$species)),
      length(unique(frequency_data$region_type)),
      length(unique(frequency_data$var_type)),
      length(unique(frequency_data$gene)),
      analysis_entries - nrow(variant_data),
      config$normalization$frequency_unit
    ),
    description = c(
      "Total analysis entries in normalized matrix",
      "Entries with observed variants (variant_count > 0)",
      "Entries with zero variants (matrix placeholders)",
      "Total observed variants preserved in analysis",
      "Original variants in input data",
      "Number of species analyzed",
      "Number of region types (CDS, IGS, intron)",
      "Number of variant types analyzed",
      "Number of genes/regions analyzed",
      "Additional entries from matrix expansion",
      "Normalization unit (per X bases)"
    ),
    stringsAsFactors = FALSE
  )
  
  # Statistics by species with improved labeling
  species_stats <- frequency_data %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
      analysis_entries = dplyr::n(),
      entries_with_variants = sum(variant_count > 0, na.rm = TRUE),
      entries_with_zero = sum(variant_count == 0, na.rm = TRUE),
      genes_regions_analyzed = length(unique(gene)),
      variant_types = length(unique(var_type)),
      observed_variants = sum(variant_count, na.rm = TRUE),
      mean_frequency = mean(frequency_per_kb, na.rm = TRUE),
      median_frequency = stats::median(frequency_per_kb, na.rm = TRUE),
      max_frequency = max(frequency_per_kb, na.rm = TRUE),
      min_frequency = min(frequency_per_kb, na.rm = TRUE),
      sd_frequency = stats::sd(frequency_per_kb, na.rm = TRUE),
      na_count = sum(is.na(frequency_per_kb)),
      .groups = "drop"
    )
  
  # Statistics by region type with improved labeling
  region_type_stats <- frequency_data %>%
    dplyr::group_by(region_type) %>%
    dplyr::summarise(
      analysis_entries = dplyr::n(),
      entries_with_variants = sum(variant_count > 0, na.rm = TRUE),
      entries_with_zero = sum(variant_count == 0, na.rm = TRUE),
      species_count = length(unique(species)),
      genes_regions_analyzed = length(unique(gene)),
      observed_variants = sum(variant_count, na.rm = TRUE),
      mean_frequency = mean(frequency_per_kb, na.rm = TRUE),
      median_frequency = stats::median(frequency_per_kb, na.rm = TRUE),
      max_frequency = max(frequency_per_kb, na.rm = TRUE),
      min_frequency = min(frequency_per_kb, na.rm = TRUE),
      sd_frequency = stats::sd(frequency_per_kb, na.rm = TRUE),
      na_count = sum(is.na(frequency_per_kb)),
      .groups = "drop"
    )
  
  # Statistics by variant type with improved labeling
  var_type_stats <- frequency_data %>%
    dplyr::group_by(var_type) %>%
    dplyr::summarise(
      analysis_entries = dplyr::n(),
      entries_with_variants = sum(variant_count > 0, na.rm = TRUE),
      entries_with_zero = sum(variant_count == 0, na.rm = TRUE),
      species_count = length(unique(species)),
      genes_regions_analyzed = length(unique(gene)),
      observed_variants = sum(variant_count, na.rm = TRUE),
      mean_frequency = mean(frequency_per_kb, na.rm = TRUE),
      median_frequency = stats::median(frequency_per_kb, na.rm = TRUE),
      max_frequency = max(frequency_per_kb, na.rm = TRUE),
      min_frequency = min(frequency_per_kb, na.rm = TRUE),
      sd_frequency = stats::sd(frequency_per_kb, na.rm = TRUE),
      na_count = sum(is.na(frequency_per_kb)),
      .groups = "drop"
    )
  
  # Top variant hotspots for Upset plot (CDS only)
  top_hotspots_per_species_cds <- frequency_data %>%
    dplyr::filter(region_type == "CDS") %>%
    dplyr::group_by(species) %>%
    dplyr::arrange(dplyr::desc(frequency_per_kb)) %>%
    dplyr::filter(frequency_per_kb > 0) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::summarise(
      top_genes = list(gene),
      .groups = 'drop'
    )
  
  # Top variant hotspots (genes/regions with highest frequencies and variant counts)
  top_hotspots <- frequency_data %>%
    dplyr::filter(!is.na(frequency_per_kb) & frequency_per_kb > 0) %>%
    dplyr::group_by(species, region_type, gene) %>%
    dplyr::summarise(
      observed_variants = sum(variant_count, na.rm = TRUE),
      analysis_entries = dplyr::n(),
      entries_with_variants = sum(variant_count > 0, na.rm = TRUE),
      mean_frequency = mean(frequency_per_kb, na.rm = TRUE),
      max_frequency = max(frequency_per_kb, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(max_frequency)) %>%
    dplyr::slice_head(n = 20)
  
  # Frequency distribution
  freq_summary <- calculate_summary_stats(frequency_data$frequency_per_kb)
  
  log_message("Normalization statistics generation completed")
  
  return(list(
    overall = overall_stats,
    by_species = species_stats,
    by_region_type = region_type_stats,
    by_var_type = var_type_stats,
    top_hotspots_per_species_cds = top_hotspots_per_species_cds,
    top_hotspots = top_hotspots,
    frequency_distribution = freq_summary
  ))
}

#' Save normalization results to the structured P03 directory
#'
#' @param norm_results Results from normalize_frequencies
#' @param config Configuration list
#' @return Invisible NULL
save_normalization_results <- function(norm_results, config) {
  session_id <- config$session_info$session_id
  if (is.null(session_id) || session_id == "") {
    stop("FATAL: session_id is missing, cannot save normalization results.")
  }

  # Use the new path manager to get the correct P03 output directory
  output_dir <- get_processing_stage_path(session_id, "P03_normalized", create_dir = TRUE)
  log_message(sprintf("Saving normalization results to: %s", output_dir))

  # Save RDS object if requested (V31 fix: add defensive programming)
  if (!is.null(config$output$save_rds) && config$output$save_rds) {
    rds_file <- file.path(output_dir, "normalization_results.rds")
    saveRDS(norm_results, rds_file)
    log_message("RDS object saved successfully.")
  }

  # Save CSV files
  if (is.null(config$output$save_csv) || config$output$save_csv) {
    # Helper to safely write CSV files (V31 fix: enhanced error processing)
    write_csv_safe <- function(data, file_name) {
      if (!is.null(data) && nrow(data) > 0) {
        file_path <- file.path(output_dir, file_name)
        if (nchar(file_path) > 0 && nchar(file_name) > 0) {
          result <- safe_write_file(data, file_path, row.names = FALSE)
          if (result) {
            log_message(sprintf("[SUCCESS] Saved: %s (%d rows)", file_name, nrow(data)))
          } else {
            log_message(sprintf("[ERROR] Failed to save: %s", file_name), level = "error")
          }
        } else {
          log_message(sprintf("[ERROR] Invalid file path for: %s", file_name), level = "error")
        }
      } else {
        log_message(sprintf("[WARNING] Skipping empty data: %s", file_name), level = "warning")
      }
    }

    write_csv_safe(norm_results$frequencies, "normalized_frequencies.csv")
    write_csv_safe(norm_results$frequencies_summary, "normalized_frequencies_summary.csv")
    write_csv_safe(norm_results$statistics$overall, "normalization_statistics.csv")
    write_csv_safe(norm_results$statistics$by_species, "normalization_species_stats.csv")
    write_csv_safe(norm_results$statistics$by_region_type, "normalization_region_stats.csv")
    write_csv_safe(norm_results$statistics$by_var_type, "normalization_variant_type_stats.csv")
    write_csv_safe(norm_results$statistics$top_hotspots, "normalization_top_hotspots.csv")
  }

  log_message("Normalization results saved successfully.")
  return(invisible(NULL))
}

#' Validate normalization results
#' 
#' @param norm_results Results from normalize_frequencies
#' @param config Configuration list
#' @return Logical, TRUE if validation passes
validate_normalization_results <- function(norm_results, config) {
  
  # Check that we have the main results object
  if (is.null(norm_results)) {
    log_message("Normalization results object is NULL", level = "error")
    return(FALSE)
  }
  
  # Check for required top-level elements
  required_elements <- c("frequencies", "frequencies_summary", "statistics", "metadata")
  if (!all(required_elements %in% names(norm_results))) {
    log_message("Missing required elements in normalization results", level = "error")
    return(FALSE)
  }
  
  # Validate detailed frequencies
  frequency_data <- norm_results$frequencies
  if (nrow(frequency_data) == 0) {
    log_message("No detailed frequency data generated", level = "error")
    return(FALSE)
  }
  
  # Check for required columns in detailed frequencies
  required_cols <- c("species", "var_type", "region_type", "gene", "frequency_per_kb")
  if (!all(required_cols %in% names(frequency_data))) {
    log_message("Missing required columns in detailed frequency data", level = "error")
    return(FALSE)
  }
  
  # Validate summary frequencies
  summary_data <- norm_results$frequencies_summary
  if (nrow(summary_data) == 0) {
    log_message("No summary frequency data generated", level = "warning")
  } else {
    # Check for required columns in summary frequencies
    if (!all(required_cols %in% names(summary_data))) {
      log_message("Missing required columns in summary frequency data", level = "error")
      return(FALSE)
    }
    
    # Validate that all summary entries have var_type = "all_types"
    if (!all(summary_data$var_type == "all_types")) {
      log_message("Summary frequency data should have var_type = 'all_types'", level = "error")
      return(FALSE)
    }
  }
  
  # Check for reasonable frequency values
  finite_frequencies <- frequency_data$frequency_per_kb[is.finite(frequency_data$frequency_per_kb)]
  
  if (length(finite_frequencies) == 0) {
    log_message("No finite frequency values found", level = "warning")
  } else {
    # Check for negative frequencies
    if (any(finite_frequencies < 0)) {
      log_message("Negative frequencies detected", level = "warning")
    }
    
    # Check for extremely high frequencies
    if (any(finite_frequencies > 1000)) {
      log_message("Extremely high frequencies detected (>1000)", level = "warning")
    }
  }
  
  log_message("Normalization validation completed")
  return(TRUE)
}

#' Validate data flow and transformation integrity
#' 
#' @param variant_data Original variant data
#' @param norm_results Results from normalize_frequencies
#' @param config Configuration list
#' @return List containing validation results and detailed explanations
validate_data_flow <- function(variant_data, norm_results, config) {
  
  log_message("=== DATA FLOW VALIDATION ===", level = "debug")
  
  # Initialize validation results
  validation_results <- list(
    passed = TRUE,
    errors = c(),
    warnings = c(),
    explanations = c(),
    statistics = list()
  )
  
  # Original data statistics
  original_variants <- nrow(variant_data)
  original_species <- length(unique(variant_data$species))
  original_var_types <- unique(variant_data$var_type)
  
  # Normalized data statistics
  normalized_entries <- nrow(norm_results$frequencies)
  normalized_variants <- sum(norm_results$frequencies$variant_count)
  normalized_species <- length(unique(norm_results$frequencies$species))
  normalized_var_types <- unique(norm_results$frequencies$var_type)
  
  # Summary data statistics
  summary_entries <- nrow(norm_results$frequencies_summary)
  summary_variants <- sum(norm_results$frequencies_summary$variant_count)
  
  # Store statistics
  validation_results$statistics <- list(
    original_variants = original_variants,
    original_species = original_species,
    original_var_types = length(original_var_types),
    normalized_entries = normalized_entries,
    normalized_variants = normalized_variants,
    normalized_species = normalized_species,
    normalized_var_types = length(normalized_var_types),
    summary_entries = summary_entries,
    summary_variants = summary_variants,
    matrix_expansion = normalized_entries - original_variants
  )
  
  log_message(sprintf("VALIDATION STATISTICS:"), level = "debug")
  log_message(sprintf("  Original: %d variants, %d species, %d variant types", 
                     original_variants, original_species, length(original_var_types)), level = "debug")
  log_message(sprintf("  Normalized: %d entries, %d variants, %d species, %d variant types", 
                     normalized_entries, normalized_variants, normalized_species, length(normalized_var_types)), level = "debug")
  log_message(sprintf("  Summary: %d entries, %d variants", summary_entries, summary_variants), level = "debug")
  log_message(sprintf("  Matrix expansion: +%d entries", normalized_entries - original_variants), level = "debug")
  
  # Validation 1: Species consistency
  if (original_species != normalized_species) {
    validation_results$errors <- c(validation_results$errors,
                                  sprintf("Species count mismatch: %d original vs %d normalized", 
                                         original_species, normalized_species))
    validation_results$passed <- FALSE
  }
  
  # Validation 2: Variant type consistency
  if (length(setdiff(original_var_types, normalized_var_types)) > 0) {
    validation_results$errors <- c(validation_results$errors,
                                  "Some original variant types missing in normalized data")
    validation_results$passed <- FALSE
  }
  
  # Validation 3: Variant count preservation
  if (normalized_variants != summary_variants) {
    validation_results$errors <- c(validation_results$errors,
                                  sprintf("Variant count mismatch between detailed (%d) and summary (%d)", 
                                         normalized_variants, summary_variants))
    validation_results$passed <- FALSE
  }
  
  # Validation 4: Matrix expansion reasonableness
  expansion_ratio <- normalized_entries / original_variants
  if (expansion_ratio > 10) {
    validation_results$warnings <- c(validation_results$warnings,
                                    sprintf("Large matrix expansion ratio: %.2f (may indicate data issues)", 
                                           expansion_ratio))
  }
  
  # Validation 5: Zero entries check
  zero_entries <- sum(norm_results$frequencies$variant_count == 0)
  zero_percentage <- (zero_entries / normalized_entries) * 100
  
  if (zero_percentage > 90) {
    validation_results$warnings <- c(validation_results$warnings,
                                    sprintf("High percentage of zero entries: %.1f%%", zero_percentage))
  }
  
  # Explanations for the transformation
  validation_results$explanations <- c(
    sprintf("MATRIX EXPANSION: %d original variants expanded to %d analysis entries", 
            original_variants, normalized_entries),
    sprintf("REASON: Complete analytical matrix includes all possible species x gene x variant_type combinations"),
    sprintf("ZERO ENTRIES: %d entries (%.1f%%) have zero variants (analytical placeholders)", 
            zero_entries, zero_percentage),
    sprintf("VARIANT PRESERVATION: %d variants preserved in analysis (%d in summary)", 
            normalized_variants, summary_variants),
    sprintf("SCIENTIFIC RATIONALE: Matrix expansion enables proper normalization and comparative analysis")
  )
  
  # Log validation results
  if (validation_results$passed) {
    log_message("DATA FLOW VALIDATION: PASSED", level = "debug")
  } else {
    log_message("DATA FLOW VALIDATION: FAILED", level = "debug")
    for (error in validation_results$errors) {
      log_message(sprintf("  ERROR: %s", error), level = "debug")
    }
  }
  
  if (length(validation_results$warnings) > 0) {
    for (warning in validation_results$warnings) {
      log_message(sprintf("  WARNING: %s", warning), level = "debug")
    }
  }
  
  for (explanation in validation_results$explanations) {
    log_message(sprintf("  EXPLANATION: %s", explanation), level = "debug")
  }
  
  return(validation_results)
}

# log_message("Frequency normalization functions loaded successfully")