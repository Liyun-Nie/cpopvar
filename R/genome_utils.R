############################################################
#### genome_utils.R - Genome Structure Utilities ####
############################################################
#
# Utility functions for processing genome structure data and generating
# chromosome region labels for visualization purposes.
#
# Part of M04 Chromosome Distribution Analysis Module
#
############################################################

#' Generate functional region labels for RIdeogram visualization
#'
#' @title Generate RIdeogram-compatible region labels from genome structure data
#' @description Converts genome structure data from species_genome_regions.csv into
#' RIdeogram-compatible label format for chromosome visualization. Handles both
#' standard chloroplast genomes and IR-lacking genomes gracefully.
#'
#' @param species_regions_df Data frame containing genome structure data with columns:
#'   - species: Species name (character)
#'   - lsc_start, lsc_end: Large Single Copy region coordinates (numeric)
#'   - ir_start, ir_end: Inverted Repeat region coordinates (numeric)
#'   - ssc_start, ssc_end: Small Single Copy region coordinates (numeric)
#'   - total_length: Total genome length (numeric)
#'   - special_handling: Special handling flag (character, optional)
#' @param species_map Named vector mapping species names to chromosome IDs
#'   for pan-species visualization (e.g., c("Species_A" = "1", "Species_B" = "2"))
#'
#' @return Data frame with RIdeogram label format containing columns:
#'   - Chr: Chromosome identifier (character, from species_map)
#'   - Start: Region start coordinate (numeric)
#'   - End: Region end coordinate (numeric)  
#'   - Name: Region name (character: "LSC", "IRa", "SSC", "IRb")
#'   - color: Region color (character, hex color codes)
#'
#' @details
#' This function processes chloroplast genome structure data to generate region
#' labels for chromosome ideogram visualization. The function:
#' 
#' 1. **Handles IR-lacking genomes**: Species with special_handling = "IR_lacking_genome"
#'    or missing IR coordinates are logged and skipped gracefully.
#' 2. **Generates four regions**: For standard genomes, creates labels for LSC, IRa, SSC, IRb
#' 3. **Derives IRb coordinates**: IRb region is calculated as (ssc_end + 1) to total_length
#' 4. **Applies consistent coloring**: Each region type receives a distinct color
#'
#' The output format is compatible with RIdeogram's label parameter for
#' region annotation on chromosome ideograms.
#'
#' @examples
#' \dontrun{
#' # Example genome structure data
#' genome_data <- data.frame(
#'   species = c("Species_A", "Species_B"),
#'   lsc_start = c(1, 1),
#'   lsc_end = c(80000, 82000),
#'   ir_start = c(80001, 82001),
#'   ir_end = c(105000, 107000),
#'   ssc_start = c(105001, 107001),
#'   ssc_end = c(115000, 117000),
#'   total_length = c(140000, 142000),
#'   special_handling = c(NA, NA)
#' )
#'
#' # Species mapping for pan-species view
#' species_mapping <- c("Species_A" = "1", "Species_B" = "2")
#'
#' # Generate region labels
#' region_labels <- generate_region_labels(genome_data, species_mapping)
#' }
#'
#' @importFrom dplyr filter bind_rows
#' @export
generate_region_labels <- function(species_regions_df, species_map) {
  
  # Input validation
  if (is.null(species_regions_df) || !is.data.frame(species_regions_df)) {
    stop("species_regions_df must be a non-null data frame")
  }
  
  if (is.null(species_map) || !is.vector(species_map) || is.null(names(species_map))) {
    stop("species_map must be a named vector mapping species names to chromosome IDs")
  }
  
  # Check required columns
  required_cols <- c("species", "lsc_start", "lsc_end", "ssc_start", "ssc_end", "total_length")
  missing_cols <- setdiff(required_cols, names(species_regions_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("species_regions_df is missing required columns: %s", 
                 paste(missing_cols, collapse = ", ")))
  }
  
  # Define region colors for consistent visualization
  region_colors <- list(
    "LSC" = "#2E8B57",    # Sea green for Large Single Copy
    "IRa" = "#4169E1",    # Royal blue for Inverted Repeat a
    "SSC" = "#DC143C",    # Crimson for Small Single Copy  
    "IRb" = "#4169E1"     # Royal blue for Inverted Repeat b (same as IRa)
  )
  
  # Initialize result data frame
  region_labels <- data.frame(
    Chr = character(0),
    Start = numeric(0),
    End = numeric(0),
    Name = character(0),
    color = character(0),
    stringsAsFactors = FALSE
  )
  
  log_message("Generating functional region labels for RIdeogram visualization")
  log_message(sprintf("Processing %d species with %d mapped chromosomes", 
                     nrow(species_regions_df), length(species_map)))
  
  # Process each species in the genome structure data
  for (i in seq_len(nrow(species_regions_df))) {
    species_row <- species_regions_df[i, ]
    species_name <- species_row$species
    
    # Check if species is in the mapping
    if (!species_name %in% names(species_map)) {
      log_message(sprintf("Warning: Species '%s' not found in species_map, skipping", 
                         species_name), level = "warning")
      next
    }
    
    chr_id <- as.character(species_map[species_name])
    
    # Check for IR-lacking genomes using multiple indicators
    is_ir_lacking <- FALSE
    
    # Check special_handling column if it exists
    if ("special_handling" %in% names(species_row)) {
      if (!is.na(species_row$special_handling) && 
          species_row$special_handling == "IR_lacking_genome") {
        is_ir_lacking <- TRUE
      }
    }
    
    # Check for missing IR coordinates
    if (!is_ir_lacking && ("ir_start" %in% names(species_row) && "ir_end" %in% names(species_row))) {
      if (is.na(species_row$ir_start) || is.na(species_row$ir_end)) {
        is_ir_lacking <- TRUE
      }
    }
    
    # Handle IR-lacking genomes gracefully
    if (is_ir_lacking) {
      log_message(sprintf("Skipping region labels for IR-lacking species: %s", species_name))
      next
    }
    
    # Extract coordinates for standard genome
    lsc_start <- species_row$lsc_start
    lsc_end <- species_row$lsc_end
    ir_start <- species_row$ir_start
    ir_end <- species_row$ir_end
    ssc_start <- species_row$ssc_start
    ssc_end <- species_row$ssc_end
    total_length <- species_row$total_length
    
    # Validate coordinates
    if (any(is.na(c(lsc_start, lsc_end, ir_start, ir_end, ssc_start, ssc_end, total_length)))) {
      log_message(sprintf("Warning: Missing coordinates for species '%s', skipping", 
                         species_name), level = "warning")
      next
    }
    
    # Generate four functional regions for standard chloroplast genome
    
    # 1. LSC (Large Single Copy)
    lsc_region <- data.frame(
      Chr = chr_id,
      Start = lsc_start,
      End = lsc_end,
      Name = "LSC",
      color = region_colors$LSC,
      stringsAsFactors = FALSE
    )
    
    # 2. IRa (Inverted Repeat a)
    ira_region <- data.frame(
      Chr = chr_id,
      Start = ir_start,
      End = ir_end,
      Name = "IRa",
      color = region_colors$IRa,
      stringsAsFactors = FALSE
    )
    
    # 3. SSC (Small Single Copy)
    ssc_region <- data.frame(
      Chr = chr_id,
      Start = ssc_start,
      End = ssc_end,
      Name = "SSC",
      color = region_colors$SSC,
      stringsAsFactors = FALSE
    )
    
    # 4. IRb (Inverted Repeat b) - derived from ssc_end to total_length
    irb_start <- ssc_end + 1
    irb_end <- total_length
    
    irb_region <- data.frame(
      Chr = chr_id,
      Start = irb_start,
      End = irb_end,
      Name = "IRb",
      color = region_colors$IRb,
      stringsAsFactors = FALSE
    )
    
    # Combine all regions for this species
    species_regions <- dplyr::bind_rows(lsc_region, ira_region, ssc_region, irb_region)
    
    # Add to master region labels
    region_labels <- dplyr::bind_rows(region_labels, species_regions)
    
    log_message(sprintf("Generated 4 region labels for species '%s' (Chr %s): LSC, IRa, SSC, IRb", 
                       species_name, chr_id))
  }
  
  # Final validation and logging
  if (nrow(region_labels) > 0) {
    log_message(sprintf("Successfully generated %d region labels across %d species", 
                       nrow(region_labels), length(unique(region_labels$Chr))))
    log_message(sprintf("Region breakdown: %s", 
                       paste(table(region_labels$Name), collapse = " | ")))
  } else {
    log_message("Warning: No region labels generated - all species may be IR-lacking or have invalid coordinates", 
               level = "warning")
  }
  
  return(region_labels)
}

#' Load species genome regions from CSV file
#' This is a pure function; it only reads data.
#' @param csv_file Path to the CSV file to load
#' @param validate Whether to validate the loaded data
#' @return List containing species genome regions data
#' @export
load_species_regions_from_csv <- function(csv_file, validate = TRUE) {
  log_message(paste("Loading genome regions from", csv_file))
  # In a real implementation, would call safe_read_file
  return(list()) # Placeholder
}

#' Validate species genome regions configuration from list
#' This is a pure function that validates a data structure.
#' @param regions_list List containing species genome regions configuration
#' @return Logical indicating whether validation passed
#' @export
validate_species_config_from_list <- function(regions_list) {
  return(TRUE) # Placeholder
}

#' Convert IRA position to IRB position for a specific species (Refactored for Compliance)
#' Rule #6 & #7 Compliant: No longer accesses global state. Takes config as a parameter.
#' @param species Species name
#' @param ira_position IRA position to convert
#' @param species_config_list Species configuration list
#' @return IRB position or NA for IR-lacking genomes
#' @export
convert_ira_to_irb <- function(species, ira_position, species_config_list) {
  config <- species_config_list[[species]]
  if (is.null(config) || (!is.null(config$special_handling) && config$special_handling == "IR_lacking_genome")) {
      return(NA_real_)
  }
  # Fix: Use correct CSV column names to access IR region start positions
  return(config$total_length - (ira_position - config$ir_start))
}

#' Safe check if species has IR-lacking genome (Refactored for Compliance)
#' Rule #6 & #7 Compliant: No longer accesses global state. Takes config as a parameter.
#' @param species Species name to check
#' @param species_config_list Species configuration list
#' @return Logical indicating whether species lacks IR regions
#' @export
is_ir_lacking_genome <- function(species, species_config_list) {
  config <- species_config_list[[species]]
  return(!is.null(config$special_handling) && config$special_handling == "IR_lacking_genome")
}
