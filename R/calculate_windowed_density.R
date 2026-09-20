############################################################
#### calculate_windowed_density.R - Windowed Density Calculation ####
############################################################
#
# Core algorithm for calculating variant density in genomic windows
# This function implements the standard windowed density approach used
# in population genomics for spatial distribution analysis.
#
# References RIdeogram.R script lines 40-150 for algorithm implementation
#
############################################################

#' Calculate variant density in genomic windows for multiple species
#'
#' This function calculates variant density across species using a sliding
#' window approach. It processes multi-species variant data and maps each
#' species to a unique chromosome ID for pan-species visualization.
#'
#' @param variant_df Data frame containing variant information with columns:
#'   - species: Species name (character)
#'   - Chr: Chromosome identifier (character) 
#'   - position: Genomic position of variant (numeric)
#'   Additional columns are preserved but not used in calculation
#' @param karyotype_df Data frame containing pan-species karyotype with columns:
#'   - Chr: Species chromosome ID (character, from species_map)
#'   - Start: Start position (numeric, typically 0)
#'   - End: End position/chromosome length (numeric)
#' @param species_map Named vector mapping species names to chromosome IDs
#' @param window_size Size of sliding windows in base pairs (numeric, default: 2000)
#'
#' @return Data frame with windowed density calculations containing:
#'   - Chr: Species chromosome identifier (from species_map)
#'   - Start: Start position of window (renamed from Window_Start for RIdeogram compatibility)
#'   - End: End position of window (renamed from Window_End for RIdeogram compatibility)
#'   - Value: Variant density value (renamed from Variant_Density for RIdeogram compatibility)
#'   - Variant_Count: Number of variants in window (preserved for debugging)
#'
#' @details
#' The function processes each chromosome independently, creating non-overlapping
#' windows of the specified size. For each window, it counts variants whose
#' positions fall within the window boundaries and calculates density as
#' count/window_size.
#'
#' The algorithm is based on the standard approach used in population genomics
#' for spatial variant analysis, as implemented in the reference RIdeogram
#' workflow.
#'
#' @examples
#' \dontrun{
#' # Example variant data
#' variants <- data.frame(
#'   species = "Species_A",
#'   Chr = c("1", "1", "2"),
#'   position = c(1500, 3500, 1000)
#' )
#' 
#' # Example karyotype data
#' karyotype <- data.frame(
#'   Chr = c("1", "2"),
#'   Start = c(1, 1),
#'   End = c(10000, 8000)
#' )
#' 
#' # Calculate windowed density
#' density_data <- calculate_windowed_density(variants, karyotype, 2000)
#' }
#'
#' @importFrom dplyr filter select %>%
#' @export
calculate_windowed_density <- function(variant_df, karyotype_df, species_map, window_size = 2000) {
  
  # Input validation
  if (is.null(variant_df) || !is.data.frame(variant_df)) {
    stop("variant_df must be a non-null data frame")
  }
  
  if (is.null(karyotype_df) || !is.data.frame(karyotype_df)) {
    stop("karyotype_df must be a non-null data frame")
  }
  
  if (is.null(species_map) || !is.vector(species_map) || is.null(names(species_map))) {
    stop("species_map must be a named vector mapping species names to chromosome IDs")
  }
  
  if (!is.numeric(window_size) || window_size <= 0) {
    stop("window_size must be a positive numeric value")
  }
  
  # Check required columns in variant_df
  required_variant_cols <- c("species", "Chr", "position")
  missing_variant_cols <- setdiff(required_variant_cols, names(variant_df))
  if (length(missing_variant_cols) > 0) {
    stop(sprintf("variant_df is missing required columns: %s", 
                 paste(missing_variant_cols, collapse = ", ")))
  }
  
  # Check required columns in karyotype_df
  required_karyotype_cols <- c("Chr", "End")
  missing_karyotype_cols <- setdiff(required_karyotype_cols, names(karyotype_df))
  if (length(missing_karyotype_cols) > 0) {
    stop(sprintf("karyotype_df is missing required columns: %s", 
                 paste(missing_karyotype_cols, collapse = ", ")))
  }
  
  # Ensure chromosome identifiers are character type for consistent matching
  variant_df$Chr <- as.character(variant_df$Chr)
  karyotype_df$Chr <- as.character(karyotype_df$Chr)
  
  # Initialize results data frame with RIdeogram-compatible column names
  results <- data.frame(
    Chr = character(),
    Start = numeric(),
    End = numeric(),
    Value = numeric(),
    Variant_Count = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Get unique species to process
  species_to_process <- unique(variant_df$species)
  
  # Validate that all species in data have corresponding entries in species_map
  missing_species <- setdiff(species_to_process, names(species_map))
  if (length(missing_species) > 0) {
    warning(sprintf("Species not found in species_map: %s", 
                   paste(missing_species, collapse = ", ")))
  }
  
  # Filter to only process species that exist in species_map
  valid_species <- intersect(species_to_process, names(species_map))
  
  if (length(valid_species) == 0) {
    warning("No valid species found that exist in both variant data and species_map")
    return(results)
  }
  
  # Process each species (mapped to chromosome IDs)
  for (species_name in valid_species) {
    
    # Get species-specific variant data
    species_variants <- variant_df %>% 
      dplyr::filter(species == species_name)
    
    # Map species name to chromosome ID using species_map
    species_chr_id <- as.character(species_map[species_name])
    
    # Get genome length for this species from pan-species karyotype
    chr_length <- karyotype_df %>% 
      dplyr::filter(Chr == species_chr_id) %>% 
      dplyr::select(End) %>% 
      as.numeric()
    
    if (length(chr_length) == 0) {
      warning(sprintf("No karyotype information found for species %s (Chr ID: %s)", 
                     species_name, species_chr_id))
      next
    }
    
    # Use first entry if multiple rows exist for same chromosome
    chr_length <- chr_length[1]
    
    if (exists("log_message")) {
      log_message(sprintf("Processing species %s -> Chr %s (length: %d bp)", 
                         species_name, species_chr_id, chr_length))
      log_message(sprintf("  - Input variants for %s: %d records", 
                         species_name, nrow(species_variants)))
      if (nrow(species_variants) > 0) {
        log_message(sprintf("  - Position range: %d - %d", 
                           min(species_variants$position), max(species_variants$position)))
      }
    }
    
    # Define windows - create non-overlapping windows across chromosome
    # Start from position 0 and create windows of specified size
    windows <- seq(0, chr_length, by = window_size)
    
    # Calculate density for each window
    for (i in 1:(length(windows) - 1)) {
      window_start <- windows[i]
      window_end <- windows[i + 1] - 1
      
      # Filter variants within current window
      # Use inclusive bounds for start and end positions
      variants_in_window <- species_variants %>%
        dplyr::filter(position >= window_start & position <= window_end)
      
      # Calculate variant count and density
      variant_count <- nrow(variants_in_window)
      variant_density <- variant_count / window_size
      
      # Add results to output data frame with RIdeogram-compatible format
      results <- rbind(results, data.frame(
        Chr = species_chr_id,
        Start = window_start,
        End = window_end,
        Value = variant_density,
        Variant_Count = variant_count,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Log summary information
  if (exists("log_message")) {
    log_message(sprintf("Calculated windowed density for %d species with window size %d bp", 
                       length(valid_species), window_size))
    log_message(sprintf("Generated %d windows total with %d total variants across all species", 
                       nrow(results), sum(results$Variant_Count)))
    
    # Log species-specific summary
    for (species_name in valid_species) {
      species_chr_id <- as.character(species_map[species_name])
      species_windows <- sum(results$Chr == species_chr_id)
      species_variants <- sum(results$Variant_Count[results$Chr == species_chr_id])
      log_message(sprintf("  Species %s (Chr %s): %d windows, %d variants", 
                         species_name, species_chr_id, species_windows, species_variants))
    }
  }
  
  # Final data structure validation before return
  if (exists("log_message") && nrow(results) > 0) {
    log_message("Final density calculation results validation:")
    log_message(sprintf("  - Total windows generated: %d", nrow(results)))
    log_message(sprintf("  - Unique Chr IDs: %s", paste(unique(results$Chr), collapse = ", ")))
    log_message(sprintf("  - Column structure: %s", paste(names(results), collapse = ", ")))
    log_message(sprintf("  - Value range: %.6f - %.6f", min(results$Value), max(results$Value)))
    log_message("  - Sample results (first 3 rows):")
    print(head(results, 3))
  }
  
  return(results)
}