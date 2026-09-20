############################################################
#### annotation_processing.R - GFF Annotation Processing Engine ####
############################################################
#
# Three-stage GFF processing pipeline for cpopvar
# 
# Architecture designed by Core Architect (Gemini CLI)
# Implementation by Core Engineer (Claude Code)
#
# This module implements the strategic upgrade from CSV-based annotation
# dependency to industry-standard GFF3 format processing, providing:
# - Stage 1: Format Conversion (external scripts, optional)
# - Stage 2: GFF Validation & Standardization (rtracklayer)
# - Stage 3: Feature Extraction & IGS/Intron Calculation (GenomicRanges)
#
# Key Features:
# - Replaces fragile CSV dependencies with robust GFF processing
# - Uses industry gold-standard tools (rtracklayer::import, GenomicRanges)
# - Generates unified output files for downstream analysis
# - Backward compatibility with legacy CSV mode
#
############################################################

#' @importFrom magrittr %>%
NULL

# Load required libraries
# NOTE: rtracklayer and GenomicRanges are essential for this module

#' Process GFF directory and generate unified annotation files
#' 
#' This is the main entry point for the three-stage GFF processing pipeline.
#' It coordinates the validation, standardization, and feature extraction
#' processes to generate unified annotation files for downstream analysis.
#' 
#' @param gff_directory Path to directory containing GFF3 files
#' @param config Configuration list containing processing parameters
#' @param output_dir Directory to save unified annotation files
#' @param progress_callback Optional function for progress updates
#' @return List containing processing results and file paths
#' @export
process_gff_directory <- function(gff_directory, 
                                 config, 
                                 output_dir, 
                                 progress_callback = NULL) {
  
  log_message("=== STARTING THREE-STAGE GFF PROCESSING PIPELINE ===")
  log_message(sprintf("GFF Directory: %s", gff_directory))
  log_message(sprintf("Output Directory: %s", output_dir))
  
  # Initialize progress tracking
  if (!is.null(progress_callback)) {
    progress_callback("Starting GFF processing pipeline", 0.1)
  }
  
  # Validate inputs
  if (!dir.exists(gff_directory)) {
    stop(sprintf("GFF directory does not exist: %s", gff_directory))
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    log_message(sprintf("Created output directory: %s", output_dir))
  }
  
  # STAGE 2: GFF Validation & Standardization
  if (!is.null(progress_callback)) {
    progress_callback("Stage 2: Validating and standardizing GFF files", 0.3)
  }
  
  log_message("=== STAGE 2: GFF VALIDATION & STANDARDIZATION ===")
  validated_gff_objects <- validate_gff_files(gff_directory, config)
  
  if (length(validated_gff_objects) == 0) {
    stop("No valid GFF files found or all files failed validation")
  }
  
  log_message(sprintf("Successfully validated %d GFF files", length(validated_gff_objects)))
  
  # STAGE 3: Feature Extraction & Region Calculation
  if (!is.null(progress_callback)) {
    progress_callback("Stage 3: Extracting features and calculating regions", 0.6)
  }
  
  log_message("=== STAGE 3: FEATURE EXTRACTION & REGION CALCULATION ===")
  extracted_features <- extract_region_features(validated_gff_objects, config)
  
  # Generate unified annotation files
  if (!is.null(progress_callback)) {
    progress_callback("Generating unified annotation files", 0.8)
  }
  
  log_message("=== GENERATING UNIFIED ANNOTATION FILES ===")
  unified_files <- generate_unified_annotations(extracted_features, output_dir, config)
  
  # Final validation and summary
  if (!is.null(progress_callback)) {
    progress_callback("Finalizing GFF processing", 1.0)
  }
  
  log_message("=== GFF PROCESSING PIPELINE COMPLETED SUCCESSFULLY ===")
  
  # Prepare results object
  results <- list(
    validated_gff_count = length(validated_gff_objects),
    species_processed = names(validated_gff_objects),
    extracted_features = extracted_features,
    unified_files = unified_files,
    processing_time = Sys.time(),
    pipeline_version = "v4.0_three_stage"
  )
  
  # Log summary statistics
  log_summary_statistics(results)
  
  return(results)
}

#' Validate and standardize GFF files using rtracklayer
#' 
#' Stage 2 of the pipeline: Uses rtracklayer::import() as the validation
#' and standardization engine. This function serves as a powerful "firewall"
#' that ensures only properly formatted GFF files proceed to feature extraction.
#' 
#' @param gff_directory Path to directory containing GFF files
#' @param config Configuration list
#' @return Named list of validated GRanges objects (one per species)
validate_gff_files <- function(gff_directory, config) {
  
  log_message("Stage 2: Validating GFF files using rtracklayer::import()")
  
  # Check for required packages
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required but not installed. Please install it using: BiocManager::install('rtracklayer')")
  }
  
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed. Please install it using: BiocManager::install('GenomicRanges')")
  }
  
  # Find all GFF files in directory
  gff_files <- list.files(gff_directory, pattern = "\\.gff$|\\.gff3$", full.names = TRUE, ignore.case = TRUE)
  
  if (length(gff_files) == 0) {
    stop(sprintf("No GFF files found in directory: %s", gff_directory))
  }
  
  log_message(sprintf("Found %d GFF files for validation", length(gff_files)))
  
  validated_objects <- list()
  validation_failures <- character()
  
  # Process each GFF file
  for (gff_file in gff_files) {
    species_name <- tools::file_path_sans_ext(basename(gff_file))
    log_message(sprintf("Validating GFF file for species: %s", species_name))
    
    tryCatch({
      # CORE VALIDATION: Use rtracklayer::import() as validation engine
      # This is the "golden standard" validation - if this fails, the GFF file is invalid
      gff_object <- rtracklayer::import(gff_file)
      
      log_message(sprintf("[OK] Successfully imported GFF file: %s (%d features)", 
                         species_name, length(gff_object)))
      
      # Standardization steps
      standardized_gff <- standardize_gff_object(gff_object, species_name, config)
      
      validated_objects[[species_name]] <- standardized_gff
      
      log_message(sprintf("[OK] Standardized GFF object for species: %s", species_name))
      
    }, error = function(e) {
      error_msg <- sprintf("[ERROR] Failed to validate GFF file %s: %s", species_name, e$message)
      log_message(error_msg, level = "error")
      validation_failures <- c(validation_failures, error_msg)
    })
  }
  
  # Report validation results
  if (length(validation_failures) > 0) {
    log_message("=== VALIDATION FAILURES ===", level = "warning")
    for (failure in validation_failures) {
      log_message(failure, level = "warning")
    }
  }
  
  log_message(sprintf("Validation completed: %d successful, %d failed", 
                     length(validated_objects), length(validation_failures)))
  
  return(validated_objects)
}

#' Standardize GFF object for consistent processing
#' 
#' Applies standardization rules to ensure consistent feature types,
#' chromosome naming, and attribute formatting across different GFF sources.
#' 
#' @param gff_object GRanges object from rtracklayer::import()
#' @param species_name Species name for this GFF file
#' @param config Configuration list
#' @return Standardized GRanges object
standardize_gff_object <- function(gff_object, species_name, config) {
  
  log_message(sprintf("Standardizing GFF object for %s", species_name), level = "debug")
  
  # Add species information to metadata
  GenomicRanges::mcols(gff_object)$species <- species_name
  
  # Standardize chromosome/sequence names (remove "chr" prefix if present)
  seqnames_standardized <- gsub("^chr", "", GenomicRanges::seqnames(gff_object))
  GenomicRanges::seqnames(gff_object) <- seqnames_standardized
  
  # Standardize feature types to internal conventions
  # Map various gene/CDS naming conventions to our standard types
  feature_type <- GenomicRanges::mcols(gff_object)$type
  if (!is.null(feature_type)) {
    # Standard mapping for common variations
    feature_type <- gsub("protein_coding_gene", "gene", feature_type)
    feature_type <- gsub("protein_coding", "CDS", feature_type)
    feature_type <- gsub("mRNA", "gene", feature_type)
    
    GenomicRanges::mcols(gff_object)$type <- feature_type
  }
  
  # Ensure essential attributes exist
  if (is.null(GenomicRanges::mcols(gff_object)$ID)) {
    # Generate IDs if missing
    GenomicRanges::mcols(gff_object)$ID <- paste0(species_name, "_feature_", seq_along(gff_object))
    log_message("Generated missing ID attributes", level = "debug")
  }
  
  log_message(sprintf("Standardization completed for %s: %d features", 
                     species_name, length(gff_object)), level = "debug")
  
  return(gff_object)
}

#' Extract region features and calculate IGS/Intron regions
#' 
#' Stage 3 of the pipeline: Uses GenomicRanges to perform sophisticated
#' genomic interval calculations to derive IGS and Intron regions from
#' the validated GFF data.
#' 
#' @param validated_gff_objects Named list of validated GRanges objects
#' @param config Configuration list
#' @return List containing extracted features for all species
extract_region_features <- function(validated_gff_objects, config) {
  
  log_message("Stage 3: Extracting region features using GenomicRanges")
  
  all_features <- list()
  
  for (species_name in names(validated_gff_objects)) {
    log_message(sprintf("Processing features for species: %s", species_name))
    
    gff_object <- validated_gff_objects[[species_name]]
    
    # Extract different feature types
    species_features <- list()
    
    # Extract genes (CDS features)
    species_features$cds <- extract_cds_features(gff_object, species_name)
    
    # Calculate IGS (intergenic sequences) using GenomicRanges::gaps()
    species_features$igs <- calculate_igs_regions(gff_object, species_name)
    
    # Calculate Introns using GenomicRanges::gaps() within genes
    species_features$introns <- calculate_intron_regions(gff_object, species_name)
    
    # Extract genome structure information
    species_features$genome_structure <- extract_genome_structure(gff_object, species_name)
    
    all_features[[species_name]] <- species_features
    
    log_message(sprintf("Feature extraction completed for %s: CDS=%d, IGS=%d, Introns=%d", 
                       species_name, 
                       nrow(species_features$cds),
                       nrow(species_features$igs), 
                       nrow(species_features$introns)))
  }
  
  return(all_features)
}

#' Extract CDS features from GFF object
#' 
#' @param gff_object Validated GRanges object
#' @param species_name Species name
#' @return Data frame with CDS features
extract_cds_features <- function(gff_object, species_name) {
  
  # Filter for CDS features
  cds_features <- gff_object[GenomicRanges::mcols(gff_object)$type == "CDS"]
  
  if (length(cds_features) == 0) {
    log_message(sprintf("No CDS features found for %s", species_name), level = "warning")
    return(data.frame())
  }
  
  # Extract information and calculate lengths
  cds_df <- data.frame(
    species = species_name,
    region_name = GenomicRanges::mcols(cds_features)$Name %||% 
                  GenomicRanges::mcols(cds_features)$ID %||% 
                  paste0("CDS_", seq_along(cds_features)),
    region_type = "CDS",
    region_start = GenomicRanges::start(cds_features),
    region_end = GenomicRanges::end(cds_features),
    region_length = GenomicRanges::width(cds_features),
    chromosome = as.character(GenomicRanges::seqnames(cds_features)),
    strand = as.character(GenomicRanges::strand(cds_features)),
    stringsAsFactors = FALSE
  )
  
  # Calculate total CDS length per gene (for genes with multiple CDS)
  cds_df$gene_name <- gsub(" CDS$", "", cds_df$region_name)
  
  total_cds_lengths <- cds_df %>%
    dplyr::group_by(species, gene_name) %>%
    dplyr::summarise(total_cds_length = sum(region_length), .groups = "drop")
  
  # Add total_cds_length to the dataframe
  cds_df <- cds_df %>%
    dplyr::left_join(total_cds_lengths, by = c("species", "gene_name"))
  
  return(cds_df)
}

#' Calculate IGS (intergenic sequences) using GenomicRanges::gaps()
#' 
#' @param gff_object Validated GRanges object
#' @param species_name Species name
#' @return Data frame with IGS regions
calculate_igs_regions <- function(gff_object, species_name) {
  
  # Filter for gene features to define genic regions
  gene_features <- gff_object[GenomicRanges::mcols(gff_object)$type == "gene"]
  
  if (length(gene_features) == 0) {
    log_message(sprintf("No gene features found for IGS calculation in %s", species_name), level = "warning")
    return(data.frame())
  }
  
  # Use GenomicRanges::gaps() to find intergenic regions
  # This is the gold standard method for calculating gaps between features
  igs_ranges <- GenomicRanges::gaps(gene_features)
  
  if (length(igs_ranges) == 0) {
    log_message(sprintf("No IGS regions calculated for %s", species_name), level = "warning")
    return(data.frame())
  }
  
  # Convert to data frame
  igs_df <- data.frame(
    species = species_name,
    region_name = paste0("IGS_", seq_along(igs_ranges)),
    region_type = "IGS",
    region_start = GenomicRanges::start(igs_ranges),
    region_end = GenomicRanges::end(igs_ranges),
    region_length = GenomicRanges::width(igs_ranges),
    chromosome = as.character(GenomicRanges::seqnames(igs_ranges)),
    strand = "*",  # IGS regions are non-stranded
    gene_name = NA_character_,  # IGS regions don't belong to specific genes
    total_cds_length = NA_real_,  # Not applicable for IGS
    stringsAsFactors = FALSE
  )
  
  return(igs_df)
}

#' Calculate Intron regions using GenomicRanges::gaps() within genes
#' 
#' @param gff_object Validated GRanges object
#' @param species_name Species name
#' @return Data frame with Intron regions
calculate_intron_regions <- function(gff_object, species_name) {
  
  # Filter for exon features
  exon_features <- gff_object[GenomicRanges::mcols(gff_object)$type == "exon"]
  
  if (length(exon_features) == 0) {
    log_message(sprintf("No exon features found for intron calculation in %s", species_name), level = "warning")
    return(data.frame())
  }
  
  # Group exons by gene and calculate gaps (introns) within each gene
  intron_list <- list()
  gene_ids <- unique(GenomicRanges::mcols(exon_features)$Parent %||% 
                     GenomicRanges::mcols(exon_features)$gene_id)
  
  if (is.null(gene_ids) || length(gene_ids) == 0) {
    log_message(sprintf("No gene grouping information found for intron calculation in %s", species_name), level = "warning")
    return(data.frame())
  }
  
  for (gene_id in gene_ids) {
    if (is.na(gene_id)) next
    
    # Get exons for this gene
    gene_exons <- exon_features[GenomicRanges::mcols(exon_features)$Parent == gene_id | 
                               GenomicRanges::mcols(exon_features)$gene_id == gene_id]
    
    if (length(gene_exons) < 2) next  # Need at least 2 exons to have introns
    
    # Calculate introns (gaps between exons)
    gene_introns <- GenomicRanges::gaps(gene_exons)
    
    if (length(gene_introns) > 0) {
      intron_list[[gene_id]] <- gene_introns
    }
  }
  
  if (length(intron_list) == 0) {
    log_message(sprintf("No introns calculated for %s", species_name), level = "info")
    return(data.frame())
  }
  
  # Combine all introns and convert to data frame
  all_introns <- do.call(c, intron_list)
  
  intron_df <- data.frame(
    species = species_name,
    region_name = paste0("Intron_", seq_along(all_introns)),
    region_type = "intron",
    region_start = GenomicRanges::start(all_introns),
    region_end = GenomicRanges::end(all_introns),
    region_length = GenomicRanges::width(all_introns),
    chromosome = as.character(GenomicRanges::seqnames(all_introns)),
    strand = as.character(GenomicRanges::strand(all_introns)),
    gene_name = rep(names(intron_list), lengths(intron_list)),
    total_cds_length = NA_real_,  # Not applicable for introns
    stringsAsFactors = FALSE
  )
  
  return(intron_df)
}

#' Extract genome structure information (LSC/IR/SSC regions)
#' 
#' @param gff_object Validated GRanges object
#' @param species_name Species name
#' @return Data frame with genome structure information
extract_genome_structure <- function(gff_object, species_name) {
  
  # Look for region features that define genome structure
  region_features <- gff_object[GenomicRanges::mcols(gff_object)$type == "region" |
                               GenomicRanges::mcols(gff_object)$type == "misc_feature"]
  
  if (length(region_features) == 0) {
    # Create a basic genome structure based on the total sequence length
    total_length <- max(GenomicRanges::end(gff_object))
    
    genome_structure <- data.frame(
      species = species_name,
      lsc_start = 1,
      lsc_end = total_length,
      ir_start = NA_integer_,
      ir_end = NA_integer_,
      ssc_start = NA_integer_,
      ssc_end = NA_integer_,
      total_length = total_length,
      special_handling = FALSE,
      stringsAsFactors = FALSE
    )
    
    log_message(sprintf("Generated basic genome structure for %s (total length: %d)", 
                       species_name, total_length), level = "info")
    
    return(genome_structure)
  }
  
  # Extract region information from GFF annotations
  # This is a simplified extraction - may need refinement based on actual GFF content
  total_length <- max(GenomicRanges::end(gff_object))
  
  genome_structure <- data.frame(
    species = species_name,
    lsc_start = 1,
    lsc_end = total_length,
    ir_start = NA_integer_,
    ir_end = NA_integer_,
    ssc_start = NA_integer_,
    ssc_end = NA_integer_,
    total_length = total_length,
    special_handling = FALSE,
    stringsAsFactors = FALSE
  )
  
  return(genome_structure)
}

#' Generate unified annotation files
#' 
#' Creates the final unified CSV files that replace the legacy CSV dependencies.
#' According to Core Architect instruction: merge region info and CDS lengths
#' into a single unified_region_annotations.csv file.
#' 
#' @param extracted_features List of extracted features for all species
#' @param output_dir Output directory
#' @param config Configuration list
#' @return List of created file paths
generate_unified_annotations <- function(extracted_features, output_dir, config) {
  
  log_message("Generating unified annotation files")
  
  # Combine all region features into a single data frame
  all_regions <- list()
  all_genome_structures <- list()
  
  for (species_name in names(extracted_features)) {
    species_features <- extracted_features[[species_name]]
    
    # Combine CDS, IGS, and Intron features
    combined_regions <- rbind(
      species_features$cds,
      species_features$igs,
      species_features$introns
    )
    
    all_regions[[species_name]] <- combined_regions
    all_genome_structures[[species_name]] <- species_features$genome_structure
  }
  
  # Create unified region annotations file
  unified_regions <- do.call(rbind, all_regions)
  unified_regions_file <- file.path(output_dir, "unified_region_annotations.csv")
  
  # Ensure column order matches legacy expectations
  expected_columns <- c("species", "region_name", "region_type", "region_start", 
                       "region_end", "region_length", "chromosome", "strand", 
                       "gene_name", "total_cds_length")
  
  # Add missing columns if needed
  for (col in expected_columns) {
    if (!col %in% names(unified_regions)) {
      unified_regions[[col]] <- NA
    }
  }
  
  # Reorder columns
  unified_regions <- unified_regions[, expected_columns]
  
  # Write unified region annotations
  write.csv(unified_regions, unified_regions_file, row.names = FALSE)
  log_message(sprintf("Created unified region annotations: %s (%d entries)", 
                     unified_regions_file, nrow(unified_regions)))
  
  # Create unified genome structures file
  unified_genome_structures <- do.call(rbind, all_genome_structures)
  unified_genome_file <- file.path(output_dir, "unified_genome_structures.csv")
  
  write.csv(unified_genome_structures, unified_genome_file, row.names = FALSE)
  log_message(sprintf("Created unified genome structures: %s (%d entries)", 
                     unified_genome_file, nrow(unified_genome_structures)))
  
  # Return list of created files
  created_files <- list(
    unified_region_annotations = unified_regions_file,
    unified_genome_structures = unified_genome_file
  )
  
  return(created_files)
}

#' Log summary statistics for GFF processing results
#' 
#' @param results Processing results object
log_summary_statistics <- function(results) {
  
  log_message("=== GFF PROCESSING SUMMARY STATISTICS ===")
  log_message(sprintf("Pipeline Version: %s", results$pipeline_version))
  log_message(sprintf("Species Processed: %d (%s)", 
                     length(results$species_processed),
                     paste(results$species_processed, collapse = ", ")))
  log_message(sprintf("Processing Time: %s", results$processing_time))
  
  # Count features by type
  total_cds <- 0
  total_igs <- 0
  total_introns <- 0
  
  for (species_features in results$extracted_features) {
    total_cds <- total_cds + nrow(species_features$cds)
    total_igs <- total_igs + nrow(species_features$igs)
    total_introns <- total_introns + nrow(species_features$introns)
  }
  
  log_message(sprintf("Features Extracted: CDS=%d, IGS=%d, Introns=%d (Total=%d)", 
                     total_cds, total_igs, total_introns, 
                     total_cds + total_igs + total_introns))
  
  log_message("Generated Files:")
  for (file_type in names(results$unified_files)) {
    log_message(sprintf("  %s: %s", file_type, results$unified_files[[file_type]]))
  }
  
  log_message("=== GFF PROCESSING COMPLETED SUCCESSFULLY ===")
}

#' Null coalescing operator
#' 
#' @description Returns the left-hand side if not NULL, otherwise returns the right-hand side.
#' This is a standard null-coalescing operator commonly used in R packages.
#' 
#' @param x First value to test (left-hand side)
#' @param y Alternative value if \code{x} is NULL
#' @return Either \code{x} (if not NULL) or \code{y}
#' 
#' @examples
#' # Basic usage
#' NULL %||% "default"        # Returns "default"
#' "value" %||% "default"     # Returns "value"
#' 
#' # Common use case
#' user_input <- NULL
#' final_value <- user_input %||% "default_value"
#' 
#' @keywords internal
#' @name null.coalesce
#' @export
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}