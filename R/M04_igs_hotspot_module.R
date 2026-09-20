############################################################
#### M04_igs_hotspot_module.R - M04 IGS Hotspot Module ####
############################################################
#
# Lightweight IGS-specific hotspot analysis module
# Focused exclusively on generating IGS variation frequency heatmaps
# Excludes PCA, Upset, correlation, and enrichment analysis
#
# Part of the M04 IGS Hotspot Analysis Module
#
############################################################

#' M04 IGS Hotspot Analysis Module Orchestrator
#' 
#' @title Run M04 IGS-specific hotspot heatmap analysis
#' @description Generates a focused heatmap visualization of IGS (Intergenic Spacer)
#' variation frequencies across species. This module provides a streamlined,
#' IGS-specific alternative to the comprehensive M03 hotspot analysis.
#' 
#' @param normalized_data Normalized frequency data frame with columns:
#'   - species: Species identifier (character)
#'   - region_type: Genomic region type (character, must include "IGS")
#'   - gene: Gene/region identifier (character)
#'   - frequency_per_kb: Normalized frequency values (numeric)
#'   - var_type: Variant type (character)
#'   Additional columns are preserved
#' @param config Configuration object containing analysis parameters
#' @param output_dir Output directory for generated plots and data
#' @param task_name Task identifier for output organization
#' @param task_params Task-specific parameters including data_scoping configuration
#' 
#' @return List containing analysis results:
#'   - success: Boolean indicating successful completion
#'   - plots_generated: Number of plots created (typically 1)
#'   - igs_regions_analyzed: Number of IGS regions included
#'   - species_count: Number of species analyzed
#'   - output_files: List of generated file paths
#'   - error_message: Error details if analysis failed
#' 
#' @details
#' The M04 module implements a focused IGS analysis pipeline:
#' 1. **Strict IGS filtering**: Applies data_scoping to retain only region_type == "IGS"
#' 2. **Data reshaping**: Transforms data into species × IGS_region frequency matrix
#' 3. **Heatmap generation**: Calls generate_global_frequency_heatmap with IGS-specific title
#' 4. **Quality logging**: Records analysis statistics and file locations
#' 
#' This module maintains consistency with M03's heatmap styling while providing
#' a dedicated focus on intergenic spacer regions. It explicitly excludes all
#' gene-centric analyses (PCA, correlation, enrichment) to maintain functional clarity.
#' 
#' @examples
#' \dontrun{
#' # Run M04 IGS analysis
#' result <- run_m04_igs_hotspot(
#'   normalized_data = frequency_data,
#'   config = analysis_config,
#'   output_dir = session_paths$plots,
#'   task_name = "M04_igs_standard",
#'   task_params = igs_task_params
#' )
#' 
#' if (result$success) {
#'   cat("IGS heatmap generated:", result$output_files$heatmap, "\n")
#' }
#' }
#' 
#' @importFrom dplyr filter select
#' @importFrom tidyr pivot_wider
#' @export
run_m04_igs_hotspot <- function(normalized_data,
                                config,
                                output_dir,
                                task_name,
                                task_params) {
  
  log_message(sprintf("=== M04 IGS HOTSPOT ANALYSIS: %s ===", task_name))
  log_message("Focused IGS variation frequency heatmap generation")
  
  # Initialize results structure
  analysis_results <- list(
    success = FALSE,
    plots_generated = 0,
    igs_regions_analyzed = 0,
    species_count = 0,
    output_files = list(),
    error_message = NULL
  )
  
  tryCatch({
    
    # Validate input data
    if (is.null(normalized_data) || nrow(normalized_data) == 0) {
      stop("M04 analysis cannot proceed with empty normalized_data")
    }
    
    log_message(sprintf("Input data: %d rows across %d species", 
                       nrow(normalized_data), 
                       length(unique(normalized_data$species))))
    
    # =============================================
    # STEP 1: Apply strict IGS data_scoping filters
    # =============================================
    
    log_message("STEP 1: Applying IGS-specific data_scoping filters")
    
    # Apply data scoping to filter for IGS regions only
    scoped_data <- apply_data_scoping(normalized_data, task_params)
    
    if (is.null(scoped_data) || nrow(scoped_data) == 0) {
      stop("No IGS data remaining after applying data_scoping filters")
    }
    
    # Verify IGS filtering was successful
    region_types_found <- unique(scoped_data$region_type)
    if (length(region_types_found) > 1 || !"IGS" %in% region_types_found) {
      log_message(sprintf("Warning: Expected only IGS regions, found: %s", 
                         paste(region_types_found, collapse = ", ")), level = "warning")
    }
    
    log_message(sprintf("IGS data after scoping: %d rows across %d species", 
                       nrow(scoped_data), 
                       length(unique(scoped_data$species))))
    
    log_message(sprintf("Region types in filtered data: %s", 
                       paste(region_types_found, collapse = ", ")))
    
    # =============================================
    # STEP 2: Clean and prepare IGS data
    # =============================================
    
    log_message("STEP 2: Cleaning and preparing IGS frequency data")
    
    # Filter out NA frequencies but keep zero frequencies for comprehensive analysis
    # Validate IGS rows with region_name rather than gene
    cleaned_data <- scoped_data %>% 
      dplyr::filter(!is.na(frequency_per_kb)) %>%
      dplyr::filter(!is.na(region_name) & !is.na(region_type)) %>%  # Use region_name for IGS
      dplyr::filter(!is.na(species))  # Ensure species information is complete
    
    if (nrow(cleaned_data) == 0) {
      stop("No valid IGS data after cleaning - all frequencies may be NA")
    }
    
    log_message(sprintf("IGS data after cleaning: %d -> %d valid observations", 
                       nrow(scoped_data), nrow(cleaned_data)))
    
    # Count unique IGS regions and species
    unique_igs_regions <- unique(cleaned_data$region_name)  # Use region_name for IGS identification
    unique_species <- unique(cleaned_data$species)
    
    analysis_results$igs_regions_analyzed <- length(unique_igs_regions)
    analysis_results$species_count <- length(unique_species)
    
    log_message(sprintf("IGS analysis scope: %d unique IGS regions across %d species",
                       length(unique_igs_regions), length(unique_species)))
    
    # =============================================
    # STEP 3: Apply species ordering
    # =============================================
    
    log_message("STEP 3: Applying species ordering for consistent visualization")
    
    # Apply species ordering for biological consistency
    species_ordering_result <- load_and_apply_species_ordering(
      data = cleaned_data,
      session_id = config$session_info$session_id,
      species_col = "species"
    )
    
    cleaned_data <- species_ordering_result$data
    log_message(species_ordering_result$message)
    
    # =============================================
    # STEP 4: Create output directory and save intermediate data
    # =============================================
    
    log_message("STEP 4: Setting up output directory and saving intermediate data")
    
    # Create M04-specific output directory
    m04_output_dir <- file.path(output_dir, "M04_igs_hotspot", task_name)
    dir.create(m04_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save cleaned IGS data for transparency and debugging
    cleaned_data_file <- file.path(m04_output_dir, "M04_igs_cleaned_data.csv")
    write.csv(cleaned_data, cleaned_data_file, row.names = FALSE)
    log_message(sprintf("Saved cleaned IGS data: %s", cleaned_data_file))
    
    # Generate IGS analysis quality report
    quality_report <- data.frame(
      analysis_type = "M04_IGS_Hotspot",
      total_observations = nrow(cleaned_data),
      unique_species = length(unique_species),
      unique_igs_regions = length(unique_igs_regions),
      region_types = paste(unique(cleaned_data$region_type), collapse = ", "),
      var_types = paste(unique(cleaned_data$var_type), collapse = ", "),
      mean_frequency = mean(cleaned_data$frequency_per_kb, na.rm = TRUE),
      median_frequency = median(cleaned_data$frequency_per_kb, na.rm = TRUE),
      frequency_range = paste(range(cleaned_data$frequency_per_kb, na.rm = TRUE), collapse = " - "),
      analysis_timestamp = Sys.time(),
      stringsAsFactors = FALSE
    )
    
    quality_report_file <- file.path(m04_output_dir, "M04_igs_quality_report.csv")
    write.csv(quality_report, quality_report_file, row.names = FALSE)
    log_message(sprintf("Saved IGS quality report: %s", quality_report_file))
    
    # =============================================
    # STEP 5: Generate IGS frequency heatmap
    # =============================================
    
    log_message("STEP 5: Generating IGS variation frequency heatmap")
    
    # Generate file path for IGS heatmap
    heatmap_filename <- sprintf("M04_igs_frequency_heatmap_%s.pdf", 
                               format(Sys.time(), "%Y%m%d_%H%M%S"))
    heatmap_file_path <- file.path(m04_output_dir, heatmap_filename)
    
    # Call the global frequency heatmap function with IGS-specific configuration and direct file save
    # This reuses M03's proven heatmap logic while maintaining IGS focus
    # The filename parameter is passed directly to pheatmap() for automatic saving
    heatmap_result <- generate_global_frequency_heatmap(
      normalized_data = cleaned_data,
      config = config,
      task_params = task_params,
      filename = heatmap_file_path
    )
    
    if (is.null(heatmap_result)) {
      stop("Failed to generate IGS frequency heatmap - function returned NULL")
    }
    
    # Track successful plot generation
    analysis_results$plots_generated <- 1
    analysis_results$output_files$heatmap_file <- heatmap_file_path
    
    log_message(sprintf("IGS heatmap successfully saved to: %s", heatmap_file_path))
    
    log_message("IGS frequency heatmap generated successfully")
    log_message(sprintf("Heatmap includes %d IGS regions across %d species", 
                       length(unique_igs_regions), length(unique_species)))
    
    # =============================================
    # STEP 6: Finalize results and create manifest
    # =============================================
    
    log_message("STEP 6: Finalizing M04 analysis results")
    
    # Record output files
    analysis_results$output_files <- list(
      cleaned_data = cleaned_data_file,
      quality_report = quality_report_file,
      output_directory = m04_output_dir,
      heatmap_result = heatmap_result  # Contains plot object and metadata
    )
    
    # Mark analysis as successful
    analysis_results$success <- TRUE
    
    log_message(sprintf("M04 IGS analysis completed successfully:"))
    log_message(sprintf("  - IGS regions analyzed: %d", analysis_results$igs_regions_analyzed))
    log_message(sprintf("  - Species included: %d", analysis_results$species_count))
    log_message(sprintf("  - Plots generated: %d", analysis_results$plots_generated))
    log_message(sprintf("  - Output directory: %s", m04_output_dir))
    
  }, error = function(e) {
    
    # Log error and record failure
    analysis_results$success <- FALSE
    analysis_results$error_message <- e$message
    log_message(sprintf("M04 IGS analysis failed: %s", e$message), level = "error")
    
  })
  
  return(analysis_results)
}