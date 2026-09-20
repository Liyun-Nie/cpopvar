############################################################
#### M05_hotspot_ideogram_module.R - M05 Hotspot Ideogram Visualization ####
############################################################
#
# M05: Chromosome-level Hotspot Visualization Module
# Integrates CDS hotspots (from M03) and poiGS hotspots (from M04)
# into unified chromosome ideogram plots
#
# Design Philosophy:
# - Single responsibility: Hotspot spatial visualization
# - Clear dependencies: Requires M03 and M04 outputs
# - Thin orchestrator: Delegates to generate_chromosome_hotspot_ideogram()
#
# Module independence:
# - M01-M04: Independent analysis modules
# - M05: Integrative visualization module
#
############################################################

#' M05 Hotspot Ideogram Module Orchestrator
#' 
#' @title Run M05 Chromosome Hotspot Ideogram Visualization
#' @description Generates chromosome ideogram plots showing spatial distribution
#' of CDS and poiGS hotspots. Integrates results from M03 (core shared CDS hotspots)
#' and M04 (core shared poiGS hotspots) to reveal co-localization patterns.
#' 
#' @param normalized_data Normalized frequency data (not used, for API consistency)
#' @param config Configuration object containing session info and input file paths
#' @param output_dir Output directory for plots and data files
#' @param task_name Task name for output organization (default: "M05_hotspot_ideogram")
#' @param task_params Task-specific parameters from configuration
#' 
#' @return List containing:
#'   \item{success}{Logical indicating if analysis completed successfully}
#'   \item{plots_generated}{Number of plots generated}
#'   \item{manifest}{List of output files for report generation}
#'   \item{error_message}{Error message if analysis failed}
#' 
#' @details
#' This module performs the following steps:
#' 1. Validates M03 and M04 output files exist
#' 2. Loads region_info_complete.csv for coordinate information
#' 3. Extracts hotspot coordinates for CDS and poiGS regions
#' 4. Creates chromosome karyotype data
#' 5. Generates RIdeogram visualization with hotspot markers
#' 
#' Dependencies:
#' - M03_hotspot: Must complete successfully to generate M03_core_shared_hotspots.csv
#' - M04_igs_hotspot: Must complete successfully to generate M04_core_shared_poigs_hotspots.csv
#' 
#' @export
run_M05_hotspot_ideogram <- function(normalized_data = NULL,
                                     config,
                                     output_dir,
                                     task_name = "M05_hotspot_ideogram",
                                     task_params = list()) {
  
  log_message("=== M05 HOTSPOT IDEOGRAM VISUALIZATION ===")
  log_message("Integrating M03 CDS hotspots and M04 poiGS hotspots for chromosome-level visualization")
  
  # Initialize results structure
  results <- list(
    success = FALSE,
    plots_generated = 0,
    manifest = list(),
    error_message = NULL
  )
  
  tryCatch({
    
    # ================================================================
    # STEP 1: Get session paths and validate dependencies
    # ================================================================
    
    log_message("STEP 1: Validating M03 and M04 dependencies")
    
    # Get session info
    session_id <- config$session_info$session_id %||% "unknown"
    session_paths <- get_session_paths(session_id)
    plots_dir <- session_paths$plots
    
    # Auto-detect M03 and M04 hotspot files
    m03_hotspots_path <- file.path(plots_dir, "M03_hotspot", "M03_hotspot", 
                                   "M03_core_shared_hotspots.csv")
    m04_hotspots_path <- file.path(plots_dir, "M04_poigs_hotspot_engine", 
                                   "M04_igs_hotspot", "M04_core_shared_poigs_hotspots.csv")
    
    log_message(sprintf("M05: Checking M03 hotspots at: %s", m03_hotspots_path))
    log_message(sprintf("M05: Checking M04 hotspots at: %s", m04_hotspots_path))
    
    # Validate input files exist
    m03_exists <- file.exists(m03_hotspots_path)
    m04_exists <- file.exists(m04_hotspots_path)
    
    log_message(sprintf("M05: M03 file exists: %s", m03_exists))
    log_message(sprintf("M05: M04 file exists: %s", m04_exists))
    
    if (!m03_exists) {
      stop(sprintf(
        "M03 hotspots file not found: %s\n\nPlease ensure:\n1. M03_hotspot task is enabled in configuration\n2. M03_hotspot has completed successfully\n3. M05 is executed AFTER M03 in the task sequence", 
        m03_hotspots_path
      ))
    }
    
    if (!m04_exists) {
      stop(sprintf(
        "M04 hotspots file not found: %s\n\nPlease ensure:\n1. M04_igs_hotspot task is enabled in configuration\n2. M04_igs_hotspot has completed successfully\n3. M05 is executed AFTER M04 in the task sequence", 
        m04_hotspots_path
      ))
    }
    
    log_message("M05: Dependency validation successful - both M03 and M04 outputs found")
    
    # ================================================================
    # STEP 2: Create output directory
    # ================================================================
    
    log_message("STEP 2: Creating output directory")
    
    ideogram_output_dir <- file.path(output_dir, "M05_hotspot_ideogram")
    dir.create(ideogram_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    log_message(sprintf("M05: Output directory: %s", ideogram_output_dir))
    
    # ================================================================
    # STEP 3: Call chromosome hotspot ideogram generation
    # ================================================================
    
    log_message("STEP 3: Generating chromosome hotspot ideogram")
    log_message("M05: Delegating to generate_chromosome_hotspot_ideogram()")
    
    # Call the existing comprehensive function
    ideogram_results <- generate_chromosome_hotspot_ideogram(
      m03_hotspots_path = m03_hotspots_path,
      m04_hotspots_path = m04_hotspots_path,
      config = config,
      output_dir = ideogram_output_dir,
      session_id = session_id,
      task_params = task_params  # V20 FIX: Pass task_params for custom_label support
    )
    
    # ================================================================
    # STEP 4: Process results
    # ================================================================
    
    if (ideogram_results$success) {
      log_message("=== M05 HOTSPOT IDEOGRAM COMPLETED SUCCESSFULLY ===")
      log_message(sprintf("M05: Generated %d plots", ideogram_results$plots_generated))
      
      # Update results
      results$success <- TRUE
      results$plots_generated <- ideogram_results$plots_generated
      results$manifest <- ideogram_results$manifest
      
      # Log manifest entries
      if (length(results$manifest) > 0) {
        log_message(sprintf("M05: Created %d manifest entries for report generation", 
                           length(results$manifest)))
      }
      
    } else {
      warning_msg <- sprintf("M05 encountered errors: %s", 
                           ideogram_results$error_message %||% "Unknown error")
      log_message(warning_msg, level = "warning")
      warning(warning_msg)
      
      results$error_message <- ideogram_results$error_message
    }
    
  }, error = function(e) {
    error_msg <- sprintf("M05 Hotspot Ideogram Module failed: %s", e$message)
    log_message(error_msg, level = "error")
    
    results$success <- FALSE
    results$error_message <- e$message
    
    # Re-throw error to ensure task dispatcher is aware of failure
    stop(error_msg)
  })
  
  return(results)
}

