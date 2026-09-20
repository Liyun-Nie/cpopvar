############################################################
#### run_m04_igs_analysis.R - M04 IGS Hotspot Analysis Module ####
############################################################
#
# Advanced M04 Position-based Orthologous IGS (poiGS) Hotspot Discovery Engine
# Transforms species-specific IGS regions into cross-species comparable units
# Implements sophisticated hotspot filtering and information-dense visualization
#
# Part of the M04 IGS Hotspot Analysis Module
#
############################################################

#' M04 Position-based Orthologous IGS Hotspot Analysis Engine
#' 
#' @title Run M04 poiGS-based IGS hotspot discovery and analysis
#' @description Implements an advanced IGS hotspot discovery engine that solves the
#' fundamental problem of IGS sequence incomparability across species. Uses position-based
#' orthologous IGS (poiGS) approach to identify evolutionarily conserved high-variation
#' IGS regions with strong biological significance.
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
#'   - plots_generated: Number of plots created
#'   - poigs_regions_analyzed: Number of poiGS regions identified
#'   - core_hotspots_identified: Number of core shared hotspots
#'   - species_count: Number of species analyzed
#'   - output_files: List of generated file paths
#'   - error_message: Error details if analysis failed
#' 
#' @details
#' The M04 engine implements a sophisticated 4-step poiGS analysis pipeline:
#' 
#' **Step 1: Position-based Orthologous IGS (poiGS) Identification**
#' - Loads IGS region position information from region_length_statistics
#' - Uses flanking gene pairs to define cross-species orthologous IGS units
#' - Generates standardized poiGS identifiers (e.g., "poiGS_psbA-trnH")
#' - Creates comprehensive poiGS mapping table
#' 
#' **Step 2: Frequency Data Aggregation by poiGS**  
#' - Links variant frequency data with poiGS mapping
#' - Aggregates frequencies by species and poiGS (handles multiple IGS per poiGS)
#' - Produces cross-species comparable frequency dataset
#' 
#' **Step 3: Two-tier Hotspot Filtering**
#' - **Candidate filtering**: Identifies high-frequency poiGS using M03 logic
#' - **Core hotspot filtering**: Retains poiGS shared across multiple species
#' - Applies configurable thresholds for biological significance
#' 
#' **Step 4: Information-dense Visualization**
#' - Generates focused heatmap containing only core shared hotspots
#' - Uses optimized visualization settings for maximum information density
#' - Produces publication-ready output with biological annotations
#' 
#' Key innovations over basic M04:
#' - Solves IGS cross-species comparability through poiGS approach
#' - Implements rigorous hotspot filtering for biological significance  
#' - Generates information-dense rather than sparse visualizations
#' - Provides comprehensive data provenance and quality metrics
#' 
#' @examples
#' \dontrun{
#' # Run advanced M04 poiGS analysis
#' result <- run_m04_igs_analysis(
#'   normalized_data = frequency_data,
#'   config = analysis_config,
#'   output_dir = session_paths$plots,
#'   task_name = "M04_poigs_hotspot_discovery",
#'   task_params = igs_task_params
#' )
#' 
#' if (result$success) {
#'   cat("Core IGS hotspots identified:", result$core_hotspots_identified, "\n")
#'   cat("Final heatmap:", result$output_files$core_heatmap, "\n")
#' }
#' }
#' 
#' @importFrom dplyr filter select left_join group_by summarise ungroup
#' @importFrom readr read_csv
#' @export
run_m04_igs_analysis <- function(normalized_data,
                                config,
                                output_dir,
                                task_name,
                                task_params) {
  
  log_message(sprintf("=== M04 POSITION-BASED ORTHOLOGOUS IGS HOTSPOT ENGINE: %s ===", task_name))
  log_message("Advanced poiGS-based IGS hotspot discovery with biological significance filtering")
  
  # Initialize comprehensive results structure
  analysis_results <- list(
    success = FALSE,
    plots_generated = 0,
    poigs_regions_analyzed = 0,
    core_hotspots_identified = 0,
    species_count = 0,
    output_files = list(),
    error_message = NULL
  )
  
  tryCatch({
    
    # =============================================
    # INPUT VALIDATION AND SETUP
    # =============================================
    
    log_message("M04 Engine: Validating inputs and initializing analysis environment")
    
    if (is.null(normalized_data) || nrow(normalized_data) == 0) {
      stop("M04 poiGS analysis requires valid normalized frequency data")
    }
    
    # Create M04-specific output directory structure
    m04_output_dir <- file.path(output_dir, "M04_poigs_hotspot_engine", task_name)
    dir.create(m04_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Read plot_types parameter (aligned with M03 architecture)
    plot_types <- get_task_parameter(task_params, config, "plot_types", 
                                     c("global_heatmap", "core_heatmap", "upset", "statistical_test"))
    
    log_message(sprintf("M04 Analysis Environment:"))
    log_message(sprintf("  - Input data: %d rows across %d species", 
                       nrow(normalized_data), length(unique(normalized_data$species))))
    log_message(sprintf("  - Output directory: %s", m04_output_dir))
    log_message(sprintf("  - Plot types enabled: %s", paste(plot_types, collapse = ", ")))
    
    # =============================================
    # STEP 1: POSITION-BASED ORTHOLOGOUS IGS (poiGS) IDENTIFICATION
    # =============================================
    
    log_message("M04 Step 1: Identifying Position-based Orthologous IGS (poiGS) regions")
    
    # V3.16.14: Use new region_info_complete.csv from P01_preprocessed
    session_id <- config$session_info$session_id
    p01_dir <- get_processing_stage_path(session_id, "P01_preprocessed")
    region_info_file <- file.path(p01_dir, "region_info_complete.csv")
    
    if (!file.exists(region_info_file)) {
      stop(sprintf("M04 Engine: Cannot find region_info_complete.csv. Expected: %s", region_info_file))
    }
    
    log_message(sprintf("Loading IGS region information from: %s", region_info_file))
    
    # Load and filter IGS region data
    # V3.16.14: Load from region_info_complete.csv and rename columns
    igs_region_info <- readr::read_csv(region_info_file, show_col_types = FALSE) %>%
      dplyr::filter(region_type == "IGS") %>%
      # Rename columns to match expected format for identify_positional_orthologs
      dplyr::rename(
        Region_Name = region_name,
        Type = region_type,
        Start_Position = region_start,
        End_Position = region_end,
        Length = region_length
      ) %>%
      dplyr::select(species, Region_Name, Type, Start_Position, End_Position, Length)
    
    if (nrow(igs_region_info) == 0) {
      stop("M04 Engine: No intergenic regions found in region statistics file")
    }
    
    log_message(sprintf("Loaded %d IGS regions across %d species for poiGS identification", 
                       nrow(igs_region_info), length(unique(igs_region_info$species))))
    
    # Generate poiGS mapping using the helper function
    log_message("Generating position-based orthologous IGS mapping...")
    poigs_mapping <- identify_positional_orthologs(igs_region_info)
    
    # Save poiGS mapping table
    poigs_mapping_file <- file.path(m04_output_dir, "M04_poigs_mapping.csv")
    write.csv(poigs_mapping, poigs_mapping_file, row.names = FALSE)
    log_message(sprintf("Saved poiGS mapping: %s", poigs_mapping_file))
    
    analysis_results$poigs_regions_analyzed <- length(unique(poigs_mapping$poiGS_ID))
    
    # =============================================
    # STEP 2: FREQUENCY DATA AGGREGATION BY poiGS
    # =============================================
    
    log_message("M04 Step 2: Aggregating variant frequencies by poiGS units")
    
    # Apply data scoping to get IGS-only frequency data
    scoped_data <- apply_data_scoping(normalized_data, task_params)
    
    if (is.null(scoped_data) || nrow(scoped_data) == 0) {
      stop("M04 Engine: No IGS frequency data remaining after data scoping")
    }
    
    log_message(sprintf("IGS frequency data after scoping: %d rows", nrow(scoped_data)))
    
    # Link frequency data with poiGS mapping
    # Join on region_name, which is the IGS identifier in this context
    frequency_with_poigs <- scoped_data %>%
      dplyr::left_join(poigs_mapping, 
                      by = c("species" = "species", "region_name" = "Region_Name")) %>%
      dplyr::filter(!is.na(poiGS_ID))  # Remove entries that couldn't be mapped
    
    if (nrow(frequency_with_poigs) == 0) {
      stop("M04 Engine: No frequency data could be linked with poiGS mapping")
    }
    
    log_message(sprintf("Successfully linked %d frequency observations with poiGS mapping", 
                       nrow(frequency_with_poigs)))
    
    # ========================================================================
    # TWO-STAGE FREQUENCY AGGREGATION FOR poiGS
    # ========================================================================
    # Stage 1: Aggregate multiple original IGS regions into single poiGS unit
    # (Handle rare cases where one poiGS maps to multiple IGS regions)
    poigs_frequencies_stage1 <- frequency_with_poigs %>%
      dplyr::group_by(species, poiGS_ID, var_type, region_type) %>%
      dplyr::summarise(
        frequency_per_kb = mean(frequency_per_kb, na.rm = TRUE),
        variant_count = sum(variant_count, na.rm = TRUE),
        region_count = dplyr::n(),
        .groups = "drop"
      )
    
    log_message(sprintf("Stage 1 aggregation complete: %d observations", nrow(poigs_frequencies_stage1)))
    
    # Stage 2: Merge all variant types for each poiGS
    # This is the KEY FIX - combine SNP, INDEL, complex, etc. into total poiGS frequency
    poigs_frequencies <- poigs_frequencies_stage1 %>%
      dplyr::group_by(species, poiGS_ID, region_type) %>%
      dplyr::summarise(
        frequency_per_kb = sum(frequency_per_kb, na.rm = TRUE),  # SUM all var_type frequencies
        variant_count = sum(variant_count, na.rm = TRUE),
        region_count = max(region_count),  # Preserve original IGS count
        var_type = paste(sort(unique(var_type)), collapse = "+"),  # Record which types present (e.g., "INDEL+snp")
        .groups = "drop"
      )
    # NOTE: No longer adding "gene" column - M04 uses poiGS_ID directly via feature_col parameter
    
    # Save aggregated poiGS frequencies
    poigs_freq_file <- file.path(m04_output_dir, "M04_poigs_frequencies.csv")
    write.csv(poigs_frequencies, poigs_freq_file, row.names = FALSE)
    log_message(sprintf("Saved poiGS aggregated frequencies: %s", poigs_freq_file))
    
    log_message(sprintf("Stage 2 aggregation complete (var_type merged): %d observations", nrow(poigs_frequencies)))
    log_message(sprintf("  - Unique poiGS units: %d", length(unique(poigs_frequencies$poiGS_ID))))
    log_message(sprintf("  - Species coverage: %d", length(unique(poigs_frequencies$species))))
    log_message(sprintf("  - Frequency range after merging: %.3f to %.3f per kb", 
                       min(poigs_frequencies$frequency_per_kb), 
                       max(poigs_frequencies$frequency_per_kb)))
    
    analysis_results$species_count <- length(unique(poigs_frequencies$species))
    
    # =============================================
    # STEP 3: TWO-TIER HOTSPOT FILTERING
    # =============================================
    
    log_message("M04 Step 3a: Identifying candidate poiGS hotspots using M03 logic")
    
    # [CRITICAL FIX] Extract and flatten the nested hotspot_parameters block
    # to ensure compatibility with the reused M03 function
    hotspot_params <- task_params$hotspot_parameters
    if (is.null(hotspot_params)) {
      stop("M04 Engine: `hotspot_parameters` block is missing in the YAML configuration for this task.")
    }
    
    # Apply M03's candidate hotspot identification with poiGS_ID as feature column
    candidate_hotspots <- identify_candidate_hotspots(
      normalized_data = poigs_frequencies,
      config = config,
      task_params = hotspot_params,  # Pass flattened parameters
      feature_col = "poiGS_ID"  # Use poiGS_ID instead of gene for biological clarity
    )
    
    # [CRITICAL FIX] identify_candidate_hotspots returns a data frame, not a list with $candidates
    if (is.null(candidate_hotspots) || nrow(candidate_hotspots) == 0) {
      stop("M04 Engine: Failed to identify candidate poiGS hotspots")
    }
    
    # Save COMPLETE threshold-annotated data (including non-hotspots)
    # This allows users to see ALL poiGS with their threshold values
    all_poigs_annotated_file <- file.path(m04_output_dir, "M04_all_poigs_with_thresholds.csv")
    write.csv(candidate_hotspots, all_poigs_annotated_file, row.names = FALSE)
    log_message(sprintf("Saved all poiGS with threshold annotations: %s", basename(all_poigs_annotated_file)))
    log_message(sprintf("  - Total poiGS observations: %d", nrow(candidate_hotspots)))
    log_message(sprintf("  - Hotspots: %d (%.1f%%)", 
                       sum(candidate_hotspots$is_potential_hotspot), 
                       100 * sum(candidate_hotspots$is_potential_hotspot) / nrow(candidate_hotspots)))
    
    # Filter for actual candidate hotspots (is_potential_hotspot == TRUE)
    candidate_hotspot_data <- candidate_hotspots %>%
      dplyr::filter(is_potential_hotspot == TRUE)
    
    log_message(sprintf("Candidate poiGS hotspots identified: %d", 
                       length(unique(candidate_hotspot_data$poiGS_ID))))
    
    # Save COMPLETE candidate hotspots data (before core filtering)
    # This allows users to see threshold-based filtering results
    candidate_hotspots_file <- file.path(m04_output_dir, "M04_candidate_poigs_hotspots.csv")
    write.csv(candidate_hotspot_data, candidate_hotspots_file, row.names = FALSE)
    log_message(sprintf("Saved candidate poiGS hotspots (threshold-filtered): %s", basename(candidate_hotspots_file)))
    log_message(sprintf("  - Total candidate entries: %d", nrow(candidate_hotspot_data)))
    log_message(sprintf("  - Unique candidate poiGS: %d", length(unique(candidate_hotspot_data$poiGS_ID))))
    
    log_message("M04 Step 3b: Filtering for core shared poiGS hotspots")
    
    # [CRITICAL FIX] Use M03-consistent dual-mode threshold calculation
    # Read threshold_mode from hotspot_parameters to determine filtering strategy
    threshold_mode <- get_task_parameter(hotspot_params, config, "threshold_mode", "absolute", log_default = TRUE)
    
    total_species <- analysis_results$species_count
    if (total_species == 0) {
      stop("M04 Engine: No species available to calculate species coverage threshold.")
    }
    
    # Calculate min_species_coverage ratio based on threshold_mode
    if (threshold_mode == "absolute") {
      # Use absolute number of species
      min_species_number <- get_task_parameter(hotspot_params, config, "min_species_number", 3, log_default = TRUE)
      min_species_coverage <- min_species_number / total_species
      log_message(sprintf("Using absolute threshold mode: min_species_number = %d (%.1f%% of %d species)", 
                         min_species_number, min_species_coverage * 100, total_species))
    } else if (threshold_mode == "ratio") {
      # Use species ratio
      min_species_coverage <- get_task_parameter(hotspot_params, config, "min_species_ratio", 0.25, log_default = TRUE)
      min_species_number <- ceiling(total_species * min_species_coverage)
      log_message(sprintf("Using ratio threshold mode: min_species_ratio = %.2f (%d species out of %d)", 
                         min_species_coverage, min_species_number, total_species))
    } else {
      stop(sprintf("M04 Engine: Invalid threshold_mode '%s'. Must be 'absolute' or 'ratio'.", threshold_mode))
    }
    
    # Apply second-tier filtering for core shared hotspots
    core_hotspots <- filter_core_hotspots(
      candidate_hotspots_data = candidate_hotspot_data,
      feature_col = "poiGS_ID",  # Use biologically meaningful column name
      species_col = "species",
      min_species_coverage = min_species_coverage
    )
    
    analysis_results$core_hotspots_identified <- length(unique(core_hotspots$poiGS_ID))
    
    log_message(sprintf("Core shared poiGS hotspots: %d (coverage >= %d species, %.1f%%)", 
                       analysis_results$core_hotspots_identified, 
                       min_species_number,
                       min_species_coverage * 100))
    
    # Modify core_hotspots to include filtering stage information
    # Calculate species coverage percentage for each core hotspot
    core_hotspots_with_metadata <- core_hotspots %>%
      dplyr::group_by(poiGS_ID) %>%
      dplyr::mutate(
        filtering_stage = "core_shared",
        min_species_threshold = min_species_number,
        species_count_for_hotspot = dplyr::n_distinct(species),
        species_coverage_pct = (species_count_for_hotspot / total_species) * 100
      ) %>%
      dplyr::ungroup()
    
    # Save with clearer filename indicating this is FINAL filtered result
    core_hotspots_file <- file.path(m04_output_dir, "M04_core_shared_hotspots_final.csv")
    write.csv(core_hotspots_with_metadata, core_hotspots_file, row.names = FALSE)
    log_message(sprintf("Saved FINAL core shared hotspots: %s", basename(core_hotspots_file)))
    log_message(sprintf("  - Core hotspots: %d poiGS", analysis_results$core_hotspots_identified))
    log_message(sprintf("  - Min species coverage: %d (%.1f%%)", min_species_number, min_species_coverage * 100))
    
    # =============================================
    # GENERATE CORE SHARED poiGS SUMMARY TABLES (aligned with M03 architecture)
    # =============================================
    
    log_message("Generating core shared poiGS summary tables (aligned with M03 output format)")
    
    # Summary Table 1: poiGS-level aggregation with species list and average frequency
    # (Equivalent to M03_core_shared_genes_summary.csv)
    core_poigs_summary <- core_hotspots_with_metadata %>%
      dplyr::group_by(poiGS_ID) %>%
      dplyr::summarise(
        n_species = dplyr::n_distinct(species),
        species_list = paste(sort(unique(species)), collapse = ", "),
        avg_frequency = mean(frequency_per_kb, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(n_species), dplyr::desc(avg_frequency))
    
    core_poigs_summary_file <- file.path(m04_output_dir, "M04_core_shared_poigs_summary.csv")
    write.csv(core_poigs_summary, core_poigs_summary_file, row.names = FALSE)
    log_message(sprintf("Saved core shared poiGS summary: %s", basename(core_poigs_summary_file)))
    log_message(sprintf("  - Summary contains %d unique poiGS", nrow(core_poigs_summary)))
    
    # Summary Table 2: poiGS-level aggregation with detailed frequency statistics
    # (Equivalent to M03_core_shared_hotspots.csv)
    core_poigs_hotspots <- core_hotspots_with_metadata %>%
      dplyr::group_by(poiGS_ID) %>%
      dplyr::summarise(
        num_species = dplyr::n_distinct(species),
        mean_frequency = mean(frequency_per_kb, na.rm = TRUE),
        max_frequency = max(frequency_per_kb, na.rm = TRUE),
        species_list = paste(sort(unique(species)), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(num_species), dplyr::desc(mean_frequency))
    
    core_poigs_hotspots_file <- file.path(m04_output_dir, "M04_core_shared_poigs_hotspots.csv")
    write.csv(core_poigs_hotspots, core_poigs_hotspots_file, row.names = FALSE)
    log_message(sprintf("Saved core shared poiGS hotspots: %s", basename(core_poigs_hotspots_file)))
    log_message(sprintf("  - Hotspots table contains %d unique poiGS", nrow(core_poigs_hotspots)))
    
    # Record new output files
    analysis_results$output_files$core_poigs_summary <- core_poigs_summary_file
    analysis_results$output_files$core_poigs_hotspots <- core_poigs_hotspots_file
    
    # =============================================
    # STEP 4A: GLOBAL poiGS HEATMAP GENERATION (ALL 128 poiGS)
    # =============================================
    
    log_message("M04 Step 4a: Generating global poiGS heatmap (all poiGS units for reference)")
    
    # Generate global heatmap with all poiGS units for user comparison
    # Use consistent naming convention (no timestamp) to align with other modules
    global_heatmap_filename <- "M04_global_poigs_heatmap.pdf"
    global_heatmap_path <- file.path(m04_output_dir, global_heatmap_filename)
    
    log_message(sprintf("Generating global heatmap with %d poiGS regions across %d species",
                       length(unique(poigs_frequencies$poiGS_ID)),
                       length(unique(poigs_frequencies$species))))
    
    global_heatmap_result <- generate_global_frequency_heatmap(
      normalized_data = poigs_frequencies,
      config = config,
      task_params = task_params,
      filename = global_heatmap_path,
      show_rownames = TRUE,
      main = "Global Position-based Orthologous IGS Frequency Distribution",
      width = 14,   # System architecture standard for heatmaps
      height = 10,  # System architecture standard for heatmaps
      feature_col = "poiGS_ID"  # Use poiGS_ID instead of gene for M04
    )
    
    if (!is.null(global_heatmap_result)) {
      log_message(sprintf("Global poiGS heatmap saved: %s", global_heatmap_path))
      analysis_results$output_files$global_heatmap <- global_heatmap_path
    } else {
      log_message("Global poiGS heatmap generation failed", level = "warning")
    }
    
    # =============================================
    # STEP 4B: INFORMATION-DENSE CORE HOTSPOT HEATMAP GENERATION
    # =============================================
    
    log_message("M04 Step 4b: Generating information-dense core poiGS hotspot heatmap")
    
    # Filter poiGS frequency data to include only core hotspots
    core_hotspot_ids <- unique(core_hotspots_with_metadata$poiGS_ID)
    core_hotspot_frequencies <- poigs_frequencies %>%
      dplyr::filter(poiGS_ID %in% core_hotspot_ids)
    
    if (nrow(core_hotspot_frequencies) == 0) {
      warning("M04 Engine: No core hotspot frequency data available for visualization")
      analysis_results$plots_generated <- 0
    } else {
      
      # Generate information-dense heatmap
      # Use consistent naming convention (no timestamp) to align with other modules
      heatmap_filename <- "M04_core_igs_heatmap.pdf"
      heatmap_file_path <- file.path(m04_output_dir, heatmap_filename)
      
      log_message(sprintf("Generating core hotspot heatmap with %d poiGS regions across %d species",
                         length(unique(core_hotspot_frequencies$poiGS_ID)),
                         length(unique(core_hotspot_frequencies$species))))
      
      # Generate heatmap with enhanced settings for information density
      heatmap_result <- generate_global_frequency_heatmap(
        normalized_data = core_hotspot_frequencies,
        config = config,
        task_params = task_params,
        filename = heatmap_file_path,
        show_rownames = TRUE,  # Show poiGS names for biological interpretation
        main = "Core Position-based Orthologous IGS Hotspots",
        width = 14,   # System architecture standard for heatmaps
        height = 10,  # System architecture standard for heatmaps
        feature_col = "poiGS_ID"  # Use poiGS_ID instead of gene for M04
      )
      
      if (!is.null(heatmap_result)) {
        analysis_results$plots_generated <- 2  # Global + Core heatmaps
        log_message(sprintf("Core hotspot heatmap saved: %s", heatmap_file_path))
      } else {
        warning("M04 Engine: Core heatmap generation failed")
        analysis_results$plots_generated <- 1  # Only global heatmap succeeded
      }
    }
    
    # =============================================
    # STEP 4C: UPSET PLOT GENERATION (poiGS SHARING PATTERNS)
    # =============================================
    
    if ("upset" %in% plot_types) {
      log_message("M04 Step 4c: Generating UpSet plot for poiGS hotspot sharing patterns")
      
      upset_plot_result <- tryCatch({
        generate_hotspot_upset_plot(
          candidate_hotspots = candidate_hotspot_data,
          config = config,
          task_params = hotspot_params,
          feature_col = "poiGS_ID"  # Use poiGS_ID for M04
        )
      }, error = function(e) {
        log_message(sprintf("Failed to generate UpSet plot: %s", e$message), level = "warning")
        NULL
      })
      
      if (!is.null(upset_plot_result)) {
        # Save UpSet plot
        upset_plot_file <- file.path(m04_output_dir, "M04_candidate_upset.pdf")
        pdf(upset_plot_file, width = 14, height = 10)
        print(upset_plot_result$plot)
        dev.off()
        
        # Save intersection data
        upset_csv_file <- file.path(m04_output_dir, "M04_upset_intersections.csv")
        write.csv(upset_plot_result$intersection_data, upset_csv_file, row.names = FALSE)
        
        analysis_results$plots_generated <- analysis_results$plots_generated + 1
        analysis_results$output_files$upset_plot <- upset_plot_file
        analysis_results$output_files$upset_intersections <- upset_csv_file
        
        log_message(sprintf("UpSet plot saved: %s", upset_plot_file))
        log_message(sprintf("UpSet intersection data saved: %s", upset_csv_file))
      } else {
        log_message("UpSet plot generation skipped or failed", level = "warning")
      }
    } else {
      log_message("M04 Step 4c: UpSet plot generation disabled by plot_types configuration")
    }
    
    # =============================================
    # STEP 4D: STATISTICAL TESTING (TWO-SPEED EVOLUTION MODEL)
    # =============================================
    
    if ("statistical_test" %in% plot_types) {
      log_message("M04 Step 4d: Performing statistical tests for poiGS sharing patterns")
      
      # CRITICAL FIX: Use ALL candidate hotspots (not just core shared) to align with M03 logic
      # This ensures the statistical test pool includes private, intermediate, and core poiGS
      # for proper Two-Speed Evolution Model testing
      poigs_sharing <- candidate_hotspot_data %>%
        dplyr::group_by(poiGS_ID) %>%
        dplyr::summarise(
          num_species = dplyr::n_distinct(species),
          mean_frequency = mean(frequency_per_kb, na.rm = TRUE),
          max_frequency = max(frequency_per_kb, na.rm = TRUE),
          species_list = paste(unique(species), collapse = ", "),
          .groups = "drop"
        )
      
      hotspot_analysis_for_stats <- list(
        all_gene_sharing = poigs_sharing,
        parameters = list(
          total_species = total_species,
          threshold_mode = threshold_mode,
          min_species_threshold = min_species_number
        )
      )
      
      statistical_test_result <- tryCatch({
        perform_sharing_statistical_test(
          hotspot_analysis_results = hotspot_analysis_for_stats,
          config = config,
          task_params = hotspot_params,
          feature_col = "poiGS_ID"  # Use poiGS_ID for M04
        )
      }, error = function(e) {
        log_message(sprintf("Failed to perform statistical tests: %s", e$message), level = "warning")
        NULL
      })
      
      if (!is.null(statistical_test_result)) {
        # Save statistical test report
        stats_report_file <- file.path(m04_output_dir, "M04_core_shared_statistical_tests.txt")
        writeLines(statistical_test_result$report, stats_report_file)
        
        analysis_results$output_files$statistical_tests <- stats_report_file
        
        log_message(sprintf("Statistical test report saved: %s", stats_report_file))
        log_message(sprintf("Two-Speed Evolution Model: %s", 
                           ifelse(statistical_test_result$tests$two_speed_model$model_supported, 
                                 "SUPPORTED", "NOT SUPPORTED")))
        
        # =============================================
        # SAVE DETAILED DATA TABLES (ALIGN WITH M03)
        # =============================================
        
        # 1. Save Chi-square test details (observed vs expected distribution)
        if (!is.null(statistical_test_result$comparison_table)) {
          chisq_csv_filename <- "M04_core_shared_chisquare_details.csv"
          chisq_csv_path <- file.path(m04_output_dir, chisq_csv_filename)
          write.csv(statistical_test_result$comparison_table, chisq_csv_path, row.names = FALSE)
          log_message(sprintf("M04: Saved chi-square test details: %s", chisq_csv_path))
        }
        
        # 2. Save Hypergeometric test details
        if (!is.null(statistical_test_result$tests$hypergeometric)) {
          
          # 2a. Save complete poiGS-level hypergeometric test results (ALIGN WITH M03)
          if (!is.null(statistical_test_result$tests$hypergeometric$all_gene_results) && 
              nrow(statistical_test_result$tests$hypergeometric$all_gene_results) > 0) {
            
            poigs_level_data <- statistical_test_result$tests$hypergeometric$all_gene_results
            
            log_message(sprintf("M04: poiGS-level hypergeometric data ready (%d poiGS, %d cols)", 
                               nrow(poigs_level_data), ncol(poigs_level_data)))
            
            # Save complete poiGS-level results
            poigs_level_csv_filename <- "M04_core_shared_hypergeometric_tests.csv"
            poigs_level_csv_path <- file.path(m04_output_dir, poigs_level_csv_filename)
            
            tryCatch({
              write.csv(poigs_level_data, poigs_level_csv_path, row.names = FALSE)
              log_message(sprintf("M04: [SUCCESS] Successfully saved poiGS-level hypergeometric tests: %s", poigs_level_csv_path))
            }, error = function(e) {
              log_message(sprintf("M04: [ERROR] Failed to save poiGS-level hypergeometric tests: %s", e$message), level = "error")
            })
            
            # 2b. Save driver poiGS summary (significant contributors)
            driver_poigs_summary <- data.frame(
              category = c(rep("Private", nrow(statistical_test_result$tests$hypergeometric$private_drivers %||% data.frame())),
                          rep("Core", nrow(statistical_test_result$tests$hypergeometric$core_drivers %||% data.frame()))),
              rbind(
                statistical_test_result$tests$hypergeometric$private_drivers %||% data.frame()[0,],
                statistical_test_result$tests$hypergeometric$core_drivers %||% data.frame()[0,]
              ),
              stringsAsFactors = FALSE
            )
            
            if (nrow(driver_poigs_summary) > 0) {
              driver_csv_filename <- "M04_core_shared_driver_poigs_summary.csv"
              driver_csv_path <- file.path(m04_output_dir, driver_csv_filename)
              
              tryCatch({
                write.csv(driver_poigs_summary, driver_csv_path, row.names = FALSE)
                log_message(sprintf("M04: [SUCCESS] Successfully saved driver poiGS summary: %s", driver_csv_path))
              }, error = function(e) {
                log_message(sprintf("M04: [ERROR] Failed to save driver poiGS summary: %s", e$message), level = "error")
              })
            }
          }
          
          # 2c. Save hypergeometric summary (legacy format for backward compatibility)
          hypergeom_data <- data.frame(
            test_method = "poiGS-by-poiGS Hypergeometric Analysis",
            total_poigs_tested = statistical_test_result$tests$hypergeometric$total_genes_tested %||% 0,
            private_poigs_tested = statistical_test_result$tests$hypergeometric$private_genes_tested %||% 0,
            core_poigs_tested = statistical_test_result$tests$hypergeometric$core_genes_tested %||% 0,
            significant_private_drivers = statistical_test_result$tests$hypergeometric$significant_private_drivers %||% 0,
            significant_core_drivers = statistical_test_result$tests$hypergeometric$significant_core_drivers %||% 0,
            alpha_level = statistical_test_result$tests$hypergeometric$alpha_level %||% 0.05,
            method_description = statistical_test_result$tests$hypergeometric$method_description %||% "poiGS-by-poiGS analysis",
            stringsAsFactors = FALSE
          )
          
          hypergeom_csv_filename <- "M04_core_shared_hypergeometric_details.csv"
          hypergeom_csv_path <- file.path(m04_output_dir, hypergeom_csv_filename)
          write.csv(hypergeom_data, hypergeom_csv_path, row.names = FALSE)
          log_message(sprintf("M04: Saved hypergeometric test summary: %s", hypergeom_csv_path))
        }
        
        # 3. Save sharing statistics summary
        if (!is.null(statistical_test_result$sharing_summary)) {
          summary_data <- data.frame(
            statistic = c("total_poigs", "total_species", "mean_sharing", "median_sharing", 
                         "max_sharing", "poigs_in_multiple_species", "poigs_in_all_species"),
            value = c(
              statistical_test_result$sharing_summary$total_genes,
              statistical_test_result$sharing_summary$total_species,
              statistical_test_result$sharing_summary$mean_sharing,
              statistical_test_result$sharing_summary$median_sharing,
              statistical_test_result$sharing_summary$max_sharing,
              statistical_test_result$sharing_summary$genes_in_multiple_species,
              statistical_test_result$sharing_summary$genes_in_all_species
            ),
            stringsAsFactors = FALSE
          )
          
          summary_csv_filename <- "M04_core_shared_sharing_statistics.csv"
          summary_csv_path <- file.path(m04_output_dir, summary_csv_filename)
          write.csv(summary_data, summary_csv_path, row.names = FALSE)
          log_message(sprintf("M04: Saved sharing statistics summary: %s", summary_csv_path))
        }
        
      } else {
        log_message("Statistical testing skipped or failed", level = "warning")
      }
    } else {
      log_message("M04 Step 4d: Statistical testing disabled by plot_types configuration")
    }
    
    # =============================================
    # STEP 5: FINALIZE RESULTS AND QUALITY REPORTING
    # =============================================
    
    log_message("M04 Step 5: Finalizing analysis results and generating quality report")
    
    # Generate comprehensive quality report
    quality_report <- data.frame(
      analysis_type = "M04_poiGS_Hotspot_Engine",
      analysis_timestamp = Sys.time(),
      
      # Input data metrics
      input_frequency_observations = nrow(scoped_data),
      input_species_count = length(unique(scoped_data$species)),
      input_igs_regions = length(unique(scoped_data$gene)),
      
      # poiGS mapping metrics
      total_igs_regions_mapped = nrow(poigs_mapping),
      unique_poigs_identified = length(unique(poigs_mapping$poiGS_ID)),
      poigs_mapping_efficiency = round(length(unique(poigs_mapping$poiGS_ID)) / nrow(poigs_mapping) * 100, 2),
      
      # Frequency aggregation metrics
      poigs_frequency_observations = nrow(poigs_frequencies),
      mean_poigs_frequency = round(mean(poigs_frequencies$frequency_per_kb, na.rm = TRUE), 4),
      median_poigs_frequency = round(median(poigs_frequencies$frequency_per_kb, na.rm = TRUE), 4),
      var_type_aggregation_method = "sum_all_var_types",
      
      # Hotspot filtering metrics
      total_poigs_with_thresholds = nrow(candidate_hotspots),
      candidate_hotspots_count = length(unique(candidate_hotspot_data$poiGS_ID)),
      candidate_hotspot_rate_pct = round(length(unique(candidate_hotspot_data$poiGS_ID)) / 
                                         length(unique(candidate_hotspots$poiGS_ID)) * 100, 2),
      core_shared_hotspots_count = analysis_results$core_hotspots_identified,
      core_to_candidate_ratio_pct = round(analysis_results$core_hotspots_identified / 
                                          length(unique(candidate_hotspot_data$poiGS_ID)) * 100, 2),
      
      # Analysis parameters (aligned with M03)
      threshold_mode = threshold_mode,
      min_species_coverage_threshold = min_species_coverage,
      min_species_number = min_species_number,
      frequency_percentile_threshold = get_task_parameter(hotspot_params, config, "frequency_percentile", 0.75),
      
      # Output metrics
      plots_generated = analysis_results$plots_generated,
      
      stringsAsFactors = FALSE
    )
    
    quality_report_file <- file.path(m04_output_dir, "M04_poigs_quality_report.csv")
    write.csv(quality_report, quality_report_file, row.names = FALSE)
    log_message(sprintf("Saved comprehensive quality report: %s", quality_report_file))
    
    # Record all output files (including new transparency-enhancing files)
    analysis_results$output_files <- list(
      poigs_mapping = poigs_mapping_file,
      poigs_frequencies = poigs_freq_file,
      all_poigs_with_thresholds = all_poigs_annotated_file,
      candidate_hotspots = candidate_hotspots_file,
      core_shared_hotspots_final = core_hotspots_file,
      quality_report = quality_report_file,
      output_directory = m04_output_dir
    )
    
    if (analysis_results$plots_generated > 0) {
      analysis_results$output_files$global_heatmap <- global_heatmap_path
      analysis_results$output_files$core_heatmap <- heatmap_file_path
      analysis_results$output_files$heatmap_result <- heatmap_result
    }
    
    # Mark analysis as successful and record output directory
    analysis_results$success <- TRUE
    analysis_results$output_dir <- m04_output_dir  # CRITICAL: Required for M05 integration to locate M04 results
    
    log_message("=== M04 POSITION-BASED ORTHOLOGOUS IGS ANALYSIS COMPLETE ===")
    log_message(sprintf("Analysis Results Summary:"))
    log_message(sprintf("  - Position-based orthologous IGS regions: %d", analysis_results$poigs_regions_analyzed))
    log_message(sprintf("  - Core shared hotspots identified: %d", analysis_results$core_hotspots_identified))
    log_message(sprintf("  - Species analyzed: %d", analysis_results$species_count))
    log_message(sprintf("  - Information-dense visualizations: %d", analysis_results$plots_generated))
    log_message(sprintf("  - Output directory: %s", m04_output_dir))
    
    if (analysis_results$core_hotspots_identified > 0) {
      log_message("SUCCESS: M04 engine identified biologically significant cross-species IGS hotspots")
    } else {
      log_message("WARNING: No core shared IGS hotspots identified - consider adjusting thresholds", level = "warning")
    }
    
  }, error = function(e) {
    
    # Comprehensive error logging and reporting
    analysis_results$success <- FALSE
    analysis_results$error_message <- e$message
    log_message(sprintf("M04 poiGS Hotspot Engine FAILED: %s", e$message), level = "error")
    log_message("Check input data quality, file paths, and configuration parameters", level = "error")
    
  })
  
  return(analysis_results)
}