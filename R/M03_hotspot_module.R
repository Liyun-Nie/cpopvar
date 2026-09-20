############################################################
#### M03_hotspot_module.R - M03 Module Orchestrator ####
############################################################
#
# Thin orchestrator for M03 Hotspot Analysis
# This is a Phase 2 "thin wrapper" that loads dependencies and delegates
# to the existing generate_M03_hotspot_plots function
#
# Phase 2 Rule #3 Compliance: Does NOT modify existing function libraries,
# only provides a clean interface for task_dispatcher.R
#
############################################################

#' M03 Hotspot Analysis Module Orchestrator
#' 
#' @param normalized_data Normalized frequency data
#' @param config Configuration object
#' @param output_dir Output directory for plots
#' @param task_name Task name for output organization
#' @param task_params Task-specific parameters
#' @return M03 analysis results
#' 
#' ARCHITECTURE: Implements M01/M02-style standardized data cleaning and output generation
#' following the established golden standard for module orchestrators.
#' run_m03_hotspot
#' @export
run_m03_hotspot <- function(normalized_data,
                           config,
                           output_dir,
                           task_name,
                           task_params) {
  
  # Function dependencies are automatically loaded in R package context
  
  log_message(sprintf("=== M03 HOTSPOT ANALYSIS: %s ===", task_name))
  
  # Validate input data
  if (nrow(normalized_data) == 0) {
    log_message("ERROR: normalized_data is empty", level = "error")
    stop("M03 analysis cannot proceed with empty normalized_data")
  }
  
  # Apply data_scoping filters (PROJECT_CHARTER.md Principle 3)
  log_message("=== APPLYING DATA_SCOPING FILTERS ===")
  scoped_data <- apply_data_scoping(normalized_data, task_params)
  
  if (nrow(scoped_data) == 0) {
    log_message("ERROR: data_scoping resulted in empty dataset for M03", level = "error")
    stop("M03 analysis cannot proceed with empty data after data_scoping")
  }
  
  log_message(sprintf("M03 data scoping: %d -> %d rows (%.2f%% retained)", 
                      nrow(normalized_data), nrow(scoped_data), 
                      ifelse(nrow(normalized_data) > 0, round(nrow(scoped_data)/nrow(normalized_data)*100, 2), 0)))
  
  # Aggregate region-level data to gene-level for hotspot analysis
  # Architecture principle: M03 hotspot detection works at gene level
  # With P03 optimization, we now have both region_name and gene columns
  # Simply aggregate by gene column (no regex parsing needed!)
  log_message("=== M03: AGGREGATING TO GENE-LEVEL FOR HOTSPOT ANALYSIS ===")
  
  nrow_before_agg <- nrow(scoped_data)
  
  # Check if data has region_name column (new P03 format)
  if ("region_name" %in% names(scoped_data)) {
    log_message("P03 optimized format detected: using 'gene' column directly")
    
    # Aggregate by gene for CDS/intron regions (gene column is already at gene-level)
    # For IGS, gene is NA, so they will be grouped separately
    scoped_data <- scoped_data %>%
      dplyr::filter(!is.na(gene)) %>%  # Remove IGS (gene=NA) as they don't have gene-level aggregation
      dplyr::group_by(species, var_type, region_type, genome_region, gene) %>%
      dplyr::summarise(
        variant_count = sum(variant_count, na.rm = TRUE),
        region_length = sum(region_length, na.rm = TRUE),  # Sum all region lengths
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        # Recalculate frequency at gene level
        frequency_per_kb = (variant_count / region_length) * 1000
      )
    
    log_message(sprintf("Gene-level aggregation: %d region-level records -> %d gene-level records (%.1f%% reduction)",
                       nrow_before_agg, nrow(scoped_data), 
                       (1 - nrow(scoped_data)/nrow_before_agg) * 100))
    log_message(sprintf("  - Unique genes after aggregation: %d", length(unique(scoped_data$gene))))
  } else {
    log_message("Legacy format detected: gene column already at appropriate level")
  }
  
  # Filter out NA frequencies but keep zero frequencies for hotspot analysis
  # Zero frequencies are important background for hotspot identification
  cleaned_data <- scoped_data %>% 
    dplyr::filter(!is.na(frequency_per_kb)) %>%
    dplyr::filter(!is.na(gene) & !is.na(region_type))  # Ensure gene and region information is complete
  
  log_message(sprintf("Data prepared for hotspot analysis: %d -> %d valid observations", 
                     nrow(normalized_data), nrow(cleaned_data)))
  log_message(sprintf("Hotspot data includes %d unique genes across %d species",
                     length(unique(cleaned_data$gene)), 
                     length(unique(cleaned_data$species))))
  
  # Apply species ordering
  species_ordering_result <- load_and_apply_species_ordering(
    data = cleaned_data,
    session_id = config$session_info$session_id,
    species_col = "species"
  )
  
  cleaned_data <- species_ordering_result$data
  log_message(species_ordering_result$message)
  
  # Create output directory
  data_output_dir <- file.path(output_dir, "M03_hotspot", task_name)
  dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Validate cleaned data for hotspot analysis requirements
  if (nrow(cleaned_data) == 0) {
    log_message("ERROR: No valid data after cleaning", level = "error")
    stop("Cannot proceed with M03 analysis - no valid hotspot data")
  }
  
  # Save cleaned hotspot data for debugging and validation
  cleaned_data_file <- file.path(data_output_dir, "M03_hotspot_data_cleaned.csv")
  write.csv(cleaned_data, cleaned_data_file, row.names = FALSE)
  log_message(sprintf("Saved cleaned hotspot data: %s", cleaned_data_file))
  
  # ========================================================================
  # VARIANT TYPE FREQUENCY AGGREGATION (aligned with M05 architecture)
  # ========================================================================
  # Check if multiple variant types are present - if so, aggregate them
  unique_var_types <- unique(cleaned_data$var_type)
  
  if (length(unique_var_types) > 1) {
    log_message(sprintf("M03: Multiple variant types detected (%s) - performing frequency aggregation", 
                       paste(unique_var_types, collapse = ", ")))
    
    # Two-stage aggregation (same logic as M05)
    # Stage 1: Handle potential duplicates for same gene-species-var_type combination
    aggregated_stage1 <- cleaned_data %>%
      dplyr::group_by(species, gene, var_type, region_type) %>%
      dplyr::summarise(
        frequency_per_kb = mean(frequency_per_kb, na.rm = TRUE),
        variant_count = sum(variant_count, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Stage 2: Merge all variant types for each gene
    aggregated_data <- aggregated_stage1 %>%
      dplyr::group_by(species, gene, region_type) %>%
      dplyr::summarise(
        frequency_per_kb = sum(frequency_per_kb, na.rm = TRUE),  # SUM all var_type frequencies
        variant_count = sum(variant_count, na.rm = TRUE),
        var_type = paste(sort(unique(var_type)), collapse = "+"),  # Record which types present
        .groups = "drop"
      )
    
    log_message(sprintf("M03: Frequency aggregation complete - %d observations (from %d)", 
                       nrow(aggregated_data), nrow(cleaned_data)))
    log_message(sprintf("M03: Aggregated frequency range: %.3f to %.3f per kb", 
                       min(aggregated_data$frequency_per_kb), 
                       max(aggregated_data$frequency_per_kb)))
    
    # Use aggregated data for hotspot analysis
    analysis_data <- aggregated_data
  } else {
    log_message(sprintf("M03: Single variant type detected (%s) - no aggregation needed", 
                       unique_var_types[1]))
    analysis_data <- cleaned_data
  }
  
  # Generate hotspot analysis quality report
  quality_report <- data.frame(
    total_observations = nrow(analysis_data),
    unique_species = length(unique(analysis_data$species)),
    unique_genes = length(unique(analysis_data$gene)),
    region_types = paste(unique(analysis_data$region_type), collapse = ", "),
    var_types = paste(unique(analysis_data$var_type), collapse = ", "),
    var_type_aggregation = ifelse(length(unique_var_types) > 1, "aggregated", "single_type"),
    mean_frequency = mean(analysis_data$frequency_per_kb, na.rm = TRUE),
    median_frequency = median(analysis_data$frequency_per_kb, na.rm = TRUE),
    max_frequency = max(analysis_data$frequency_per_kb, na.rm = TRUE),
    species_ordering_applied = !is.null(species_ordering_result$species_order),
    timestamp = Sys.time()
  )
  
  quality_report_file <- file.path(data_output_dir, "M03_hotspot_quality_report.csv")
  write.csv(quality_report, quality_report_file, row.names = FALSE)
  log_message(sprintf("Generated hotspot quality report: %s", quality_report_file))
  
  # Call the main M03 hotspot analysis function
  log_message("Delegating to M03 hotspot analysis engine...")
  
  tryCatch({
    results <- generate_M03_hotspot_plots(
      normalized_data = analysis_data,  # Use aggregated data if multiple var_types present
      config = config,
      output_dir = data_output_dir,
      task_name = task_name,
      task_params = task_params
    )
    
    if (is.null(results) || !results$success) {
      log_message("M03 hotspot analysis failed", level = "error")
      return(list(success = FALSE, error = "Hotspot analysis failed"))
    }
    
    # Save additional CSV outputs if available
    if (!is.null(results$plots$upset) && !is.null(results$plots$upset$intersection_data)) {
      upset_csv_path <- file.path(data_output_dir, sprintf("%s_upset_intersections.csv", task_name))
      write.csv(results$plots$upset$intersection_data, upset_csv_path, row.names = FALSE)
      log_message(sprintf("Saved UpSet intersection data: %s", upset_csv_path))
    }
    
    # ================================================================
    # Phase 1 Pilot: Generate Manifest for Scientific Report System  
    # ================================================================
    # ===============================================================
    # V4 ARCHITECTURE FIX: Robust Manifest Generation
    # ===============================================================
    log_message("Generating output manifest for scientific report system (V4 Robust Version)")

    manifest_entries <- list()

    # Process all single-file-output results first
    for (result_name in names(results$plots)) {
      # Skip the multi-file pca result, it will be handled separately
      if (result_name == 'global_pca') next

      plot_info <- results$plots[[result_name]]

      if (!is.null(plot_info$file_path) && file.exists(plot_info$file_path)) {
        # Determine role and priority
        role <- "individual_plot" # Default
        if (result_name %in% c("global_heatmap")) role <- "composite_plot"
        if (result_name %in% c("upset", "hotspot_correlation")) role <- "key_plot"

        priority <- switch(role, "composite_plot" = 1, "key_plot" = 2, 3)

        entry <- create_manifest_entry(
          path = plot_info$file_path,
          role = role,
          title = get_standard_title("M03", role, result_name),
          caption = get_standard_caption("M03", role, result_name),
          module = "M03",
          file_type = "plot",
          priority = priority,
          parameters = task_params  # V7.1 Enhancement: Pass task parameters
        )
        manifest_entries[[result_name]] <- entry
        log_message(sprintf("[SUCCESS] Created manifest entry for %s (%s)", result_name, role))
      }
    }

    # Now, specifically handle the multi-file PCA results
    if (!is.null(results$plots$global_pca) && !is.null(results$plots$global_pca$file_paths)) {
      log_message(sprintf("Processing %d files from global_pca result...", length(results$plots$global_pca$file_paths)))

      for (pca_file_path in results$plots$global_pca$file_paths) {
        if (file.exists(pca_file_path)) {
          file_basename <- basename(pca_file_path)
          # Infer role and title from filename
          role <- "individual_plot"
          if (grepl("composite", file_basename)) role <- "composite_plot"
          if (grepl("scree|loadings", file_basename)) role <- "key_plot"

          priority <- switch(role, "composite_plot" = 1, "key_plot" = 2, 3)

          entry_name <- tools::file_path_sans_ext(file_basename)

          entry <- create_manifest_entry(
            path = pca_file_path,
            role = role,
            title = get_standard_title("M03", role, entry_name),
            caption = get_standard_caption("M03", role, entry_name),
            module = "M03",
            file_type = "plot",
            priority = priority,
            parameters = task_params  # V7.1 Enhancement: Pass task parameters
          )
          manifest_entries[[entry_name]] <- entry
          log_message(sprintf("[SUCCESS] Created manifest entry for PCA file %s (%s)", file_basename, role))
        }
      }
    }
    
    # Process M03 data files (CSV outputs)
    if (!is.null(results$data_files)) {
      for (data_name in names(results$data_files)) {
        data_path <- results$data_files[[data_name]]
        
        if (file.exists(data_path)) {
          
          # Classify data files
          role <- if (grepl("summary|core_shared", data_name)) {
            "summary_table"
          } else if (grepl("statistical|test", data_name)) {
            "statistical_report"  
          } else {
            "raw_data"
          }
          
          title <- get_standard_title("M03", role, data_name)
          caption <- get_standard_caption("M03", role)
          
          entry <- create_manifest_entry(
            path = data_path,
            role = role,
            title = title,
            caption = caption,
            module = "M03",
            file_type = "table",
            priority = if (role == "summary_table") 4 else 6,
            parameters = task_params  # V7.1 Enhancement: Pass task parameters
          )
          
          manifest_entries[[data_name]] <- entry
          log_message(sprintf("[SUCCESS] Created manifest entry for %s data file", data_name))
        }
      }
    }
    
    # Check for UpSet intersection data that was saved
    upset_csv_path <- file.path(data_output_dir, sprintf("%s_upset_intersections.csv", task_name))
    if (file.exists(upset_csv_path)) {
      title <- get_standard_title("M03", "summary_table", "upset_intersections")
      caption <- get_standard_caption("M03", "summary_table")
      
      entry <- create_manifest_entry(
        path = upset_csv_path,
        role = "summary_table",
        title = title,
        caption = caption,
        module = "M03",
        file_type = "table",
        priority = 4,
        parameters = task_params  # V7.1 Enhancement: Pass task parameters
      )
      
      manifest_entries[["upset_intersections"]] <- entry
      log_message("[SUCCESS] Created manifest entry for UpSet intersection data")
    }
    
    # Create manifest collection
    session_id <- config$session_info$session_id %||% "unknown"
    manifest <- create_manifest_collection(manifest_entries, "M03", session_id)
    
    log_message(sprintf("=== M03 HOTSPOT ANALYSIS COMPLETED: %s ===", task_name))
    # Report plot generation results with proper validation
    plot_count <- if (!is.null(results$plots_generated)) {
      length(results$plots_generated)
    } else if (!is.null(results$plots)) {
      length(names(results$plots))
    } else {
      0
    }
    log_message(sprintf("Generated %d plot types", plot_count))
    log_message(sprintf("Generated manifest with %d entries", length(manifest_entries)))
    
    return(list(
      success = TRUE,
      plots = results$plots,
      data = cleaned_data,
      quality_report = quality_report,
      plots_generated = results$plots_generated,
      task_name = task_name,
      output_dir = data_output_dir,
      manifest = manifest
    ))
    
  }, error = function(e) {
    log_message(sprintf("M03 hotspot analysis error: %s", e$message), level = "error")
    return(list(success = FALSE, error = e$message))
  })
}