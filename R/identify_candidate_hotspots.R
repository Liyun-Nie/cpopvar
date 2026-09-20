#' Identify candidate hotspot regions
#' 
#' @title Identify candidate hotspot regions
#' @description Identifies genomic regions with high variant frequencies as candidate hotspots
#' @param normalized_data Normalized frequency dataset
#' @param config Configuration object
#' @param task_params Task-specific parameters
#' @param feature_col Name of the feature column to rank
#' @return List containing identified hotspot candidates and analysis data
#' @export
identify_candidate_hotspots <- function(normalized_data, config, task_params, 
                                        feature_col = "gene") {
  
  # Get parameters from config, with a final fallback default.
  frequency_percentile <- get_task_parameter(task_params, config, "frequency_percentile", 0.75, log_default = TRUE)

  tryCatch({
    log_message(sprintf("M03 Step 1: Ranking all %s features and calculating hotspot thresholds (percentile: %.2f)", 
                       feature_col, frequency_percentile))
    
    # Debug input data
    log_message(sprintf("Input data for feature ranking: %d rows, %d columns", 
                       nrow(normalized_data), ncol(normalized_data)))
    log_message(sprintf("Column names: %s", paste(names(normalized_data), collapse = ", ")))
    
    # Use internal standard column names (boundary principle)
    species_col <- "species"
    gene_col <- feature_col  # Use parameterized feature column
    freq_col <- "frequency_per_kb"
    region_col <- "region_type"
    var_type_col <- "var_type"
    
    # Check if required columns exist
    # var_type is now OPTIONAL for M05 compatibility
    required_cols <- c(species_col, gene_col, freq_col, region_col)
    optional_cols <- c(var_type_col)
    
    missing_required <- setdiff(required_cols, names(normalized_data))
    if (length(missing_required) > 0) {
      log_message(sprintf("CRITICAL: Missing required columns for feature ranking: %s", 
                         paste(missing_required, collapse = ", ")), level = "error")
      return(NULL)
    }
    
    # Check optional columns
    has_var_type <- var_type_col %in% names(normalized_data)
    
    # Trust the input data - it's already been scoped and filtered by the orchestrator (Rule #11)
    # The orchestrator ensures NA values and zero frequencies are already removed
    filtered_data <- normalized_data
    
    log_message(sprintf("Hotspot threshold calculation data:"))
    log_message(sprintf("  - Input data (pre-filtered by orchestrator): %d observations", nrow(filtered_data)))
    log_message(sprintf("  - Unique genes for threshold calculation: %d", length(unique(filtered_data[[gene_col]]))))
    
    log_message(sprintf("Data validation summary:"))
    log_message(sprintf("  - Region types: %s", paste(unique(normalized_data[[region_col]]), collapse = ", ")))
    if (has_var_type) {
      log_message(sprintf("  - Variant types: %s", paste(unique(normalized_data[[var_type_col]]), collapse = ", ")))
    } else {
      log_message(sprintf("  - Variant types: [aggregated/not present]"))
    }
    log_message(sprintf("  - Frequency range: %.3f to %.3f per kb", 
                       min(normalized_data[[freq_col]], na.rm = TRUE),
                       max(normalized_data[[freq_col]], na.rm = TRUE)))
    
    if (nrow(filtered_data) == 0) {
      log_message("No data found after filtering for gene ranking", level = "warning")
      log_message("Checking available values in key columns:", level = "warning")
      log_message(sprintf("  - region_type values: %s", 
                         paste(unique(normalized_data[[region_col]]), collapse = ", ")), level = "warning")
      if (has_var_type) {
        log_message(sprintf("  - var_type values: %s", 
                           paste(unique(normalized_data[[var_type_col]]), collapse = ", ")), level = "warning")
      }
      log_message(sprintf("  - frequency_per_kb range: %.3f to %.3f", 
                         min(normalized_data[[freq_col]], na.rm = TRUE),
                         max(normalized_data[[freq_col]], na.rm = TRUE)), level = "warning")
      return(NULL)
    }
    
    log_message(sprintf("Filtered data for ranking: %d rows across %d species", 
                       nrow(filtered_data), length(unique(filtered_data[[species_col]]))))
    
    # --- Correct Threshold Calculation Logic ---

    # 1. Create a temporary dataset that only includes non-zero frequencies, specifically for calculating thresholds
    nonzero_data <- filtered_data %>%dplyr::filter(.data[[freq_col]] > 0)

    # 2. Calculate threshold for each species based only on non-zero data
    species_thresholds <- nonzero_data %>%
      dplyr::group_by(.data[[species_col]]) %>%
      dplyr::summarise(
        threshold = quantile(.data[[freq_col]], frequency_percentile, na.rm = TRUE),
        .groups = "drop"
      )

    log_message(sprintf("Calculated hotspot thresholds for %d species using non-zero frequency data.",
                       nrow(species_thresholds)))

    # 3. Apply calculated thresholds to the complete dataset including zero values
    all_genes_ranked <- filtered_data %>%
      dplyr::left_join(species_thresholds, by = species_col) %>%
      # If a species has no non-zero genes, its threshold will be NA, replace with Inf to ensure they are not selected as hotspots
      dplyr::mutate(threshold = ifelse(is.na(threshold), Inf, threshold)) %>%
      dplyr::mutate(
        is_potential_hotspot = (.data[[freq_col]] >= threshold)
      ) %>%
      dplyr::arrange(.data[[species_col]],dplyr::desc(.data[[freq_col]]))

    # 
      
    # Add validation and results logging
    log_message(sprintf("Applied thresholds to complete dataset: %d total entries", nrow(all_genes_ranked)))
    
    # Return the ranked genes with threshold information
    # Rule #11 compliance: Trust orchestrator-provided data completely  
    # The orchestrator has already filtered the data appropriately for hotspot analysis
    all_genes_with_threshold <- all_genes_ranked
      
    log_message(sprintf("Completed gene ranking and threshold annotation for %d gene entries", 
                       nrow(all_genes_with_threshold)))
    
    # Final validation summary
    total_genes_in_result <- nrow(all_genes_with_threshold)
    hotspot_genes_count <- sum(all_genes_with_threshold$is_potential_hotspot, na.rm = TRUE)
    unique_hotspot_genes <- length(unique(all_genes_with_threshold[[gene_col]][all_genes_with_threshold$is_potential_hotspot]))
    
    log_message(sprintf("THRESHOLD CALCULATION - Final Summary:"))
    log_message(sprintf("  - Total analyzed genes: %d (orchestrator-filtered data)", total_genes_in_result))
    log_message(sprintf("  - Genes meeting hotspot criteria: %d entries from %d unique genes", 
                       hotspot_genes_count, unique_hotspot_genes))
    log_message(sprintf("  - Hotspot rate: %.1f%%", 
                       100 * hotspot_genes_count / total_genes_in_result))
    
    # Validate the results
    if (unique_hotspot_genes > 0 && unique_hotspot_genes < total_genes_in_result) {
      log_message("[SUCCESS] THRESHOLD VALIDATION: Success - reasonable number of hotspot genes identified")
    } else if (unique_hotspot_genes == total_genes_in_result) {
      log_message("[WARNING] THRESHOLD VALIDATION: Warning - all genes identified as hotspots", level = "warning")
    } else {
      log_message("[ERROR] THRESHOLD VALIDATION: Error - no hotspot genes identified", level = "error")
    }
    
    # IMPORTANT: This function returns ranked genes with threshold information
    # Rule #11 compliance: Works with orchestrator-filtered data only
    # No additional filtering - trusts the data provided by the orchestrator
    return(all_genes_with_threshold)
    
  }, error = function(e) {
    log_message(sprintf("Failed to rank genes and calculate thresholds: %s", e$message), level = "error")
    return(NULL)
  })
}
