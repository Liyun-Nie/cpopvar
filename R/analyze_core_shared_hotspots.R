#' Analyze core shared hotspots
#'
#' @param candidate_hotspots Data frame of candidate hotspots
#' @param config Configuration object
#' @param task_params Task parameters
#' @return List containing core shared genes analysis results
#' @importFrom magrittr %>%
#' @export
analyze_core_shared_hotspots <- function(candidate_hotspots, 
                                       config,
                                       task_params = NULL) {
  
  # Get parameters from config, with a final fallback default.
  threshold_mode <- get_task_parameter(task_params, config, "threshold_mode", "ratio", log_default = TRUE)
  min_species_ratio <- get_task_parameter(task_params, config, "min_species_ratio", 0.5, log_default = TRUE)
  min_species_number <- get_task_parameter(task_params, config, "min_species_number", 3, log_default = TRUE)
  
  tryCatch({
    # Determine threshold mode and value
    if (threshold_mode == "absolute" || !is.null(min_species_number)) {
      threshold_mode <- "absolute"
      threshold_value <- min_species_number
    } else if (threshold_mode == "ratio" || !is.null(min_species_ratio)) {
      threshold_mode <- "ratio"
      threshold_value <- min_species_ratio
    } else {
      # Default to ratio mode
      threshold_mode <- "ratio"
      threshold_value <- 0.5
    }
    
    log_message(sprintf("M03 Step 3: Analyzing core shared hotspots (%s mode: %s)", 
                       threshold_mode, 
                       ifelse(threshold_mode == "ratio", sprintf("%.2f", threshold_value), as.character(threshold_value))))
    
    # Use internal standard column names (boundary principle)
    gene_col <- "gene"
    species_col <- "species"
    freq_col <- "frequency_per_kb"
    
    # Calculate how many species each gene appears in as a hotspot
    gene_sharing <- candidate_hotspots %>%
      dplyr::group_by(.data[[gene_col]]) %>%
      dplyr::summarise(
        num_species = dplyr::n(),
        mean_frequency = mean(.data[[freq_col]], na.rm = TRUE),
        max_frequency = max(.data[[freq_col]], na.rm = TRUE),
        species_list = paste(.data[[species_col]], collapse = ", "),
        .groups = "drop"
      )
    
    total_species <- length(unique(candidate_hotspots[[species_col]]))
    
    # Calculate threshold based on mode
    if (threshold_mode == "ratio") {
      min_species_threshold <- ceiling(total_species * threshold_value)
    } else {
      min_species_threshold <- threshold_value
    }
    
    log_message(sprintf("Total species: %d, minimum species threshold: %d", total_species, min_species_threshold))
    
    # Identify core shared hotspot genes (genes in >= threshold species)
    core_shared_genes <- gene_sharing %>%
      dplyr::filter(num_species >= min_species_threshold) %>%
      dplyr::arrange(dplyr::desc(num_species), dplyr::desc(mean_frequency))
    
    # Identify private genes (genes in exactly 1 species)
    private_genes <- gene_sharing %>%
      dplyr::filter(num_species == 1) %>%
      dplyr::arrange(dplyr::desc(mean_frequency))
    
    # Identify intermediate genes (genes in 2 to threshold-1 species)
    intermediate_genes <- gene_sharing %>%
      dplyr::filter(num_species > 1, num_species < min_species_threshold) %>%
      dplyr::arrange(dplyr::desc(num_species), dplyr::desc(mean_frequency))
    
    # Log detailed statistics
    log_message(sprintf("Gene sharing analysis results:"))
    log_message(sprintf("  - Core shared genes (>=%d species): %d", min_species_threshold, nrow(core_shared_genes)))
    log_message(sprintf("  - Private genes (1 species): %d", nrow(private_genes)))
    log_message(sprintf("  - Intermediate genes (2-%d species): %d", min_species_threshold-1, nrow(intermediate_genes)))
    log_message(sprintf("  - Total candidate hotspot genes: %d", nrow(gene_sharing)))
    
    if (nrow(core_shared_genes) == 0) {
      log_message(sprintf("No core shared hotspot genes found with threshold=%d species", min_species_threshold), level = "warning")
      # Show distribution for debugging
      sharing_dist <- gene_sharing %>% 
        dplyr::count(num_species) %>% 
        dplyr::arrange(dplyr::desc(num_species))
      log_message(sprintf("Gene sharing distribution: %s", 
                         paste(sprintf("%d genes in %d species", sharing_dist$n, sharing_dist$num_species), 
                               collapse = ", ")), level = "info")
    } else {
      log_message(sprintf("Top core shared genes: %s", 
                         paste(core_shared_genes[[gene_col]][1:min(5, nrow(core_shared_genes))], collapse = ", ")))
    }
    
    if (nrow(private_genes) > 0) {
      log_message(sprintf("Example private genes: %s", 
                         paste(private_genes[[gene_col]][1:min(5, nrow(private_genes))], collapse = ", ")))
    }
    
    # NOTE: Functional enrichment analysis has been moved to post-statistical analysis phase
    # This ensures statistical tests are completed before enrichment analysis per user requirements
    log_message("Functional enrichment analysis will be performed AFTER statistical testing")
    
    return(list(
      core_shared_genes = core_shared_genes,
      private_genes = private_genes,
      intermediate_genes = intermediate_genes,
      all_gene_sharing = gene_sharing,
      parameters = list(
        threshold_mode = threshold_mode,
        threshold_value = threshold_value,
        min_species_threshold = min_species_threshold,
        total_species = total_species
      )
    ))
    
  }, error = function(e) {
    log_message(sprintf("Failed to analyze core shared hotspots: %s", e$message), level = "error")
    return(list(genes = data.frame(), enrichment = NULL))
  })
}
