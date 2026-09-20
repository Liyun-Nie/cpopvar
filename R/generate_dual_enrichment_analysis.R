#' generate_dual_enrichment_analysis
#' @param hotspot_analysis_results Results from hotspot analysis
#' @param gene_function_data Gene function mapping data
#' @param background_genes Background genes for analysis
#' @param config Configuration list
#' @param task_params Task-specific parameters
#' @return List containing dual enrichment analysis results
#' @export
generate_dual_enrichment_analysis <- function(hotspot_analysis_results,
                                             gene_function_data,
                                             background_genes = NULL,
                                             config,
                                             task_params = NULL) {
  
  tryCatch({
    log_message("M03 Step 3b: Performing dual enrichment analysis (core shared + private genes independently)")
    
    results <- list()
    
    # Use internal standard column name (boundary principle)
    gene_col <- "gene"
    
    # Extract both core shared and private genes for independent analysis
    core_shared_genes <- hotspot_analysis_results$core_shared_genes[[gene_col]]
    private_genes <- hotspot_analysis_results$private_genes[[gene_col]]
    
    # Log the counts to verify we have both gene sets
    log_message(sprintf("Dual enrichment analysis: %d core shared genes (>=3 species), %d private genes (1 species)", 
                       length(core_shared_genes), length(private_genes)))
    
    # Core shared genes enrichment (main analysis)
    if (length(core_shared_genes) > 0) {
      log_message("Analyzing functional enrichment for core shared hotspot genes")
      # Get expected categories from configuration
      expected_categories <- get_task_parameter(task_params, config, "expected_gene_categories", 
                                               c("Other Functional Genes", "Unknown Function Genes", 
                                                 "Photosynthesis Related", "Self-Replication Related"))
      
      core_enrichment <- perform_functional_enrichment(
        target_genes = core_shared_genes,
        gene_function_data = gene_function_data,
        background_genes = background_genes,
        expected_categories = expected_categories
      )
      
      if (nrow(core_enrichment) > 0) {
        core_enrichment$gene_set_type <- "Core Shared Hotspots"
        results$core_shared_enrichment <- core_enrichment
        
        significant_core <- sum(core_enrichment$p_adjusted < 0.05, na.rm = TRUE)
        log_message(sprintf("Core shared genes enrichment: %d gene categories tested, %d significant (p<0.05)", 
                           nrow(core_enrichment), significant_core))
        
        # Report top significant functions
        if (significant_core > 0) {
          top_functions <- core_enrichment %>% 
dplyr::filter(p_adjusted < 0.05) %>% 
dplyr::arrange(p_adjusted) %>% 
            head(3)
          log_message(sprintf("Top enriched gene categories: %s", 
                             paste(top_functions$functional_category, collapse = ", ")))
        }
      } else {
        log_message("No enrichment results for core shared genes", level = "warning")
        results$core_shared_enrichment <- data.frame()
      }
    } else {
      log_message("No core shared genes available for enrichment analysis", level = "warning")
      results$core_shared_enrichment <- data.frame()
    }
    
    # Private genes enrichment (independent analysis)
    if (length(private_genes) > 0) {
      log_message("Analyzing functional enrichment for private hotspot genes")
      private_enrichment <- perform_functional_enrichment(
        target_genes = private_genes,
        gene_function_data = gene_function_data,
        background_genes = background_genes,
        expected_categories = expected_categories
      )
      
      if (nrow(private_enrichment) > 0) {
        private_enrichment$gene_set_type <- "Private Hotspots"
        results$private_enrichment <- private_enrichment
        
        significant_private <- sum(private_enrichment$p_adjusted < 0.05, na.rm = TRUE)
        log_message(sprintf("Private genes enrichment: %d gene categories tested, %d significant (p<0.05)", 
                           nrow(private_enrichment), significant_private))
        
        # Report top significant functions
        if (significant_private > 0) {
          top_functions <- private_enrichment %>% 
dplyr::filter(p_adjusted < 0.05) %>% 
dplyr::arrange(p_adjusted) %>% 
            head(3)
          log_message(sprintf("Top enriched gene categories (private): %s", 
                             paste(top_functions$functional_category, collapse = ", ")))
        }
      } else {
        log_message("No enrichment results for private genes", level = "warning")
        results$private_enrichment <- data.frame()
      }
    } else {
      log_message("No private genes available for enrichment analysis", level = "warning")
      results$private_enrichment <- data.frame()
    }
    
    # Create combined enrichment by merging both analyses (for compatibility)
    combined_enrichment <- data.frame()
    if ("core_shared_enrichment" %in% names(results) && nrow(results$core_shared_enrichment) > 0) {
      combined_enrichment <- rbind(combined_enrichment, results$core_shared_enrichment)
    }
    if ("private_enrichment" %in% names(results) && nrow(results$private_enrichment) > 0) {
      combined_enrichment <- rbind(combined_enrichment, results$private_enrichment)
    }
    results$combined_enrichment <- combined_enrichment
    
    metadata <- list(
      core_shared_genes_count = length(core_shared_genes),
      private_genes_count = length(private_genes),
      core_enrichment_groups = ifelse("core_shared_enrichment" %in% names(results), 
                                    nrow(results$core_shared_enrichment), 0),
      private_enrichment_groups = ifelse("private_enrichment" %in% names(results), 
                                       nrow(results$private_enrichment), 0),
      has_significant_core = ifelse("core_shared_enrichment" %in% names(results), 
                                  any(results$core_shared_enrichment$p_adjusted < 0.05, na.rm = TRUE), FALSE),
      has_significant_private = ifelse("private_enrichment" %in% names(results), 
                                     any(results$private_enrichment$p_adjusted < 0.05, na.rm = TRUE), FALSE),
      analysis_focus = "dual_independent"  # Flag to indicate dual independent analysis
    )
    
    results$metadata <- metadata
    
    log_message("Dual enrichment analysis completed")
    return(results)
    
  }, error = function(e) {
    log_message(sprintf("Failed to perform dual enrichment analysis: %s", e$message), level = "error")
    return(list())
  })
}
