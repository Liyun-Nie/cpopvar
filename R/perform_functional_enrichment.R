#' Perform functional enrichment analysis
#' @param target_genes Vector of target genes for enrichment analysis
#' @param gene_function_data Data frame containing gene function information
#' @param background_genes Background genes for comparison (default: NULL)
#' @param expected_categories Expected functional categories for analysis (default: NULL)
#' @return Data frame with enrichment results
#' @export
perform_functional_enrichment <- function(target_genes, gene_function_data, background_genes = NULL, expected_categories = NULL) {
  
  tryCatch({
    log_message("Performing functional enrichment analysis with Fisher's exact test")
    
    # Use all genes in function mapping as background if not specified
    if (is.null(background_genes)) {
      background_genes <- unique(gene_function_data$gene)
    }
    
    # Filter gene function data to background genes
    background_function_data <- gene_function_data %>%
      dplyr::filter(gene %in% background_genes)
    
    # CRITICAL FIX: Force use of gene_category (4 categories) instead of functional_group (102 categories)
    # This addresses the user requirement for focused functional enrichment analysis
    if (!"gene_category" %in% names(background_function_data)) {
      log_message("CRITICAL ERROR: gene_category column missing from gene function data. Enrichment analysis requires gene_category (4 categories), not functional_group (102 categories).", level = "error")
      return(data.frame())
    }
    
    # Use gene_category for enrichment analysis with configurable categories
    # Default to standard 4 categories if not specified
    if (is.null(expected_categories)) {
      expected_categories <- c(
        "Other Functional Genes",
        "Unknown Function Genes", 
        "Photosynthesis Related",
        "Self-Replication Related"
      )
    }
    
    # Clean and standardize gene categories to avoid duplicates
    raw_categories <- background_function_data$gene_category
    log_message(sprintf("Raw gene_category data: %d total entries", length(raw_categories)))
    
    # Step 1: Basic cleaning
    functional_groups <- trimws(raw_categories)  # Remove leading/trailing whitespace
    functional_groups <- functional_groups[!is.na(functional_groups) & functional_groups != ""]
    
    # Step 2: Get unique categories before normalization
    unique_raw <- unique(functional_groups)
    log_message(sprintf("Unique categories before normalization: %d (%s)", 
                       length(unique_raw), paste(unique_raw, collapse = ", ")))
    
    # Step 3: Aggressive normalization and standardization
    functional_groups_clean <- character()
    for (fg in unique_raw) {
      # Normalize spaces, convert to title case, remove special characters
      fg_normalized <- gsub("\\s+", " ", trimws(fg))
      fg_normalized <- tools::toTitleCase(tolower(fg_normalized))
      
      # Try to match with expected categories (case-insensitive, flexible matching)
      matched_category <- NULL
      for (expected in expected_categories) {
        if (grepl(gsub(" ", ".*", tolower(expected)), tolower(fg_normalized), ignore.case = TRUE) ||
            grepl(gsub(" ", ".*", tolower(fg_normalized)), tolower(expected), ignore.case = TRUE)) {
          matched_category <- expected
          break
        }
      }
      
      # Use matched category or keep normalized version
      final_category <- if (!is.null(matched_category)) matched_category else fg_normalized
      
      if (!final_category %in% functional_groups_clean) {
        functional_groups_clean <- c(functional_groups_clean, final_category)
      }
    }
    
    functional_groups <- functional_groups_clean
    
    # Step 4: Validation - ensure we have exactly 4 expected categories
    if (length(functional_groups) != 4) {
      log_message(sprintf("WARNING: Expected 4 gene categories, found %d: (%s)", 
                         length(functional_groups), paste(functional_groups, collapse = ", ")), level = "warning")
      log_message("Expected categories: Other Functional Genes, Unknown Function Genes, Photosynthesis Related, Self-Replication Related", level = "warning")
      
      # Filter to only expected categories if possible
      valid_categories <- functional_groups[functional_groups %in% expected_categories]
      if (length(valid_categories) >= 3) {  # Accept if at least 3 of 4 expected categories found
        functional_groups <- valid_categories
        log_message(sprintf("Using %d valid categories: (%s)", 
                           length(functional_groups), paste(functional_groups, collapse = ", ")))
      }
    }
    
    log_message(sprintf("Using gene_category for enrichment analysis: %d unique categories (%s)", 
                       length(functional_groups), paste(functional_groups, collapse = ", ")))
    
    enrichment_results <- data.frame()
    
    for (func_group in functional_groups) {
      # Create contingency table
      # Rows: in functional group vs not in functional group
      # Cols: in target genes vs not in target genes
      
      # Use standardized matching with normalized gene categories
      background_categories_normalized <- gsub("\\s+", " ", trimws(background_function_data$gene_category))
      genes_in_func_group <- background_function_data$gene[background_categories_normalized == func_group]
      
      # Count genes in each category
      target_in_func <- sum(target_genes %in% genes_in_func_group)
      target_not_in_func <- length(target_genes) - target_in_func
      background_in_func <- length(genes_in_func_group)
      background_not_in_func <- length(background_genes) - background_in_func
      
      # Create contingency table
      contingency_table <- matrix(
        c(target_in_func, target_not_in_func, 
          background_in_func - target_in_func, background_not_in_func - target_not_in_func),
        nrow = 2, ncol = 2,
        dimnames = list(
          c("In_Function", "Not_In_Function"),
          c("Target", "Background")
        )
      )
      
      # Perform Fisher's exact test (only if we have enough data)
      if (all(contingency_table >= 0) && target_in_func > 0) {
        fisher_result <- fisher.test(contingency_table, alternative = "greater")
        
        # Calculate enrichment ratio
        expected_ratio <- background_in_func / length(background_genes)
        observed_ratio <- target_in_func / length(target_genes)
        enrichment_ratio <- observed_ratio / expected_ratio
        
        # Store results
        enrichment_results <- rbind(enrichment_results, data.frame(
          functional_category = func_group,  # Changed from functional_group to functional_category
          target_count = target_in_func,
          target_total = length(target_genes),
          background_count = background_in_func,
          background_total = length(background_genes),
          observed_ratio = observed_ratio,
          expected_ratio = expected_ratio,
          enrichment_ratio = enrichment_ratio,
          p_value = fisher_result$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    if (nrow(enrichment_results) > 0) {
      # Adjust p-values for multiple testing
      enrichment_results$p_adjusted <- p.adjust(enrichment_results$p_value, method = "fdr")
      
      # Sort by significance
      enrichment_results <- enrichment_results %>%
        dplyr::arrange(p_adjusted,dplyr::desc(enrichment_ratio))
      
      log_message(sprintf("Functional enrichment completed: %d gene categories tested", nrow(enrichment_results)))
      
      # Log significant results
      significant_results <- enrichment_results[enrichment_results$p_adjusted < 0.05, ]
      if (nrow(significant_results) > 0) {
        log_message(sprintf("Found %d significantly enriched gene categories", nrow(significant_results)))
      } else {
        log_message("No significantly enriched gene categories found")
      }
    }
    
    return(enrichment_results)
    
  }, error = function(e) {
    log_message(sprintf("Functional enrichment analysis failed: %s", e$message), level = "error")
    return(data.frame())
  })
}
