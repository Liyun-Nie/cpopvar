#' Perform statistical comparisons between groups
#' 
#' @title Perform statistical comparisons between groups
#' @description Conducts statistical tests to compare values between different groups
#' @param data Dataset containing group and value columns
#' @param group_col Name of the grouping column
#' @param value_col Name of the value column to compare
#' @return List containing statistical test results
#' @export
perform_statistical_comparisons <- function(data, group_col, value_col) {
  tryCatch({
    # Check if ggpubr is available
    if (!requireNamespace("ggpubr", quietly = TRUE)) {
      return(list(error = "ggpubr package not available"))
    }
    
    group_levels <- unique(data[[group_col]])
    results <- list()
    
    # Overall test (Kruskal-Wallis if more than 2 groups, Wilcox if 2 groups)
    if (length(group_levels) > 2) {
      kruskal_result <- kruskal.test(formula(paste(value_col, "~", group_col)), data = data)
      results$overall_test <- list(
        method = "Kruskal-Wallis",
        statistic = kruskal_result$statistic,
        p.value = kruskal_result$p.value,
        significant = kruskal_result$p.value < 0.05
      )
    } else if (length(group_levels) == 2) {
      wilcox_result <- wilcox.test(formula(paste(value_col, "~", group_col)), data = data)
      results$overall_test <- list(
        method = "Wilcoxon",
        statistic = wilcox_result$statistic,
        p.value = wilcox_result$p.value,
        significant = wilcox_result$p.value < 0.05
      )
    }
    
    # Pairwise comparisons (if more than 1 group and reasonable number)
    if (length(group_levels) >= 2 && length(group_levels) <= 8) {
      pairwise_results <- list()
      comparisons <- utils::combn(group_levels, 2, simplify = FALSE)
      
      for (i in seq_along(comparisons)) {
        pair <- comparisons[[i]]
        subset_data <- data[data[[group_col]] %in% pair, ]
        
        if (nrow(subset_data) > 3) {  # Need minimum data for test
          pair_test <- wilcox.test(formula(paste(value_col, "~", group_col)), data = subset_data)
          
          comparison_name <- paste(pair[1], "vs", pair[2])
          pairwise_results[[comparison_name]] <- list(
            groups = pair,
            method = "Wilcoxon",
            statistic = pair_test$statistic,
            p.value = pair_test$p.value,
            significant = pair_test$p.value < 0.05,
            p.signif = dplyr::case_when(
              pair_test$p.value < 0.001 ~ "***",
              pair_test$p.value < 0.01 ~ "**",
              pair_test$p.value < 0.05 ~ "*",
              pair_test$p.value < 0.1 ~ ".",
              TRUE ~ "ns"
            )
          )
        }
      }
      
      # Apply multiple testing correction
      if (length(pairwise_results) > 1) {
        p_values <- sapply(pairwise_results, function(x) x$p.value)
        adjusted_p <- p.adjust(p_values, method = "BH")  # Benjamini-Hochberg
        
        for (i in seq_along(pairwise_results)) {
          pairwise_results[[i]]$p.adjusted <- adjusted_p[i]
          pairwise_results[[i]]$significant_adjusted <- adjusted_p[i] < 0.05
          pairwise_results[[i]]$p.signif_adjusted <- dplyr::case_when(
            adjusted_p[i] < 0.001 ~ "***",
            adjusted_p[i] < 0.01 ~ "**", 
            adjusted_p[i] < 0.05 ~ "*",
            adjusted_p[i] < 0.1 ~ ".",
            TRUE ~ "ns"
          )
        }
      }
      
      results$pairwise_tests <- pairwise_results
    }
    
    # Add summary information
    results$summary <- list(
      total_groups = length(group_levels),
      group_names = group_levels,
      total_comparisons = if(length(group_levels) >= 2) choose(length(group_levels), 2) else 0,
      significant_comparisons = if(!is.null(results$pairwise_tests)) {
        sum(sapply(results$pairwise_tests, function(x) x$significant))
      } else 0
    )
    
    return(results)
    
  }, error = function(e) {
    return(list(error = paste("Statistical comparison failed:", e$message)))
  })
}
