#' Perform comprehensive statistical analysis
#' 
#' @title Perform comprehensive statistical analysis
#' @description Conducts multi-factor statistical analysis with multiple comparison corrections
#' @param data Dataset for analysis
#' @param response_var Name of the response variable (default: "frequency_per_kb")
#' @param grouping_vars Vector of grouping variable names
#' @param correction_method Multiple comparison correction method (default: "bonferroni")
#' @param alpha_level Significance threshold (default: 0.05)
#' @return List containing comprehensive statistical results
#' @export
perform_comprehensive_statistical_analysis <- function(data,
                                                      response_var = "frequency_per_kb",
                                                      grouping_vars,
                                                      correction_method = "bonferroni",
                                                      alpha_level = 0.05) {
  
  tryCatch({
    log_message("Performing comprehensive statistical analysis")
    
    results <- list()
    
    # Single factor analyses
    for (var in grouping_vars) {
      if (var %in% names(data)) {
        groups <- unique(data[[var]])
        if (length(groups) > 1) {
          
          # Kruskal-Wallis test
          kw_result <- kruskal.test(data[[response_var]], data[[var]])
          results[[paste0(var, "_kruskal_wallis")]] <- kw_result
          
          # Pairwise comparisons if significant
          if (kw_result$p.value < alpha_level && length(groups) > 2) {
            pairwise_result <- pairwise.wilcox.test(
              data[[response_var]], 
              data[[var]], 
              p.adjust.method = correction_method
            )
            results[[paste0(var, "_pairwise")]] <- pairwise_result
          }
        }
      }
    }
    
    # Two-factor interaction analysis
    if (length(grouping_vars) >= 2) {
      var1 <- grouping_vars[1]
      var2 <- grouping_vars[2]
      
      if (var1 %in% names(data) && var2 %in% names(data)) {
        # Two-way ANOVA
        formula_str <- sprintf("%s ~ %s * %s", response_var, var1, var2)
        anova_result <- tryCatch({
          aov_result <- aov(as.formula(formula_str), data = data)
          list(
            anova = aov_result,
            summary = summary(aov_result),
            tukey = TukeyHSD(aov_result)
          )
        }, error = function(e) {
          log_message("Two-way ANOVA failed", level = "warning")
          NULL
        })
        results$two_way_anova <- anova_result
      }
    }
    
    # Effect size calculations
    if (requireNamespace("effsize", quietly = TRUE)) {
      for (var in grouping_vars) {
        if (var %in% names(data)) {
          groups <- unique(data[[var]])
          if (length(groups) == 2) {
            # Cohen's d for two groups
            group1_data <- data[data[[var]] == groups[1], response_var]
            group2_data <- data[data[[var]] == groups[2], response_var]
            
            cohen_d <- tryCatch({
              effsize::cohen.d(group1_data, group2_data)
            }, error = function(e) NULL)
            
            if (!is.null(cohen_d)) {
              results[[paste0(var, "_cohens_d")]] <- cohen_d
            }
          }
        }
      }
    }
    
    return(results)
    
  }, error = function(e) {
    log_message(sprintf("Statistical analysis failed: %s", e$message), level = "error")
    return(list(error = e$message))
  })
}
