#' Generate significance star notation
#' 
#' @title Generate significance star notation
#' @description Converts p-value matrix to significance star notation
#' @param pvalue_matrix Matrix of p-values to convert
#' @param significance_levels Named list defining star levels and thresholds
#' @param show_nonsignificant Whether to show non-significant values (default: FALSE)
#' @return Matrix of significance stars
#' @export
generate_significance_stars <- function(pvalue_matrix, 
                                     significance_levels = list(
                                       "***" = 0.001,
                                       "**" = 0.01, 
                                       "*" = 0.05,
                                       "." = 0.1,
                                       "ns" = 1.0
                                     ),
                                     show_nonsignificant = FALSE) {
  
  stars_matrix <- matrix("", nrow = nrow(pvalue_matrix), ncol = ncol(pvalue_matrix))
  rownames(stars_matrix) <- rownames(pvalue_matrix)
  colnames(stars_matrix) <- colnames(pvalue_matrix)
  
  # Apply significance levels in order (most stringent first)
  sorted_levels <- significance_levels[order(unlist(significance_levels))]
  
  for (symbol in names(sorted_levels)) {
    threshold <- sorted_levels[[symbol]]
    
    if (symbol == "ns" && !show_nonsignificant) {
      # Skip non-significant marker if not requested
      next
    }
    
    if (symbol == "ns") {
      # Non-significant: everything above the previous threshold
      prev_threshold <- if (length(sorted_levels) > 1) {
        max(unlist(sorted_levels[-length(sorted_levels)]))
      } else {
        0
      }
      stars_matrix[pvalue_matrix > prev_threshold] <- symbol
    } else {
      # Significant levels: at or below threshold (only mark if not already marked)
      stars_matrix[pvalue_matrix <= threshold & stars_matrix == ""] <- symbol
    }
  }
  
  # Clear diagonal (self-correlations)
  diag(stars_matrix) <- ""
  
  return(stars_matrix)
}
