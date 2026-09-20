#' Create annotated correlation text matrix
#' 
#' @title Create annotated correlation text matrix  
#' @description Creates text annotations combining correlation values with significance stars
#' @param cor_matrix Matrix of correlation coefficients
#' @param pvalue_matrix Matrix of p-values corresponding to correlations
#' @return Matrix of formatted text combining correlations and significance
#' @export
create_annotated_correlation_text <- function(cor_matrix, pvalue_matrix) {
  if (is.null(pvalue_matrix)) {
    log_message("Warning: P-value matrix is NULL, showing correlations without significance", level = "warning")
    return(sprintf("%.2f", cor_matrix))
  }
  
  stars <- generate_significance_stars(pvalue_matrix)
  
  # Debug: count significant correlations
  n_significant <- sum(pvalue_matrix <= 0.05, na.rm = TRUE)
  log_message(sprintf("Correlation significance: %d out of %d pairs significant (p <= 0.05)", 
                     n_significant, length(pvalue_matrix)))
  
  annotated_matrix <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  rownames(annotated_matrix) <- rownames(cor_matrix)
  colnames(annotated_matrix) <- colnames(cor_matrix)
  
  for (i in 1:nrow(cor_matrix)) {
    for (j in 1:ncol(cor_matrix)) {
      if (i == j) {
        annotated_matrix[i, j] <- "1.00"  # Diagonal elements (perfect self-correlation)
      } else {
        annotated_matrix[i, j] <- sprintf("%.2f%s", cor_matrix[i, j], stars[i, j])
      }
    }
  }
  
  return(annotated_matrix)
}
