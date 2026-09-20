#' Calculate correlation p-values for data matrix
#' 
#' @title Calculate correlation p-values for data matrix
#' @description Calculates pairwise correlation p-values between rows of a data matrix
#' @param data_matrix Numeric matrix for correlation analysis
#' @param method Correlation method to use (default: "spearman")
#' @return Matrix of p-values for pairwise correlations
#' @export
calculate_correlation_pvalues <- function(data_matrix, method = "spearman") {
  tryCatch({
    n_vars <- nrow(data_matrix)
    pvalue_matrix <- matrix(1, nrow = n_vars, ncol = n_vars)
    rownames(pvalue_matrix) <- rownames(data_matrix)
    colnames(pvalue_matrix) <- rownames(data_matrix)
    
    # Calculate pairwise correlation p-values
    for (i in 1:(n_vars-1)) {
      for (j in (i+1):n_vars) {
        x <- data_matrix[i, ]
        y <- data_matrix[j, ]
        
        # Remove NA pairs
        valid_pairs <- !is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y)
        
        if (sum(valid_pairs) >= 3) {  # Need at least 3 points for meaningful correlation
          test_result <- tryCatch({
            stats::cor.test(x[valid_pairs], y[valid_pairs], method = method)
          }, error = function(e) {
            list(p.value = 1)  # Non-significant if test fails
          })
          
          pvalue_matrix[i, j] <- test_result$p.value
          pvalue_matrix[j, i] <- test_result$p.value
        }
      }
    }
    
    # Diagonal elements are always 0 (self-correlation is always significant)
    diag(pvalue_matrix) <- 0
    
    return(pvalue_matrix)
    
  }, error = function(e) {
    log_message(sprintf("Failed to calculate correlation p-values: %s", e$message), level = "warning")
    return(NULL)
  })
}
