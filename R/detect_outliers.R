############################################################
#### detect_outliers.R - Unified Outlier Detection System ####
############################################################
#
# Comprehensive outlier detection using multiple statistical methods
# This function unifies and replaces:
# - M01_functions/enhanced_outlier_detection.R (archived)
# - common_helpers/detect_outliers_by_species.R (archived)  
# - src/utils/outlier_detection.R (reference implementation)
#
# Design: Pure utility function with complete dependency injection
############################################################

#' Unified outlier detection function
#' 
#' Detects outliers using multiple statistical methods with full configurability
#' through dependency injection. Supports both single-group and grouped analysis.
#' 
#' @param data Data frame containing the data to analyze
#' @param value_col Column name for values to analyze (default: "frequency_per_kb")
#' @param group_col Column name for grouping (NULL for global analysis, default: NULL)
#' @param method Detection method ("IQR", "Z_score", "modified_Z", default: "IQR")
#' @param iqr_multiplier Multiplier for IQR method (default: 1.5)
#' @param z_threshold Threshold for Z-score method (default: 3.0)
#' @param modified_z_threshold Threshold for modified Z-score method (default: 3.5)
#' @param min_group_size Minimum group size for outlier detection (default: 4)
#' @param return_type Return format ("outliers_only", "full_analysis", default: "outliers_only")
#' @return Data frame with outliers or list with full analysis results
#' @importFrom magrittr %>%
#' @export
detect_outliers <- function(data,
                           value_col = "frequency_per_kb",
                           group_col = NULL,
                           method = "IQR",
                           iqr_multiplier = 1.5,
                           z_threshold = 3.0,
                           modified_z_threshold = 3.5,
                           min_group_size = 4,
                           return_type = "outliers_only") {
  
  tryCatch({
    # Input validation
    if (is.null(data) || nrow(data) == 0) {
      log_message("No data provided for outlier detection", level = "warning")
      return(if (return_type == "outliers_only") data.frame() else create_empty_outlier_result())
    }
    
    if (!value_col %in% names(data)) {
      stop(sprintf("Value column '%s' not found in data", value_col))
    }
    
    if (!is.null(group_col) && !group_col %in% names(data)) {
      stop(sprintf("Group column '%s' not found in data", group_col))
    }
    
    # Remove NA values for analysis
    clean_data <- data[!is.na(data[[value_col]]), ]
    
    if (nrow(clean_data) == 0) {
      log_message("No non-NA values found for outlier detection", level = "warning")
      return(if (return_type == "outliers_only") data.frame() else create_empty_outlier_result())
    }
    
    log_message(sprintf("Starting outlier detection: method=%s, groups=%s, n=%d", 
                       method, ifelse(is.null(group_col), "none", "grouped"), nrow(clean_data)))
    
    # Perform outlier detection
    if (is.null(group_col)) {
      # Global outlier detection
      outlier_results <- detect_outliers_single_group(
        clean_data, value_col, method, iqr_multiplier, z_threshold, modified_z_threshold
      )
    } else {
      # Group-wise outlier detection
      outlier_results <- detect_outliers_by_group(
        clean_data, value_col, group_col, method, iqr_multiplier, 
        z_threshold, modified_z_threshold, min_group_size
      )
    }
    
    # Return results based on requested format
    if (return_type == "outliers_only") {
      # Extract just the outlier rows
      if (is.null(group_col)) {
        outlier_mask <- outlier_results$is_outlier
        outliers <- clean_data[outlier_mask, ]
      } else {
        outliers <- outlier_results$data[outlier_results$data$is_outlier, ]
        outliers$is_outlier <- NULL  # Remove the temporary column
      }
      
      # Add metadata columns
      if (nrow(outliers) > 0) {
        outliers$outlier_method <- method
        outliers$detection_timestamp <- Sys.time()
      }
      
      log_message(sprintf("Outlier detection completed: %d outliers detected using %s method", 
                         nrow(outliers), method))
      return(outliers)
      
    } else {
      # Return full analysis results
      log_message(sprintf("Outlier detection completed: %d outliers from %d total points", 
                         sum(outlier_results$is_outlier, na.rm = TRUE), nrow(clean_data)))
      return(outlier_results)
    }
    
  }, error = function(e) {
    log_message(sprintf("Outlier detection failed: %s", e$message), level = "error")
    return(if (return_type == "outliers_only") data.frame() else list(error = e$message))
  })
}

#' Detect outliers for single group (global analysis)
#' 
#' @param data Data frame
#' @param value_col Value column name
#' @param method Detection method
#' @param iqr_multiplier IQR multiplier
#' @param z_threshold Z-score threshold
#' @param modified_z_threshold Modified Z-score threshold
#' @return List with outlier analysis results
detect_outliers_single_group <- function(data, value_col, method, iqr_multiplier, 
                                        z_threshold, modified_z_threshold) {
  
  values <- data[[value_col]]
  n_total <- length(values)
  
  # Calculate outliers by method
  outlier_info <- calculate_outliers_by_method(values, method, iqr_multiplier, 
                                              z_threshold, modified_z_threshold)
  
  result <- list(
    data = data,
    is_outlier = outlier_info$is_outlier,
    outlier_info = outlier_info,
    method = method,
    summary = create_outlier_summary(values, outlier_info$is_outlier),
    n_total = n_total,
    n_outliers = sum(outlier_info$is_outlier, na.rm = TRUE)
  )
  
  return(result)
}

#' Detect outliers by group
#' 
#' @param data Data frame
#' @param value_col Value column name
#' @param group_col Group column name
#' @param method Detection method
#' @param iqr_multiplier IQR multiplier
#' @param z_threshold Z-score threshold
#' @param modified_z_threshold Modified Z-score threshold
#' @param min_group_size Minimum group size for analysis
#' @return List with outlier analysis results
detect_outliers_by_group <- function(data, value_col, group_col, method, iqr_multiplier,
                                    z_threshold, modified_z_threshold, min_group_size) {
  
  # Initialize empty result for outliers
  all_outliers <- data.frame()
  
  # Process each group separately
  unique_groups <- unique(data[[group_col]])
  
  for (group in unique_groups) {
    group_data <- data[data[[group_col]] == group, ]
    
    if (nrow(group_data) < min_group_size) {
      # For small groups, apply simplified detection if at least 2 points exist
      if (nrow(group_data) >= 2) {
        values <- group_data[[value_col]]
        mean_val <- mean(values, na.rm = TRUE)
        sd_val <- stats::sd(values, na.rm = TRUE)
        
        if (!is.na(sd_val) && sd_val > 0) {
          # Use 2*SD for very small groups
          outlier_mask <- abs(values - mean_val) > 2 * sd_val
          if (any(outlier_mask, na.rm = TRUE)) {
            outliers <- group_data[outlier_mask, ]
            all_outliers <- rbind(all_outliers, outliers)
          }
        }
      }
      next  # Skip to next group
    }
    
    # Calculate outliers for this group using specified method
    values <- group_data[[value_col]]
    outlier_info <- calculate_outliers_by_method(values, method, iqr_multiplier,
                                                z_threshold, modified_z_threshold)
    
    # Add outliers to results
    if (any(outlier_info$is_outlier, na.rm = TRUE)) {
      outliers <- group_data[outlier_info$is_outlier, ]
      all_outliers <- rbind(all_outliers, outliers)
    }
  }
  
  # Create summary by group
  group_summary <- data %>%
    dplyr::group_by(!!rlang::sym(group_col)) %>%
    dplyr::summarise(
      n_total = dplyr::n(),
      n_outliers = sum(.data[[value_col]] %in% all_outliers[[value_col]], na.rm = TRUE),
      outlier_percentage = round(100 * n_outliers / n_total, 2),
      mean_value = mean(.data[[value_col]], na.rm = TRUE),
      median_value = stats::median(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create combined result with outlier flags
  data_with_flags <- data
  data_with_flags$is_outlier <- data[[value_col]] %in% all_outliers[[value_col]]
  
  result <- list(
    data = data_with_flags,
    is_outlier = data_with_flags$is_outlier,
    method = method,
    group_summary = group_summary,
    overall_summary = list(
      n_total = nrow(data),
      n_outliers = nrow(all_outliers),
      outlier_percentage = round(100 * nrow(all_outliers) / nrow(data), 2)
    )
  )
  
  return(result)
}

#' Calculate outliers using specific method
#' 
#' @param values Numeric vector
#' @param method Detection method
#' @param iqr_multiplier IQR multiplier
#' @param z_threshold Z-score threshold
#' @param modified_z_threshold Modified Z-score threshold
#' @return List with outlier information
calculate_outliers_by_method <- function(values, method, iqr_multiplier, 
                                       z_threshold, modified_z_threshold) {
  
  n <- length(values)
  is_outlier <- rep(FALSE, n)
  outlier_info <- list()
  
  if (method == "IQR") {
    # Interquartile range method
    q1 <- stats::quantile(values, 0.25, na.rm = TRUE)
    q3 <- stats::quantile(values, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    
    lower_bound <- q1 - iqr_multiplier * iqr
    upper_bound <- q3 + iqr_multiplier * iqr
    
    is_outlier <- values < lower_bound | values > upper_bound
    
    outlier_info <- list(
      method = "IQR",
      q1 = q1,
      q3 = q3,
      iqr = iqr,
      multiplier = iqr_multiplier,
      lower_bound = lower_bound,
      upper_bound = upper_bound,
      is_outlier = is_outlier
    )
    
  } else if (method == "Z_score") {
    # Z-score method
    mean_val <- mean(values, na.rm = TRUE)
    sd_val <- stats::sd(values, na.rm = TRUE)
    
    if (sd_val == 0) {
      z_scores <- rep(0, n)
    } else {
      z_scores <- abs(values - mean_val) / sd_val
    }
    
    is_outlier <- z_scores > z_threshold
    
    outlier_info <- list(
      method = "Z_score",
      mean = mean_val,
      sd = sd_val,
      threshold = z_threshold,
      z_scores = z_scores,
      is_outlier = is_outlier
    )
    
  } else if (method == "modified_Z") {
    # Modified Z-score method (using median)
    median_val <- stats::median(values, na.rm = TRUE)
    mad_val <- stats::mad(values, na.rm = TRUE)
    
    if (mad_val == 0) {
      modified_z_scores <- rep(0, n)
    } else {
      modified_z_scores <- 0.6745 * (values - median_val) / mad_val
    }
    
    is_outlier <- abs(modified_z_scores) > modified_z_threshold
    
    outlier_info <- list(
      method = "modified_Z",
      median = median_val,
      mad = mad_val,
      threshold = modified_z_threshold,
      modified_z_scores = modified_z_scores,
      is_outlier = is_outlier
    )
    
  } else {
    stop("Unknown outlier detection method: ", method)
  }
  
  return(outlier_info)
}

#' Create outlier summary statistics
#' 
#' @param values Numeric vector
#' @param is_outlier Logical vector indicating outliers
#' @return List with summary statistics
create_outlier_summary <- function(values, is_outlier) {
  
  n_total <- length(values)
  n_outliers <- sum(is_outlier, na.rm = TRUE)
  
  summary_stats <- list(
    n_total = n_total,
    n_outliers = n_outliers,
    outlier_percentage = round(100 * n_outliers / n_total, 2),
    
    # Overall statistics
    overall = list(
      mean = mean(values, na.rm = TRUE),
      median = stats::median(values, na.rm = TRUE),
      sd = stats::sd(values, na.rm = TRUE),
      min = min(values, na.rm = TRUE),
      max = max(values, na.rm = TRUE),
      q1 = stats::quantile(values, 0.25, na.rm = TRUE),
      q3 = stats::quantile(values, 0.75, na.rm = TRUE)
    ),
    
    # Non-outlier statistics
    non_outliers = if (n_outliers < n_total) {
      clean_values <- values[!is_outlier]
      list(
        mean = mean(clean_values, na.rm = TRUE),
        median = stats::median(clean_values, na.rm = TRUE),
        sd = stats::sd(clean_values, na.rm = TRUE),
        min = min(clean_values, na.rm = TRUE),
        max = max(clean_values, na.rm = TRUE)
      )
    } else {
      NULL
    },
    
    # Outlier statistics
    outliers = if (n_outliers > 0) {
      outlier_values <- values[is_outlier]
      list(
        mean = mean(outlier_values, na.rm = TRUE),
        median = stats::median(outlier_values, na.rm = TRUE),
        sd = stats::sd(outlier_values, na.rm = TRUE),
        min = min(outlier_values, na.rm = TRUE),
        max = max(outlier_values, na.rm = TRUE),
        values = outlier_values
      )
    } else {
      NULL
    }
  )
  
  return(summary_stats)
}

#' Create empty outlier result for edge cases
#' 
#' @return Empty outlier result list
create_empty_outlier_result <- function() {
  return(list(
    data = data.frame(),
    is_outlier = logical(0),
    method = "none",
    summary = list(n_total = 0, n_outliers = 0, outlier_percentage = 0),
    n_total = 0,
    n_outliers = 0
  ))
}

# Log successful loading
# log_message("Unified outlier detection system loaded successfully")