#' Generate Dual Factor Composite Plot
#' 
#' Creates a composite plot combining all dual-factor analysis plots
#' using intelligent grid layout with configurable maximum columns.
#' Handles larger plots with intelligent sizing for dual-factor visualizations.
#'
#' @param dual_factor_results List of dual factor analysis results
#' @param config Configuration object containing layout parameters
#' @param task_params Task-specific parameters
#' @return List containing composite plot, metadata, and summary information
#'
#' @examples
#' \dontrun{
#' # Generate composite plot from multiple dual factor analyses
#' composite_result <- generate_dual_factor_composite_plot(
#'   dual_factor_results = list(
#'     "region_type_var_type" = dual_result1,
#'     "phylogeny_region_type" = dual_result2,
#'     "life_form_var_type" = dual_result3
#'   ),
#'   config = config,
#'   task_params = task_params
#' )
#' }
#' @export
generate_dual_factor_composite_plot <- function(dual_factor_results,
                                              config,
                                              task_params = NULL) {
  
  tryCatch({
    log_message("Generating dual-factor composite plot")
    
    # Validate inputs
    if (is.null(dual_factor_results) || length(dual_factor_results) == 0) {
      warning("No dual factor results provided for composite plot")
      return(list(
        composite_plot =ggplot2::ggplot() + 
ggplot2::ggtitle("No Dual Factor Analyses Available") + 
ggplot2::labs(subtitle = "Cannot create composite plot without analysis results"),
        metadata = list(
          plot_count = 0,
          factor_pairs = character(0),
          total_observations = 0
        ),
        summary = data.frame()
      ))
    }
    
    # Extract all individual plots
    individual_plots <- list()
    factor_pairs <- character()
    total_observations <- 0
    factorial_designs <- character()
    
    for (pair_name in names(dual_factor_results)) {
      result <- dual_factor_results[[pair_name]]
      
      # Validate individual result structure
      if (!is.null(result) && !is.null(result$plot)) {
        individual_plots[[pair_name]] <- result$plot
        factor_pairs <- c(factor_pairs, pair_name)
        
        # Add observation count if available
        if (!is.null(result$metadata$total_observations)) {
          total_observations <- total_observations + result$metadata$total_observations
        }
        
        # Collect factorial design information
        if (!is.null(result$metadata$factorial_design)) {
          factorial_designs <- c(factorial_designs, result$metadata$factorial_design)
        }
        
        log_message(sprintf("  -> Added plot for factor pair: %s", pair_name))
      } else {
        log_message(sprintf("  -> Skipped invalid result for factor pair: %s", pair_name), level = "warning")
      }
    }
    
    # Check if we have any valid plots
    if (length(individual_plots) == 0) {
      warning("No valid plots found in dual factor results")
      return(list(
        composite_plot =ggplot2::ggplot() + 
ggplot2::ggtitle("No Valid Dual Factor Plots") + 
ggplot2::labs(subtitle = "All dual factor analyses failed or returned invalid plots"),
        metadata = list(
          plot_count = 0,
          factor_pairs = factor_pairs,
          total_observations = 0
        ),
        summary = data.frame()
      ))
    }
    
    # Get layout configuration - dual factor plots may need fewer columns due to complexity
    max_columns <- get_task_parameter(task_params, config, "max_columns", 3)
    if (is.null(max_columns) || !is.numeric(max_columns) || max_columns < 1) {
      max_columns <- 3  # Fallback to default
    }
    
    # Intelligent column adjustment for dual factor plots
    # Since dual factor plots are typically more complex, consider reducing max columns
    if (length(individual_plots) >= 4 && max_columns >= 3) {
      adjusted_max_cols <- max(2, max_columns - 1)  # Reduce by 1 but minimum of 2
      log_message(sprintf("Adjusting max columns for dual-factor plots: %d -> %d", 
                         max_columns, adjusted_max_cols))
      max_columns <- adjusted_max_cols
    }
    
    # Dependencies automatically loaded in R package context
    
    # Create composite plot title
    composite_title <- "Dual-Factor Analysis"
    
    # Generate composite plot using smart grid layout
    composite_plot <- arrange_plots_smart_grid(
      plot_list = individual_plots,
      max_cols = max_columns,
      auto_arrange = TRUE,
      common_title = composite_title,
      title_size = 14
    )
    
    log_message(sprintf("[SUCCESS] Created dual-factor composite plot with %d individual plots (max_cols=%d)", 
                       length(individual_plots), max_columns))
    
    # Create summary statistics
    summary_data <- data.frame(
      factor_pair = factor_pairs,
      plot_included = TRUE,
      stringsAsFactors = FALSE
    )
    
    # Add metadata from individual results if available
    for (i in seq_along(factor_pairs)) {
      pair_name <- factor_pairs[i]
      result <- dual_factor_results[[pair_name]]
      
      if (!is.null(result$metadata)) {
        # Extract factor names from metadata
        factors <- result$metadata$factors_analyzed
        if (!is.null(factors) && length(factors) >= 2) {
          summary_data$factor1[i] <- factors[1]
          summary_data$factor2[i] <- factors[2]
        }
        
        summary_data$factorial_design[i] <- result$metadata$factorial_design %||% "Unknown"
        summary_data$total_observations[i] <- result$metadata$total_observations %||% NA
        summary_data$statistical_test[i] <- result$metadata$statistical_test %||% "Unknown"
        summary_data$plot_type[i] <- result$metadata$plot_type %||% "Unknown"
        
        # Add factor-specific information
        summary_data$n_factor1_categories[i] <- result$metadata$n_factor1_categories %||% NA
        summary_data$n_factor2_categories[i] <- result$metadata$n_factor2_categories %||% NA
      }
    }
    
    # Create comprehensive metadata
    composite_metadata <- list(
      plot_count = length(individual_plots),
      factor_pairs = factor_pairs,
      max_columns_used = max_columns,
      total_unique_observations = total_observations,
      factorial_designs = factorial_designs,
      layout_type = "smart_grid_adjusted",
      composite_title = composite_title,
      created_timestamp = Sys.time(),
      complexity_adjustment = ifelse(length(individual_plots) >= 4, "column_reduction_applied", "no_adjustment")
    )
    
    # Return comprehensive results
    return(list(
      composite_plot = composite_plot,
      individual_plots = individual_plots,
      metadata = composite_metadata,
      summary = summary_data
    ))
    
  }, error = function(e) {
    log_message(sprintf("Failed to generate dual-factor composite plot: %s", e$message), level = "error")
    
    # Return error plot
    error_plot <-ggplot2::ggplot() + 
ggplot2::ggtitle("Dual-Factor Composite Plot Generation Failed") + 
ggplot2::labs(
        subtitle = paste("Error:", e$message),
        caption = "Check individual dual-factor analysis results and configuration"
      ) +
      theme_minimal() +
ggplot2::theme(
        plot.title =ggplot2::element_text(color = "red", size = 12, hjust = 0.5),
        plot.subtitle =ggplot2::element_text(color = "darkred", size = 10, hjust = 0.5),
        plot.caption =ggplot2::element_text(size = 8, hjust = 0.5)
      )
    
    return(list(
      composite_plot = error_plot,
      individual_plots = list(),
      metadata = list(
        plot_count = 0,
        factor_pairs = character(0),
        error = e$message,
        total_observations = 0
      ),
      summary = data.frame()
    ))
  })
}

# Helper function for null coalescing (if not already defined)
if (!exists("%||%")) {
  `%||%` <- function(x, y) {
    if (is.null(x) || (is.atomic(x) && length(x) == 1 && is.na(x))) y else x
  }
}