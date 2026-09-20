#' Generate Single Factor Composite Plot
#' 
#' Creates a composite plot combining all single-factor analysis plots
#' using intelligent grid layout with configurable maximum columns.
#'
#' @param single_factor_results List of single factor analysis results
#' @param config Configuration object containing layout parameters
#' @param task_params Task-specific parameters
#' @return List containing composite plot, metadata, and summary information
#'
#' @examples
#' \dontrun{
#' # Generate composite plot from multiple single factor analyses
#' composite_result <- generate_single_factor_composite_plot(
#'   single_factor_results = list(
#'     region_type = single_region_result,
#'     var_type = single_var_result,
#'     phylogeny = single_phylogeny_result
#'   ),
#'   config = config,
#'   task_params = task_params
#' )
#' }
#' @export
generate_single_factor_composite_plot <- function(single_factor_results,
                                                config,
                                                task_params = NULL) {
  
  tryCatch({
    log_message("Generating single-factor composite plot")
    
    # Validate inputs
    if (is.null(single_factor_results) || length(single_factor_results) == 0) {
      warning("No single factor results provided for composite plot")
      return(list(
        composite_plot =ggplot2::ggplot() + 
ggplot2::ggtitle("No Single Factor Analyses Available") + 
ggplot2::labs(subtitle = "Cannot create composite plot without analysis results"),
        metadata = list(
          plot_count = 0,
          factor_names = character(0),
          total_observations = 0
        ),
        summary = data.frame()
      ))
    }
    
    # Extract all individual plots
    individual_plots <- list()
    factor_names <- character()
    total_observations <- 0
    
    for (factor_name in names(single_factor_results)) {
      result <- single_factor_results[[factor_name]]
      
      # Validate individual result structure
      if (!is.null(result) && !is.null(result$plot)) {
        individual_plots[[factor_name]] <- result$plot
        factor_names <- c(factor_names, factor_name)
        
        # Add observation count if available
        if (!is.null(result$metadata$total_observations)) {
          total_observations <- total_observations + result$metadata$total_observations
        }
        
        log_message(sprintf("  -> Added plot for factor: %s", factor_name))
      } else {
        log_message(sprintf("  -> Skipped invalid result for factor: %s", factor_name), level = "warning")
      }
    }
    
    # Check if we have any valid plots
    if (length(individual_plots) == 0) {
      warning("No valid plots found in single factor results")
      return(list(
        composite_plot =ggplot2::ggplot() + 
ggplot2::ggtitle("No Valid Single Factor Plots") + 
ggplot2::labs(subtitle = "All single factor analyses failed or returned invalid plots"),
        metadata = list(
          plot_count = 0,
          factor_names = factor_names,
          total_observations = 0
        ),
        summary = data.frame()
      ))
    }
    
    # Get layout configuration
    max_columns <- get_task_parameter(task_params, config, "max_columns", 3)
    if (is.null(max_columns) || !is.numeric(max_columns) || max_columns < 1) {
      max_columns <- 3  # Fallback to default
    }
    
    # Dependencies automatically loaded in R package context
    
    # Create composite plot title
    composite_title <- "Single-Factor Analysis"
    
    # Generate composite plot using smart grid layout
    composite_plot <- arrange_plots_smart_grid(
      plot_list = individual_plots,
      max_cols = max_columns,
      auto_arrange = TRUE,
      common_title = composite_title,
      title_size = 14
    )
    
    log_message(sprintf("[SUCCESS] Created single-factor composite plot with %d individual plots (max_cols=%d)", 
                       length(individual_plots), max_columns))
    
    # Create summary statistics
    summary_data <- data.frame(
      factor_name = factor_names,
      plot_included = TRUE,
      stringsAsFactors = FALSE
    )
    
    # Add metadata from individual results if available
    for (i in seq_along(factor_names)) {
      factor_name <- factor_names[i]
      result <- single_factor_results[[factor_name]]
      
      if (!is.null(result$metadata)) {
        summary_data$n_categories[i] <- result$metadata$n_categories %||% NA
        summary_data$total_observations[i] <- result$metadata$total_observations %||% NA
        summary_data$statistical_test[i] <- result$metadata$statistical_test %||% "Unknown"
        summary_data$plot_type[i] <- result$metadata$plot_type %||% "Unknown"
      }
    }
    
    # Create comprehensive metadata
    composite_metadata <- list(
      plot_count = length(individual_plots),
      factor_names = factor_names,
      max_columns_used = max_columns,
      total_unique_observations = total_observations,
      layout_type = "smart_grid",
      composite_title = composite_title,
      created_timestamp = Sys.time()
    )
    
    # Return comprehensive results
    return(list(
      composite_plot = composite_plot,
      individual_plots = individual_plots,
      metadata = composite_metadata,
      summary = summary_data
    ))
    
  }, error = function(e) {
    log_message(sprintf("Failed to generate single-factor composite plot: %s", e$message), level = "error")
    
    # Return error plot
    error_plot <-ggplot2::ggplot() + 
ggplot2::ggtitle("Single-Factor Composite Plot Generation Failed") + 
ggplot2::labs(
        subtitle = paste("Error:", e$message),
        caption = "Check individual single-factor analysis results and configuration"
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
        factor_names = character(0),
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