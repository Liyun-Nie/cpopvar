#' Generate M01 data distribution plots
#' 
#' @title Generate M01 data distribution plots
#' @description Main orchestrator for generating data distribution visualizations in M01 module
#' @param pie_data Dataset for pie chart matrix generation
#' @param freq_data Dataset for frequency-based plots
#' @param config Configuration object
#' @param task_params Task-specific parameters
#' @param output_dir Output directory for plots
#' @param task_name Name of the current task
#' @return List containing plot results and metadata
#' @export
generate_M01_data_distribution_plots <- function(pie_data,
                                                freq_data,
                                                config,
                                                task_params,
                                                output_dir = NULL,
                                                task_name = NULL) {
  
  # Dependencies automatically loaded in R package context
  
  log_message("=== M01 DISTRIBUTION PLOTTING ===")
  log_message(sprintf("Generating M01 distribution plots for task: %s", task_name %||% "unnamed"))
  
  # Validate input data structures
  if (is.null(pie_data) || !is.data.frame(pie_data) || nrow(pie_data) == 0) {
    stop("ERROR: pie_data must be a non-empty data frame")
  }
  
  if (is.null(freq_data) || !is.data.frame(freq_data) || nrow(freq_data) == 0) {
    stop("ERROR: freq_data must be a non-empty data frame")
  }
  
  log_message(sprintf("Received cleaned data - pie_data: %d rows, freq_data: %d rows", 
                      nrow(pie_data), nrow(freq_data)))
  
  # V3.16.14: Unified directory structure - use M01_distribution_standard for all outputs
  if (!is.null(task_name) && !is.null(output_dir)) {
    task_output_dir <- file.path(output_dir, "M01_distribution", "M01_distribution_standard")
    dir.create(task_output_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    task_output_dir <- output_dir %||% "plots"
  }
  
  # Extract task-specific parameters
  plot_types <- get_task_parameter(task_params, config, "plot_types", 
                                   default = c("region_density", "region_boxplot_composite", "histogram", "enhanced_pie_matrix"))
  comprehensive_overview <- get_task_parameter(task_params, config, "comprehensive_overview", default = FALSE)
  outlier_detection <- get_task_parameter(task_params, config, "outlier_detection", default = TRUE)
  color_palette <- get_task_parameter(task_params, config, "color_palette", default = "default")
  
  results <- list()
  plots_generated <- character(0)
  
  # Direct routing with cleaned data
  plot_result <- NULL
  
  # Generate plots based on specified types
  for (plot_type in plot_types) {
    log_message(sprintf("Generating %s plot for task %s", plot_type, task_name))
    
    # Direct data routing - pie_matrix gets pie_data, all others get freq_data
    plot_data <- NULL
    
    if (plot_type == "enhanced_pie_matrix") {
        plot_data <- pie_data
    } else {
        # histogram, region_density, and region_boxplot_composite all use freq_data
        plot_data <- freq_data
    }
    
    log_message(sprintf("Using %d rows of cleaned data for %s plot", nrow(plot_data), plot_type))
    
    if (plot_type == "histogram") {
      plot_result <- generate_distribution_histogram(
        plot_data,
        config = config,
        task_params = task_params
      )
      
    } else if (plot_type == "region_density") {
      # Generate region-specific density distribution plot
      plot_result <- generate_region_density_plot(
        plot_data,
        config = config,
        task_params = task_params
      )
      
    } else if (plot_type == "region_boxplot_composite") {
      # Generate composite boxplot + scatter plot
      log_message(sprintf("Data for boxplot: %d rows, %d cols", nrow(plot_data), ncol(plot_data)))
      plot_result <- generate_region_boxplot_composite(
        plot_data,
        config = config,
        task_params = task_params
      )
      
    } else if (plot_type == "enhanced_pie_matrix") {
      # Generate enhanced pie chart matrix
      log_message(sprintf("Data for pie_matrix: %d rows, %d cols", nrow(plot_data), ncol(plot_data)))
      plot_result <- generate_enhanced_pie_matrix(
        plot_data,
        config = config,
        task_params = task_params
      )
      
    } else {
      warning(sprintf("Unknown plot type: %s", plot_type))
      next
    }
    
    if (!is.null(plot_result)) {
      # Save plot with task-specific naming
      plot_filename <- sprintf("%s_%s.png", task_name %||% "M01_distribution", plot_type)
      plot_path <- file.path(task_output_dir, plot_filename)
      
      # Use appropriate dimensions for different plot types
      if (plot_type == "boxplot" && !is.null(plot_result$metadata$panels_generated) && plot_result$metadata$panels_generated == 4) {
        # 4-panel distribution overview
        plot_base_path <- tools::file_path_sans_ext(plot_path)
        saved_path <- save_plot(plot_result$plot, plot_base_path, config, width = 16, height = 12)
      } else if (plot_type == "region_boxplot_composite") {
        # 1×3 composite plot
        plot_base_path <- tools::file_path_sans_ext(plot_path)
        saved_path <- save_plot(plot_result$plot, plot_base_path, config, width = 18, height = 6)  
      } else if (plot_type == "enhanced_pie_matrix") {
        # Large pie matrix
        plot_base_path <- tools::file_path_sans_ext(plot_path)
        saved_path <- save_plot(plot_result$plot, plot_base_path, config, width = 20, height = 12)
      } else {
        # Standard size for other plots
        plot_base_path <- tools::file_path_sans_ext(plot_path)
        saved_path <- save_plot(plot_result$plot, plot_base_path, config)
      }
      
      results[[plot_type]] <- list(
        plot = plot_result$plot,
        data = plot_result$data,
        file_path = saved_path %||% plot_path,
        metadata = plot_result$metadata
      )
      plots_generated <- c(plots_generated, plot_type)
      
      log_message(sprintf("Saved %s plot: %s", plot_type, plot_path))
    }
  }
  
  # Generate task summary
  task_summary <- list(
    task_name = task_name,
    module = "M01_distribution",
    plots_generated = plots_generated,
    output_directory = task_output_dir,
    data_rows_pie = nrow(pie_data),
    data_rows_freq = nrow(freq_data),
    parameters_used = list(
      plot_types = plot_types,
      comprehensive_overview = comprehensive_overview,
      outlier_detection = outlier_detection,
      color_palette = color_palette
    )
  )
  
  log_message(sprintf("M01 distribution analysis completed: %d plots generated", length(plots_generated)))
  
  return(list(
    plots = results,
    summary = task_summary,
    success = length(plots_generated) > 0
  ))
}
