#' @importFrom stats as.dist binom.test density end hclust median sd quantile mad pnorm pchisq phyper start TukeyHSD aov kruskal.test shapiro.test as.formula cor cor.test dbinom fisher.test formula frequency p.adjust pairwise.wilcox.test prcomp reorder runif setNames var wilcox.test
#' @importFrom grDevices col2rgb colorRampPalette rainbow rgb
#' @importFrom utils capture.output head installed.packages read.csv read.table str tail write.csv write.table
NULL

#' Run cpopvar Analysis Pipeline
#'
#' Executes the complete back-end analysis pipeline using a configuration file.
#'
#' @param config_file Path to the YAML configuration file.
#' @param session_id Optional session identifier for tracking.
#' @param output_dir Working directory for output. Defaults to current directory.
#' @return Invisibly returns the path to the session output directory.
#' @export
#' @importFrom yaml read_yaml write_yaml
run_analysis <- function(config_file, session_id = NULL, output_dir = getwd()) {

  # --- Step 1: Input validation and configuration loading ---
  if (is.null(config_file) || !file.exists(config_file)) {
    stop(paste("Configuration file not found:", config_file), call. = FALSE)
  }

  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  previous_output_dir <- getOption("cpopvar.output_dir", NULL)
  options(cpopvar.output_dir = output_dir)
  on.exit(options(cpopvar.output_dir = previous_output_dir), add = TRUE)

  # Note: Using null coalescing operator (%||%) exported from annotation_processing.R

  # Load the raw config data to get the session_id
  config_data <- yaml::read_yaml(config_file)

  # V3.13: Inject config file path into config object for downstream access
  config_file_absolute <- normalizePath(config_file, mustWork = TRUE)
  if (is.null(config_data$session_info)) {
    config_data$session_info <- list()
  }
  config_data$session_info$config_file_path <- config_file_absolute

  # Determine session_id: function parameter > config file > auto-generate
  session_id <- session_id %||% config_data$session_info$session_id %||% paste0("session_", format(Sys.time(), "%Y%m%d_%H%M%S"))

  # --- Step 2: Session initialization and logging setup ---
  
  # Initialize session directories and logging with dynamic output directory
  session_paths <- get_session_paths(session_id, create_dirs = TRUE, output_dir = output_dir)
  initialize_session_logging(session_id, output_dir = output_dir)

  log_message("=== cpopvar R Package Analysis Run Initialized ===", session_id = session_id)
  log_message(paste("Config file:", config_file), session_id = session_id)

  # --- Step 3: Task configuration parsing ---

  log_message("Parsing task configuration...", session_id = session_id)
  # V3.14: Pass config_data object to preserve injected config_file_path
  # Note: Removed override_mode since in package mode, this should be fully controlled by config file
  parsed_task_config <- parse_task_configuration(config_data = config_data)

  # Override session_id in the parsed config if it was provided via function parameter
  if (!is.null(session_id)) {
      parsed_task_config$config_data$session_info$session_id <- session_id
  }

  # Display the execution plan to the user
  log_message("--- Execution Plan Summary ---", session_id = session_id)
  log_message(get_execution_plan_summary(parsed_task_config$execution_plan), session_id = session_id)
  log_message("----------------------------", session_id = session_id)

  # --- Step 4: Core analysis execution ---

  log_message("Dispatching analysis tasks...", session_id = session_id)
  # Remove progress_callback or adapt to be more suitable for R package logging
  analysis_results <- dispatch_analysis_tasks(
    parsed_config = parsed_task_config,
    session_id = session_id,
    output_dir = output_dir
  )

  # --- Step 5: Results summary and export ---

  log_message("Analysis complete. Saving summary...", session_id = session_id)
  summary_data <- list(
    session_id = session_id,
    config_file = config_file,
    execution_mode = "full_analysis",  # V3.16: Simplified configuration
    tasks_executed = names(analysis_results$task_results),
    end_time = Sys.time(),
    session_paths = session_paths
  )

  summary_file <- file.path(session_paths$results, "run_summary.rds")
  saveRDS(summary_data, summary_file)

  # Save the exact configuration that was used for this run
  config_used_file <- file.path(session_paths$results, "config_used.yml")
  yaml::write_yaml(parsed_task_config$config_data, config_used_file)

  log_message(paste("All results are available in:", session_paths$base), session_id = session_id)
  log_message("=== cpopvar R Package Analysis Finished ===", session_id = session_id)

  # --- Step 6: Return session path ---
  invisible(session_paths$base)
}