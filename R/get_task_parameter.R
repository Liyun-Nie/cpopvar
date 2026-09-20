#' Get task parameter value with hierarchical lookup
#' @param task_params Task parameters
#' @param config Configuration object
#' @param param_name Parameter name
#' @param default Default value if not found
#' @param log_default Whether to log when using default
#' @return Parameter value
#' @export
get_task_parameter <- function(task_params, config, param_name, default = NULL, log_default = FALSE) {
  found_value <- NULL
  found_source <- NULL
  
  # First priority: task-specific parameters (direct level)
  if (!is.null(task_params) && param_name %in% names(task_params)) {
    found_value <- task_params[[param_name]]
    found_source <- "task_params (direct)"
  }
  
  # Second priority: task-specific parameters (nested under 'parameters')
  else if (!is.null(task_params) && "parameters" %in% names(task_params) && 
           !is.null(task_params$parameters) && param_name %in% names(task_params$parameters)) {
    found_value <- task_params$parameters[[param_name]]
    found_source <- "task_params (nested)"
  }
  
  # Third priority: config-specific parameters
  else if (!is.null(config) && !is.null(config$visualization) && param_name %in% names(config$visualization)) {
    found_value <- config$visualization[[param_name]]
    found_source <- "config$visualization"
  }
  
  # Fallback to default with optional warning
  else {
    found_value <- default
    if (log_default && !is.null(default)) {
      log_message(
        sprintf("WARNING: Parameter '%s' not found in config. Using default value: %s", 
                param_name, 
                paste(deparse(default), collapse = "")),
        level = "WARN"
      )
    }
  }
  
  return(found_value)
}
