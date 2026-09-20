#' Get task parameter (alias function)
#' @param task_params Task parameters
#' @param config Configuration object
#' @param param_name Parameter name
#' @param default Default value if not found
#' @param log_default Whether to log when using default
#' @return Parameter value
#' @export
get_task_param <- function(task_params, config, param_name, default = NULL, log_default = FALSE) {
  return(get_task_parameter(task_params, config, param_name, default, log_default))
}
