#' Get ordered factor levels
#' 
#' @title Get ordered factor levels
#' @description Determines factor level ordering based on preferences or data characteristics
#' @param factor_values Vector of factor values to order
#' @param factor_name Name of the factor being processed
#' @param preferred_order Preferred ordering of factor levels (optional)
#' @return Ordered vector of factor levels
#' @export
get_factor_levels <- function(factor_values, factor_name, preferred_order = NULL) {
  # Pure utility function - uses dependency injection for ordering (Rule #11)
  if (!is.null(preferred_order)) {
    ordered_values <- intersect(preferred_order, factor_values)
    remaining <- setdiff(factor_values, ordered_values)
    return(c(ordered_values, sort(remaining)))
  } else {
    # Fallback to alphabetical ordering if no preferred order provided
    return(sort(factor_values))
  }
}
