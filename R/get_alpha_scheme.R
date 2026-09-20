#' Get alpha transparency scheme
#' 
#' @title Get alpha transparency scheme
#' @description Returns predefined alpha transparency values for plot elements
#' @param scheme_name Name of the alpha scheme (default: "standard")
#' @return List of alpha values for different plot elements
#' @export
get_alpha_scheme <- function(scheme_name = "standard") {
  
  # Define the standard alpha values for plot elements
  standard_scheme <- list(
    # Main fill for larger geoms like boxplots, violins, densities
    fill_main = 0.7,
    
    # Alpha for primary data points in scatter/jitter plots
    point_main = 0.4,
    
    # Alpha for highlighted outlier points to make them stand out
    point_outlier = 0.9,
    
    # Fill for histograms
    histogram_fill = 0.7
  )
  
  # Return the selected scheme (currently only one is defined)
  return(standard_scheme)
}
