############################################################
#### apply_species_ordering.R - Species Ordering Utility ####
############################################################
#
# Provides a pure function to apply a specific ordering to
# a species column in a dataframe. Refactored from visualization_themes.R
# to be a compliant, stateless utility.
#
############################################################

#' Reorder a species factor based on a provided order vector (Refactored for Compliance)
#' Rule #6 & #9 Compliant: No longer reads files or uses config.
#' @param data Data frame containing species column
#' @param species_order_vector A character vector defining the desired order of species.
#' @param species_col The name of the species column in the data frame.
#' @return Data frame with the species column as an ordered factor.
#' apply_species_ordering
#' @export
apply_species_ordering <- function(data, species_order_vector, species_col = "species") {
  if (is.null(species_order_vector) || !species_col %in% names(data)) {
    return(data)
  }

  data_species <- unique(as.character(data[[species_col]]))
  missing_species <- setdiff(data_species, species_order_vector)
  final_order <- c(intersect(species_order_vector, data_species), missing_species)

  data[[species_col]] <- factor(data[[species_col]], levels = final_order)

  return(data)
}
