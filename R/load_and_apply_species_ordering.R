############################################################
#### load_and_apply_species_ordering.R - Unified Species Ordering ####
############################################################
#
# Centralized species ordering logic for module orchestrators
# Eliminates code duplication between M02 and M03 modules
#
# ARCHITECTURE: Pure utility function following Rule #13 compliance
# - No config dependencies, uses dependency injection
# - Unified error handling and logging
# - Consistent behavior across all modules
#
############################################################

#' Load and Apply Species Ordering from Session Data
#' 
#' This function centralizes the logic for loading user-provided species ordering
#' from session files and applying it to data. It eliminates code duplication
#' between M02 and M03 module orchestrators.
#' 
#' @param data Data frame to apply species ordering to
#' @param session_id Session ID to locate species ordering file
#' @param species_col Name of the species column (default: "species")
#' @return List containing:
#'   - data: Data frame with species ordering applied (if successful)
#'   - species_order: Character vector of species order used (or NULL)
#'   - custom_label_mapping: Named vector mapping species to custom labels (or NULL)
#'   - success: Logical indicating if ordering was applied
#'   - message: Status message for logging
#' load_and_apply_species_ordering
#' @export
load_and_apply_species_ordering <- function(data, session_id, species_col = "species") {
  
  # Input validation
  if (is.null(data) || nrow(data) == 0) {
    return(list(
      data = data,
      species_order = NULL,
      custom_label_mapping = NULL,
      success = FALSE,
      message = "Input data is empty or NULL"
    ))
  }
  
  if (is.null(session_id)) {
    return(list(
      data = data,
      species_order = NULL,
      custom_label_mapping = NULL,
      success = FALSE,
      message = "Session ID is NULL"
    ))
  }
  
  # Check if species column exists
  if (!species_col %in% names(data)) {
    return(list(
      data = data,
      species_order = NULL,
      custom_label_mapping = NULL,
      success = FALSE,
      message = sprintf("Species column '%s' not found in data", species_col)
    ))
  }
  
  tryCatch({
    # Get session paths for species ordering file
    session_paths <- get_session_paths(session_id)
    species_order_file <- file.path(session_paths$raw, "species_label_order.csv")
    
    # Check if species ordering file exists
    if (!file.exists(species_order_file)) {
      log_message("Species ordering file not found, using default alphabetical order", level = "warning")
      return(list(
        data = data,
        species_order = NULL,
        custom_label_mapping = NULL,
        success = FALSE,
        message = sprintf("Species ordering file not found: %s", species_order_file)
      ))
    }
    
    # V19 FIX: Make file reading more robust by explicitly setting header=TRUE
    # and using check.names=FALSE to handle potentially non-standard column names.
    log_message(sprintf("Loading species ordering from: %s", species_order_file))
    species_order_data <- read.csv(species_order_file, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
    
    # V19 FIX: Unconditionally use the first column as the species order vector.
    species_order <- species_order_data[[1]]
    log_message(sprintf("Loaded species ordering using column '%s': %d species", 
                        names(species_order_data)[1], 
                        length(species_order)))
    
    # V20 ENHANCEMENT: Check for custom_label column
    custom_label_mapping <- NULL
    if ("custom_label" %in% names(species_order_data)) {
      custom_label_mapping <- setNames(species_order_data$custom_label, species_order)
      log_message(sprintf("Found custom_label column with %d labels", length(custom_label_mapping)))
      
      # Validate custom labels (check for duplicates)
      duplicates <- custom_label_mapping[duplicated(custom_label_mapping)]
      if (length(duplicates) > 0) {
        log_message(sprintf("WARNING: Duplicate custom labels found: %s. This may cause visualization issues.",
                           paste(unique(duplicates), collapse = ", ")), level = "warning")
      }
    } else {
      log_message("No custom_label column found - will use default labeling", level = "info")
    }
    
    # Apply species ordering to data
    ordered_data <- apply_species_ordering(data, 
                                         species_order_vector = species_order, 
                                         species_col = species_col)
    
    return(list(
      data = ordered_data,
      species_order = species_order,
      custom_label_mapping = custom_label_mapping,  # V20: New field
      success = TRUE,
      message = sprintf("Successfully applied species ordering with %d species", length(species_order))
    ))
    
  }, error = function(e) {
    error_message <- sprintf("Failed to load/apply species ordering: %s", e$message)
    log_message(error_message, level = "error")
    
    return(list(
      data = data,
      species_order = NULL,
      custom_label_mapping = NULL,
      success = FALSE,
      message = error_message
    ))
  })
}