#' Safely load and merge group information
#' 
#' @title Safely load and merge group information
#' @description Loads group_info.csv and merges it with main data using first column principle
#' @param main_data Main dataset to merge with
#' @param session_paths Session paths object containing raw data directory
#' @return Merged dataset with group information
#' @export
merge_with_group_info <- function(main_data, session_paths) {
  group_info_path <- file.path(session_paths$raw, "group_info.csv")
  if (!file.exists(group_info_path)) {
    log_message("Group info file not found. Skipping merge.", level = "warning")
    return(main_data)
  }
  tryCatch({
    group_data <- read.csv(group_info_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(group_data) < 1) {
      log_message("group_info.csv is empty. Merge failed.", level = "error")
      return(main_data)
    }
    original_species_col_name <- names(group_data)[1]
    log_message(sprintf("Assuming first column '%s' as species identifier.", original_species_col_name))
    names(group_data)[1] <- "species"

    main_data$species <- as.character(main_data$species)
    group_data$species <- as.character(group_data$species)

  # Merge data
  merged_data <- dplyr::left_join(main_data, group_data, by = "species")
  
  # Do not filter NA values globally here
  # Architecture principle: NA filtering should be per-analysis, not global
  # Reason: A species with Phylogeny=NA should still be included in region_type analysis
  # Each analysis module (single/dual factor) will filter NA for its specific factors
  
  # Log group columns for transparency
  group_columns <- setdiff(names(merged_data), names(main_data))
  if (length(group_columns) > 0) {
    log_message(sprintf("Merged group columns: %s", paste(group_columns, collapse = ", ")))
    
    # Log NA counts per group column (for user awareness)
    for (col in group_columns) {
      na_count <- sum(is.na(merged_data[[col]]))
      if (na_count > 0) {
        log_message(sprintf("  WARNING: %d rows have NA values in '%s' (will be excluded from %s-specific analyses)", 
                           na_count, col, col))
      }
    }
  }
  
  return(merged_data)
    
  }, error = function(e) {
    log_message(sprintf("Critical error in merge_with_group_info: %s", e$message), level = "error")
    return(main_data) # Return original data on error
  })
}