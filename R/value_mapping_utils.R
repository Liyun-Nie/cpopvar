############################################################
#### value_mapping_utils.R - Universal Value Mapping System ####
############################################################
#
# Universal value mapping utilities for user-defined data transformations
# Supports rename, merge, and delete operations on data frame columns
# Created for configuration-driven frontend/backend value mapping
#
############################################################

#' Apply user-defined value mapping to data column
#' 
#' @title Apply user-defined value mapping to data column
#' @description Applies user-defined mapping rules to transform values in a specific column.
#'              Supports three operations: rename (1:1), merge (N:1), and delete (to NA).
#' @param data Data frame to process
#' @param target_column Character string specifying the column name to transform
#' @param mapping_rules Named list or vector defining the mapping rules.
#'        Format: list("source_value1" = "target_value", "source_value2" = "target_value", "delete_me" = NA)
#'        - Rename: "old_name" = "new_name"  
#'        - Merge: "type1" = "COMBINED", "type2" = "COMBINED"
#'        - Delete: "unwanted" = NA
#' @return Data frame with transformed values in the target column.
#'         Rows with NA values (from delete operations) are removed from the result.
#' @examples
#' # Example 1: Rename and merge variant types
#' data <- data.frame(
#'   species = c("A", "B", "C"),
#'   var_type = c("del", "ins", "snp"),
#'   region_type = c("CDS", "rRNA", "tRNA")
#' )
#' rules <- list("del" = "INDEL", "ins" = "INDEL")  # Merge del and ins to INDEL
#' result <- apply_value_mapping(data, "var_type", rules)
#' 
#' # Example 2: Delete unwanted categories
#' rules <- list("rRNA" = NA, "tRNA" = NA)  # Delete rRNA and tRNA entries
#' result <- apply_value_mapping(data, "region_type", rules)
#' 
#' @export
apply_value_mapping <- function(data, target_column, mapping_rules) {
  
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame")
  }
  
  if (!is.character(target_column) || length(target_column) != 1) {
    stop("Input 'target_column' must be a single character string")
  }
  
  if (!target_column %in% names(data)) {
    stop(paste("Column", target_column, "not found in data frame"))
  }
  
  if (is.null(mapping_rules) || length(mapping_rules) == 0) {
    log_message(paste("No mapping rules provided for", target_column, "- returning original data"))
    return(data)
  }
  
  # Convert to named list if needed
  if (!is.list(mapping_rules)) {
    mapping_rules <- as.list(mapping_rules)
  }
  
  log_message(sprintf("Applying value mapping to column '%s' with %d rules", 
                     target_column, length(mapping_rules)))
  
  # Track transformation statistics
  original_values <- unique(data[[target_column]])
  transformation_log <- list()
  
  # Create a copy to avoid modifying the original data
  result_data <- data
  
  # Convert target column to character to support gsub operations safely
  if (is.factor(result_data[[target_column]])) {
    result_data[[target_column]] <- as.character(result_data[[target_column]])
  }
  
  # Apply each mapping rule
  for (source_value in names(mapping_rules)) {
    target_value <- mapping_rules[[source_value]]
    
    # Skip if source value doesn't exist in data
    if (!source_value %in% result_data[[target_column]]) {
      transformation_log[[source_value]] <- "source_not_found"
      next
    }
    
    # Count occurrences before transformation
    occurrence_count <- sum(result_data[[target_column]] == source_value, na.rm = TRUE)
    
    if (is.na(target_value) || is.null(target_value)) {
      # Delete operation: mark for removal
      result_data[[target_column]][result_data[[target_column]] == source_value] <- NA
      transformation_log[[source_value]] <- sprintf("deleted_%d_rows", occurrence_count)
      log_message(sprintf("  DELETE: %s (%d occurrences) -> removed", source_value, occurrence_count))
      
    } else {
      # Rename/Merge operation
      result_data[[target_column]][result_data[[target_column]] == source_value] <- target_value
      transformation_log[[source_value]] <- sprintf("mapped_to_%s_%d_rows", target_value, occurrence_count)
      
      if (source_value == target_value) {
        log_message(sprintf("  KEEP: %s (%d occurrences) -> unchanged", source_value, occurrence_count))
      } else {
        log_message(sprintf("  MAP: %s (%d occurrences) -> %s", source_value, occurrence_count, target_value))
      }
    }
  }
  
  # Remove rows where target column is NA (result of delete operations)
  initial_rows <- nrow(result_data)
  result_data <- result_data[!is.na(result_data[[target_column]]), ]
  final_rows <- nrow(result_data)
  deleted_rows <- initial_rows - final_rows
  
  if (deleted_rows > 0) {
    log_message(sprintf("Removed %d rows with deleted values from %s column", 
                       deleted_rows, target_column))
  }
  
  # Log summary statistics
  final_values <- unique(result_data[[target_column]])
  log_message(sprintf("Value mapping completed for '%s': %d -> %d unique values, %d -> %d total rows", 
                     target_column, length(original_values), length(final_values), 
                     initial_rows, final_rows))
  
  # Store transformation metadata as attribute (for debugging)
  attr(result_data, "value_mapping_log") <- list(
    target_column = target_column,
    original_values = original_values,
    final_values = final_values,
    transformation_log = transformation_log,
    rows_deleted = deleted_rows
  )
  
  return(result_data)
}

#' Convert UI preprocessing rules to backend mapping format
#' 
#' @title Convert UI preprocessing rules to backend mapping format  
#' @description Converts the YAML-style preprocessing rules from config file
#'              to the named list format expected by apply_value_mapping()
#' @param preprocessing_rules List from config$preprocessing$value_mapping_rules
#' @param rule_type Character string: "variant_types" or "region_types"
#' @return Named list suitable for apply_value_mapping() or NULL if no rules
#' @examples
#' # YAML config structure:
#' # preprocessing:
#' #   value_mapping_rules:
#' #     variant_types:
#' #       - source: del
#' #         target: INDEL
#' #       - source: ins  
#' #         target: INDEL
#' rules <- list(variant_types = list(
#'   list(source = "del", target = "INDEL"),
#'   list(source = "ins", target = "INDEL")
#' ))
#' mapping <- convert_preprocessing_rules_to_mapping(rules, "variant_types")
#' # Result: list("del" = "INDEL", "ins" = "INDEL")
#' @export
convert_preprocessing_rules_to_mapping <- function(preprocessing_rules, rule_type) {
  
  if (is.null(preprocessing_rules) || !is.list(preprocessing_rules)) {
    return(NULL)
  }
  
  rules_for_type <- preprocessing_rules[[rule_type]]
  
  if (is.null(rules_for_type) || length(rules_for_type) == 0) {
    return(NULL)
  }
  
  # Convert list of source/target pairs to named list
  mapping_list <- list()
  
  for (rule in rules_for_type) {
    if (is.list(rule) && "source" %in% names(rule)) {
      source_value <- rule$source
      target_value <- rule$target
      
      # Handle YAML null values (represented as NULL in R)
      if (is.null(target_value)) {
        target_value <- NA
      }
      
      mapping_list[[source_value]] <- target_value
    }
  }
  
  if (length(mapping_list) == 0) {
    return(NULL)
  }
  
  log_message(sprintf("Converted %d preprocessing rules for %s: %s", 
                     length(mapping_list), rule_type,
                     paste(names(mapping_list), "->", unlist(mapping_list), collapse = ", ")))
  
  return(mapping_list)
}

#' Extract mapping rules from current UI state
#' 
#' @title Extract mapping rules from current UI state
#' @description Helper function to collect user-defined mapping rules from Shiny UI inputs.
#'              Used by build_ui_configuration() to capture preprocessing settings.
#' @param input Shiny input object containing all UI input values
#' @param rv Reactive values object containing ui_choices data
#' @return List with variant_types and region_types mapping rules in YAML format
#' @examples
#' # This function is called internally by build_ui_configuration()
#' # Returns structure like:
#' # list(
#' #   variant_types = list(
#' #     list(source = "del", target = "INDEL"),
#' #     list(source = "ins", target = "INDEL")
#' #   ),
#' #   region_types = list(
#' #     list(source = "intergenic", target = "IGS")
#' #   )
#' # )
collect_ui_preprocessing_rules <- function(input, rv) {
  
  tryCatch({
    # Initialize result structure
    preprocessing_rules <- list(
      variant_types = list(),
      region_types = list()
    )
    
    # Collect variant type mapping rules
    variant_types <- rv$ui_choices$variant_types
    if (!is.null(variant_types) && length(variant_types) > 0) {
      
      for (source_type in variant_types) {
        input_id <- paste0("variant_target_", gsub("[^A-Za-z0-9]", "_", source_type))
        target_value <- input[[input_id]]
        
        # Only include if user has modified the value (target != source)
        if (!is.null(target_value) && !is.na(target_value)) {
          target_value <- trimws(target_value)
          
          # Check if this is a meaningful change
          if (target_value != source_type) {
            rule <- list(source = source_type)
            
            # Handle empty target (delete operation)
            if (target_value == "" || target_value == "DELETE") {
              rule$target <- NULL  # Will become YAML null
            } else {
              rule$target <- target_value
            }
            
            preprocessing_rules$variant_types <- append(preprocessing_rules$variant_types, list(rule))
          }
        }
      }
    }
    
    # Collect region type mapping rules  
    region_types <- rv$ui_choices$region_types
    if (!is.null(region_types) && length(region_types) > 0) {
      
      for (source_type in region_types) {
        input_id <- paste0("region_target_", gsub("[^A-Za-z0-9]", "_", source_type))
        target_value <- input[[input_id]]
        
        # Only include if user has modified the value (target != source)
        if (!is.null(target_value) && !is.na(target_value)) {
          target_value <- trimws(target_value)
          
          # Check if this is a meaningful change
          if (target_value != source_type) {
            rule <- list(source = source_type)
            
            # Handle empty target (delete operation)
            if (target_value == "" || target_value == "DELETE") {
              rule$target <- NULL  # Will become YAML null
            } else {
              rule$target <- target_value
            }
            
            preprocessing_rules$region_types <- append(preprocessing_rules$region_types, list(rule))
          }
        }
      }
    }
    
    # Return NULL if no rules were collected
    if (length(preprocessing_rules$variant_types) == 0 && 
        length(preprocessing_rules$region_types) == 0) {
      log_message("No preprocessing rules found in UI state")
      return(NULL)
    }
    
    log_message(sprintf("Collected %d variant type rules and %d region type rules from UI",
                       length(preprocessing_rules$variant_types),
                       length(preprocessing_rules$region_types)))
    
    return(preprocessing_rules)
    
  }, error = function(e) {
    log_message(paste("Error collecting UI preprocessing rules:", e$message), level = "error")
    return(NULL)
  })
}

log_message("Universal value mapping utilities loaded successfully")