############################################################
#### data_utils.R - Data Processing Utilities ####
############################################################
#
# Provides common utilities for data validation, cleaning,
# and transformation. Migrated from src/utils and made compliant.
#
############################################################

# Load required libraries

#' Validate data frame structure
#' Validate data frame structure
#' @param data Data frame to validate
#' @param required_columns Required column names
#' @param min_rows Minimum number of rows
#' @param data_name Name of data for error messages
#' @return Logical indicating if validation passed
#' @importFrom rlang sym
#' @export
validate_data_frame <- function(data, required_columns = NULL, min_rows = 1, data_name = "data") {
  if (!is.data.frame(data)) return(FALSE)
  if (nrow(data) < min_rows) return(FALSE)
  if (!is.null(required_columns)) {
    if (length(setdiff(required_columns, names(data))) > 0) return(FALSE)
  }
  return(TRUE)
}

#' Clean and standardize species names
#' @param species_names Character vector of species names to clean
clean_species_names <- function(species_names) {
  species_names <- trimws(species_names)
  species_names <- gsub(" ", "_", species_names)
  species_names <- gsub("[^A-Za-z0-9_]", "", species_names)
  return(species_names)
}

#' Standardize variant type names
#' @param var_types Character vector of variant types to standardize
standardize_variant_types <- function(var_types) {
  var_types <- tolower(var_types)
  dplyr::case_when(
    var_types %in% c("snp", "substitution") ~ "snp",
    var_types %in% c("del", "deletion") ~ "del",
    var_types %in% c("ins", "insertion") ~ "ins",
    var_types %in% c("indel", "del", "ins") ~ "indel",
    var_types %in% c("complex", "complex_indel") ~ "complex",
    var_types %in% c("mnp", "mnv") ~ "mnp",
    TRUE ~ var_types
  )
}

#' Standardize species names for consistent matching
#' @param species_names Character vector of species names to standardize
standardize_species_names <- function(species_names) {
  species_names %>%
    trimws() %>%
    gsub("\\s+", "_", .) %>%
    gsub("_+$", "_", .) %>%
    gsub("^[_]|[_]$", "", .)
}

#' Calculate summary statistics for a numeric vector
#' 
#' @title Calculate summary statistics for a numeric vector
#' @description Calculates basic descriptive statistics for a numeric vector
#' @param x Numeric vector to analyze
#' @param na.rm Logical, whether to remove NA values (default: TRUE)
#' @return List containing summary statistics
#' @export
calculate_summary_stats <- function(x, na.rm = TRUE) {
  if (!is.numeric(x)) return(NULL)
  list(
    n = length(x),
    n_missing = sum(is.na(x)),
    mean = mean(x, na.rm = na.rm),
    median = median(x, na.rm = na.rm),
    sd = sd(x, na.rm = na.rm),
    min = min(x, na.rm = na.rm),
    max = max(x, na.rm = na.rm)
  )
}

#' Convert data to long format for analysis
#' Pivot data to long format
#' @param data Data frame to pivot
#' @param id_cols ID columns
#' @param measure_cols Measure columns
#' @param measure_name Name for measure column
#' @param value_name Name for value column
#' @return Long format data frame
#' @importFrom tidyr pivot_longer
#' @export
pivot_to_long <- function(data, id_cols, measure_cols, measure_name = "variable", value_name = "value") {
tidyr::pivot_longer(data, cols =dplyr::all_of(measure_cols), names_to = measure_name, values_to = value_name)
}

#' Arrange genes by functional category (Refactored for Compliance)
#' Rule #6 & #9 Compliant: No longer reads files or uses config.
#' Takes a pre-loaded mapping data frame.
#' Apply gene functional ordering
#' @param data Data frame with gene data
#' @param gene_function_map Gene function mapping
#' @param gene_col Gene column name
#' @return Data frame with functional ordering applied
#' @importFrom magrittr %>%
#' @export
apply_gene_functional_ordering <- function(data, gene_function_map, gene_col = "gene") {
    if (is.null(gene_function_map) || !"gene" %in% names(gene_function_map) || !"gene_category" %in% names(gene_function_map)) {
        return(data) # Return original data if map is invalid
    }
    # Dynamically determine category order from the map itself
    category_order <- unique(gene_function_map$gene_category)

    data_with_function <- merge(data, gene_function_map, by = gene_col, all.x = TRUE)
    
    data_with_function$gene_category <- factor(data_with_function$gene_category, levels = category_order)
    
    data_sorted <- data_with_function %>%
dplyr::arrange(gene_category, functional_group, !!sym(gene_col)) %>%
dplyr::select(-dplyr::any_of(c("gene_category", "functional_group")))
    return(data_sorted)
}
