############################################################
#### manifest_utils.R - Analysis Output Manifest System ####
############################################################
#
# Scientific report generation system upgrade - Manifest architecture
# Provides structured metadata for analysis outputs to enable
# high-quality, academic-grade report generation.
#
# Key Features:
# - Standardized manifest data structure
# - Role-based content classification
# - Academic title and caption support
# - Path conversion utilities
#
############################################################

# Load required libraries

# Note: Using null coalescing operator (%||%) exported from annotation_processing.R

#' Create a standardized output manifest entry
#' 
#' @param path Absolute path to the output file
#' @param role Role of the file (composite_plot, individual_plot, summary_table, raw_data, etc.)
#' @param title Academic title for the output (e.g., "Figure 1: Distribution Analysis")
#' @param caption Detailed caption/description of the output
#' @param module Source module that generated this output (M01, M02, M03, etc.)
#' @param file_type Type of file (plot, table, data, report)
#' @param priority Display priority (1=highest, used for report ordering)
#' @param parameters Analysis parameters used for this output (V7.1 Enhancement)
#' @return Standardized manifest entry (list)
#' create_manifest_entry
#' @export
create_manifest_entry <- function(path, role, title, caption, 
                                 module = NULL, file_type = NULL, priority = 5, parameters = NULL) {
  
  # Validate required parameters
  if (missing(path) || missing(role) || missing(title)) {
    stop("path, role, and title are required parameters")
  }
  
  # Auto-detect file type if not provided
  if (is.null(file_type)) {
    file_ext <- tools::file_ext(path)
    file_type <- dplyr::case_when(
      tolower(file_ext) %in% c("png", "pdf", "jpg", "jpeg") ~ "plot",
      tolower(file_ext) %in% c("csv", "tsv") ~ "table", 
      tolower(file_ext) %in% c("txt", "md") ~ "report",
      tolower(file_ext) %in% c("rds", "rdata") ~ "data",
      TRUE ~ "unknown"
    )
  }
  
  # Auto-detect module if not provided (from path)
  if (is.null(module)) {
    if (grepl("M01", path, ignore.case = TRUE)) {
      module <- "M01"
    } else if (grepl("M02", path, ignore.case = TRUE)) {
      module <- "M02"
    } else if (grepl("M03", path, ignore.case = TRUE)) {
      module <- "M03"
    } else {
      module <- "Unknown"
    }
  }
  
  # Create standardized entry
  entry <- list(
    path = path,
    role = role,
    title = title,
    caption = caption %||% "",
    module = module,
    file_type = file_type,
    priority = as.integer(priority),
    parameters = parameters,  # V7.1 Enhancement: Store analysis parameters
    file_size = ifelse(file.exists(path), file.size(path), 0),
    created_time = Sys.time()
  )
  
  class(entry) <- c("manifest_entry", "list")
  return(entry)
}

#' Create a manifest collection for a module
#' 
#' @param entries List of manifest entries
#' @param module_name Name of the module (M01, M02, M03)
#' @param session_id Session identifier
#' @return Manifest collection (list)
create_manifest_collection <- function(entries, module_name, session_id) {
  
  # Validate entries
  if (!is.list(entries) || length(entries) == 0) {
    warning(sprintf("Empty or invalid entries for module %s", module_name))
    entries <- list()
  }
  
  # Sort entries by priority and role
  if (length(entries) > 0) {
    priorities <- sapply(entries, function(x) x$priority %||% 5)
    entries <- entries[order(priorities)]
  }
  
  collection <- list(
    module = module_name,
    session_id = session_id,
    entries = entries,
    count = length(entries),
    created_time = Sys.time()
  )
  
  class(collection) <- c("manifest_collection", "list")
  return(collection)
}

#' Define standard role types for output classification
#' 
#' @return Named list of role definitions with descriptions
get_standard_roles <- function() {
  list(
    # High-priority roles (for final report)
    composite_plot = "Composite figure combining multiple analyses",
    summary_table = "Summary table with key results",
    key_plot = "Key individual plot of primary importance",
    
    # Medium-priority roles  
    individual_plot = "Individual analysis plot",
    detailed_table = "Detailed data table",
    statistical_report = "Statistical analysis report",
    
    # Low-priority roles (typically excluded from final report)
    raw_data = "Raw processed data file",
    intermediate_data = "Intermediate processing data",
    debug_data = "Debug or quality control data",
    process_log = "Processing log file"
  )
}

#' Define standard academic titles for each module and role
#' 
#' @param module Module name (M01, M02, M03)
#' @param role Output role
#' @param file_name Optional file name for context
#' @return Standard academic title
get_standard_title <- function(module, role, file_name = NULL) {
  
  # Define title templates
  title_templates <- list(
    M01 = list(
      composite_plot = "Figure %d: Data Distribution Analysis",
      individual_plot = "Figure %d: %s Distribution",
      summary_table = "Table %d: Distribution Summary Statistics"
    ),
    M02 = list(
      composite_plot = "Figure %d: Comparative Analysis Results", 
      individual_plot = "Figure %d: %s Comparison",
      statistical_report = "Supplementary File %d: Statistical Analysis Report",
      summary_table = "Table %d: Comparative Analysis Summary"
    ),
    M03 = list(
      composite_plot = "Figure %d: Hotspot Analysis Overview",
      individual_plot = "Figure %d: %s Hotspot Analysis", 
      summary_table = "Table %d: Hotspot Gene Summary",
      statistical_report = "Supplementary File %d: Hotspot Statistical Tests"
    )
  )
  
  # Get template
  template <- title_templates[[module]][[role]]
  if (is.null(template)) {
    template <- "Output %d: %s"
  }
  
  # Generate title (numbering will be handled at report generation)
  if (grepl("%s", template) && !is.null(file_name)) {
    # Extract meaningful name from file name
    clean_name <- gsub("^.*_", "", tools::file_path_sans_ext(file_name))
    clean_name <- gsub("_", " ", clean_name)
    clean_name <- tools::toTitleCase(clean_name)
    return(sprintf(template, 1, clean_name))
  } else {
    return(sprintf(template, 1))
  }
}

#' Generate standard captions for common plot types
#' 
#' @param module Module name
#' @param role Output role
#' @param plot_type Specific plot type (optional)
#' @return Standard caption text
get_standard_caption <- function(module, role, plot_type = NULL) {
  
  captions <- list(
    M01 = list(
      composite_plot = "Comprehensive data distribution analysis showing variant frequency patterns across different genomic regions and species. Each panel represents a different analytical perspective of the same underlying data.",
      individual_plot = "Distribution analysis focusing on specific genomic regions or variant types, providing detailed insights into frequency patterns.",
      summary_table = "Summary statistics of variant frequency distributions including mean, median, and variance measures across all analyzed categories."
    ),
    M02 = list(
      composite_plot = "Comparative analysis results showing statistical differences between groups. Statistical significance is indicated where applicable using standard notation (* p<0.05, ** p<0.01, *** p<0.001).",
      individual_plot = "Individual factor comparison showing group differences with appropriate statistical tests applied based on data distribution properties.",
      statistical_report = "Comprehensive statistical analysis report including normality tests, variance homogeneity tests, and appropriate post-hoc comparisons."
    ),
    M03 = list(
      composite_plot = "Hotspot analysis overview combining multiple analytical approaches to identify and characterize variant hotspot regions across species.",
      individual_plot = "Specific aspect of hotspot analysis providing detailed insights into hotspot gene characteristics and distributions.",
      statistical_report = "Statistical validation of hotspot sharing patterns including tests for randomness and evolutionary significance."
    )
  )
  
  caption <- captions[[module]][[role]]
  return(caption %||% "Analysis output generated by the cpopvar pipeline.")
}

# log_message("Manifest utilities loaded successfully") # Will be available when loaded in proper context