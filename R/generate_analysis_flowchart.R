############################################################
#### generate_analysis_flowchart.R - Analysis Pipeline Flowchart Generator ####
############################################################
#
# V7.1 Enhancement: Pure text ASCII art flowchart generator for the entire pipeline
# Creates a clear, readable visual representation of the complete analysis workflow
# using simple text characters for maximum compatibility.
#
# Key Features:
# - ASCII art style flowchart using |, v, +, - characters
# - Configuration-driven content adaptation
# - Compatible with all Markdown renderers
# - No external dependencies required
#
# Note: This is separate from generate_flowchart_text.R which is specifically 
# for M02 statistical engine workflow visualization.
#
############################################################

#' Generate complete analysis pipeline flowchart as pure text ASCII art
#' 
#' @param config Configuration object containing pipeline settings
#' @return Character string containing formatted ASCII flowchart
#' generate_analysis_flowchart
#' @export
generate_analysis_flowchart <- function(config = NULL) {
  
  # Basic flowchart structure - always included
  flowchart_lines <- c(
    "```",
    "                    [ Raw Data Input ]",
    "                           |",
    "                           v",
    "                +--------------------+",
    "                |   Data Preprocessing     |",
    "                | (P01 Preprocessing) |",
    "                +--------------------+",
    "                           |",
    "                           v",
    "                +--------------------+",
    "                |   Variant Filtering      |",
    "                |  (P02 Filtering)   |",
    "                +--------------------+",
    "                           |",
    "                           v",
    "                +--------------------+",
    "                |   Frequency Normalization |",
    "                | (P03 Normalization) |",
    "                +--------------------+",
    "                           |",
    "                           v",
    "                +--------------------+",
    "                |    Visualization Analysis |",
    "                | M01: Data Distribution    |",
    "                | M02: compareanalysis       |",
    "                | M03: Hotspot Analysis    |",
    "                +--------------------+",
    "                           |",
    "                           v",
    "                    [ Final Report Generation ]",
    "```"
  )
  
  # V3.16: Simplified configuration - always full_analysis mode
  # Add execution mode information
  mode_annotation <- c(
    "",
    "**Execution Mode**: Full Analysis Mode (V3.16 Simplified Config)",
    "- Execute complete data processing pipeline",
    "- Complete analysis from raw data to final report"
  )
  
  # Insert annotation after flowchart
  flowchart_lines <- c(flowchart_lines, mode_annotation)
  
  # Add detailed stage descriptions
  stage_details <- c(
    "",
    "### Processing Phase Description",
    "",
    "- **P01 Data Preprocessing**: Region annotation, IR coordinate conversion, data deduplication",
    "- **P02 Variant Filtering**: Filter valid variants based on frequency thresholds",
    "- **P03 Frequency Normalization**: Calculate variant frequency per kilobase",
    "- **M01 Data Distribution**: Frequency distribution pattern visualization analysis",
    "- **M02 Comparative Analysis**: Between-group statistical testing and visualization",
    "- **M03 Hotspot Analysis**: Hotspot gene identification and functional enrichment analysis"
  )
  
  # Combine all parts
  complete_flowchart <- c(flowchart_lines, stage_details)
  
  # Return as single string
  return(paste(complete_flowchart, collapse = "\n"))
}

#' Generate detailed parameter summary for pipeline context
#' 
#' @param config Configuration object
#' @return Character string with parameter summary
generate_pipeline_parameter_summary <- function(config = NULL) {
  
  if (is.null(config)) {
    return("Configuration parameter information not available")
  }
  
  summary_lines <- c("### Key Configuration Parameters")
  
  # V3.16: Simplified configuration - always full_analysis mode
  summary_lines <- c(summary_lines, 
                    "- **Execution Mode**: full_analysis (V3.16 simplified config)")
  
  # Session info
  if (!is.null(config$session_info$session_id)) {
    summary_lines <- c(summary_lines,
                      sprintf("- **Session ID**: %s", config$session_info$session_id))
  }
  
  # Analysis tasks summary
  if (!is.null(config$analysis_tasks)) {
    enabled_tasks <- sum(sapply(config$analysis_tasks, function(t) t$enabled %||% FALSE))
    total_tasks <- length(config$analysis_tasks)
    summary_lines <- c(summary_lines,
                      sprintf("- **Analysis Tasks**: %d/%d tasks enabled", enabled_tasks, total_tasks))
  }
  
  # Filtering parameters
  if (!is.null(config$filtering)) {
    snp_threshold <- config$filtering$snp_threshold %||% "default"
    indel_threshold <- config$filtering$indel_threshold %||% "default"
    summary_lines <- c(summary_lines,
                      sprintf("- **Filter Thresholds**: SNP >= %s, INDEL >= %s", snp_threshold, indel_threshold))
  }
  
  return(paste(summary_lines, collapse = "\n"))
}

# Log successful loading
# if (exists("log_message")) {
#   log_message("Analysis pipeline flowchart generator loaded successfully")
# } else {
#   cat("Analysis pipeline flowchart generator loaded successfully\n")
# }