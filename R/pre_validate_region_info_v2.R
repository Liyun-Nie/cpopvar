#' Pre-validate Region Info Data Quality (Version 2)
#'
#' @description
#' Pre-validation module V2 for Geneious Prime exports, with improved
#' detection of independent gene copies. Issues are reported only; no
#' automatic repairs are applied.
#'
#' Features:
#' \itemize{
#'   \item Geneious Prime CSV export support
#'   \item automatic tRNA filtering
#'   \item truncated-gene detection (coordinates containing `<` or `>`)
#'   \item multi-copy checks grouped by physical location
#'   \item missing trans-spliced rps12 introns reported as INFO
#' }
#'
#' @param input_path Path to annotation CSV file
#' @param gene_function_map_path Path to gene_function_map.csv (optional)
#' @param strict_mode Fail on warnings (default: FALSE)
#' @param report_path Path to save quality report (default: auto-generated)
#'
#' @return List with status and validation results
#' @export
pre_validate_region_info_v2 <- function(input_path,
                                        gene_function_map_path = NULL,
                                        strict_mode = FALSE,
                                        report_path = NULL) {
  
  message("=======================================================")
  message("  REGION INFO PRE-VALIDATION MODULE V2                 ")
  message("  WARNING: Supports Geneious Prime export format       ")
  message("  WARNING: Only checks and reports issues              ")
  message("=======================================================\n")
  
  # Load and standardize data
  message(sprintf("[1/8] Loading data from: %s", input_path))
  raw_data <- load_and_standardize_annotation_data(input_path)
  
  message(sprintf("  Loaded %d records from %d species", 
                 nrow(raw_data), 
                 length(unique(raw_data$Species))))
  message(sprintf("  Filtered out %d tRNA genes\n", 
                 attr(raw_data, "trna_filtered")))
  
  # Run all checks
  all_issues <- list()
  
  message("[2/8] Checking format...")
  all_issues <- c(all_issues, validate_format_v2(raw_data))
  
  message("[3/8] Detecting truncated genes...")
  all_issues <- c(all_issues, check_truncated_genes(raw_data))
  
  message("[4/8] Checking completeness (multi-copy aware)...")
  all_issues <- c(all_issues, check_gene_completeness_v2(raw_data))
  
  message("[5/8] Detecting duplicates...")
  all_issues <- c(all_issues, check_duplicates_v2(raw_data))
  
  message("[6/8] Validating gene names...")
  all_issues <- c(all_issues, check_gene_names_v2(raw_data, gene_function_map_path))
  
  message("[7/8] Checking negative coordinates...")
  all_issues <- c(all_issues, check_negative_coordinates_v2(raw_data))
  
  message("[8/8] Validating lengths...\n")
  all_issues <- c(all_issues, check_lengths_v2(raw_data))
  
  # Generate report
  if (is.null(report_path)) {
    report_path <- sprintf("region_info_pre_validation_report_%s.json", 
                          format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  report <- generate_quality_report_v2(all_issues, report_path)
  
  # Determine final status
  critical_count <- sum(sapply(all_issues, function(x) x$severity == "CRITICAL"))
  warning_count <- sum(sapply(all_issues, function(x) x$severity == "WARNING"))
  
  status <- "PASS"
  if (critical_count > 0) {
    status <- "FAIL"
  } else if (strict_mode && warning_count > 0) {
    status <- "FAIL"
  }
  
  # Print summary
  if (status == "FAIL") {
    message("\n")
    message("=======================================================")
    message("  VALIDATION FAILED")
    message("=======================================================")
    message(sprintf("  Critical issues: %d", critical_count))
    message(sprintf("  Warnings: %d", warning_count))
    message(sprintf("\n  Full report: %s", report_path))
    message("\n  WARNING: Please fix critical issues and resubmit.\n")
  } else {
    message("\n")
    message("=======================================================")
    message("  VALIDATION PASSED")
    message("=======================================================")
    if (warning_count > 0) {
      message(sprintf("  Warnings: %d (review recommended)", warning_count))
    }
    message(sprintf("\n  Full report: %s\n", report_path))
  }
  
  return(invisible(list(
    status = status,
    report = report$report,
    path = report$path,
    summary = list(
      total_issues = length(all_issues),
      critical = critical_count,
      warnings = warning_count
    )
  )))
}


#' Load and Standardize Annotation Data for Pre-validation
#'
#' @description
#' Load annotation data and standardise column names. Supported formats:
#' - Legacy: Species, Name, Gene, Type, Minimum, Maximum, Length, Intervals
#' - Geneious: Sequence Name, Name, gene, Type, Minimum, Maximum, Length, # Intervals
#'
#' tRNA genes are filtered automatically.
#'
#' @keywords internal
load_and_standardize_annotation_data <- function(input_path) {
  
  # Detect file format
  if (grepl("\\.txt$", input_path)) {
    raw_data <- readr::read_tsv(input_path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  } else {
    raw_data <- readr::read_csv(input_path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  }
  
  # Count original rows
  original_rows <- nrow(raw_data)
  
  # Detect format and standardize column names
  if ("Sequence Name" %in% colnames(raw_data)) {
    # Geneious format
    raw_data <- raw_data %>%
      dplyr::rename(
        Species = `Sequence Name`,
        Gene = gene
      )
    
    if ("# Intervals" %in% colnames(raw_data)) {
      raw_data <- raw_data %>%
        dplyr::rename(Intervals = `# Intervals`)
    }
  }
  
  # Filter out tRNA genes
  trna_count <- sum(grepl("^trn", raw_data$Gene, ignore.case = TRUE) | 
                   grepl("tRNA", raw_data$Type, ignore.case = TRUE))
  
  raw_data <- raw_data %>%
    dplyr::filter(!grepl("^trn", Gene, ignore.case = TRUE)) %>%
    dplyr::filter(!grepl("tRNA", Type, ignore.case = TRUE))
  
  # Store metadata
  attr(raw_data, "trna_filtered") <- trna_count
  
  # Convert numeric columns (but keep original for truncation detection)
  raw_data <- raw_data %>%
    dplyr::mutate(
      Minimum_raw = Minimum,
      Maximum_raw = Maximum,
      Length_raw = Length,
      Minimum = as.numeric(gsub("[<>]", "", Minimum)),
      Maximum = as.numeric(gsub("[<>]", "", Maximum)),
      Length = as.numeric(gsub("[<>]", "", Length)),
      Intervals = as.numeric(Intervals)
    )
  
  return(raw_data)
}


#' Validate Format V2
#' @keywords internal
validate_format_v2 <- function(raw_data) {
  
  issues <- list()
  
  # Check required columns
  required_cols <- c("Species", "Name", "Gene", "Type", "Minimum", "Maximum", "Length", "Intervals")
  missing_cols <- setdiff(required_cols, colnames(raw_data))
  
  if (length(missing_cols) > 0) {
    issues[[length(issues) + 1]] <- list(
      severity = "CRITICAL",
      category = "FORMAT",
      issue = "Missing required columns",
      details = sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")),
      suggestion = "Ensure input file has all required columns"
    )
    return(issues)
  }
  
  # Check data types
  if (!is.numeric(raw_data$Minimum) || !is.numeric(raw_data$Maximum)) {
    issues[[length(issues) + 1]] <- list(
      severity = "CRITICAL",
      category = "FORMAT",
      issue = "Invalid coordinate data types",
      details = "Minimum/Maximum must be numeric",
      suggestion = "Check for non-numeric values in coordinate columns"
    )
  }
  
  # Check for empty values
  empty_species <- sum(is.na(raw_data$Species) | raw_data$Species == "")
  empty_gene <- sum(is.na(raw_data$Gene) | raw_data$Gene == "")
  
  if (empty_species > 0 || empty_gene > 0) {
    issues[[length(issues) + 1]] <- list(
      severity = "CRITICAL",
      category = "FORMAT",
      issue = "Empty required fields",
      details = sprintf("Empty Species: %d, Empty Gene: %d", empty_species, empty_gene),
      suggestion = "Fill in all required fields"
    )
  }
  
  # Check for (modified) suffix in species names
  modified_species <- raw_data %>%
    dplyr::filter(grepl("\\(modified\\)", Species, ignore.case = TRUE)) %>%
    dplyr::select(Species) %>%
    dplyr::distinct()
  
  if (nrow(modified_species) > 0) {
    issues[[length(issues) + 1]] <- list(
      severity = "WARNING",
      category = "FORMAT",
      issue = "Species names contain '(modified)' suffix",
      details = sprintf("%d species with '(modified)' suffix: %s", 
                       nrow(modified_species),
                       paste(unique(modified_species$Species), collapse = ", ")),
      suggestion = paste(
        "The '(modified)' suffix is added by Geneious Prime when annotations are edited.",
        "This suffix should be removed from species names to ensure proper species matching in downstream analysis.",
        "Please remove '(modified)' from all species names in the source file.",
        sep = "\n"
      )
    )
  }
  
  return(issues)
}


#' Check Truncated Genes
#' @keywords internal
check_truncated_genes <- function(raw_data) {
  
  issues <- list()
  
  if (!"Minimum_raw" %in% colnames(raw_data)) {
    return(issues)
  }
  
  # Detect < or > in coordinates
  truncated <- raw_data %>%
    dplyr::filter(grepl("[<>]", Minimum_raw) | grepl("[<>]", Maximum_raw)) %>%
    dplyr::select(Species, Gene, Name, Minimum_raw, Maximum_raw) %>%
    dplyr::distinct()
  
  if (nrow(truncated) > 0) {
    # Group by species
    for (species in unique(truncated$Species)) {
      species_truncated <- truncated %>%
        dplyr::filter(Species == species)
      
      issues[[length(issues) + 1]] <- list(
        severity = "INFO",
        category = "TRUNCATED",
        issue = "Truncated gene detected",
        details = sprintf("%s: %d genes with truncated coordinates (%s)", 
                         species,
                         nrow(species_truncated),
                         paste(unique(species_truncated$Gene), collapse = ", ")),
        suggestion = "Review these genes - coordinates contain < or > indicating truncation",
        affected_genes = paste(species_truncated$Gene, collapse = "|")
      )
    }
  }
  
  return(issues)
}


#' Check Gene Completeness V2 (Multi-copy Aware)
#'
#' @description
#' Completeness check that groups gene copies by physical location and
#' validates the exon/intron structure of each copy independently.
#'
#' @keywords internal
check_gene_completeness_v2 <- function(raw_data) {
  
  issues <- list()
  
  # Get multi-interval genes (CDS with Intervals > 1)
  multi_interval_genes <- raw_data %>%
    dplyr::filter(Type == "CDS", Intervals > 1) %>%
    dplyr::select(Species, Gene, Intervals, Minimum, Maximum) %>%
    dplyr::distinct()
  
  if (nrow(multi_interval_genes) == 0) {
    return(issues)
  }
  
  # For each multi-interval gene, group by physical location
  for (i in 1:nrow(multi_interval_genes)) {
    species <- multi_interval_genes$Species[i]
    gene <- multi_interval_genes$Gene[i]
    expected_intervals <- multi_interval_genes$Intervals[i]
    cds_min <- multi_interval_genes$Minimum[i]
    cds_max <- multi_interval_genes$Maximum[i]
    
    # Get all CDS records for this gene (may have multiple copies)
    gene_cds_records <- raw_data %>%
      dplyr::filter(Species == species, Gene == gene, Type == "CDS")
    
    # For each CDS record, find its exons/introns within its coordinate range
    for (j in 1:nrow(gene_cds_records)) {
      cds_record <- gene_cds_records[j, ]
      cds_start <- cds_record$Minimum
      cds_end <- cds_record$Maximum
      cds_intervals <- cds_record$Intervals
      
      # Find exons within this CDS range (with 10bp tolerance for boundary differences)
      exons <- raw_data %>%
        dplyr::filter(
          Species == species,
          Gene == gene,
          Type == "exon",
          Minimum >= (cds_start - 10),
          Maximum <= (cds_end + 10)
        )
      
      # Find introns within this CDS range
      introns <- raw_data %>%
        dplyr::filter(
          Species == species,
          Gene == gene,
          Type == "intron",
          Minimum >= (cds_start - 10),
          Maximum <= (cds_end + 10)
        )
      
      exon_count <- nrow(exons)
      intron_count <- nrow(introns)
      expected_exons <- cds_intervals
      expected_introns <- cds_intervals - 1
      
      # Check exon count
      if (exon_count != expected_exons) {
        issues[[length(issues) + 1]] <- list(
          severity = "CRITICAL",
          category = "COMPLETENESS",
          issue = "Missing or extra exon annotations",
          details = sprintf("%s - %s (copy at %d-%d): Expected %d exons, found %d", 
                           species, gene, cds_start, cds_end, expected_exons, exon_count),
          suggestion = sprintf("Check exon annotations for %s in this region", gene),
          affected_gene = sprintf("%s|%s|%d-%d", species, gene, cds_start, cds_end)
        )
      }
      
      # Check intron count (special handling for rps12)
      if (intron_count != expected_introns) {
        severity_level <- if (gene == "rps12") "INFO" else "CRITICAL"
        
        issues[[length(issues) + 1]] <- list(
          severity = severity_level,
          category = "COMPLETENESS",
          issue = if (gene == "rps12") "rps12 trans-spliced gene (expected)" else "Missing or extra intron annotations",
          details = sprintf("%s - %s (copy at %d-%d): Expected %d introns, found %d%s", 
                           species, gene, cds_start, cds_end, expected_introns, intron_count,
                           if (gene == "rps12") " (trans-spliced gene, introns may span large regions)" else ""),
          suggestion = if (gene == "rps12") {
            "rps12 is a trans-spliced gene with introns spanning large genomic regions - this is normal"
          } else {
            sprintf("Check intron annotations for %s in this region", gene)
          },
          affected_gene = sprintf("%s|%s|%d-%d", species, gene, cds_start, cds_end)
        )
      }
    }
  }
  
  return(issues)
}


#' Check Duplicates V2
#' @keywords internal
check_duplicates_v2 <- function(raw_data) {
  
  issues <- list()
  
  # Check exact duplicates
  duplicates <- raw_data %>%
    dplyr::group_by(Species, Name, Type, Minimum, Maximum) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()
  
  if (nrow(duplicates) > 0) {
    # Group duplicates by species and gene for reporting
    duplicate_summary <- duplicates %>%
      dplyr::mutate(
        Gene = ifelse(is.na(Gene) | Gene == "", 
                     stringr::str_extract(Name, "^[^_]+"), 
                     Gene)
      ) %>%
      dplyr::group_by(Species, Gene, Name, Type) %>%
      dplyr::summarise(
        count = dplyr::n(),
        coordinates = sprintf("%d-%d", min(Minimum), max(Maximum)),
        .groups = "drop"
      )
    
    # Create detailed report for each duplicate
    for (i in seq_len(nrow(duplicate_summary))) {
      row <- duplicate_summary[i, ]
      issues[[length(issues) + 1]] <- list(
        severity = "CRITICAL",
        category = "DUPLICATE",
        issue = "Exact duplicate records",
        details = sprintf("%s - %s (%s, %s): %d duplicate records at %s", 
                         row$Species, row$Gene, row$Name, row$Type, 
                         row$count, row$coordinates),
        suggestion = sprintf("Remove duplicate %s records for %s in %s", 
                           row$Type, row$Gene, row$Species),
        affected_gene = sprintf("%s|%s", row$Species, row$Gene)
      )
    }
  }
  
  # Check coordinate overlaps (within same gene, same type, but different Name)
  # This catches cases like slightly different coordinates for same annotation
  overlap_check <- raw_data %>%
    dplyr::filter(Type %in% c("exon", "intron")) %>%
    dplyr::group_by(Species, Gene, Type) %>%
    dplyr::arrange(Minimum) %>%
    dplyr::mutate(
      next_name = dplyr::lead(Name),
      next_min = dplyr::lead(Minimum),
      next_max = dplyr::lead(Maximum),
      coord_diff = pmin(abs(Minimum - next_min), abs(Maximum - next_max))
    ) %>%
    dplyr::filter(
      Name != next_name,  # Different names
      coord_diff < 10     # But coordinates very close
    ) %>%
    dplyr::ungroup()
  
  if (nrow(overlap_check) > 0) {
    # Group by species and gene
    for (species in unique(overlap_check$Species)) {
      for (gene in unique(overlap_check$Gene[overlap_check$Species == species])) {
        overlaps <- overlap_check %>%
          dplyr::filter(Species == species, Gene == gene)
        
        issues[[length(issues) + 1]] <- list(
          severity = "WARNING",
          category = "DUPLICATE",
          issue = "Overlapping coordinates",
          details = sprintf("%s - %s: %d records with overlapping coordinates (%d-%d)", 
                           species, gene, nrow(overlaps),
                           min(overlaps$Minimum), max(overlaps$Maximum)),
          suggestion = "Keep only one record per region (remove duplicates from source data)",
          affected_region = sprintf("%s|%s", species, gene)
        )
      }
    }
  }
  
  return(issues)
}


#' Check Gene Names V2
#' @keywords internal
check_gene_names_v2 <- function(raw_data, gene_function_map_path = NULL) {
  
  issues <- list()
  
  # Load gene function map
  if (is.null(gene_function_map_path)) {
    gene_function_map_path <- system.file("inst/config", "gene_function_map.csv", 
                                         package = "cpopvar")
    if (!file.exists(gene_function_map_path)) {
      gene_function_map_path <- "cpopvar/inst/config/gene_function_map.csv"
    }
  }
  
  if (!file.exists(gene_function_map_path)) {
    issues[[length(issues) + 1]] <- list(
      severity = "WARNING",
      category = "GENE_NAME",
      issue = "Gene function map not found",
      details = sprintf("Cannot validate gene names: %s", gene_function_map_path),
      suggestion = "Ensure gene_function_map.csv exists"
    )
    return(issues)
  }
  
  gene_function_map <- readr::read_csv(gene_function_map_path, show_col_types = FALSE)
  valid_gene_names <- unique(gene_function_map$gene_id)
  
  # Check CDS gene names
  cds_genes <- raw_data %>%
    dplyr::filter(Type == "CDS") %>%
    dplyr::select(Species, Gene) %>%
    dplyr::distinct()
  
  invalid_genes <- cds_genes %>%
    dplyr::filter(!Gene %in% valid_gene_names)
  
  if (nrow(invalid_genes) > 0) {
    for (species in unique(invalid_genes$Species)) {
      species_invalid <- invalid_genes %>%
        dplyr::filter(Species == species)
      
      issues[[length(issues) + 1]] <- list(
        severity = "CRITICAL",
        category = "GENE_NAME",
        issue = "Gene name not in reference list",
        details = sprintf("%s: %d gene names not found in reference list: %s", 
                         species,
                         nrow(species_invalid),
                         paste(species_invalid$Gene, collapse = ", ")),
        suggestion = paste(
          "These gene names are not in the built-in gene_function_map.csv reference list.",
          "Please verify:",
          "  1. Check if the gene name spelling is correct",
          "  2. Confirm this is a valid chloroplast gene",
          "  3. If correct, you may need to add it to the reference list",
          sprintf("Reference file: %s", gene_function_map_path),
          sep = "\n"
        ),
        affected_genes = paste(species_invalid$Gene, collapse = "|")
      )
    }
  }
  
  return(issues)
}


#' Check Negative Coordinates V2
#' @keywords internal
check_negative_coordinates_v2 <- function(raw_data) {
  
  issues <- list()
  
  # Check for truly negative coordinates (data errors)
  negative_coords <- raw_data %>%
    dplyr::filter(Minimum < 0 | Maximum < 0)
  
  if (nrow(negative_coords) > 0) {
    issues[[length(issues) + 1]] <- list(
      severity = "CRITICAL",
      category = "COORDINATE",
      issue = "Negative coordinates detected",
      details = sprintf("%d records with negative coordinates", nrow(negative_coords)),
      suggestion = "Fix negative coordinate values in source data"
    )
  }
  
  return(issues)
}


#' Check Lengths V2
#' @keywords internal
check_lengths_v2 <- function(raw_data) {
  
  issues <- list()
  
  # Check CDS length
  unusual_cds <- raw_data %>%
    dplyr::filter(Type == "CDS", Length < 100 | Length > 10000)
  
  if (nrow(unusual_cds) > 0) {
    issues[[length(issues) + 1]] <- list(
      severity = "INFO",
      category = "LENGTH",
      issue = "Unusual CDS length",
      details = sprintf("%d CDS records with unusual length (<100bp or >10000bp)", 
                       nrow(unusual_cds)),
      suggestion = "Review these records for potential annotation errors"
    )
  }
  
  return(issues)
}


#' Generate Quality Report V2
#' @keywords internal
generate_quality_report_v2 <- function(all_issues, output_path) {
  
  # Categorize issues
  critical <- Filter(function(x) x$severity == "CRITICAL", all_issues)
  warnings <- Filter(function(x) x$severity == "WARNING", all_issues)
  info <- Filter(function(x) x$severity == "INFO", all_issues)
  
  # Create report
  report <- list(
    summary = list(
      total_issues = length(all_issues),
      critical = length(critical),
      warnings = length(warnings),
      info = length(info),
      status = if (length(critical) > 0) "FAIL" else "PASS"
    ),
    issues = list(
      critical = critical,
      warnings = warnings,
      info = info
    ),
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  # Save JSON report
  jsonlite::write_json(report, output_path, pretty = TRUE, auto_unbox = TRUE)
  
  # Also save TXT report
  txt_path <- sub("\\.json$", ".txt", output_path)
  generate_txt_report_v2(report, txt_path)
  
  return(list(
    report = report,
    path = output_path
  ))
}


#' Generate TXT Report V2
#' @keywords internal
generate_txt_report_v2 <- function(report, output_path) {
  
  con <- file(output_path, "w", encoding = "UTF-8")
  on.exit(close(con))
  
  # Header
  writeLines("=======================================================", con)
  writeLines("       REGION INFO PRE-VALIDATION QUALITY REPORT V2    ", con)
  writeLines("=======================================================", con)
  writeLines("", con)
  
  # Summary
  writeLines(sprintf("Validation Status: %s", report$summary$status), con)
  writeLines(sprintf("Total Issues: %d", report$summary$total_issues), con)
  writeLines(sprintf("  - Critical: %d", report$summary$critical), con)
  writeLines(sprintf("  - Warnings: %d", report$summary$warnings), con)
  writeLines(sprintf("  - Info: %d", report$summary$info), con)
  writeLines(sprintf("Timestamp: %s", report$timestamp), con)
  writeLines("", con)
  
  # Critical issues
  if (length(report$issues$critical) > 0) {
    writeLines("-------------------------------------------------------", con)
    writeLines("CRITICAL ISSUES (Must Fix)", con)
    writeLines("-------------------------------------------------------", con)
    for (i in seq_along(report$issues$critical)) {
      issue <- report$issues$critical[[i]]
      writeLines(sprintf("\n[%d] %s", i, issue$issue), con)
      writeLines(sprintf("    Category: %s", issue$category), con)
      writeLines(sprintf("    Details: %s", issue$details), con)
      writeLines(sprintf("    Suggestion: %s", issue$suggestion), con)
    }
    writeLines("", con)
  }
  
  # Warnings
  if (length(report$issues$warnings) > 0) {
    writeLines("-------------------------------------------------------", con)
    writeLines("WARNINGS (Review Recommended)", con)
    writeLines("-------------------------------------------------------", con)
    for (i in seq_along(report$issues$warnings)) {
      issue <- report$issues$warnings[[i]]
      writeLines(sprintf("\n[%d] %s", i, issue$issue), con)
      writeLines(sprintf("    Category: %s", issue$category), con)
      writeLines(sprintf("    Details: %s", issue$details), con)
      writeLines(sprintf("    Suggestion: %s", issue$suggestion), con)
    }
    writeLines("", con)
  }
  
  # Info
  if (length(report$issues$info) > 0) {
    writeLines("-------------------------------------------------------", con)
    writeLines("INFORMATIONAL MESSAGES", con)
    writeLines("-------------------------------------------------------", con)
    for (i in seq_along(report$issues$info)) {
      issue <- report$issues$info[[i]]
      writeLines(sprintf("\n[%d] %s", i, issue$issue), con)
      writeLines(sprintf("    Category: %s", issue$category), con)
      writeLines(sprintf("    Details: %s", issue$details), con)
      if (!is.null(issue$suggestion)) {
        writeLines(sprintf("    Note: %s", issue$suggestion), con)
      }
    }
    writeLines("", con)
  }
  
  writeLines("=======================================================", con)
  writeLines("                    END OF REPORT                       ", con)
  writeLines("=======================================================", con)
}

