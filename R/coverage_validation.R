############################################################
#### coverage_validation.R - Region Coverage Validation ####
############################################################

#' @importFrom magrittr %>%
NULL

#' Validate region coverage and identify missing variants
#' 
#' @description Perform coverage validation to ensure
#' all genome regions are properly covered and no variants are lost between P02 and P03.
#' This validation consists of two core steps:
#' 1. Region backfill: Map region_info_complete.csv back to processed_genome_regions.csv
#' 2. Variant backfill: Check if P02 variants fall within uncovered regions
#' 
#' @param region_info_df Data frame from region_info_complete.csv
#' @param genome_regions_df Data frame from processed_genome_regions.csv
#' @param filtered_variants_df Data frame from variants_filtered.csv (P02 output)
#' @param output_dir Directory to save coverage validation report
#' @param log_file Path to log file
#' @return List containing validation results and statistics
#' @export
validate_region_coverage <- function(region_info_df,
                                    genome_regions_df,
                                    filtered_variants_df,
                                    output_dir,
                                    log_file = NULL) {
  
  log_message("=== GATE 4.3 S3: REGION COVERAGE VALIDATION ===")
  log_message("Starting comprehensive coverage validation")
  
  # Initialize results
  validation_results <- list(
    region_backfill = list(),
    variant_backfill = list(),
    summary_stats = list()
  )
  
  # ========== STEP 1: REGION BACKFILL ==========
  log_message("STEP 1: REGION BACKFILL - Mapping constructed regions to genome coordinates")
  
  region_backfill_results <- perform_region_backfill(
    region_info_df,
    genome_regions_df
  )
  
  validation_results$region_backfill <- region_backfill_results
  
  # ========== STEP 2: VARIANT BACKFILL (Exclude RNA) ==========
  log_message("STEP 2: VARIANT BACKFILL - Checking if P02 variants fall in uncovered regions")
  log_message("NOTE: RNA variants are excluded from this validation")
  
  # Filter out RNA variants
  non_rna_variants <- filtered_variants_df %>%
    dplyr::filter(region_type != "RNA")
  
  log_message(sprintf("Total P02 variants: %d", nrow(filtered_variants_df)))
  log_message(sprintf("Non-RNA variants for validation: %d", nrow(non_rna_variants)))
  log_message(sprintf("RNA variants excluded: %d", nrow(filtered_variants_df) - nrow(non_rna_variants)))
  
  variant_backfill_results <- perform_variant_backfill(
    non_rna_variants,
    region_backfill_results$uncovered_gaps
  )
  
  validation_results$variant_backfill <- variant_backfill_results
  
  # ========== STEP 3: GENERATE SUMMARY STATISTICS ==========
  log_message("STEP 3: GENERATING SUMMARY STATISTICS")
  
  summary_stats <- generate_coverage_summary(
    region_backfill_results,
    variant_backfill_results,
    genome_regions_df
  )
  
  validation_results$summary_stats <- summary_stats
  
  # ========== STEP 4: EXPORT VALIDATION REPORT ==========
  log_message("STEP 4: EXPORTING VALIDATION REPORT")
  
  report_path <- file.path(output_dir, "coverage_validation_report.txt")
  export_coverage_report(validation_results, report_path)
  
  log_message(sprintf("Coverage validation report saved: %s", report_path))
  log_message("=== COVERAGE VALIDATION COMPLETED ===")
  
  return(validation_results)
}


#' Perform region backfill validation
#' 
#' @description Map constructed regions back to genome coordinates to identify gaps
#' @param region_info_df Region information data frame
#' @param genome_regions_df Genome regions configuration data frame
#' @return List containing coverage results and uncovered gaps
#' @keywords internal
perform_region_backfill <- function(region_info_df, genome_regions_df) {
  
  log_message("=== REGION BACKFILL VALIDATION (Robust Version) ===")
  
  # Initialize results
  all_coverage_stats <- list()
  all_uncovered_gaps <- list()
  
  # Process each species
  species_list <- unique(genome_regions_df$species)
  
  for (sp in species_list) {
    log_message(sprintf("--- Processing species: %s ---", sp))
    
    # Get species-specific data
    sp_regions <- region_info_df %>%
      dplyr::filter(species == sp) %>%
      dplyr::arrange(region_start)
    
    sp_genome <- genome_regions_df %>%
      dplyr::filter(species == sp)
    
    log_message(sprintf("[%s] Found %d regions in region_info_df.", sp, nrow(sp_regions)))
    log_message(sprintf("[%s] Found %d rows in genome_regions_df.", sp, nrow(sp_genome)))
    
    if (nrow(sp_regions) == 0) {
      log_message(sprintf("WARNING: No regions found for %s. Skipping species.", sp), level = "warning")
      next
    }
    
    if (nrow(sp_genome) == 0) {
      log_message(sprintf("WARNING: No genome configuration found for %s. Skipping species.", sp), level = "warning")
      next
    }
    
    # Extract single row values safely
    special_handling <- if ("special_handling" %in% names(sp_genome)) sp_genome$special_handling[1] else NA
    log_message(sprintf("[%s] special_handling: '%s'", sp, ifelse(is.na(special_handling), "NA", special_handling)))

    # Identify genome structure
    genome_segments <- list()
    is_ir_lacking <- !is.na(special_handling) && special_handling == "IR_lacking_genome"
    
    if (is_ir_lacking) {
      log_message(sprintf("[%s] Detected IR-lacking genome.", sp))
      if (!is.na(sp_genome$lsc_start[1]) && !is.na(sp_genome$lsc_end[1])) {
        genome_segments <- c(genome_segments, list(list(name = "LSC", start = sp_genome$lsc_start[1], end = sp_genome$lsc_end[1])))
      }
      if (!is.na(sp_genome$ssc_start[1]) && !is.na(sp_genome$ssc_end[1])) {
        genome_segments <- c(genome_segments, list(list(name = "SSC", start = sp_genome$ssc_start[1], end = sp_genome$ssc_end[1])))
      }
      if (length(genome_segments) == 0 && !is.na(sp_genome$total_length[1])) {
        genome_segments <- c(genome_segments, list(list(name = "whole_genome", start = 1, end = sp_genome$total_length[1])))
        log_message(sprintf("[%s] Using 'whole_genome' as fallback.", sp))
      }
    } else {
      log_message(sprintf("[%s] Detected standard genome.", sp))
      if (!is.na(sp_genome$lsc_start[1]) && !is.na(sp_genome$lsc_end[1])) {
        genome_segments <- c(genome_segments, list(list(name = "LSC", start = sp_genome$lsc_start[1], end = sp_genome$lsc_end[1])))
      }
      if (!is.na(sp_genome$ir_start[1]) && !is.na(sp_genome$ir_end[1])) {
        genome_segments <- c(genome_segments, list(list(name = "IRA", start = sp_genome$ir_start[1], end = sp_genome$ir_end[1])))
      }
      if (!is.na(sp_genome$ssc_start[1]) && !is.na(sp_genome$ssc_end[1])) {
        genome_segments <- c(genome_segments, list(list(name = "SSC", start = sp_genome$ssc_start[1], end = sp_genome$ssc_end[1])))
      }
    }

    if (length(genome_segments) == 0) {
      log_message(sprintf("WARNING: No valid genome segments defined for %s. Skipping species.", sp), level = "warning")
      next
    }
    
    log_message(sprintf("[%s] Defined %d genome segments for validation.", sp, length(genome_segments)))

    # Calculate coverage for each segment
    for (seg in genome_segments) {
      seg_name <- seg$name
      seg_start <- seg$start
      seg_end <- seg$end
      
      log_message(sprintf("[%s] Analyzing segment: %s (%d - %d)", sp, seg_name, seg_start, seg_end))

      # --- CRITICAL FIX: Input Validation ---
      if (is.na(seg_start) || is.na(seg_end)) {
        log_message(sprintf("WARNING: [%s] Segment '%s' has NA coordinates. Skipping segment.", sp, seg_name), level = "warning")
        next
      }

      seg_length <- seg_end - seg_start + 1
      if (seg_length <= 0) {
          log_message(sprintf("WARNING: [%s] Segment '%s' has non-positive length (%d). Skipping segment.", sp, seg_name, seg_length), level = "warning")
          next
      }

      # Get regions in this segment
      seg_regions <- sp_regions %>%
        dplyr::filter(region_start >= seg_start, region_end <= seg_end) %>%
        dplyr::arrange(region_start)
      
      log_message(sprintf("[%s] Segment %s has %d regions.", sp, seg_name, nrow(seg_regions)))

      gaps <- list()
      
      # Use a variable to track the covered position
      current_pos <- seg_start
      
      if (nrow(seg_regions) > 0) {
        # Check for a gap before the first region
        first_region_start <- seg_regions$region_start[1]
        if (first_region_start > current_pos) {
          gap_len <- first_region_start - current_pos
          gaps[[length(gaps) + 1]] <- list(
            start = current_pos, end = first_region_start - 1, len = gap_len,
            reason = if (gap_len < 3) "Below IGS threshold" else "Missing IGS"
          )
        }
        
        # Iterate through regions to find inter-region gaps
        for (i in 1:nrow(seg_regions)) {
          region <- seg_regions[i, ]
          # Gap between last position and current region start
          if (region$region_start > current_pos) {
             gap_len <- region$region_start - current_pos
             gaps[[length(gaps) + 1]] <- list(
                start = current_pos, end = region$region_start - 1, len = gap_len,
                reason = if (gap_len < 3) "Below IGS threshold" else "Missing IGS"
             )
          }
          # Move current position to the end of the current region
          current_pos <- max(current_pos, region$region_end + 1)
        }
      }
      
      # Check for a gap after the last region
      if (current_pos <= seg_end) {
        gap_len <- seg_end - current_pos + 1
        gaps[[length(gaps) + 1]] <- list(
          start = current_pos, end = seg_end, len = gap_len,
          reason = if (nrow(seg_regions) == 0) "Entire segment uncovered" else "Missing IGS at end"
        )
      }

      # Calculate covered bases using a robust method to handle overlaps
      if (nrow(seg_regions) > 0) {
        # Create a logical vector representing the segment
        coverage_map <- rep(FALSE, seg_length)
        for (i in 1:nrow(seg_regions)) {
          start_in_map <- max(1, seg_regions$region_start[i] - seg_start + 1)
          end_in_map <- min(seg_length, seg_regions$region_end[i] - seg_start + 1)
          if (start_in_map <= end_in_map) {
            coverage_map[start_in_map:end_in_map] <- TRUE
          }
        }
        covered_bases <- sum(coverage_map)
      } else {
        covered_bases <- 0
      }
      
      coverage_rate <- if (seg_length > 0) (covered_bases / seg_length) * 100 else 0
      
      all_coverage_stats[[length(all_coverage_stats) + 1]] <- data.frame(
        species = sp,
        genome_segment = seg_name,
        segment_length = seg_length,
        covered_bases = covered_bases,
        coverage_rate = coverage_rate,
        gaps_count = length(gaps),
        stringsAsFactors = FALSE
      )

      # Convert gaps to data frame
      if (length(gaps) > 0) {
        gaps_df <- data.frame(
          species = sp,
          genome_region = seg_name,
          gap_start = sapply(gaps, `[[`, "start"),
          gap_end = sapply(gaps, `[[`, "end"),
          gap_length = sapply(gaps, `[[`, "len"),
          reason = sapply(gaps, `[[`, "reason"),
          stringsAsFactors = FALSE
        )
        all_uncovered_gaps[[length(all_uncovered_gaps) + 1]] <- gaps_df
      }
    }
  }
  
  # Combine results into final data frames
  final_coverage_stats <- if (length(all_coverage_stats) > 0) dplyr::bind_rows(all_coverage_stats) else data.frame()
  final_uncovered_gaps <- if (length(all_uncovered_gaps) > 0) dplyr::bind_rows(all_uncovered_gaps) else data.frame()

  log_message(sprintf("Region backfill complete. Found %d gaps across all species.", nrow(final_uncovered_gaps)))
  
  return(list(
    coverage_by_species = final_coverage_stats,
    uncovered_gaps = final_uncovered_gaps
  ))
}


#' Perform variant backfill validation
#' 
#' @description Check if P02 filtered variants fall within uncovered gaps
#' @param filtered_variants_df Filtered variants from P02 (non-RNA only)
#' @param uncovered_gaps_df Data frame of uncovered gaps
#' @return List containing variant mapping results
#' @keywords internal
perform_variant_backfill <- function(filtered_variants_df, uncovered_gaps_df) {
  
  log_message("=== VARIANT BACKFILL VALIDATION ===")
  
  if (nrow(uncovered_gaps_df) == 0) {
    log_message("No uncovered gaps detected - all regions are properly covered")
    return(list(
      variants_in_gaps = data.frame(),
      variants_in_gaps_count = 0,
      coverage_complete = TRUE
    ))
  }
  
  # Initialize results
  variants_in_gaps <- data.frame()
  
  # Process each species
  species_list <- unique(uncovered_gaps_df$species)
  
  for (sp in species_list) {
    sp_variants <- filtered_variants_df %>%
      dplyr::filter(species == sp) %>%
      dplyr::select(species, position, var_type, region_type, genome_region, gene)
    
    sp_gaps <- uncovered_gaps_df %>%
      dplyr::filter(species == sp)
    
    if (nrow(sp_variants) == 0) {
      next
    }
    
    # Check each gap
    for (i in 1:nrow(sp_gaps)) {
      gap <- sp_gaps[i, ]
      
      # Safely extract gap coordinates
      gap_start_val <- if (length(gap$gap_start) > 0 && !is.na(gap$gap_start[1])) gap$gap_start[1] else NA
      gap_end_val <- if (length(gap$gap_end) > 0 && !is.na(gap$gap_end[1])) gap$gap_end[1] else NA
      gap_length_val <- if (length(gap$gap_length) > 0 && !is.na(gap$gap_length[1])) gap$gap_length[1] else NA
      gap_reason_val <- if (length(gap$reason) > 0 && !is.na(gap$reason[1])) gap$reason[1] else "Unknown"
      gap_region_val <- if (length(gap$genome_region) > 0 && !is.na(gap$genome_region[1])) gap$genome_region[1] else "Unknown"
      
      # Skip if gap coordinates are invalid
      if (is.na(gap_start_val) || is.na(gap_end_val)) {
        log_message(sprintf("WARNING: %s: Invalid gap coordinates at row %d, skipping", sp, i), level = "warning")
        next
      }
      
      # Find variants within this gap
      gap_variants <- sp_variants %>%
        dplyr::filter(position >= gap_start_val & position <= gap_end_val)
      
      if (nrow(gap_variants) > 0) {
        # Add gap information to variants
        gap_variants$gap_start <- gap_start_val
        gap_variants$gap_end <- gap_end_val
        gap_variants$gap_length <- gap_length_val
        gap_variants$gap_reason <- gap_reason_val
        gap_variants$gap_genome_region <- gap_region_val
        
        variants_in_gaps <- dplyr::bind_rows(variants_in_gaps, gap_variants)
      }
    }
    
    if (nrow(variants_in_gaps) > 0) {
      log_message(sprintf("%s: Found %d variants in uncovered gaps", sp, nrow(variants_in_gaps)))
    }
  }
  
  # Summary
  if (nrow(variants_in_gaps) > 0) {
    log_message(sprintf("WARNING: %d variants found in uncovered gaps", nrow(variants_in_gaps)), level = "warning")
    
    # Breakdown by reason
    reason_summary <- variants_in_gaps %>%
      dplyr::group_by(gap_reason) %>%
      dplyr::summarise(variant_count = dplyr::n(), .groups = "drop")
    
    log_message("Breakdown by gap reason:")
    for (i in 1:nrow(reason_summary)) {
      log_message(sprintf("  %s: %d variants", reason_summary$gap_reason[i], reason_summary$variant_count[i]))
    }
  } else {
    log_message("VALIDATION PASSED: No variants found in uncovered gaps")
  }
  
  return(list(
    variants_in_gaps = variants_in_gaps,
    variants_in_gaps_count = nrow(variants_in_gaps),
    coverage_complete = nrow(variants_in_gaps) == 0
  ))
}


#' Generate coverage validation summary statistics
#' 
#' @description Compile comprehensive summary statistics from validation results
#' @param region_backfill_results Results from region backfill
#' @param variant_backfill_results Results from variant backfill
#' @param genome_regions_df Genome regions configuration
#' @return List containing summary statistics
#' @keywords internal
generate_coverage_summary <- function(region_backfill_results,
                                     variant_backfill_results,
                                     genome_regions_df) {
  
  log_message("=== GENERATING COVERAGE SUMMARY ===")
  
  # Extract data frames from results
  coverage_stats_df <- region_backfill_results$coverage_by_species
  uncovered_gaps_df <- region_backfill_results$uncovered_gaps
  
  # Overall statistics
  total_species <- length(unique(coverage_stats_df$species))
  total_gaps <- nrow(uncovered_gaps_df)
  total_variants_in_gaps <- variant_backfill_results$variants_in_gaps_count
  
  log_message(sprintf("Processing %d species, %d gaps, %d variants in gaps",
                     total_species, total_gaps, total_variants_in_gaps))
  
  # Gap statistics by reason
  if (total_gaps > 0) {
    gaps_by_reason <- uncovered_gaps_df %>%
      dplyr::group_by(reason) %>%
      dplyr::summarise(
        gap_count = dplyr::n(),
        total_bases = sum(gap_length, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    gaps_by_reason <- data.frame(
      reason = character(),
      gap_count = integer(),
      total_bases = numeric(),
      stringsAsFactors = FALSE
    )
  }
  
  # Species-level statistics
  # Aggregate coverage stats by species
  species_summary <- coverage_stats_df %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
      total_genome_length = sum(segment_length, na.rm = TRUE),
      covered_bases = sum(covered_bases, na.rm = TRUE),
      total_segments = dplyr::n(),
      avg_coverage_rate = mean(coverage_rate, na.rm = TRUE),
      total_gaps_count = sum(gaps_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      uncovered_bases = total_genome_length - covered_bases,
      coverage_rate = dplyr::if_else(total_genome_length > 0, 
                                     (covered_bases / total_genome_length) * 100, 
                                     0)
    )
  
  # Add variant counts per species
  if (nrow(variant_backfill_results$variants_in_gaps) > 0) {
    variants_per_species <- variant_backfill_results$variants_in_gaps %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(variants_in_gaps = dplyr::n(), .groups = "drop")
    
    species_summary <- species_summary %>%
      dplyr::left_join(variants_per_species, by = "species") %>%
      dplyr::mutate(variants_in_gaps = tidyr::replace_na(variants_in_gaps, 0))
  } else {
    species_summary <- species_summary %>%
      dplyr::mutate(variants_in_gaps = 0)
  }
  
  log_message(sprintf("Summary complete: %d species analyzed", nrow(species_summary)))
  
  return(list(
    total_species = total_species,
    total_gaps = total_gaps,
    total_variants_in_gaps = total_variants_in_gaps,
    gaps_by_reason = gaps_by_reason,
    species_summary = species_summary,
    validation_passed = total_variants_in_gaps == 0
  ))
}


#' Export coverage validation report
#' 
#' @description Generate human-readable text report of coverage validation
#' @param validation_results Complete validation results
#' @param report_path Path to save report
#' @keywords internal
export_coverage_report <- function(validation_results, report_path) {
  
  # Open file connection
  report_conn <- file(report_path, "w")
  
  tryCatch({
    # Header
    writeLines("=============================================================", report_conn)
    writeLines("   GATE 4.3 S3: REGION COVERAGE VALIDATION REPORT", report_conn)
    writeLines("=============================================================", report_conn)
    writeLines(sprintf("Report generated: %s", Sys.time()), report_conn)
    writeLines("", report_conn)
    
    # Summary
    writeLines("--- SUMMARY STATISTICS ---", report_conn)
    summary <- validation_results$summary_stats
    writeLines(sprintf("Total species analyzed: %d", summary$total_species), report_conn)
    writeLines(sprintf("Total uncovered gaps: %d", summary$total_gaps), report_conn)
    writeLines(sprintf("Variants in uncovered gaps: %d", summary$total_variants_in_gaps), report_conn)
    writeLines(sprintf("Validation status: %s",
                      if (summary$validation_passed) "PASSED" else "FAILED"), report_conn)
    writeLines("", report_conn)
    
    # Gaps by reason
    if (nrow(summary$gaps_by_reason) > 0) {
      writeLines("--- UNCOVERED GAPS BY REASON ---", report_conn)
      for (i in 1:nrow(summary$gaps_by_reason)) {
        writeLines(sprintf("  %s: %d gaps, %d bp",
                          summary$gaps_by_reason$reason[i],
                          summary$gaps_by_reason$gap_count[i],
                          summary$gaps_by_reason$total_bases[i]), report_conn)
      }
      writeLines("", report_conn)
    }
    
    # Species-level coverage
    writeLines("--- SPECIES-LEVEL COVERAGE ---", report_conn)
    for (i in 1:nrow(summary$species_summary)) {
      sp_row <- summary$species_summary[i, ]
      writeLines(sprintf("%s:", sp_row$species), report_conn)
      writeLines(sprintf("  Genome length: %d bp", sp_row$total_genome_length), report_conn)
      writeLines(sprintf("  Covered: %d bp (%.2f%%)", sp_row$covered_bases, sp_row$coverage_rate), report_conn)
      writeLines(sprintf("  Uncovered: %d bp", sp_row$uncovered_bases), report_conn)
      writeLines(sprintf("  Gaps: %d", sp_row$gaps_count), report_conn)
      writeLines(sprintf("  Variants in gaps: %d", sp_row$variants_in_gaps), report_conn)
      writeLines("", report_conn)
    }
    
    # Detailed gap listing
    if (nrow(validation_results$region_backfill$uncovered_gaps) > 0) {
      writeLines("--- DETAILED GAP LISTING ---", report_conn)
      gaps <- validation_results$region_backfill$uncovered_gaps
      for (i in 1:nrow(gaps)) {
        gap <- gaps[i, ]
        writeLines(sprintf("%s [%s]: %d-%d (%d bp) - %s",
                          gap$species, gap$genome_region,
                          gap$gap_start, gap$gap_end, gap$gap_length,
                          gap$reason), report_conn)
      }
      writeLines("", report_conn)
    }
    
    # Variants in gaps
    if (nrow(validation_results$variant_backfill$variants_in_gaps) > 0) {
      writeLines("--- VARIANTS IN UNCOVERED GAPS (CRITICAL) ---", report_conn)
      vars <- validation_results$variant_backfill$variants_in_gaps
      for (i in 1:nrow(vars)) {
        var <- vars[i, ]
        writeLines(sprintf("%s: pos %d (%s, %s) in gap %d-%d (%s)",
                          var$species, var$position, var$var_type, var$region_type,
                          var$gap_start, var$gap_end, var$gap_reason), report_conn)
      }
      writeLines("", report_conn)
    }
    
    # Footer
    writeLines("=============================================================", report_conn)
    writeLines("END OF REPORT", report_conn)
    writeLines("=============================================================", report_conn)
    
  }, finally = {
    close(report_conn)
  })
  
  log_message("Coverage validation report exported successfully")
}

