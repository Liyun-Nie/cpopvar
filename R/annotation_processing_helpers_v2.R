#' Annotation Processing Helpers (Python-Compatible Version)
#'
#' This module provides functions to process chloroplast genome annotation data
#' from PCG_all_cds format, fully replicating the Python scripts logic:
#' - species_length_statistics_new.py
#' - calculate_cds_length.py
#' - gene_intergenic_intron_length.py
#'
#' Key Design Decisions:
#' 1. For multi-exon genes: Remove CDS row, keep only exon rows
#' 2. For single-exon genes: Keep CDS row (no exon rows exist)
#' 3. Add Gene column to all records for easy aggregation
#' 4. IGS calculation uses Python statistical method
#' 5. IGS naming format: "IGS_geneA-geneB" (no type suffixes)
#'
#' @name annotation_processing_helpers_v2
#' @keywords internal
NULL

#' Find Overlapping CDS Intervals
#'
#' Replicates the Python function `find_overlapping_intervals()`.
#' Identifies overlapping CDS regions and calculates total overlap length.
#'
#' @param intervals Data frame with columns: start, end, name, intervals_count
#'
#' @return List with:
#'   - overlaps: Data frame of overlapping intervals
#'   - total_overlap_length: Total length of overlaps
#'
#' @keywords internal
find_overlapping_intervals_r <- function(intervals) {
  
  if (nrow(intervals) == 0) {
    return(list(overlaps = data.frame(), total_overlap_length = 0))
  }
  
  # Ensure numeric types for coordinates
  intervals <- intervals %>%
    dplyr::mutate(
      start = as.numeric(start),
      end = as.numeric(end),
      intervals_count = as.numeric(intervals_count)
    )
  
  # Sort by start position
  intervals <- intervals %>% dplyr::arrange(start)
  
  overlaps <- list()
  total_overlap_length <- 0
  
  if (nrow(intervals) < 2) {
    return(list(overlaps = data.frame(), total_overlap_length = 0))
  }
  
  current_start <- intervals$start[1]
  current_end <- intervals$end[1]
  current_name <- intervals$name[1]
  current_intervals <- intervals$intervals_count[1]
  
  for (i in 2:nrow(intervals)) {
    start_i <- intervals$start[i]
    end_i <- intervals$end[i]
    name_i <- intervals$name[i]
    intervals_i <- intervals$intervals_count[i]
    
    if (start_i <= current_end) {
      # Overlap detected
      overlap_start <- start_i
      overlap_end <- min(current_end, end_i)
      overlap_length <- overlap_end - overlap_start + 1
      
      total_overlap_length <- total_overlap_length + overlap_length
      
      overlaps[[length(overlaps) + 1]] <- data.frame(
        cds1_name = current_name,
        intervals1 = current_intervals,
        start1 = current_start,
        end1 = current_end,
        cds2_name = name_i,
        intervals2 = intervals_i,
        start2 = start_i,
        end2 = end_i,
        overlap_length = overlap_length,
        stringsAsFactors = FALSE
      )
      
      current_end <- max(current_end, end_i)
    } else {
      current_start <- start_i
      current_end <- end_i
      current_name <- name_i
      current_intervals <- intervals_i
    }
  }
  
  if (length(overlaps) > 0) {
    overlaps_df <- dplyr::bind_rows(overlaps)
  } else {
    overlaps_df <- data.frame()
  }
  
  return(list(overlaps = overlaps_df, total_overlap_length = total_overlap_length))
}

#' Calculate IGS Length Using Python Statistical Method
#'
#' Replicates the Python function `calculate_statistics()`.
#' Formula: intergenic_length = genome_length - overlaps - total_CDS - intron_total
#'
#' @param raw_data Raw annotation data from PCG_all_cds
#'
#' @return Data frame with columns: species, intergenic_length
#'
#' @keywords internal
calculate_igs_total_length <- function(raw_data) {
  
  species_list <- unique(raw_data$Species)
  results <- list()
  
  for (sp in species_list) {
    sp_data <- raw_data %>% dplyr::filter(Species == sp)
    
    # 1. Genome length
    genome_length <- max(sp_data$Maximum, na.rm = TRUE)
    
    # 2. Find overlapping CDS intervals
    cds_data <- sp_data %>% dplyr::filter(Type == "CDS")
    intervals_df <- cds_data %>%
      dplyr::select(start = Minimum, end = Maximum, name = Name, intervals_count = Intervals)
    
    overlap_result <- find_overlapping_intervals_r(intervals_df)
    total_overlapping_length <- overlap_result$total_overlap_length
    
    # Log overlaps if found
    if (nrow(overlap_result$overlaps) > 0) {
      message(sprintf("  Species: %s has %d overlapping CDS regions (total overlap: %d bp)",
                      sp, nrow(overlap_result$overlaps), total_overlapping_length))
    }
    
    # 3. CDS without intron length
    cds_without_intron_length <- cds_data %>%
      dplyr::filter(Intervals == 1) %>%
      dplyr::pull(Length) %>%
      sum(na.rm = TRUE)
    
    # 4. CDS with intron length
    selected_cds_length <- cds_data %>%
      dplyr::filter(Intervals > 1) %>%
      dplyr::pull(Length) %>%
      sum(na.rm = TRUE)
    
    # 5. Exon length
    exon_length <- sp_data %>%
      dplyr::filter(Type == "exon") %>%
      dplyr::pull(Length) %>%
      sum(na.rm = TRUE)
    
    # 6. rps12 special handling
    exon_rps12_length <- sp_data %>%
      dplyr::filter(stringr::str_detect(Name, "rps12"), Type == "exon") %>%
      dplyr::pull(Length) %>%
      sum(na.rm = TRUE)
    
    intron_rps12_length <- sp_data %>%
      dplyr::filter(stringr::str_detect(Name, "rps12"), Type == "intron") %>%
      dplyr::pull(Length) %>%
      sum(na.rm = TRUE)
    
    # 7. Calculate intron total length
    intron_total_length <- selected_cds_length - (exon_length - exon_rps12_length) + intron_rps12_length
    
    # 8. Calculate total CDS length
    total_cds_length <- cds_without_intron_length + exon_length - total_overlapping_length
    
    # 9. Calculate intergenic length (Python formula)
    intergenic_length <- genome_length - total_overlapping_length - total_cds_length - intron_total_length
    
    results[[length(results) + 1]] <- data.frame(
      species = sp,
      genome_length = genome_length,
      total_cds_length = total_cds_length,
      intron_total_length = intron_total_length,
      intergenic_length = intergenic_length,
      total_overlapping_length = total_overlapping_length,
      stringsAsFactors = FALSE
    )
  }
  
  dplyr::bind_rows(results)
}

#' Generate IGS Records from Coordinate Gaps
#'
#' Generates individual IGS records based on gaps between adjacent features.
#' IGS naming follows the format: "IGS_geneA-geneB" (no type suffixes).
#'
#' @param combined_records Combined CDS/exon and intron records
#'
#' @return Data frame of IGS records
#'
#' @keywords internal
generate_igs_from_gaps <- function(combined_records) {
  
  species_list <- unique(combined_records$species)
  igs_records <- list()
  
  for (sp in species_list) {
    sp_features <- combined_records %>%
      dplyr::filter(species == sp) %>%
      dplyr::arrange(Start_Position)
    
    if (nrow(sp_features) < 2) next
    
    for (i in 1:(nrow(sp_features) - 1)) {
      current_feature <- sp_features[i, ]
      next_feature <- sp_features[i + 1, ]
      
      # Skip if any critical values are NA
      if (is.na(current_feature$End_Position) || is.na(next_feature$Start_Position)) {
        next
      }
      
      # Added check: only generate IGS between different genes.
      if (current_feature$Gene == next_feature$Gene) {
        next
      }
      
      # Calculate gap between features
      gap_start <- current_feature$End_Position + 1
      gap_end <- next_feature$Start_Position - 1
      gap_length <- gap_end - gap_start + 1
      
      # Check if gap_length is valid (not NA and positive)
      if (!is.na(gap_length) && gap_length > 0) {
        # Extract gene names (remove ALL type suffixes: _CDS, _exon, _intron)
        current_gene <- stringr::str_remove(current_feature$Region_Name, "_(CDS|exon|intron).*")
        next_gene <- stringr::str_remove(next_feature$Region_Name, "_(CDS|exon|intron).*")
        
        igs_name <- sprintf("IGS_%s-%s", current_gene, next_gene)
        
        igs_records[[length(igs_records) + 1]] <- data.frame(
          species = sp,
          Region_Name = igs_name,
          Type = "IGS",
          Start_Position = gap_start,
          End_Position = gap_end,
          Length = gap_length,
          Gene = NA_character_,  # IGS is intergenic and should NOT be assigned to any gene
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(igs_records) > 0) {
    dplyr::bind_rows(igs_records)
  } else {
    data.frame(
      species = character(),
      Region_Name = character(),
      Type = character(),
      Start_Position = integer(),
      End_Position = integer(),
      Length = integer(),
      Gene = character(),
      stringsAsFactors = FALSE
    )
  }
}

#' Process and Validate Region Info from PCG_all_cds (Python-Compatible)
#'
#' Main entry point that fully replicates Python script logic.
#'
#' **Processing Steps**:
#' 1. Load and clean raw data
#' 2. Extract CDS/exon/intron records with Gene column
#' 3. Remove CDS rows for multi-exon genes
#' 4. Calculate IGS total length using Python statistical method
#' 5. Generate individual IGS records from coordinate gaps
#' 6. Combine all records
#'
#' @param pcg_all_cds_path Path to the PCG_all_cds.csv file
#'
#' @return A data frame with columns: species, Region_Name, Type,
#'   Start_Position, End_Position, Length, Gene
#'
#' @export
process_and_validate_region_info_v2 <- function(pcg_all_cds_path) {
  
  message("=== Processing Region Info (Python-Compatible Mode) ===")
  message(sprintf("Input file: %s", pcg_all_cds_path))
  
  # Step 1: Load raw data
  message("\n[Step 1/6] Loading raw annotation data...")
  raw_data <- readr::read_tsv(pcg_all_cds_path, show_col_types = FALSE)
  
  # Validate required columns
  required_cols <- c("Species", "Name", "Gene", "Type", "Minimum", "Maximum", "Length", "Intervals")
  missing_cols <- setdiff(required_cols, names(raw_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Remove empty rows
  raw_data <- raw_data %>%
    dplyr::filter(!is.na(Species) & !is.na(Name) & !is.na(Type))
  
  # Ensure numeric types
  raw_data <- raw_data %>%
    dplyr::mutate(
      Minimum = as.numeric(Minimum),
      Maximum = as.numeric(Maximum),
      Length = as.numeric(Length),
      Intervals = as.numeric(Intervals)
    )
  
  # Remove rows with NA in critical columns
  raw_data <- raw_data %>%
    dplyr::filter(!is.na(Minimum) & !is.na(Maximum) & !is.na(Length) & !is.na(Intervals))
  
  message(sprintf("  Loaded %d records from %d species", nrow(raw_data), length(unique(raw_data$Species))))
  
  # Step 2: Extract CDS/exon/intron records with Gene column
  message("\n[Step 2/6] Extracting CDS, exon, and intron records...")
  
  cds_exon_intron <- raw_data %>%
    dplyr::filter(Type %in% c("CDS", "exon", "intron")) %>%
    dplyr::mutate(
      Region_Name = Name,
      Start_Position = Minimum,
      End_Position = Maximum
    ) %>%
    dplyr::select(species = Species, Region_Name, Type, Start_Position, End_Position, Length, Gene) %>%
    # Force recalculate length to ensure consistency, overriding source data issues
    dplyr::mutate(Length = End_Position - Start_Position + 1)
  
  message(sprintf("  Extracted %d CDS, %d exon, %d intron records",
                  sum(cds_exon_intron$Type == "CDS"),
                  sum(cds_exon_intron$Type == "exon"),
                  sum(cds_exon_intron$Type == "intron")))
  
  # Step 3: Remove CDS rows for multi-exon genes
  message("\n[Step 3/6] Removing CDS rows for multi-exon genes...")
  
  # Identify multi-exon genes
  multi_exon_genes <- cds_exon_intron %>%
    dplyr::filter(Type == "exon") %>%
    dplyr::group_by(species, Gene) %>%
    dplyr::summarise(exon_count = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(exon_count > 0) %>%
    dplyr::select(species, Gene)
  
  # Remove CDS rows for these genes
  cds_exon_intron_filtered <- cds_exon_intron %>%
    dplyr::anti_join(
      multi_exon_genes %>% dplyr::mutate(Type = "CDS"),
      by = c("species", "Gene", "Type")
    )
  
  removed_cds_count <- nrow(cds_exon_intron) - nrow(cds_exon_intron_filtered)
  message(sprintf("  Removed %d CDS rows for multi-exon genes", removed_cds_count))
  
  # Remove duplicate records (e.g., duplicate intron annotations in source data)
  before_dedup <- nrow(cds_exon_intron_filtered)
  cds_exon_intron_filtered <- cds_exon_intron_filtered %>%
    dplyr::distinct(species, Region_Name, Type, Start_Position, End_Position, .keep_all = TRUE)
  dedup_count <- before_dedup - nrow(cds_exon_intron_filtered)
  
  if (dedup_count > 0) {
    message(sprintf("  Removed %d duplicate records", dedup_count))
  }
  
  message(sprintf("  Retained %d records (CDS for single-exon genes, all exon/intron)",
                  nrow(cds_exon_intron_filtered)))
  
  # Step 4: Calculate IGS total length using Python method
  message("\n[Step 4/6] Calculating IGS total length (Python statistical method)...")
  
  igs_summary <- calculate_igs_total_length(raw_data)
  
  message(sprintf("  Calculated IGS totals for %d species", nrow(igs_summary)))
  message("  Formula: intergenic_length = genome_length - overlaps - total_CDS - intron_total")
  
  # Step 5: Generate individual IGS records
  message("\n[Step 5/6] Generating individual IGS coordinate records...")
  
  igs_records <- generate_igs_from_gaps(cds_exon_intron_filtered)
  
  message(sprintf("  Generated %d IGS records", nrow(igs_records)))
  
  # Step 6: Combine all records
  message("\n[Step 6/6] Combining all records...")
  
  final_records <- dplyr::bind_rows(cds_exon_intron_filtered, igs_records) %>%
    dplyr::arrange(species, Start_Position) %>%
    # Final length recalculation for absolute consistency across all record types
    dplyr::mutate(Length = End_Position - Start_Position + 1)
  
  message(sprintf("[OK] Successfully processed %d total records", nrow(final_records)))
  message(sprintf("  - CDS/exon: %d", sum(final_records$Type %in% c("CDS", "exon"))))
  message(sprintf("  - intron: %d", sum(final_records$Type == "intron")))
  message(sprintf("  - IGS: %d", sum(final_records$Type == "IGS")))
  
  message("\n=== Processing Complete ===\n")
  
  return(final_records)
}

