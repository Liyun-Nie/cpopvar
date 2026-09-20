# ============================================================
# M01 Relational Analysis Helper Functions
# ============================================================
#
# This file contains helper functions for M01-2: Structural & Relational Analysis
# Specifically for poiGS-Gene correlation analysis.

#' Build poiGS Mapping from IGS Region Statistics
#' 
#' @description Constructs poiGS mapping by reading IGS region statistics file
#' @param session_id Session identifier
#' @param config Configuration object
#' @return Data frame with poiGS mapping
#' @keywords internal
build_poigs_mapping <- function(session_id, config) {
  
  log_message("M01: Building poiGS mapping from region_info_complete.csv")
  
  session_paths <- get_session_paths(session_id)
  
  # V3.16.14: Use new region_info_complete.csv from P01_preprocessed
  p01_dir <- get_processing_stage_path(session_id, "P01_preprocessed")
  region_info_file <- file.path(p01_dir, "region_info_complete.csv")
  
  if (!file.exists(region_info_file)) {
    stop(sprintf("M01: Cannot find region_info_complete.csv. Expected: %s", region_info_file))
  }
  
  log_message(sprintf("M01: Loading IGS region information from: %s", region_info_file))
  
  # Load and filter for IGS regions
  igs_region_info <- readr::read_csv(region_info_file, show_col_types = FALSE) %>%
    dplyr::filter(region_type == "IGS") %>%
    # Rename columns to match expected format for identify_positional_orthologs
    dplyr::rename(
      Region_Name = region_name,
      Type = region_type,
      Start_Position = region_start,
      End_Position = region_end,
      Length = region_length
    ) %>%
    dplyr::select(species, Region_Name, Type, Start_Position, End_Position, Length)
  
  if (nrow(igs_region_info) == 0) {
    stop("M01: No IGS regions found in region_info_complete.csv")
  }
  
  log_message(sprintf("M01: Loaded %d IGS regions across %d species", 
                     nrow(igs_region_info), length(unique(igs_region_info$species))))
  
  poigs_mapping <- identify_positional_orthologs(igs_region_info)
  
  log_message(sprintf("M01: Built poiGS mapping: %d IGS regions mapped to %d unique poiGS units",
                     nrow(poigs_mapping),
                     length(unique(poigs_mapping$poiGS_ID))))
  
  return(poigs_mapping)
}


#' Build poiGS Frequencies from Normalized Data
#' 
#' @description Aggregates variant frequencies by poiGS units
#' @param normalized_data Normalized frequency data (after data_scoping)
#' @param poigs_mapping poiGS mapping table
#' @param config Configuration object
#' @return List with frequencies and mapping
#' @keywords internal
build_poigs_frequencies <- function(normalized_data, poigs_mapping, config) {
  
  log_message("M01: Aggregating variant frequencies by poiGS units")
  
  igs_data <- normalized_data %>%
    dplyr::filter(region_type == "IGS")
  
  if (nrow(igs_data) == 0) {
    stop("M01: No IGS data found in normalized_data after filtering")
  }
  
  log_message(sprintf("M01: Filtered %d IGS frequency records", nrow(igs_data)))
  
  # V3.16.15: Debug - check column names before join
  log_message(sprintf("[DEBUG] igs_data columns: %s", paste(names(igs_data), collapse = ", ")))
  log_message(sprintf("[DEBUG] poigs_mapping columns: %s", paste(names(poigs_mapping), collapse = ", ")))
  
  # For IGS data, region_name holds the IGS identifier (e.g., "IGS_gene1-gene2")
  # The 'gene' column is NA for IGS regions. Match 'region_name' to 'Region_Name' in poigs_mapping
  frequency_with_poigs <- igs_data %>%
    dplyr::left_join(poigs_mapping, 
                    by = c("species" = "species", "region_name" = "Region_Name")) %>%
    dplyr::filter(!is.na(poiGS_ID))
  
  if (nrow(frequency_with_poigs) == 0) {
    stop("M01: No frequency data could be linked with poiGS mapping")
  }
  
  log_message(sprintf("M01: Successfully linked %d frequency observations", 
                     nrow(frequency_with_poigs)))
  
  # Stage 1: Aggregate by poiGS and var_type
  poigs_frequencies_stage1 <- frequency_with_poigs %>%
    dplyr::group_by(species, poiGS_ID, var_type, region_type) %>%
    dplyr::summarise(
      frequency_per_kb = mean(frequency_per_kb, na.rm = TRUE),
      variant_count = sum(variant_count, na.rm = TRUE),
      region_count = dplyr::n(),
      .groups = "drop"
    )
  
  # Stage 2: Merge all variant types
  poigs_frequencies <- poigs_frequencies_stage1 %>%
    dplyr::group_by(species, poiGS_ID, region_type) %>%
    dplyr::summarise(
      frequency_per_kb = sum(frequency_per_kb, na.rm = TRUE),
      variant_count = sum(variant_count, na.rm = TRUE),
      region_count = max(region_count),
      var_type = paste(sort(unique(var_type)), collapse = "+"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      var_type = ifelse(frequency_per_kb == 0, NA_character_, var_type)
    )
  
  log_message(sprintf("M01: Built poiGS frequencies: %d observations", nrow(poigs_frequencies)))
  
  return(list(
    frequencies = poigs_frequencies,
    mapping = poigs_mapping
  ))
}


#' Extract Flanking Genes from poiGS Identifier
#' 
#' @param poigs_id Character vector of poiGS identifiers
#' @return Data frame with gene1, gene2, gene_pair
#' @keywords internal
extract_flanking_genes_from_poigs <- function(poigs_id) {
  
  if (is.null(poigs_id) || length(poigs_id) == 0) {
    stop("Input poigs_id is NULL or empty")
  }
  
  gene_pairs <- stringr::str_replace(poigs_id, "^poiGS_", "")
  gene_list <- strsplit(gene_pairs, "-")
  
  result <- data.frame(
    poiGS_ID = poigs_id,
    gene1 = sapply(gene_list, function(x) if(length(x) >= 1) x[1] else NA_character_),
    gene2 = sapply(gene_list, function(x) if(length(x) >= 2) x[2] else NA_character_),
    gene_pair = gene_pairs,
    stringsAsFactors = FALSE
  )
  
  return(result)
}


#' Merge poiGS and Gene Frequency Data
#' 
#' @description Links each poiGS with average frequency of its flanking genes
#' @param poigs_freq Data frame with poiGS frequency data
#' @param gene_freq Data frame with gene frequency data
#' @param poigs_mapping Data frame with poiGS mapping
#' @param keep_all_zeros Whether to keep pairs where both values are zero.
#' @return Merged data frame
#' @keywords internal
merge_poigs_gene_frequencies <- function(poigs_freq, gene_freq, poigs_mapping, keep_all_zeros = TRUE) {
  
  if (is.null(poigs_freq) || nrow(poigs_freq) == 0) {
    stop("poigs_freq is NULL or empty")
  }
  
  if (is.null(gene_freq) || nrow(gene_freq) == 0) {
    stop("gene_freq is NULL or empty")
  }
  
  log_message("M01: Merging poiGS and gene frequency data")
  log_message(sprintf("  - Input poiGS records: %d", nrow(poigs_freq)))
  log_message(sprintf("  - Input gene records: %d", nrow(gene_freq)))
  
  # Extract gene pairs
  poigs_with_genes <- poigs_freq %>%
    dplyr::left_join(
      poigs_mapping %>% dplyr::select(species, poiGS_ID, gene_pair),
      by = c("species", "poiGS_ID")
    ) %>%
    tidyr::separate(gene_pair, into = c("gene1", "gene2"), sep = "-", remove = FALSE)
  
  # Prepare gene frequency lookup
  gene_freq_lookup <- gene_freq %>%
    dplyr::select(species, gene, gene_frequency = frequency_per_kb)
  
  # Merge with gene frequencies
  merged_with_gene1 <- poigs_with_genes %>%
    dplyr::left_join(gene_freq_lookup, by = c("species" = "species", "gene1" = "gene")) %>%
    dplyr::rename(gene1_freq = gene_frequency)
  
  merged_with_both_genes <- merged_with_gene1 %>%
    dplyr::left_join(gene_freq_lookup, by = c("species" = "species", "gene2" = "gene")) %>%
    dplyr::rename(gene2_freq = gene_frequency)
  
  # Calculate averages and categorize zero patterns
  merged_data <- merged_with_both_genes %>%
    dplyr::mutate(
      gene1_freq = ifelse(is.na(gene1_freq), 0, gene1_freq),
      gene2_freq = ifelse(is.na(gene2_freq), 0, gene2_freq),
      gene_freq_avg = (gene1_freq + gene2_freq) / 2
    ) %>%
    dplyr::rename(poiGS_freq = frequency_per_kb) %>%
    dplyr::mutate(
      zero_pattern = dplyr::case_when(
        poiGS_freq > 0 & gene_freq_avg > 0 ~ "neither_zero",
        poiGS_freq == 0 & gene_freq_avg == 0 ~ "both_zero",
        poiGS_freq > 0 & gene_freq_avg == 0 ~ "poiGS_only",
        poiGS_freq == 0 & gene_freq_avg > 0 ~ "gene_only",
        TRUE ~ "unknown"
      )
    ) %>%
    dplyr::select(species, poiGS_ID, poiGS_freq, gene1, gene2, gene1_freq, gene2_freq, 
                  gene_freq_avg, gene_pair, zero_pattern)
  
  # Filter based on keep_all_zeros
  if (!keep_all_zeros) {
    merged_data <- merged_data %>% dplyr::filter(zero_pattern != "both_zero")
    log_message("  - Filtered out [0,0] pairs (both_zero pattern)")
  } else {
    log_message("  - Kept all zero patterns including [0,0]")
  }
  
  log_message(sprintf("  - Final merged records: %d", nrow(merged_data)))
  log_message(sprintf("  - Zero pattern distribution:"))
  pattern_summary <- merged_data %>% dplyr::count(zero_pattern)
  for (i in 1:nrow(pattern_summary)) {
    log_message(sprintf("    * %s: %d (%.1f%%)", 
                       pattern_summary$zero_pattern[i], 
                       pattern_summary$n[i],
                       pattern_summary$n[i] / nrow(merged_data) * 100))
  }
  
  return(merged_data)
}


#' Calculate Within-Species Correlation
#' 
#' @description Calculates Pearson and Spearman correlation between poiGS and gene frequencies per species
#' @param merged_data Merged poiGS-Gene data frame from merge_poigs_gene_frequencies()
#' @param min_pairs_per_species Minimum number of pairs required per species (default: 5)
#' @return List containing correlation_results, merged_data, and summary_stats
#' @keywords internal
calculate_within_species_correlation <- function(merged_data, 
                                                min_pairs_per_species = 5) {
  
  # Validate input
  if (is.null(merged_data) || nrow(merged_data) == 0) {
    stop("calculate_within_species_correlation: merged_data is NULL or empty")
  }
  
  required_cols <- c("species", "poiGS_freq", "gene_freq_avg")
  missing_cols <- setdiff(required_cols, names(merged_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("calculate_within_species_correlation: Missing columns: %s", 
                paste(missing_cols, collapse = ", ")))
  }
  
  log_message("M01 Correlation: Calculating within-species correlations")
  log_message(sprintf("  - Minimum pairs per species threshold: %d", min_pairs_per_species))
  
  # Calculate correlation for each species
  correlation_results <- merged_data %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      n_unique_poigs = length(unique(poiGS_ID)),
      n_unique_gene1 = length(unique(gene1)),
      n_unique_gene2 = length(unique(gene2)),
      mean_poigs_freq = mean(poiGS_freq, na.rm = TRUE),
      mean_gene_freq = mean(gene_freq_avg, na.rm = TRUE),
      sd_poigs_freq = sd(poiGS_freq, na.rm = TRUE),
      sd_gene_freq = sd(gene_freq_avg, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Calculate correlation coefficients for species with sufficient data
  correlation_stats <- lapply(unique(merged_data$species), function(sp) {
    sp_data <- merged_data %>% dplyr::filter(species == sp)
    
    if (nrow(sp_data) < min_pairs_per_species) {
      return(data.frame(
        species = sp,
        pearson_r = NA_real_,
        pearson_p = NA_real_,
        spearman_r = NA_real_,
        spearman_p = NA_real_,
        sufficient_data = FALSE
      ))
    }
    
    # Pearson correlation
    pearson_test <- tryCatch({
      cor.test(sp_data$poiGS_freq, sp_data$gene_freq_avg, method = "pearson")
    }, error = function(e) {
      list(estimate = NA_real_, p.value = NA_real_)
    })
    
    # Spearman correlation
    spearman_test <- tryCatch({
      cor.test(sp_data$poiGS_freq, sp_data$gene_freq_avg, method = "spearman")
    }, error = function(e) {
      list(estimate = NA_real_, p.value = NA_real_)
    })
    
    data.frame(
      species = sp,
      pearson_r = as.numeric(pearson_test$estimate),
      pearson_p = pearson_test$p.value,
      spearman_r = as.numeric(spearman_test$estimate),
      spearman_p = spearman_test$p.value,
      sufficient_data = TRUE
    )
  })
  
  correlation_stats_df <- dplyr::bind_rows(correlation_stats)
  
  # Merge correlation stats with summary stats
  final_results <- correlation_results %>%
    dplyr::left_join(correlation_stats_df, by = "species")
  
  # Calculate overall summary
  valid_species <- final_results %>% dplyr::filter(sufficient_data == TRUE)
  
  summary_stats <- list(
    total_species = nrow(final_results),
    species_with_sufficient_data = nrow(valid_species),
    species_excluded = nrow(final_results) - nrow(valid_species),
    total_pairs = sum(final_results$n_pairs),
    mean_pairs_per_species = mean(final_results$n_pairs),
    mean_pearson_r = mean(valid_species$pearson_r, na.rm = TRUE),
    mean_spearman_r = mean(valid_species$spearman_r, na.rm = TRUE),
    significant_pearson_count = sum(valid_species$pearson_p < 0.05, na.rm = TRUE),
    significant_spearman_count = sum(valid_species$spearman_p < 0.05, na.rm = TRUE)
  )
  
  log_message(sprintf("  - Total species analyzed: %d", summary_stats$total_species))
  log_message(sprintf("  - Species with sufficient data (n >= %d): %d", 
                     min_pairs_per_species, summary_stats$species_with_sufficient_data))
  log_message(sprintf("  - Species excluded (insufficient data): %d", 
                     summary_stats$species_excluded))
  log_message(sprintf("  - Mean Pearson r: %.3f", summary_stats$mean_pearson_r))
  log_message(sprintf("  - Mean Spearman r: %.3f", summary_stats$mean_spearman_r))
  log_message(sprintf("  - Significant correlations (p < 0.05): Pearson=%d, Spearman=%d",
                     summary_stats$significant_pearson_count,
                     summary_stats$significant_spearman_count))
  
  return(list(
    correlation_results = final_results,
    merged_data = merged_data,
    summary_stats = summary_stats
  ))
}


#' Apply FDR Correction to Correlation P-values
#' 
#' @description Applies False Discovery Rate (FDR) correction using Benjamini-Hochberg
#' method to correlation p-values
#' @param correlation_results Data frame with correlation results
#' @param fdr_method FDR correction method (default: "BH")
#' @return Data frame with FDR-corrected columns
#' @keywords internal
apply_fdr_correction <- function(correlation_results, fdr_method = "BH") {
  
  if (is.null(correlation_results) || nrow(correlation_results) == 0) {
    warning("apply_fdr_correction: correlation_results is NULL or empty")
    return(correlation_results)
  }
  
  log_message(sprintf("M01 Correlation: Applying FDR correction (method: %s)", fdr_method))
  
  # Apply FDR correction to Pearson p-values
  if ("pearson_p" %in% names(correlation_results)) {
    valid_pearson <- !is.na(correlation_results$pearson_p)
    correlation_results$pearson_p_fdr <- NA_real_
    
    if (sum(valid_pearson) > 0) {
      correlation_results$pearson_p_fdr[valid_pearson] <- 
        p.adjust(correlation_results$pearson_p[valid_pearson], method = fdr_method)
      
      correlation_results$pearson_significant <- 
        !is.na(correlation_results$pearson_p_fdr) & 
        correlation_results$pearson_p_fdr < 0.05
    }
  }
  
  # Apply FDR correction to Spearman p-values
  if ("spearman_p" %in% names(correlation_results)) {
    valid_spearman <- !is.na(correlation_results$spearman_p)
    correlation_results$spearman_p_fdr <- NA_real_
    
    if (sum(valid_spearman) > 0) {
      correlation_results$spearman_p_fdr[valid_spearman] <- 
        p.adjust(correlation_results$spearman_p[valid_spearman], method = fdr_method)
      
      correlation_results$spearman_significant <- 
        !is.na(correlation_results$spearman_p_fdr) & 
        correlation_results$spearman_p_fdr < 0.05
    }
  }
  
  # Log FDR correction results
  if ("pearson_significant" %in% names(correlation_results)) {
    n_sig_pearson <- sum(correlation_results$pearson_significant, na.rm = TRUE)
    log_message(sprintf("  - Significant Pearson correlations after FDR: %d", n_sig_pearson))
  }
  
  if ("spearman_significant" %in% names(correlation_results)) {
    n_sig_spearman <- sum(correlation_results$spearman_significant, na.rm = TRUE)
    log_message(sprintf("  - Significant Spearman correlations after FDR: %d", n_sig_spearman))
  }
  
  return(correlation_results)
}


#' Generate Correlation Summary Statistics
#' 
#' @description Creates summary statistics for correlation results
#' @param correlation_results Data frame with correlation results
#' @param merged_data Original merged data
#' @return List with summary statistics
#' @keywords internal
generate_correlation_summary <- function(correlation_results, merged_data) {
  
  if (nrow(correlation_results) == 0) {
    return(list(
      overall = "No correlation results available",
      by_method = data.frame(),
      zero_patterns = data.frame()
    ))
  }
  
  # Overall summary
  overall_summary <- correlation_results %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      n_species = dplyr::n(),
      mean_correlation = mean(correlation, na.rm = TRUE),
      median_correlation = median(correlation, na.rm = TRUE),
      n_significant = sum(significant, na.rm = TRUE),
      n_sig_fdr = sum(sig_fdr, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Zero pattern summary
  zero_summary <- merged_data %>%
    dplyr::group_by(species, zero_pattern) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = zero_pattern, values_from = count, values_fill = 0)
  
  return(list(
    overall = overall_summary,
    by_species = correlation_results,
    zero_patterns = zero_summary
  ))
}

