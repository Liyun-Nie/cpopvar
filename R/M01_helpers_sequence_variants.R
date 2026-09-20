# ============================================================
# M01 Sequence Variant Characterization Helper Functions
# ============================================================
#
# This file contains helper functions for M01-3: Sequence Variant Characterization.
# Migrated from M05c (CDS SNP annotation) and M05d (IGS variant features).

# ============================================================
# Part 1: M05c - CDS SNP Annotation Functions
# ============================================================

#' Extract Filtered CDS SNP Positions
#' @param filtered_data_file Path to P02_filtered/variants_filtered.csv
#' @return Data frame with species, position
#' @keywords internal
extract_filtered_cds_snp_positions <- function(filtered_data_file) {
  
  if (!file.exists(filtered_data_file)) {
    stop("Filtered data file not found: ", filtered_data_file)
  }
  
  whitelist <- readr::read_csv(filtered_data_file, show_col_types = FALSE) %>%
    dplyr::filter(region_type == "CDS", var_type == "snp") %>%
    dplyr::select(species, position) %>%
    dplyr::distinct()
  
  return(whitelist)
}

#' Extract Sample Level Annotations
#' @param preprocessed_data_file Path to P01 preprocessed data
#' @param whitelist Data frame with species and position
#' @return Data frame with sample annotations
#' @keywords internal
extract_sample_level_annotations <- function(preprocessed_data_file, whitelist) {
  
  if (!file.exists(preprocessed_data_file)) {
    stop("Preprocessed data file not found: ", preprocessed_data_file)
  }
  
  if (nrow(whitelist) == 0) {
    warning("Whitelist is empty")
    return(data.frame())
  }
  
  sample_annotations <- readr::read_csv(preprocessed_data_file, show_col_types = FALSE)
  
  # Apply column mapping for backward compatibility
  if ("V6" %in% colnames(sample_annotations)) {
    sample_annotations <- sample_annotations %>%
      dplyr::rename(alt_allele = V6, functional_annotation = V12)
  }
  
  # CRITICAL FIX: Filter for var_type == "snp" ONLY to prevent contamination
  # Issue: At mixed-type positions, different samples may have different var_types
  #        (e.g., some samples have SNPs, others have complex/INDEL variants)
  # Solution: Only include samples where var_type == "snp" to ensure SNP analysis integrity
  # Example: Pisum_sativum position 6737 has 5 SNP samples + 2 complex samples
  #          Without this filter, complex variants would be wrongly included in SNP analysis
  sample_annotations <- sample_annotations %>%
    dplyr::semi_join(whitelist, by = c("species", "position")) %>%
    dplyr::filter(var_type == "snp") %>%  # ADDED: strict SNP-only filter
    dplyr::select(sample_id, species, position, var_type, region_type, gene, 
                  alt_allele, functional_annotation)
  
  return(sample_annotations)
}

#' Parse Variant Type
#' @param sample_annotations Data frame with functional_annotation
#' @return Data frame with variant_class
#' @keywords internal
parse_variant_type_from_v12 <- function(sample_annotations) {
  
  parsed_annotations <- sample_annotations %>%
    dplyr::mutate(
      variant_class = dplyr::case_when(
        grepl("synonymous_variant", functional_annotation, ignore.case = TRUE) ~ "synonymous",
        grepl("missense_variant", functional_annotation, ignore.case = TRUE) ~ "nonsynonymous",
        TRUE ~ "other"
      )
    ) %>%
    dplyr::filter(variant_class %in% c("synonymous", "nonsynonymous"))
  
  return(parsed_annotations)
}

#' Classify SNP Sites Using Majority Rule
#' @param parsed_annotations Data frame with variant_class (output from parse_variant_type_from_v12)
#' @param majority_threshold Threshold for majority classification
#' @return Data frame with site classification
#' @details
#' This function uses ALLELIC polymorphism (based on actual base variants) 
#' as the primary criterion for determining polymorphism.
#' 
#' Classification logic:
#' - Case 1: Allelic monomorphic → Directly assigned to variant_class
#' - Case 2: Allelic polymorphic but functionally monomorphic → Assigned to variant_class
#' - Case 3: Both polymorphic, MVF > threshold → Assigned to majority type
#' - Case 4: Both polymorphic, 0.5 ≤ MVF ≤ threshold → Assigned to "mixed"
#' 
#' @keywords internal
classify_snp_sites_majority_rule <- function(parsed_annotations, majority_threshold) {
  
  if (nrow(parsed_annotations) == 0) {
    warning("No annotated variants to classify")
    return(data.frame())
  }
  
  if (missing(majority_threshold)) {
    stop("majority_threshold parameter is required (must be provided via config)")
  }
  
  if (majority_threshold < 0.5 || majority_threshold > 1) {
    stop("majority_threshold must be between 0.5 and 1")
  }
  
  # Step 1: Calculate allelic polymorphism (based on actual base variants)
  allelic_info <- parsed_annotations %>%
    dplyr::group_by(species, position) %>%
    dplyr::summarise(
      n_alleles = dplyr::n_distinct(alt_allele),
      alleles = paste(sort(unique(alt_allele)), collapse = "/"),
      is_allelic_polymorphic = n_alleles > 1,
      .groups = "drop"
    )
  
  # Step 2: Extract gene information for each site
  site_gene_info <- parsed_annotations %>%
    dplyr::group_by(species, position) %>%
    dplyr::summarise(
      gene = dplyr::first(gene),
      .groups = "drop"
    )
  
  # Step 3: Calculate functional-level statistics and polymorphism
  site_classification <- parsed_annotations %>%
    # Count samples for each variant class at each site
    dplyr::group_by(species, position, variant_class) %>%
    dplyr::summarise(n_samples = dplyr::n(), .groups = "drop") %>%
    # Calculate proportions within each site
    dplyr::group_by(species, position) %>%
    dplyr::mutate(
      total_samples_at_site = sum(n_samples),
      proportion = n_samples / total_samples_at_site,
      is_functional_polymorphic = dplyr::n() > 1
    ) %>%
    # Keep only the major variant type for each site
    dplyr::filter(proportion == max(proportion)) %>%
    # Handle 50/50 ties by keeping the first one
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  # Step 4: Merge allelic polymorphism information
  site_classification <- site_classification %>%
    dplyr::left_join(allelic_info, by = c("species", "position")) %>%
    dplyr::left_join(site_gene_info, by = c("species", "position"))
  
  # Step 5: Apply classification rules (using ALLELIC polymorphism as primary criterion)
  site_classification <- site_classification %>%
    dplyr::mutate(
      site_type = dplyr::case_when(
        # Case 1: Allelic monomorphic
        !is_allelic_polymorphic ~ variant_class,
        
        # Case 2: Allelic polymorphic but functionally monomorphic
        is_allelic_polymorphic & !is_functional_polymorphic ~ variant_class,
        
        # Case 3: Both polymorphic - apply Majority Rule
        is_allelic_polymorphic & is_functional_polymorphic & proportion > majority_threshold ~ variant_class,
        
        # Case 4: Both polymorphic but no clear majority
        is_allelic_polymorphic & is_functional_polymorphic & 
          proportion >= 0.5 & proportion <= majority_threshold ~ "mixed",
        
        # Default (should not happen)
        TRUE ~ "unknown"
      )
    )
  
  return(site_classification)
}

#' Calculate S/N Ratio
#' @param site_classification Data frame with site types
#' @param pseudocount Pseudocount for ratio calculation
#' @param min_sites_per_gene Minimum sites per gene
#' @return Data frame with S/N ratios
#' @keywords internal
calculate_sn_ratio <- function(site_classification, pseudocount, min_sites_per_gene) {
  
  if (nrow(site_classification) == 0) {
    warning("No site classification data to calculate S/N ratio")
    return(list(
      species_stats = data.frame(),
      gene_stats = data.frame()
    ))
  }
  
  # Validate required parameters
  if (missing(pseudocount)) {
    stop("pseudocount parameter is required (must be provided via config)")
  }
  
  if (missing(min_sites_per_gene)) {
    stop("min_sites_per_gene parameter is required (must be provided via config)")
  }
  
  # Species-level statistics
  species_stats <- site_classification %>%
    dplyr::group_by(species, site_type) %>%
    dplyr::summarise(site_count = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = site_type,
      values_from = site_count,
      values_fill = 0
    )
  
  # Ensure all required columns exist
  if (!"synonymous" %in% colnames(species_stats)) {
    species_stats$synonymous <- 0
  }
  if (!"nonsynonymous" %in% colnames(species_stats)) {
    species_stats$nonsynonymous <- 0
  }
  if (!"mixed" %in% colnames(species_stats)) {
    species_stats$mixed <- 0
  }
  
  species_stats <- species_stats %>%
    dplyr::mutate(
      total_sites = synonymous + nonsynonymous + mixed,
      # S/N Ratio (FIXED: aligned with gene-level, only denominator uses pseudocount)
      # This ensures: synonymous=2, nonsynonymous=0 → sn_ratio = 2/1 = 2
      sn_ratio = synonymous / (nonsynonymous + pseudocount),
      # Percentages (excluding mixed sites)
      synonymous_pct = dplyr::if_else(
        (synonymous + nonsynonymous) > 0,
        synonymous / (synonymous + nonsynonymous) * 100,
        NA_real_
      ),
      nonsynonymous_pct = dplyr::if_else(
        (synonymous + nonsynonymous) > 0,
        nonsynonymous / (synonymous + nonsynonymous) * 100,
        NA_real_
      )
    )
  
  # Add statistical test for species-level data
  if (nrow(species_stats) > 0) {
    species_stats <- species_stats %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        p_value = dplyr::if_else(
          (synonymous + nonsynonymous) > 0,
          stats::binom.test(synonymous, synonymous + nonsynonymous, p = 0.5)$p.value,
          NA_real_
        ),
        selection_pressure = dplyr::case_when(
          is.na(p_value) ~ "Insufficient Data",
          p_value < 0.05 & sn_ratio > 1 ~ "Purifying Selection",
          p_value < 0.05 & sn_ratio < 1 ~ "Positive Selection",
          TRUE ~ "Neutral"
        ),
        significance = dplyr::case_when(
          is.na(p_value) ~ "",
          p_value < 0.001 ~ "***",
          p_value < 0.01  ~ "**",
          p_value < 0.05  ~ "*",
          TRUE ~ ""
        )
      ) %>%
      dplyr::ungroup()
  }
  
  # Gene-level statistics
  if ("gene" %in% colnames(site_classification)) {
    gene_stats <- site_classification %>%
      dplyr::filter(!is.na(gene) & gene != "") %>%
      dplyr::group_by(species, gene, site_type) %>%
      dplyr::summarise(site_count = dplyr::n(), .groups = "drop") %>%
      tidyr::pivot_wider(
        names_from = site_type,
        values_from = site_count,
        values_fill = 0
      )
  } else {
    gene_stats <- data.frame()
  }
  
  # Ensure all required columns exist for gene stats
  if (nrow(gene_stats) > 0) {
    if (!"synonymous" %in% colnames(gene_stats)) {
      gene_stats$synonymous <- 0
    }
    if (!"nonsynonymous" %in% colnames(gene_stats)) {
      gene_stats$nonsynonymous <- 0
    }
    if (!"mixed" %in% colnames(gene_stats)) {
      gene_stats$mixed <- 0
    }
    
    gene_stats <- gene_stats %>%
      dplyr::mutate(
        total_sites = synonymous + nonsynonymous + mixed,
        # FIXED: Gene-level sn_ratio should NOT use pseudocount for both numerator and denominator
        # Use pseudocount ONLY for denominator to avoid division by zero
        # This ensures: synonymous=1, nonsynonymous=0 → sn_ratio = 1/1 = 1 (not 2/1 = 2)
        sn_ratio = synonymous / (nonsynonymous + pseudocount)
      ) %>%
      dplyr::filter(total_sites >= min_sites_per_gene)
  }
  
  return(list(
    species_stats = species_stats,
    gene_stats = gene_stats
  ))
}

# ============================================================
# Part 2: M05d - IGS Variant Features Functions
# ============================================================

#' Calculate IGS SNP/INDEL Ratio
#' @param filtered_data_file Path to variants_filtered.csv
#' @param pseudocount Pseudocount value
#' @return Data frame with SNP/INDEL ratios
#' @keywords internal
calculate_igs_snp_indel_ratio <- function(filtered_data_file, pseudocount) {
  
  if (missing(pseudocount)) {
    stop("pseudocount parameter is required")
  }
  
  if (!file.exists(filtered_data_file)) {
    stop("Filtered data file not found: ", filtered_data_file)
  }
  
  filtered_data <- readr::read_csv(filtered_data_file, show_col_types = FALSE)
  
  igs_variants <- filtered_data %>%
    dplyr::filter(region_type == "IGS", var_type %in% c("snp", "INDEL")) %>%
    dplyr::select(species, position, var_type) %>%
    dplyr::distinct()
  
  if (nrow(igs_variants) == 0) {
    warning("No IGS variants found")
    return(data.frame())
  }
  
  species_stats <- igs_variants %>%
    dplyr::group_by(species, var_type) %>%
    dplyr::summarise(n_sites = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = var_type, values_from = n_sites, values_fill = 0) %>%
    dplyr::mutate(
      snp_indel_ratio = (snp + pseudocount) / (INDEL + pseudocount),
      total_sites = snp + INDEL,
      snp_pct = snp / total_sites * 100,
      indel_pct = INDEL / total_sites * 100
    )
  
  return(species_stats)
}

#' Extract Filtered IGS INDEL Positions
#' @param filtered_data_file Path to variants_filtered.csv
#' @return Data frame with INDEL positions
#' @keywords internal
extract_filtered_igs_indel_positions <- function(filtered_data_file) {
  
  if (!file.exists(filtered_data_file)) {
    stop("Filtered data file not found: ", filtered_data_file)
  }
  
  whitelist <- readr::read_csv(filtered_data_file, show_col_types = FALSE) %>%
    dplyr::filter(region_type == "IGS", var_type == "INDEL") %>%
    dplyr::select(species, position) %>%
    dplyr::distinct()
  
  return(whitelist)
}

#' Extract INDEL Length Information
#' @param preprocessed_data_file Path to P01 data
#' @param whitelist INDEL position whitelist
#' @return Data frame with INDEL details
#' @keywords internal
extract_indel_length_info <- function(preprocessed_data_file, whitelist) {
  
  if (!file.exists(preprocessed_data_file)) {
    stop("Preprocessed data file not found: ", preprocessed_data_file)
  }
  
  if (nrow(whitelist) == 0) {
    warning("Whitelist is empty")
    return(data.frame())
  }
  
  sample_data <- readr::read_csv(preprocessed_data_file, show_col_types = FALSE)
  
  if ("V5" %in% colnames(sample_data)) {
    sample_data <- sample_data %>% dplyr::rename(ref_allele = V5, alt_allele = V6)
  }
  
  indel_data <- sample_data %>%
    dplyr::semi_join(whitelist, by = c("species", "position")) %>%
    dplyr::select(sample_id, species, position, var_type, region_type, gene, ref_allele, alt_allele) %>%
    dplyr::mutate(
      ref_length = nchar(ref_allele),
      alt_length = nchar(alt_allele),
      indel_length = abs(alt_length - ref_length),
      indel_type = dplyr::case_when(
        alt_length > ref_length ~ "insertion",
        alt_length < ref_length ~ "deletion",
        TRUE ~ "unknown"
      ),
      indel_sequence = dplyr::case_when(
        indel_type == "insertion" ~ stringr::str_sub(alt_allele, ref_length + 1),
        indel_type == "deletion" ~ stringr::str_sub(ref_allele, alt_length + 1),
        TRUE ~ ""
      )
    ) %>%
    dplyr::filter(indel_length > 0)
  
  return(indel_data)
}

#' Calculate INDEL Length Statistics
#' @param indel_data INDEL data frame
#' @param min_indel_length Minimum length threshold
#' @return List with statistics
#' @keywords internal
calculate_indel_length_stats <- function(indel_data, min_indel_length) {
  
  if (missing(min_indel_length)) {
    stop("min_indel_length parameter is required")
  }
  
  if (nrow(indel_data) == 0) {
    warning("No INDEL data available")
    return(list(
      species_stats_all = data.frame(),
      species_stats_filtered = data.frame(),
      position_stats = data.frame()
    ))
  }
  
  # All INDELs statistics
  species_stats_all <- indel_data %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(
      n_indels = dplyr::n(),
      mean_length = mean(indel_length, na.rm = TRUE),
      median_length = median(indel_length, na.rm = TRUE),
      max_length = max(indel_length, na.rm = TRUE),
      min_length = min(indel_length, na.rm = TRUE),
      sd_length = sd(indel_length, na.rm = TRUE),
      n_length_1 = sum(indel_length == 1),
      pct_length_1 = n_length_1 / dplyr::n() * 100,
      n_insertions = sum(indel_type == "insertion"),
      n_deletions = sum(indel_type == "deletion"),
      insertion_ratio = n_insertions / (n_insertions + n_deletions),
      .groups = "drop"
    )
  
  # Filtered INDELs (length >= threshold)
  indel_data_filtered <- indel_data %>% dplyr::filter(indel_length >= min_indel_length)
  
  if (nrow(indel_data_filtered) > 0) {
    species_stats_filtered <- indel_data_filtered %>%
      dplyr::group_by(species) %>%
      dplyr::summarise(
        n_indels = dplyr::n(),
        mean_length = mean(indel_length, na.rm = TRUE),
        median_length = median(indel_length, na.rm = TRUE),
        max_length = max(indel_length, na.rm = TRUE),
        min_length = min(indel_length, na.rm = TRUE),
        sd_length = sd(indel_length, na.rm = TRUE),
        n_insertions = sum(indel_type == "insertion"),
        n_deletions = sum(indel_type == "deletion"),
        insertion_ratio = n_insertions / (n_insertions + n_deletions),
        .groups = "drop"
      )
  } else {
    species_stats_filtered <- data.frame()
  }
  
  # Position-level statistics
  # Position-level statistics (for shared INDEL sequence analysis)
  # CRITICAL: Include polymorphism metrics (aligned with M05)
  position_stats <- indel_data %>%
    dplyr::group_by(species, position) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      n_unique_sequences = dplyr::n_distinct(indel_sequence),  # ADDED: unique INDEL sequences
      n_unique_lengths = dplyr::n_distinct(indel_length),      # ADDED: unique INDEL lengths
      mean_length = mean(indel_length, na.rm = TRUE),
      length_sd = sd(indel_length, na.rm = TRUE),
      is_polymorphic = n_unique_sequences > 1 | n_unique_lengths > 1,  # ADDED: polymorphism flag
      dominant_type = names(sort(table(indel_type), decreasing = TRUE))[1],
      .groups = "drop"
    )
  
  return(list(
    species_stats_all = species_stats_all,
    species_stats_filtered = species_stats_filtered,
    position_stats = position_stats
  ))
}

#' Identify Shared INDEL Sequences
#' @param indel_data Sample-level INDEL data
#' @param min_species Minimum species threshold
#' @return Data frame with shared sequences
#' @keywords internal
identify_shared_indel_sequences <- function(indel_data, min_species) {
  
  if (missing(min_species)) {
    stop("min_species parameter is required")
  }
  
  if (nrow(indel_data) == 0 || !"indel_sequence" %in% colnames(indel_data)) {
    warning("No INDEL data or missing indel_sequence column")
    return(data.frame())
  }
  
  shared_sequences <- indel_data %>%
    dplyr::group_by(indel_sequence, indel_length) %>%
    dplyr::summarise(
      n_species = dplyr::n_distinct(species),
      n_samples = dplyr::n(),
      species_list = paste(sort(unique(species)), collapse = "; "),
      mean_length_per_sample = mean(indel_length, na.rm = TRUE),
      n_positions_total = dplyr::n_distinct(paste(species, position)),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_species >= min_species) %>%
    dplyr::arrange(dplyr::desc(n_species), dplyr::desc(n_samples))
  
  return(shared_sequences)
}

#' Prepare UpSet Matrix by Sequence
#' @param indel_data Sample-level INDEL data
#' @return Binary matrix for UpSet plot
#' @keywords internal
prepare_upset_matrix_by_sequence <- function(indel_data) {
  
  if (nrow(indel_data) == 0) {
    warning("No INDEL data for UpSet matrix")
    return(matrix(nrow = 0, ncol = 0))
  }
  
  binary_matrix <- indel_data %>%
    dplyr::select(indel_sequence, species) %>%
    dplyr::distinct() %>%
    dplyr::mutate(present = 1) %>%
    tidyr::pivot_wider(names_from = species, values_from = present, values_fill = 0) %>%
    tibble::column_to_rownames("indel_sequence") %>%
    as.matrix()
  
  return(binary_matrix)
}

#' Cluster Species by Shared INDELs (Jaccard Similarity)
#' @param upset_matrix Binary matrix of sequences x species
#' @return List with hclust object and distance matrix
#' @keywords internal
cluster_species_by_shared_indels <- function(upset_matrix) {
  
  if (nrow(upset_matrix) == 0 || ncol(upset_matrix) == 0) {
    warning("Empty upset matrix")
    return(NULL)
  }
  
  n_species <- ncol(upset_matrix)
  
  if (n_species < 3) {
    warning("Need at least 3 species for clustering")
    return(NULL)
  }
  
  if (!all(upset_matrix %in% c(0, 1))) {
    warning("Upset matrix should contain only 0/1 values")
    return(NULL)
  }
  
  # Calculate Jaccard distance matrix
  species_names <- colnames(upset_matrix)
  dist_matrix <- matrix(0, nrow = n_species, ncol = n_species)
  rownames(dist_matrix) <- species_names
  colnames(dist_matrix) <- species_names
  
  pairwise_stats <- data.frame()
  
  for (i in 1:(n_species - 1)) {
    for (j in (i + 1):n_species) {
      sp_i <- upset_matrix[, i]
      sp_j <- upset_matrix[, j]
      
      intersection <- sum(sp_i & sp_j)
      union <- sum(sp_i | sp_j)
      
      jaccard_sim <- if (union == 0) 0 else intersection / union
      jaccard_dist <- 1 - jaccard_sim
      
      dist_matrix[i, j] <- jaccard_dist
      dist_matrix[j, i] <- jaccard_dist
      
      pairwise_stats <- rbind(pairwise_stats, data.frame(
        species_A = species_names[i],
        species_B = species_names[j],
        n_shared_sequences = intersection,
        n_sequences_A = sum(sp_i),
        n_sequences_B = sum(sp_j),
        n_union = union,
        jaccard_similarity = jaccard_sim,
        jaccard_distance = jaccard_dist,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  pairwise_stats <- pairwise_stats[order(-pairwise_stats$jaccard_similarity), ]
  
  # Hierarchical clustering (UPGMA)
  dist_obj <- as.dist(dist_matrix)
  hc <- hclust(dist_obj, method = "average")
  
  return(list(
    hc = hc,
    dist_matrix = dist_matrix,
    pairwise_stats = pairwise_stats
  ))
}

# ============================================================
# Part 3: M01-3e - Genome Region Statistics Functions (NEW)
# ============================================================

#' Calculate S/N Ratio by Genome Region
#' @description Calculate synonymous/nonsynonymous ratio for CDS SNPs grouped by genome region (LSC/SSC/IRA)
#' @param preprocessed_data_file Path to P01 preprocessed data
#' @param filtered_data_file Path to P02 filtered data
#' @param majority_rule_threshold Threshold for majority rule classification
#' @param pseudocount Pseudocount for ratio calculation (added to denominator only)
#' @return Data frame with region-level S/N statistics
#' @keywords internal
calculate_sn_ratio_by_genome_region <- function(
  preprocessed_data_file,
  filtered_data_file,
  majority_rule_threshold,
  pseudocount
) {
  
  if (!file.exists(preprocessed_data_file)) {
    stop("Preprocessed data file not found: ", preprocessed_data_file)
  }
  
  if (!file.exists(filtered_data_file)) {
    stop("Filtered data file not found: ", filtered_data_file)
  }
  
  # Step 1: Extract CDS SNP whitelist from P02 (with genome_region)
  whitelist <- readr::read_csv(filtered_data_file, show_col_types = FALSE) %>%
    dplyr::filter(region_type == "CDS", var_type == "snp") %>%
    dplyr::filter(genome_region %in% c("LSC", "SSC", "IRA", "IRB")) %>%
    dplyr::select(species, position, genome_region) %>%
    dplyr::distinct()
  
  if (nrow(whitelist) == 0) {
    warning("No CDS SNP positions found with valid genome regions")
    return(data.frame())
  }
  
  # Create a mapping from position to region and a simple position whitelist
  region_mapping <- whitelist %>%
    dplyr::group_by(species, position) %>%
    dplyr::summarise(genome_region = dplyr::first(genome_region), .groups = "drop")
  
  position_whitelist <- region_mapping %>% dplyr::select(species, position)
  
  # Step 2: Extract sample annotations from P01 using only species and position
  sample_annotations <- readr::read_csv(preprocessed_data_file, show_col_types = FALSE) %>%
    dplyr::semi_join(position_whitelist, by = c("species", "position")) %>%
    dplyr::filter(var_type == "snp") %>%
    dplyr::select(sample_id, species, position, gene, functional_annotation, alt_allele)
  
  if (nrow(sample_annotations) == 0) {
    warning("No sample annotations found for whitelisted CDS SNP positions")
    return(data.frame())
  }
  
  # Step 3: Parse variant types from functional_annotation (V12)
  parsed_data <- sample_annotations %>%
    dplyr::mutate(
      variant_class = dplyr::case_when(
        grepl("synonymous_variant", functional_annotation, ignore.case = TRUE) ~ "synonymous",
        grepl("missense_variant", functional_annotation, ignore.case = TRUE) ~ "nonsynonymous",
        TRUE ~ "other"
      )
    )
  
  # Step 4: Apply the SAME classification logic as cds_snp_annotation module
  # This ensures consistency with M01_cds_snp_site_classification.csv
  site_classification_no_region <- classify_snp_sites_majority_rule(
    parsed_annotations = parsed_data,
    majority_threshold = majority_rule_threshold
  )
  
  if (nrow(site_classification_no_region) == 0) {
    warning("Site classification was empty after applying majority rule.")
    return(data.frame())
  }
  
  # Step 5: Join genome_region back to the classified sites
  site_classification <- site_classification_no_region %>%
    dplyr::left_join(region_mapping, by = c("species", "position"))
  
  # Step 6: Calculate final region-level statistics
  region_stats <- site_classification %>%
    dplyr::group_by(species, genome_region, site_type) %>%
    dplyr::summarise(n_sites = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = site_type,
      values_from = n_sites,
      values_fill = 0
    )
  
  # Ensure all columns exist with default values
  if (!"synonymous" %in% names(region_stats)) {
    region_stats$synonymous <- 0
  }
  if (!"nonsynonymous" %in% names(region_stats)) {
    region_stats$nonsynonymous <- 0
  }
  if (!"mixed" %in% names(region_stats)) {
    region_stats$mixed <- 0
  }
  
  # Calculate ratios and percentages
  region_stats <- region_stats %>%
    dplyr::mutate(
      total_sites = synonymous + nonsynonymous + mixed,
      sn_ratio = synonymous / (nonsynonymous + pseudocount),
      syn_pct = synonymous / total_sites * 100,
      nonsyn_pct = nonsynonymous / total_sites * 100,
      mixed_pct = mixed / total_sites * 100
    )
  
  # Add binomial test for each (species, genome_region) combination
  # Test if observed S/N ratio deviates from neutral expectation (0.5)
  region_stats <- region_stats %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      binom_test_p = {
        if (synonymous + nonsynonymous > 0) {
          binom.test(synonymous, synonymous + nonsynonymous, p = 0.5, alternative = "two.sided")$p.value
        } else {
          NA_real_
        }
      },
      significance = dplyr::case_when(
        is.na(binom_test_p) ~ "",
        binom_test_p < 0.001 ~ "***",
        binom_test_p < 0.01 ~ "**",
        binom_test_p < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    dplyr::ungroup()
  
  # Return both aggregated stats and detailed site classification
  return(list(
    region_stats = region_stats,
    site_classification = site_classification %>% dplyr::select(species, position, gene, site_type, genome_region)
  ))
}

#' Calculate SNP/INDEL Ratio by Genome Region
#' @description Calculate SNP/INDEL ratio for all region types grouped by genome region (LSC/SSC/IRA)
#' @param filtered_data_file Path to P02 filtered data
#' @param pseudocount Pseudocount for ratio calculation
#' @return Data frame with region-level SNP/INDEL statistics
#' @keywords internal
calculate_snp_indel_ratio_by_genome_region <- function(
  filtered_data_file,
  pseudocount
) {
  
  if (!file.exists(filtered_data_file)) {
    stop("Filtered data file not found: ", filtered_data_file)
  }
  
  filtered_data <- readr::read_csv(filtered_data_file, show_col_types = FALSE)
  
  # Count by genome region (all region types, exclude whole_genome species)
  # Per user request (2025-11-01), this ratio should include all region_types (CDS, IGS, intron etc.)
  region_stats <- filtered_data %>%
    dplyr::filter(var_type %in% c("snp", "INDEL")) %>%
    dplyr::filter(genome_region %in% c("LSC", "SSC", "IRA", "IRB")) %>%  # Exclude whole_genome
    dplyr::select(species, position, genome_region, var_type) %>%
    dplyr::distinct() %>%
    dplyr::group_by(species, genome_region, var_type) %>%
    dplyr::summarise(n_sites = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = var_type,
      values_from = n_sites,
      values_fill = 0
    )
  
  # Ensure all columns exist with default values
  if (!"snp" %in% names(region_stats)) {
    region_stats$snp <- 0
  }
  if (!"INDEL" %in% names(region_stats)) {
    region_stats$INDEL <- 0
  }
  
  # Calculate ratios and percentages
  region_stats <- region_stats %>%
    dplyr::mutate(
      total_sites = snp + INDEL,
      snp_indel_ratio = (snp + pseudocount) / (INDEL + pseudocount),
      snp_pct = snp / total_sites * 100,
      indel_pct = INDEL / total_sites * 100
    )
  
  return(region_stats)
}

