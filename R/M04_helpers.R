############################################################
#### M04_helpers.R - M04 IGS Hotspot Analysis Helpers ####
############################################################
#
# Core helper functions for M04 position-based orthologous IGS (poiGS) analysis
# These functions implement the "poiGS" approach to make IGS regions comparable
# across species by using flanking gene pairs as positional markers
#
# Part of the M04 IGS Hotspot Analysis Module
#
############################################################

#' Identify Position-based Orthologous IGS (poiGS) regions
#' 
#' @title Identify Position-based Orthologous IGS regions
#' @description Converts species-specific IGS regions into cross-species comparable 
#' position-based orthologous IGS (poiGS) units by analyzing flanking gene pairs.
#' This function solves the fundamental problem that IGS sequences vary dramatically
#' across species and cannot be directly compared.
#' 
#' @param igs_info_df Data frame containing IGS region information with columns:
#'   - species: Species identifier (character)
#'   - Region_Name: IGS region name in format "IGS_geneA-geneB" (character)
#'   - Type: Region type, should be "IGS" (character)
#'   - Other columns are preserved but not used
#' 
#' @return Data frame with poiGS mapping containing columns:
#'   - species: Original species identifier
#'   - Region_Name: Original IGS region name  
#'   - poiGS_ID: Standardized position-based orthologous IGS identifier
#'   - gene_pair: Standardized gene pair (sorted alphabetically)
#'   - original_gene_pair: Original gene pair before sorting
#' 
#' @details
#' The poiGS approach works by:
#' 1. **Extracting flanking genes**: Parse Region_Name (format: IGS_geneA-geneB) 
#'    to identify the genes flanking each IGS region
#' 2. **Standardizing gene pairs**: Sort gene names alphabetically to ensure 
#'    IGS_psbA-trnH and IGS_trnH-psbA are treated as the same orthologous unit
#' 3. **Generating poiGS IDs**: Create identifiers such as
#'    `poiGS_geneA-geneB` for cross-species comparison
#' 
#' This transforms the analysis unit from species-specific IGS regions to 
#' position-based orthologous units that can be meaningfully compared across species.
#' 
#' @examples
#' \dontrun{
#' # Load IGS region information
#' igs_data <- read.csv("output_gene_intergenic_intron_pos_length.csv") %>%
#'   filter(Type == "IGS")
#' 
#' # Generate poiGS mapping
#' poigs_mapping <- identify_positional_orthologs(igs_data)
#' 
#' # Check results
#' head(poigs_mapping)
#' length(unique(poigs_mapping$poiGS_ID))  # Number of orthologous IGS units
#' }
#' 
#' @importFrom dplyr mutate select filter
#' @importFrom stringr str_extract str_replace
#' @export
identify_positional_orthologs <- function(igs_info_df) {
  
  # Validate input data
  if (is.null(igs_info_df) || nrow(igs_info_df) == 0) {
    stop("identify_positional_orthologs: Input igs_info_df is NULL or empty")
  }
  
  # Check required columns
  required_cols <- c("species", "Region_Name", "Type")
  missing_cols <- setdiff(required_cols, names(igs_info_df))
  if (length(missing_cols) > 0) {
    stop(sprintf("identify_positional_orthologs: Missing required columns: %s", 
                paste(missing_cols, collapse = ", ")))
  }
  
  # Filter for IGS regions only
  # V3.16.12: Use only "IGS" standard term
  igs_only <- igs_info_df %>%
    dplyr::filter(Type == "IGS")
  
  if (nrow(igs_only) == 0) {
    stop("identify_positional_orthologs: No intergenic regions found in input data")
  }
  
  # Extract gene pairs from Region_Name (format: IGS_geneA-geneB)
  # Use regex to extract the gene pair part after "IGS_"
  gene_pair_extracted <- igs_only %>%
    dplyr::mutate(
      # Extract everything after "IGS_" prefix
      gene_pair_raw = stringr::str_replace(Region_Name, "^IGS_", ""),
      
      # Split the gene pair and extract individual genes
      # Handle cases where gene names might contain underscores or dashes
      gene1 = sapply(strsplit(gene_pair_raw, "-"), function(x) x[1]),
      gene2 = sapply(strsplit(gene_pair_raw, "-"), function(x) x[2]),
      
      # Store original gene pair for reference
      original_gene_pair = gene_pair_raw
    )
  
  # Validate gene extraction
  invalid_extractions <- gene_pair_extracted %>%
    dplyr::filter(is.na(gene1) | is.na(gene2) | gene1 == "" | gene2 == "")
  
  if (nrow(invalid_extractions) > 0) {
    warning(sprintf("identify_positional_orthologs: Failed to extract gene pairs from %d regions. Check Region_Name format.", 
                   nrow(invalid_extractions)))
  }
  
  # Remove invalid extractions
  valid_extractions <- gene_pair_extracted %>%
    dplyr::filter(!is.na(gene1) & !is.na(gene2) & gene1 != "" & gene2 != "")
  
  if (nrow(valid_extractions) == 0) {
    stop("identify_positional_orthologs: No valid gene pairs could be extracted from Region_Name column")
  }
  
  # Standardize gene pairs by alphabetical sorting
  # This ensures IGS_psbA-trnH and IGS_trnH-psbA are treated as the same orthologous unit
  poigs_mapping <- valid_extractions %>%
    dplyr::mutate(
      # Sort genes alphabetically within each pair
      sorted_gene1 = ifelse(gene1 <= gene2, gene1, gene2),
      sorted_gene2 = ifelse(gene1 <= gene2, gene2, gene1),
      
      # Create standardized gene pair
      gene_pair = paste(sorted_gene1, sorted_gene2, sep = "-"),
      
      # Generate unique poiGS identifier
      poiGS_ID = paste0("poiGS_", gene_pair)
    ) %>%
    
    # Select and organize output columns
    dplyr::select(
      species,
      Region_Name,
      poiGS_ID,
      gene_pair,
      original_gene_pair
    )
  
  # Log summary statistics
  num_species <- length(unique(poigs_mapping$species))
  num_original_igs <- nrow(poigs_mapping)
  num_unique_poigs <- length(unique(poigs_mapping$poiGS_ID))
  
  cat(sprintf("poiGS mapping summary:\n"))
  cat(sprintf("  - Species: %d\n", num_species))
  cat(sprintf("  - Original IGS regions: %d\n", num_original_igs))
  cat(sprintf("  - Unique poiGS units: %d\n", num_unique_poigs))
  cat(sprintf("  - Average IGS per poiGS: %.2f\n", num_original_igs / num_unique_poigs))
  
  return(poigs_mapping)
}


#' Filter Core Shared Hotspots
#' 
#' @title Filter core shared hotspots from candidate hotspots
#' @description Filters candidate hotspots to identify those that are shared 
#' across a minimum number of species, representing core evolutionary hotspots
#' 
#' @param candidate_hotspots_data Data frame containing candidate hotspot information
#' @param feature_col Name of the column containing feature identifiers (e.g., "poiGS_ID")
#' @param species_col Name of the column containing species identifiers
#' @param min_species_coverage Minimum fraction of species that must contain the hotspot (0-1)
#' 
#' @return Data frame containing only core shared hotspots
#' 
#' @details
#' This function implements the second-round filtering for M04 analysis:
#' 1. **Species counting**: For each candidate hotspot, count how many species contain it
#' 2. **Coverage filtering**: Retain only hotspots present in >= min_species_coverage * total_species
#' 3. **Core hotspot identification**: These represent evolutionarily conserved high-variation regions
#' 
#' @examples
#' \dontrun{
#' # Filter for hotspots present in at least 50% of species
#' core_hotspots <- filter_core_hotspots(
#'   candidate_hotspots_data = candidate_results,
#'   feature_col = "poiGS_ID",
#'   species_col = "species", 
#'   min_species_coverage = 0.5
#' )
#' }
#' 
#' @importFrom dplyr group_by summarise filter
#' @export
filter_core_hotspots <- function(candidate_hotspots_data, 
                                feature_col = "poiGS_ID",
                                species_col = "species",
                                min_species_coverage = 0.5) {
  
  # Validate inputs
  if (is.null(candidate_hotspots_data) || nrow(candidate_hotspots_data) == 0) {
    warning("filter_core_hotspots: No candidate hotspots provided")
    return(data.frame())
  }
  
  if (!feature_col %in% names(candidate_hotspots_data)) {
    stop(sprintf("filter_core_hotspots: Feature column '%s' not found in data", feature_col))
  }
  
  if (!species_col %in% names(candidate_hotspots_data)) {
    stop(sprintf("filter_core_hotspots: Species column '%s' not found in data", species_col))
  }
  
  # Calculate total number of species
  total_species <- length(unique(candidate_hotspots_data[[species_col]]))
  min_species_count <- ceiling(total_species * min_species_coverage)
  
  cat(sprintf("Core hotspot filtering:\n"))
  cat(sprintf("  - Total species: %d\n", total_species))
  cat(sprintf("  - Min species coverage: %.1f%% (%d species)\n", 
             min_species_coverage * 100, min_species_count))
  
  # Count species per hotspot feature
  feature_species_counts <- candidate_hotspots_data %>%
    dplyr::group_by(!!rlang::sym(feature_col)) %>%
    dplyr::summarise(
      species_count = length(unique(!!rlang::sym(species_col))),
      .groups = "drop"
    )
  
  # Filter for core shared hotspots
  core_features <- feature_species_counts %>%
    dplyr::filter(species_count >= min_species_count)
  
  # Filter original data to include only core shared hotspots
  core_hotspots <- candidate_hotspots_data %>%
    dplyr::filter(!!rlang::sym(feature_col) %in% core_features[[feature_col]])
  
  cat(sprintf("  - Candidate hotspots: %d\n", length(unique(candidate_hotspots_data[[feature_col]]))))
  cat(sprintf("  - Core shared hotspots: %d\n", length(unique(core_hotspots[[feature_col]]))))
  
  return(core_hotspots)
}