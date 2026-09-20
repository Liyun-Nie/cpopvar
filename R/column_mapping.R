############################################################
#### column_mapping.R - Enhanced Column Mapping System ####
############################################################
#
# Comprehensive column mapping utilities with multi-file support
# Provides detailed field descriptions and intelligent mapping
# Based on DATA_REQUIREMENTS_SPEC.md standards
#
############################################################

# Load required libraries

#' Get field descriptions for interactive help
#' 
#' @param file_type Type of file
#' @return List of field descriptions
#' get_field_descriptions
#' @export
get_field_descriptions <- function(file_type = "main_data") {
  
  if (file_type == "main_data") {
    descriptions <- list(
      sample_id = list(
        title = "Sample unique identifier",
        description = "Unique identifier for each sample, used for sample identification and counting",
        examples = "CB579883, DQ073553, JN897307",
        data_type = "character",
        required = TRUE,
        notes = "Must ensure uniqueness, no duplicate values allowed"
      ),
      species = list(
        title = "Species name",
        description = "Species scientific name, used for data grouping and species-specific processing",
        examples = "Glycine_max, Glycine_soja, Pisum_sativum",
        data_type = "character",
        required = TRUE,
        notes = "Must use underscore-connected scientific name format"
      ),
      position = list(
        title = "Genome position coordinate",
        description = "Precise coordinate position of variant site in genome",
        examples = "10245, 85467, 121456",
        data_type = "integer",
        required = TRUE,
        notes = "Used for IR coordinate conversion and heatmap analysis, typically ranges 1-200000"
      ),
      variant_type = list(
        title = "Variant type",
        description = "Specific type of variant, affects filtering threshold setting",
        examples = "snp, ins, del, complex",
        data_type = "character",
        required = TRUE,
        notes = "ins and del will be automatically merged as INDEL for statistical analysis"
      ),
      ref_allele = list(
        title = "Reference allele",
        description = "Nucleotide sequence at this site in reference genome",
        examples = "A, T, C, G, AT, GC",
        data_type = "character",
        required = TRUE,
        notes = "Can be single nucleotide or short sequence"
      ),
      alt_allele = list(
        title = "Variant allele",
        description = "Nucleotide sequence at this site in variant sample",
        examples = "G, C, A, T, GC, AT",
        data_type = "character",
        required = TRUE,
        notes = "Variant sequence corresponding to reference allele"
      ),
      var_location = list(
        title = "Variant location type",
        description = "Functional region type of variant site in genome",
        examples = "CDS, intergenic, intron",
        data_type = "character",
        required = FALSE,
        notes = "Used for region-specific analysis, corresponds to V8 column in raw data"
      ),
      gene_name = list(
        title = "Gene name",
        description = "Gene identifier where CDS variant is located",
        examples = "matK, rbcL, atpF, ndhF",
        data_type = "character",
        required = FALSE,
        notes = "Only CDS variants have this information, used for gene-specific analysis, corresponds to V14 column in raw data"
      )
    )
  } else if (file_type == "cds_lengths") {
    descriptions <- list(
      gene_name = list(
        title = "Gene identifier",
        description = "Standard name or ID of gene",
        examples = "matK, rbcL, ndhF",
        data_type = "character",
        required = TRUE,
        notes = "Must match gene names in main data file"
      ),
      gene_length = list(
        title = "Gene length (bp)",
        description = "Base pair length of gene, used for frequency normalization",
        examples = "1530, 1428, 2337",
        data_type = "integer",
        required = TRUE,
        notes = "Used to calculate variants per kilobase frequency = variant_count / gene_length * 1000"
      ),
      species = list(
        title = "Species name",
        description = "Species-specific gene length(optional)",
        examples = "Glycine_max, Glycine_soja",
        data_type = "character",
        required = FALSE,
        notes = "If exists, this is species-specific length data"
      )
    )
  } else if (file_type == "region_info") {
    descriptions <- list(
      region_name = list(
        title = "Region identifier",
        description = "Unique identifier for non-coding regions",
        examples = "IGS_trnK-rps16, intron_atpF",
        data_type = "character",
        required = TRUE,
        notes = "Used to identify specific intergenic regions or introns"
      ),
      start_pos = list(
        title = "Region start position",
        description = "Start coordinate of region in genome",
        examples = "2570, 15234",
        data_type = "integer",
        required = TRUE,
        notes = "Must be less than end position"
      ),
      end_pos = list(
        title = "Region end position",
        description = "End coordinate of region in genome",
        examples = "3120, 16789",
        data_type = "integer",
        required = TRUE,
        notes = "Must be greater than start position"
      ),
      region_length = list(
        title = "Region length(bp)",
        description = "Base pair length of region",
        examples = "550, 1555",
        data_type = "integer",
        required = TRUE,
        notes = "Used for region length normalization analysis"
      ),
      region_type = list(
        title = "Region type",
        description = "Functional classification of region",
        examples = "intergenic, intron",
        data_type = "character",
        required = TRUE,
        notes = "Used for region type grouping analysis"
      )
    )
  } else if (file_type == "genome_regions") {
    descriptions <- list(
      species = list(
        title = "Species name",
        description = "Scientific name identifier of species",
        examples = "Glycine_max, Pisum_sativum",
        data_type = "character",
        required = TRUE,
        notes = "Must match species name in main data file"
      ),
      LSC_start = list(
        title = "LSCRegion start position",
        description = "Large single copy region(Large Single Copy)start coordinate",
        examples = "1, 1",
        data_type = "integer",
        required = TRUE,
        notes = "Usually starts from position 1 of genome"
      ),
      LSC_end = list(
        title = "LSCRegion end position",
        description = "Large single copy regionend coordinate",
        examples = "83142, 76778",
        data_type = "integer",
        required = TRUE,
        notes = "LSC region is usually the largest genome region"
      ),
      IR_start = list(
        title = "IRRegion start position",
        description = "Inverted repeat region(Inverted Repeat)start coordinate",
        examples = "83143, 76779",
        data_type = "integer",
        required = TRUE,
        notes = "IRA and IRB have same length, system will auto-convert"
      ),
      IR_end = list(
        title = "IRRegion end position",
        description = "Inverted repeat regionend coordinate",
        examples = "108874, 101510",
        data_type = "integer",
        required = TRUE,
        notes = "Used for IRA to IRB coordinate conversion calculation"
      ),
      SSC_start = list(
        title = "SSCRegion start position",
        description = "Small single copy region(Small Single Copy)start coordinate",
        examples = "108875, 101511",
        data_type = "integer",
        required = TRUE,
        notes = "SSC region is usually the smallest genome region"
      ),
      SSC_end = list(
        title = "SSCRegion end position",
        description = "Small single copy regionend coordinate",
        examples = "121581, 114247",
        data_type = "integer",
        required = TRUE,
        notes = "SSC end position is close to total genome length"
      ),
      special_handling = list(
        title = "Special handling marker",
        description = "Mark genome types requiring special handling",
        examples = "IR_lacking_genome, normal",
        data_type = "character",
        required = FALSE,
        notes = "Species like Pisum_sativum lacking IR regions need marking"
      )
    )
  } else if (file_type == "group_info") {
    descriptions <- list(
      sample_id = list(
        title = "Sample identifier",
        description = "Sample ID corresponding to sample_id in main data file",
        examples = "CB579883, DQ073553",
        data_type = "character",
        required = TRUE,
        notes = "Must exactly match sample ID in main data file"
      ),
      group_name = list(
        title = "Group name",
        description = "Main group category of samples",
        examples = "Cultivated, Wild, Landrace",
        data_type = "character",
        required = TRUE,
        notes = "Used for group comparison analysis"
      ),
      subgroup = list(
        title = "Subgroup",
        description = "More detailed group information(optional)",
        examples = "Asia, America, Europe",
        data_type = "character",
        required = FALSE,
        notes = "Can be used for more detailed stratified analysis"
      )
    )
  }
  
  return(descriptions)
}

#' Generate smart column mapping suggestions with file type support
#' 
#' @param column_names Vector of column names from user's data
#' @param file_type Type of file ("main_data", "cds_lengths", "region_info", "genome_regions", "group_info")
#' @return List of mapping suggestions for each expected column
generate_smart_mapping <- function(column_names, file_type = "main_data") {
  
  # Define comprehensive mapping rules for different file types
  if (file_type == "main_data") {
    mapping_rules <- list(
      # Core required fields (original V1-V6)
      sample_id = c("sample_id", "sample", "id", "sample", "sample_id", "sample_name", "accession", "V1", "sampleid"),
      species = c("species", "organism", "species_name", "species", "species_name", "species_id", "spp", "V2", "taxon"),
      position = c("position", "pos", "location", "position", "coordinate", "coordinate", "chr_pos", "genome_pos", "V3"),
      variant_type = c("variant_type", "var_type", "mutation_type", "Variant type", "type", "variant", "variant", "V4", "vartype"),
      ref_allele = c("ref_allele", "ref", "reference", "Reference allele", "reference_allele", "V5", "ref_seq"),
      alt_allele = c("alt_allele", "alt", "alternative", "Variant allele", "alternative_allele", "V6", "alt_seq"),
      
      # Extended fields (original V7-V15)
      var_location = c("var_location", "location_type", "variant position", "genomic_location", "annotation", "V8", "feature"),
      gene_name = c("gene_name", "gene", "gene_name", "gene_id", "locus", "V14", "geneid")
    )
  } else if (file_type == "cds_lengths") {
    mapping_rules <- list(
      gene_name = c("gene_name", "gene", "gene_name", "gene_id", "locus", "name", "geneid"),
      gene_length = c("gene_length", "length", "Gene length", "size", "bp", "basepairs", "len"),
      species = c("species", "organism", "species_name", "species", "species_name", "taxon")
    )
  } else if (file_type == "region_info") {
    mapping_rules <- list(
      region_name = c("region_name", "region", "region_name", "name", "region_id", "feature"),
      start_pos = c("start_pos", "start", "start position", "start_position", "begin", "from"),
      end_pos = c("end_pos", "end", "end position", "end_position", "finish", "to"),
      region_length = c("region_length", "length", "Region length", "size", "len"),
      region_type = c("region_type", "type", "Region type", "category", "feature_type")
    )
  } else if (file_type == "genome_regions") {
    mapping_rules <- list(
      species = c("species", "organism", "species_name", "species", "species_name", "taxon"),
      LSC_start = c("LSC_start", "lsc_start", "LSC start", "large_single_copy_start"),
      LSC_end = c("LSC_end", "lsc_end", "LSC end", "large_single_copy_end"),
      IR_start = c("IR_start", "ir_start", "IR start", "inverted_repeat_start"),
      IR_end = c("IR_end", "ir_end", "IR end", "inverted_repeat_end"),
      SSC_start = c("SSC_start", "ssc_start", "SSC start", "small_single_copy_start"),
      SSC_end = c("SSC_end", "ssc_end", "SSC end", "small_single_copy_end"),
      special_handling = c("special_handling", "special", "special handling", "notes", "comment")
    )
  } else if (file_type == "group_info") {
    mapping_rules <- list(
      sample_id = c("sample_id", "sample", "id", "sample", "sample_id", "sample_name", "accession", "sampleid"),
      group_name = c("group_name", "group", "group", "category", "class", "type"),
      subgroup = c("subgroup", "sub_group", "Subgroup", "subcategory", "subclass")
    )
  } else {
    stop(paste("Unknown file_type:", file_type))
  }
  
  # Initialize suggested mapping
  suggested_mapping <- list()
  
  # For each expected column, find the best match
  for (standard_col in names(mapping_rules)) {
    possible_names <- mapping_rules[[standard_col]]
    
    # Normalize column names for better matching (handle underscores and dots)
    normalized_user_cols <- tolower(gsub("[._-]", "", column_names))
    normalized_possible <- tolower(gsub("[._-]", "", possible_names))
    
    # First try exact match (case-insensitive, normalized)
    exact_match_idx <- which(normalized_user_cols %in% normalized_possible)
    
    if (length(exact_match_idx) > 0) {
      suggested_mapping[[standard_col]] <- column_names[exact_match_idx[1]]
    } else {
      # Try partial match for keywords
      partial_matches <- c()
      for (keyword in possible_names) {
        normalized_keyword <- tolower(gsub("[._-]", "", keyword))
        matches_idx <- which(grepl(normalized_keyword, normalized_user_cols, fixed = TRUE))
        if (length(matches_idx) > 0) {
          partial_matches <- c(partial_matches, column_names[matches_idx[1]])
        }
      }
      
      if (length(partial_matches) > 0) {
        suggested_mapping[[standard_col]] <- partial_matches[1]
      } else {
        suggested_mapping[[standard_col]] <- ""
      }
    }
  }
  
  log_message(paste("Smart mapping suggestions generated for", file_type, "with", length(column_names), "columns"))
  
  return(suggested_mapping)
}

#' Validate column mapping completeness with file type support
#' 
#' @param mapping_list List of user-selected column mappings
#' @param file_type Type of file
#' @return List with validation results
validate_column_mapping <- function(mapping_list, file_type = "main_data") {
  
  # Define required columns by file type
  if (file_type == "main_data") {
    required_columns <- c("species", "variant_type", "position", "sample_id", "ref_allele", "alt_allele")
  } else if (file_type == "cds_lengths") {
    required_columns <- c("gene_name", "gene_length")
  } else if (file_type == "region_info") {
    required_columns <- c("region_name", "start_pos", "end_pos", "region_length", "region_type")
  } else if (file_type == "genome_regions") {
    required_columns <- c("species", "LSC_start", "LSC_end", "IR_start", "IR_end", "SSC_start", "SSC_end")
  } else if (file_type == "group_info") {
    required_columns <- c("sample_id", "group_name")
  } else {
    required_columns <- c()
  }
  
  # Check for missing required mappings
  missing_required <- c()
  for (required_col in required_columns) {
    if (is.null(mapping_list[[required_col]]) || mapping_list[[required_col]] == "") {
      missing_required <- c(missing_required, required_col)
    }
  }
  
  # Check for duplicate mappings
  mapped_columns <- unlist(mapping_list)
  mapped_columns <- mapped_columns[mapped_columns != ""]
  duplicate_mappings <- mapped_columns[duplicated(mapped_columns)]
  
  # Generate validation result
  validation_result <- list(
    is_valid = length(missing_required) == 0,
    missing_required = missing_required,
    duplicate_mappings = unique(duplicate_mappings),
    warnings = c(),
    errors = c()
  )
  
  # Generate error messages
  if (length(missing_required) > 0) {
    descriptions <- get_field_descriptions(file_type)
    missing_titles <- sapply(missing_required, function(col) {
      desc <- descriptions[[col]]
      if (!is.null(desc)) desc$title else col
    })
    validation_result$errors <- c(validation_result$errors, 
                                  paste("Please complete mapping for required fields:", paste(missing_titles, collapse = ", ")))
  }
  
  # Generate warning messages
  if (length(duplicate_mappings) > 0) {
    validation_result$warnings <- c(validation_result$warnings,
                                    paste("Following columns are mapped multiple times:", paste(duplicate_mappings, collapse = ", ")))
  }
  
  return(validation_result)
}

#' Apply column mapping to data frame
#' 
#' @param data Original data frame
#' @param mapping_list List of column mappings (standard_name -> original_name)
#' @param file_type Type of file for validation
#' @return Data frame with renamed columns
apply_column_mapping <- function(data, mapping_list, file_type = "main_data") {
  
  if (is.null(data) || nrow(data) == 0) {
    log_message("No data available for column mapping", level = "warning")
    return(data)
  }
  
  mapped_data <- data
  successful_mappings <- c()
  failed_mappings <- c()
  
  # Apply each mapping
  for (standard_name in names(mapping_list)) {
    original_name <- mapping_list[[standard_name]]
    
    # Skip empty mappings
    if (is.null(original_name) || original_name == "") {
      next
    }
    
    # Check if original column exists
    if (original_name %in% names(mapped_data)) {
      # Rename the column
      names(mapped_data)[names(mapped_data) == original_name] <- standard_name
      successful_mappings <- c(successful_mappings, paste(original_name, "->", standard_name))
    } else {
      failed_mappings <- c(failed_mappings, original_name)
    }
  }
  
  # Log results
  if (length(successful_mappings) > 0) {
    log_message(paste("Column mapping applied for", file_type, ":", paste(successful_mappings, collapse = ", ")))
  }
  
  if (length(failed_mappings) > 0) {
    log_message(paste("Column mapping failed for", file_type, ":", paste(failed_mappings, collapse = ", ")), level = "warning")
  }
  
  return(mapped_data)
}

#' Generate enhanced data quality summary
#' 
#' @param data Data frame to analyze
#' @param file_type Type of file for context-specific analysis
#' @return Character string with data quality summary
generate_data_quality_summary <- function(data, file_type = "main_data") {
  
  if (is.null(data) || nrow(data) == 0) {
    return("No data available for analysis")
  }
  
  # Basic statistics
  n_rows <- nrow(data)
  n_cols <- ncol(data)
  
  # Missing values analysis
  missing_counts <- sapply(data, function(x) sum(is.na(x) | x == ""))
  total_missing <- sum(missing_counts)
  missing_percentage <- round(total_missing / (n_rows * n_cols) * 100, 2)
  
  # Identify columns with high missing rates
  high_missing_cols <- names(missing_counts[missing_counts > n_rows * 0.5])
  
  # File-type specific analysis
  specific_summary <- ""
  if (file_type == "main_data" && "species" %in% names(data)) {
    species_count <- length(unique(data$species[!is.na(data$species) & data$species != ""]))
    specific_summary <- paste("Detected", species_count, "species")
  } else if (file_type == "cds_lengths" && "gene_name" %in% names(data)) {
    gene_count <- length(unique(data$gene_name[!is.na(data$gene_name) & data$gene_name != ""]))
    specific_summary <- paste("Contains", gene_count, "gene length information")
  }
  
  # Generate summary
  summary_lines <- c(
    paste("Data dimensions:", n_rows, "rows x", n_cols, "columns"),
    paste("Missing values:", total_missing, "items (", missing_percentage, "%)"),
    paste("Column names:", paste(names(data), collapse = ", "))
  )
  
  if (specific_summary != "") {
    summary_lines <- c(summary_lines, specific_summary)
  }
  
  if (length(high_missing_cols) > 0) {
    summary_lines <- c(summary_lines, 
                       paste("High missing rate columns (>50%):", paste(high_missing_cols, collapse = ", ")))
  }
  
  return(paste(summary_lines, collapse = "\n"))
}

#' Create sample data for testing different file types
#' 
#' @param n_rows Number of rows to generate
#' @param file_type Type of file to generate
#' @return Data frame with sample data
create_sample_data <- function(n_rows = 100, file_type = "main_data") {
  
  if (file_type == "main_data") {
    sample_data <- data.frame(
      Sample_ID = paste0("S", 1:n_rows),
      Organism = sample(c("Glycine_max", "Glycine_soja", "Pisum_sativum"), n_rows, replace = TRUE),
      Location = sample(1:100000, n_rows, replace = TRUE),
      Mutation_Type = sample(c("snp", "ins", "del", "complex"), n_rows, replace = TRUE),
      Reference = sample(c("A", "T", "C", "G"), n_rows, replace = TRUE),
      Alternative = sample(c("A", "T", "C", "G"), n_rows, replace = TRUE),
      stringsAsFactors = FALSE
    )
  } else if (file_type == "cds_lengths") {
    genes <- c("matK", "rbcL", "atpF", "ndhF", "rpoC1", "rpoC2")
    sample_data <- data.frame(
      Gene = rep(genes, ceiling(n_rows/length(genes)))[1:n_rows],
      Length = sample(500:3000, n_rows, replace = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  return(sample_data)
}

# log_message("Enhanced column mapping utilities loaded successfully", level = "info")