#' Annotation Preprocessing Helpers
#'
#' @description
#' Build a complete `region_info_complete.csv` table from `Annotations.csv`.
#' The helpers cover rps12 trans-splicing, gene-overlap detection, and IGS
#' length filtering.
#'
#' @name annotation_preprocessing_helpers
#' @keywords internal
NULL


#' Assign Genome Region Based on Coordinates
#'
#' @description
#' Assign `genome_region` (`LSC`, `IRA`, `SSC`, or `IRB`) from coordinates.
#' Exon-level intervals are supported. Features that span multiple plastome
#' partitions are assigned to the partition with the largest overlap.
#'
#' @param start Region start coordinate
#' @param end Region end coordinate
#' @param species_name Species name
#' @param species_genome_data Species genome configuration data frame
#' @return Character string: "LSC", "IRA", "IRB", "SSC", "whole_genome", "unknown"
#' @keywords internal
assign_genome_region <- function(start, end, species_name, species_genome_data) {
  
  # Validate inputs
  if (is.null(species_genome_data) || nrow(species_genome_data) == 0) {
    return("unknown")
  }
  
  if (is.na(start) || is.na(end) || is.na(species_name)) {
    return("unknown")
  }
  
  # Get species configuration
  species_config <- species_genome_data[species_genome_data$species == species_name, ]
  
  if (nrow(species_config) == 0) {
    message(sprintf("  WARNING: No genome configuration found for species: %s", species_name))
    return("unknown")
  }
  
  # Check for special handling (IR-lacking genome)
  if (!is.null(species_config$special_handling) && 
      !is.na(species_config$special_handling) && 
      species_config$special_handling == "IR_lacking_genome") {
    return("whole_genome")
  }
  
  # Calculate overlap with each plastome partition (majority-overlap rule)
  
  # Helper function to calculate overlap between two ranges
  calculate_overlap <- function(r1_start, r1_end, r2_start, r2_end) {
    overlap_start <- max(r1_start, r2_start)
    overlap_end <- min(r1_end, r2_end)
    if (overlap_start <= overlap_end) {
      return(overlap_end - overlap_start + 1)
    } else {
      return(0)
    }
  }
  
  # Calculate overlap with each genome region
  length_in_lsc <- calculate_overlap(start, end, species_config$lsc_start, species_config$lsc_end)
  length_in_ssc <- calculate_overlap(start, end, species_config$ssc_start, species_config$ssc_end)
  
  # For IR regions, need to determine IRA vs IRB
  # IRA is before SSC, IRB is after SSC
  ira_start <- species_config$ir_start
  ira_end <- min(species_config$ir_end, species_config$ssc_start - 1)
  irb_start <- max(species_config$ir_start, species_config$ssc_end + 1)
  irb_end <- species_config$ir_end
  
  length_in_ira <- calculate_overlap(start, end, ira_start, ira_end)
  length_in_irb <- calculate_overlap(start, end, irb_start, irb_end)
  
  # Find the region with maximum overlap
  max_length <- max(length_in_lsc, length_in_ira, length_in_irb, length_in_ssc)
  
  # If no overlap with any region, return unknown
  if (max_length == 0) {
    message(sprintf("  WARNING: Region [%d-%d] has no overlap with any genome region in %s", 
                   start, end, species_name))
    return("unknown")
  }
  
  # Check if region spans multiple genome regions (for logging purposes)
  regions_with_overlap <- sum(c(length_in_lsc, length_in_ira, length_in_irb, length_in_ssc) > 0)
  if (regions_with_overlap > 1) {
    total_length <- end - start + 1
    message(sprintf("  INFO: Region [%d-%d] spans multiple genome regions in %s (LSC:%d, IRA:%d, SSC:%d, IRB:%d bp), assigned to region with max overlap", 
                   start, end, species_name, length_in_lsc, length_in_ira, length_in_ssc, length_in_irb))
  }
  
  # Return the region with maximum overlap
  if (length_in_lsc == max_length) {
    return("LSC")
  } else if (length_in_ira == max_length) {
    return("IRA")
  } else if (length_in_irb == max_length) {
    return("IRB")
  } else if (length_in_ssc == max_length) {
    return("SSC")
  } else {
    # This should never happen, but as a fallback
    return("unknown")
  }
}


#' Extract Gene Boundaries with Intelligent rps12 Handling
#'
#' @description
#' Extract gene boundaries with special handling for trans-spliced rps12:
#' - standard genes use CDS coordinates
#' - rps12 exons are grouped by gap size and intron annotation
#'
#' @param annotations_data Standardized annotation data
#' @param gap_threshold Gap size to determine trans-splicing (default: 5000bp)
#' @return Data frame with gene boundaries
#' @keywords internal
extract_gene_boundaries_smart <- function(annotations_data, 
                                          gap_threshold = 5000) {
  
  message("Extracting gene boundaries with intelligent rps12 handling...")
  
  # Collect CDS records for all genes
  all_cds <- annotations_data %>%
    dplyr::filter(Type == "CDS") %>%
    dplyr::select(Species, Gene, Minimum, Maximum, Intervals) %>%
    dplyr::distinct()
  
  # Separate rps12 from other genes
  rps12_genes <- all_cds %>% dplyr::filter(Gene == "rps12")
  other_genes <- all_cds %>% dplyr::filter(Gene != "rps12")
  
  # Standard genes: use CDS coordinates as boundaries
  other_boundaries <- other_genes %>%
    dplyr::mutate(
      gene_start = Minimum,
      gene_end = Maximum
    ) %>%
    dplyr::select(Species, Gene, gene_start, gene_end)
  
  # Handle rps12 genes separately
  rps12_boundaries <- list()
  
  for (i in seq_len(nrow(rps12_genes))) {
    species <- rps12_genes$Species[i]
    
    # Collect rps12 exons for this species, ordered by position
    rps12_exons <- annotations_data %>%
      dplyr::filter(Species == species, Gene == "rps12", Type == "exon") %>%
      dplyr::arrange(Minimum) %>%
      dplyr::select(Minimum, Maximum, Length)
    
    if (nrow(rps12_exons) == 0) {
      # Fallback: use CDS coordinates
      message(sprintf("  WARNING: %s rps12 has no exon records, using CDS boundary", species))
      rps12_boundaries[[length(rps12_boundaries) + 1]] <- data.frame(
        Species = species,
        Gene = "rps12",
        gene_start = rps12_genes$Minimum[i],
        gene_end = rps12_genes$Maximum[i],
        stringsAsFactors = FALSE
      )
      next
    }
    
    if (nrow(rps12_exons) == 1) {
      # Single exon: use it as the boundary
      message(sprintf("  INFO: %s rps12 has 1 exon - single boundary", species))
      rps12_boundaries[[length(rps12_boundaries) + 1]] <- data.frame(
        Species = species,
        Gene = "rps12_part1",
        gene_start = rps12_exons$Minimum[1],
        gene_end = rps12_exons$Maximum[1],
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Multiple exons: group by gap type
    # Step 1: compute gaps between adjacent exons
    exon_gaps <- data.frame(
      gap_index = 1:(nrow(rps12_exons) - 1),
      gap_start = rps12_exons$Maximum[1:(nrow(rps12_exons) - 1)] + 1,
      gap_end = rps12_exons$Minimum[2:nrow(rps12_exons)] - 1,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(gap_size = gap_end - gap_start + 1)
    
    # Step 2: check whether intron annotations cover these gaps
    rps12_introns <- annotations_data %>%
      dplyr::filter(Species == species, Gene == "rps12", Type == "intron")
    
    exon_gaps <- exon_gaps %>%
      dplyr::mutate(
        has_intron_annotation = FALSE
      )
    
    # Mark gaps that have a matching intron annotation
    if (nrow(rps12_introns) > 0) {
      for (j in 1:nrow(exon_gaps)) {
        gap <- exon_gaps[j, ]
        # Allow a +/- 50 bp tolerance when matching intron coordinates
        has_intron <- any(
          abs(rps12_introns$Minimum - gap$gap_start) <= 50 &
          abs(rps12_introns$Maximum - gap$gap_end) <= 50
        )
        exon_gaps$has_intron_annotation[j] <- has_intron
      }
    }
    
    # Step 3: treat large gaps or gaps without intron annotation as trans-splicing
    exon_gaps <- exon_gaps %>%
      dplyr::mutate(
        is_trans_splicing = (gap_size > gap_threshold) | !has_intron_annotation
      )
    
    # Step 4: split exons at trans-splicing gaps
    # Identify trans-splicing gap indices
    trans_splicing_gaps <- which(exon_gaps$is_trans_splicing)
    
    if (length(trans_splicing_gaps) == 0) {
      # No trans-splicing: treat all exons as one contiguous fragment
      message(sprintf("  INFO: %s rps12 has %d exons - all continuous (no trans-splicing)", 
                      species, nrow(rps12_exons)))
      rps12_boundaries[[length(rps12_boundaries) + 1]] <- data.frame(
        Species = species,
        Gene = "rps12",
        gene_start = rps12_exons$Minimum[1],
        gene_end = rps12_exons$Maximum[nrow(rps12_exons)],
        stringsAsFactors = FALSE
      )
    } else {
      # Trans-splicing present: split into parts
      # Split after the exon that precedes each trans-splicing gap
      split_points <- c(0, trans_splicing_gaps, nrow(rps12_exons))
      
      message(sprintf("  INFO: %s rps12 has %d exons with %d trans-splicing gaps - creating %d parts", 
                      species, nrow(rps12_exons), 
                      length(trans_splicing_gaps),
                      length(split_points) - 1))
      
      # Create a boundary record for each exon group
      for (part_idx in 1:(length(split_points) - 1)) {
        start_exon <- split_points[part_idx] + 1
        end_exon <- split_points[part_idx + 1]
        
        rps12_boundaries[[length(rps12_boundaries) + 1]] <- data.frame(
          Species = species,
          Gene = sprintf("rps12_part%d", part_idx),
          gene_start = rps12_exons$Minimum[start_exon],
          gene_end = rps12_exons$Maximum[end_exon],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Combine rps12 boundary records
  if (length(rps12_boundaries) > 0) {
    rps12_boundaries_df <- dplyr::bind_rows(rps12_boundaries)
  } else {
    rps12_boundaries_df <- data.frame(
      Species = character(),
      Gene = character(),
      gene_start = integer(),
      gene_end = integer(),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all gene boundaries
  all_boundaries <- dplyr::bind_rows(
    other_boundaries,
    rps12_boundaries_df
  ) %>%
    dplyr::arrange(Species, gene_start)
  
  message(sprintf("  Extracted %d gene boundaries (%d standard + %d rps12 parts)",
                 nrow(all_boundaries),
                 nrow(other_boundaries),
                 nrow(rps12_boundaries_df)))
  
  return(all_boundaries)
}


#' Detect Gene Overlaps
#'
#' @description
#' Detect overlapping gene boundaries.
#' The logic matches the Python helper `find_overlapping_intervals`.
#'
#' @param gene_boundaries Gene boundaries data frame
#' @return List with overlaps data frame and total overlap length
#' @keywords internal
detect_gene_overlaps <- function(gene_boundaries) {
  
  message("Detecting gene overlaps...")
  
  overlaps_list <- list()
  total_overlap_length <- 0
  
  for (sp in unique(gene_boundaries$Species)) {
    sp_genes <- gene_boundaries %>%
      dplyr::filter(Species == sp) %>%
      dplyr::arrange(gene_start)
    
    if (nrow(sp_genes) < 2) next
    
    for (i in 1:(nrow(sp_genes) - 1)) {
      current <- sp_genes[i, ]
      next_gene <- sp_genes[i + 1, ]
      
      # Overlap when the next gene starts at or before the current gene ends
      if (next_gene$gene_start <= current$gene_end) {
        overlap_start <- next_gene$gene_start
        overlap_end <- min(current$gene_end, next_gene$gene_end)
        overlap_length <- overlap_end - overlap_start + 1
        
        total_overlap_length <- total_overlap_length + overlap_length
        
        overlaps_list[[length(overlaps_list) + 1]] <- data.frame(
          Species = sp,
          gene1 = current$Gene,
          gene1_start = current$gene_start,
          gene1_end = current$gene_end,
          gene2 = next_gene$Gene,
          gene2_start = next_gene$gene_start,
          gene2_end = next_gene$gene_end,
          overlap_start = overlap_start,
          overlap_end = overlap_end,
          overlap_length = overlap_length,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(overlaps_list) > 0) {
    overlaps_df <- dplyr::bind_rows(overlaps_list)
    message(sprintf("  WARNING: Found %d gene overlaps, total overlap length: %d bp",
                   nrow(overlaps_df), total_overlap_length))
    
    # Summarise overlaps by species
    overlap_summary <- overlaps_df %>%
      dplyr::group_by(Species) %>%
      dplyr::summarise(
        overlap_count = dplyr::n(),
        total_overlap_bp = sum(overlap_length),
        .groups = "drop"
      )
    
    message("  Overlap summary by species:")
    print(overlap_summary)
  } else {
    overlaps_df <- data.frame()
    message("  No gene overlaps detected")
  }
  
  return(list(
    overlaps = overlaps_df,
    total_overlap_length = total_overlap_length
  ))
}


#' Generate IGS Records with Overlap and Length Filtering
#'
#' @description
#' Build IGS records from gene boundaries. Overlapping gene intervals
#' produce no IGS, intervals shorter than `min_igs_length` are dropped,
#' and each retained IGS is assigned a `genome_region`.
#'
#' @param gene_boundaries Gene boundaries data frame
#' @param min_igs_length Minimum IGS length to keep (default: 3bp)
#' @param species_genome_data Species genome configuration (optional)
#' @return List with IGS records and filtering statistics
#' @keywords internal
generate_igs_with_filtering <- function(gene_boundaries, 
                                        min_igs_length = 3,
                                        species_genome_data = NULL) {
  
  message(sprintf("Generating IGS records (minimum length: %d bp)...", min_igs_length))
  
  igs_records <- list()
  filter_stats <- list(
    total_gaps = 0,
    overlaps = 0,
    too_short = 0,
    valid_igs = 0,
    boundary_igs = 0  # track boundary IGS (currently unused)
  )
  
  for (sp in unique(gene_boundaries$Species)) {
    sp_genes <- gene_boundaries %>%
      dplyr::filter(Species == sp) %>%
      dplyr::arrange(gene_start)
    
    if (nrow(sp_genes) < 1) next  # at least one gene is required
    
    # Do not emit a start-of-genome IGS. For circular genomes with truncated
    # IR regions, the interval from genome start to the first gene is not a
    # complete spacer and would create false-positive hotspot calls.
    
    # Process inter-gene IGS (only if 2+ genes exist)
    if (nrow(sp_genes) < 2) next
    
    for (i in 1:(nrow(sp_genes) - 1)) {
      current_gene <- sp_genes[i, ]
      next_gene <- sp_genes[i + 1, ]
      
      filter_stats$total_gaps <- filter_stats$total_gaps + 1
      
      # Compute the intergenic gap
      igs_start <- current_gene$gene_end + 1
      igs_end <- next_gene$gene_start - 1
      igs_length <- igs_end - igs_start + 1
      
      # Filter 1: overlapping genes (igs_length <= 0)
      if (igs_length <= 0) {
        filter_stats$overlaps <- filter_stats$overlaps + 1
        next  # skip overlapping intervals
      }
      
      # Filter 2: IGS shorter than min_igs_length
      if (igs_length < min_igs_length) {
        filter_stats$too_short <- filter_stats$too_short + 1
        message(sprintf("  WARNING: %s: Very short IGS between %s and %s (%d bp) - filtered out",
                       sp, current_gene$Gene, next_gene$Gene, igs_length))
        next  # skip short IGS intervals
      }
      
      # Passed filters: create an IGS record
      filter_stats$valid_igs <- filter_stats$valid_igs + 1
      
      # Build IGS names, stripping rps12_partN suffixes
      gene1_name <- sub("_part[0-9]+$", "", current_gene$Gene)
      gene2_name <- sub("_part[0-9]+$", "", next_gene$Gene)
      
      # Assign genome_region from coordinates
      if (!is.null(species_genome_data)) {
        genome_region <- assign_genome_region(igs_start, igs_end, sp, species_genome_data)
      } else {
        genome_region <- "unknown"
      }
      
      igs_records[[length(igs_records) + 1]] <- data.frame(
        species = sp,
        region_name = sprintf("IGS_%s-%s", gene1_name, gene2_name),
        region_type = "IGS",
        gene = NA_character_,
        region_start = igs_start,
        region_end = igs_end,
        region_length = igs_length,
        genome_region = genome_region,
        stringsAsFactors = FALSE
      )
    }
    
    # Do not emit an end-of-genome IGS. The last-gene-to-genome-end interval
    # can also be incomplete on circular genomes with truncated IR regions.
  }
  
  # Report filtering statistics
  message("\n=== IGS Filtering Statistics ===")
  message(sprintf("  Total gene gaps:        %d", filter_stats$total_gaps))
  message(sprintf("  Boundary IGS:           %d (start + end)", filter_stats$boundary_igs))
  message(sprintf("  - Gene overlaps:        %d (filtered out)", filter_stats$overlaps))
  message(sprintf("  - Too short (<%d bp):   %d (filtered out)", 
                 min_igs_length, filter_stats$too_short))
  message(sprintf("  - Valid IGS generated:  %d (including %d boundary IGS)", 
                 filter_stats$valid_igs, filter_stats$boundary_igs))
  
  if (length(igs_records) > 0) {
    igs_df <- dplyr::bind_rows(igs_records)
    return(list(
      igs = igs_df,
      stats = filter_stats
    ))
  } else {
    return(list(
      igs = data.frame(
        species = character(),
        region_name = character(),
        region_type = character(),
        gene = character(),
        region_start = integer(),
        region_end = integer(),
        region_length = integer(),
        genome_region = character(),
        stringsAsFactors = FALSE
      ),
      stats = filter_stats
    ))
  }
}


#' Preprocess Annotations to Complete Region Info
#'
#' @description
#' Complete annotation preprocessing, including:
#' - rps12 trans-splicing handling
#' - gene-overlap detection
#' - IGS length filtering
#' - a summary statistics report
#'
#' @param annotations_path Path to validated Annotations.csv
#' @param species_genome_data Optional species genome configuration data frame
#' @param output_path Output path for region_info_complete.csv (optional)
#' @param min_igs_length Minimum IGS length (default: 3bp)
#' @param gap_threshold Gap size for trans-splicing detection (default: 5000bp)
#' @return List with region_info data and statistics
#' @export
preprocess_annotations_to_region_info <- function(annotations_path, 
                                                  species_genome_data = NULL,
                                                  output_path = NULL,
                                                  min_igs_length = 3,
                                                  gap_threshold = 5000) {
  
  message("=======================================================")
  message("  PREPROCESSING: Annotations -> Region Info Complete  ")
  message("=======================================================")
  message(sprintf("Input: %s", annotations_path))
  message(sprintf("IGS minimum length: %d bp", min_igs_length))
  message(sprintf("Trans-splicing gap threshold: %d bp", gap_threshold))
  
  # Validate species_genome_data if provided
  if (!is.null(species_genome_data)) {
    message(sprintf("Genome region assignment: ENABLED (%d species configurations)", 
                   nrow(species_genome_data)))
  } else {
    message("Genome region assignment: DISABLED (no species_genome_data provided)")
  }
  message("")
  
  # Step 1: load and standardise annotation data
  message("[1/7] Loading and standardizing annotation data...")
  annotations <- load_and_standardize_annotation_data(annotations_path)
  message(sprintf("  Loaded %d records from %d species", 
                 nrow(annotations), 
                 length(unique(annotations$Species))))
  
  # Step 2: extract gene boundaries, with rps12 handling
  message("\n[2/7] Extracting gene boundaries with intelligent rps12 handling...")
  gene_boundaries <- extract_gene_boundaries_smart(annotations, gap_threshold)
  
  # Step 3: detect gene overlaps
  message("\n[3/7] Detecting gene overlaps...")
  overlap_result <- detect_gene_overlaps(gene_boundaries)
  
  # Step 4: generate filtered IGS records
  message("\n[4/7] Generating IGS records with filtering...")
  igs_result <- generate_igs_with_filtering(gene_boundaries, min_igs_length, species_genome_data)
  
  # Step 5: extract CDS records at exon resolution and assign genome_region
  message("\n[5/7] Extracting CDS records with exon-level splitting and genome_region assignment...")
  
  # CDS extraction strategy:
  # 1. Single-exon genes (Intervals=1): use CDS rows; keep independent copies
  # 2. Multi-exon genes (Intervals>1): use exon rows; keep independent copies
  # 3. Truncated coordinates: classify using the visible coordinate range
  
  # Collect CDS summary rows, keeping coordinates so duplicate copies can be split
  cds_summary <- annotations %>%
    dplyr::filter(Type == "CDS") %>%
    dplyr::select(Species, Gene, Intervals, Minimum, Maximum, Length) %>%
    dplyr::distinct()
  
  # Detect independent copies by grouping on species and gene
  cds_summary <- cds_summary %>%
    dplyr::group_by(Species, Gene) %>%
    dplyr::mutate(
      copy_count = dplyr::n(),
      copy_index = dplyr::row_number()
    ) %>%
    dplyr::ungroup()
  
  # Initialise the CDS record list
  cds_records_list <- list()
  
  # Process each gene copy
  for (i in 1:nrow(cds_summary)) {
    species_name <- cds_summary$Species[i]
    gene_name <- cds_summary$Gene[i]
    intervals <- cds_summary$Intervals[i]
    copy_count <- cds_summary$copy_count[i]
    copy_index <- cds_summary$copy_index[i]
    cds_minimum <- cds_summary$Minimum[i]
    cds_maximum <- cds_summary$Maximum[i]
    
    # Build a gene-name prefix that distinguishes independent copies
    if (copy_count > 1) {
      gene_prefix <- sprintf("%s_copy%d", gene_name, copy_index)
      if (copy_index == 1) {
        message(sprintf("  INFO: %s %s: Detected %d independent copies", 
                       species_name, gene_name, copy_count))
      }
    } else {
      gene_prefix <- gene_name
    }
    
    if (intervals == 1) {
      # Single-exon gene: use the CDS row (match coordinates to split copies)
      cds_row <- annotations %>%
        dplyr::filter(Species == species_name, Gene == gene_name, Type == "CDS",
                     Minimum == cds_minimum, Maximum == cds_maximum)
      
      if (nrow(cds_row) > 0) {
        cds_row <- cds_row[1, ]
        
        # Handle truncated coordinates (strip ">" markers)
        coord_end <- cds_row$Maximum
        coord_length <- cds_row$Length
        is_truncated <- FALSE
        
        if (is.character(coord_end) && grepl(">", coord_end)) {
          is_truncated <- TRUE
          coord_end <- as.numeric(gsub(">", "", coord_end))
          coord_length <- as.numeric(gsub(">", "", coord_length))
          message(sprintf("  INFO: %s %s: Truncated gene detected [%d->%d], using visible coordinates", 
                         species_name, gene_prefix, cds_row$Minimum, coord_end))
        }
        
        # Assign genome_region; truncated genes use visible coordinates
        if (!is.null(species_genome_data)) {
          genome_region <- assign_genome_region(cds_row$Minimum, coord_end, 
                                                species_name, species_genome_data)
        } else {
          genome_region <- "unknown"
        }
        
        cds_records_list[[length(cds_records_list) + 1]] <- data.frame(
          species = species_name,
          region_name = gene_prefix,
          region_type = "CDS",
          gene = gene_name,
          region_start = cds_row$Minimum,
          region_end = coord_end,
          region_length = coord_length,
          genome_region = genome_region,
          stringsAsFactors = FALSE
        )
      }
    } else {
      # Multi-exon gene: use exon rows (match the copy coordinate range)
      exon_rows <- annotations %>%
        dplyr::filter(Species == species_name, Gene == gene_name, Type == "exon",
                     Minimum >= cds_minimum, Maximum <= cds_maximum) %>%
        dplyr::arrange(Minimum)
      
      if (nrow(exon_rows) > 0) {
        for (j in 1:nrow(exon_rows)) {
          exon_row <- exon_rows[j, ]
          
          # Name exons and include copy information when needed
          region_name <- sprintf("%s_exon%d", gene_prefix, j)
          
          # Assign genome_region
          if (!is.null(species_genome_data)) {
            genome_region <- assign_genome_region(exon_row$Minimum, exon_row$Maximum, 
                                                  species_name, species_genome_data)
          } else {
            genome_region <- "unknown"
          }
          
          cds_records_list[[length(cds_records_list) + 1]] <- data.frame(
            species = species_name,
            region_name = region_name,
            region_type = "CDS",
            gene = gene_name,
            region_start = exon_row$Minimum,
            region_end = exon_row$Maximum,
            region_length = exon_row$Length,
            genome_region = genome_region,
            stringsAsFactors = FALSE
          )
        }
      } else {
        # Fallback: if exon rows are missing, use the CDS summary row
        message(sprintf("  WARNING: Gene %s (%s) has Intervals=%d but no exon rows, using CDS summary", 
                       gene_prefix, species_name, intervals))
        cds_row <- annotations %>%
          dplyr::filter(Species == species_name, Gene == gene_name, Type == "CDS",
                       Minimum == cds_minimum, Maximum == cds_maximum)
        
        if (nrow(cds_row) > 0) {
          cds_row <- cds_row[1, ]
          
          # Handle truncated coordinates
          coord_end <- cds_row$Maximum
          coord_length <- cds_row$Length
          
          if (is.character(coord_end) && grepl(">", coord_end)) {
            coord_end <- as.numeric(gsub(">", "", coord_end))
            coord_length <- as.numeric(gsub(">", "", coord_length))
          }
          
          if (!is.null(species_genome_data)) {
            genome_region <- assign_genome_region(cds_row$Minimum, coord_end, 
                                                  species_name, species_genome_data)
          } else {
            genome_region <- "unknown"
          }
          
          cds_records_list[[length(cds_records_list) + 1]] <- data.frame(
            species = species_name,
            region_name = gene_prefix,
            region_type = "CDS",
            gene = gene_name,
            region_start = cds_row$Minimum,
            region_end = coord_end,
            region_length = coord_length,
            genome_region = genome_region,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  # Combine all CDS records
  if (length(cds_records_list) > 0) {
    cds_records <- dplyr::bind_rows(cds_records_list)
  } else {
    cds_records <- data.frame(
      species = character(),
      region_name = character(),
      region_type = character(),
      gene = character(),
      region_start = numeric(),
      region_end = numeric(),
      region_length = numeric(),
      genome_region = character(),
      stringsAsFactors = FALSE
    )
  }
  
  message(sprintf("  Extracted %d CDS records (exon-level)", nrow(cds_records)))
  if (!is.null(species_genome_data)) {
    genome_region_counts <- table(cds_records$genome_region)
    message(sprintf("  Genome region distribution: %s", 
                   paste(names(genome_region_counts), "=", genome_region_counts, collapse = ", ")))
  }
  
  # Step 6: extract intron records and assign genome_region
  message("\n[6/7] Extracting intron records with genome_region assignment...")
  intron_records <- annotations %>%
    dplyr::filter(Type == "intron") %>%
    dplyr::group_by(Species, Gene) %>%
    dplyr::arrange(Minimum) %>%
    dplyr::mutate(
      intron_number = dplyr::row_number()
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      region_name = sprintf("%s_intron_%d", Gene, intron_number),
      region_type = "intron",
      gene = Gene,
      region_start = Minimum,
      region_end = Maximum,
      region_length = Length
    ) %>%
    dplyr::select(species = Species, region_name, region_type, gene,
                  region_start, region_end, region_length) %>%
    dplyr::distinct()
  
  # Add genome_region
  if (!is.null(species_genome_data) && nrow(intron_records) > 0) {
    intron_records$genome_region <- sapply(1:nrow(intron_records), function(i) {
      assign_genome_region(
        intron_records$region_start[i],
        intron_records$region_end[i],
        intron_records$species[i],
        species_genome_data
      )
    })
  } else {
    intron_records$genome_region <- "unknown"
  }
  
  message(sprintf("  Extracted %d intron records", nrow(intron_records)))
  if (!is.null(species_genome_data) && nrow(intron_records) > 0) {
    intron_genome_counts <- table(intron_records$genome_region)
    message(sprintf("  Genome region distribution: %s", 
                   paste(names(intron_genome_counts), "=", intron_genome_counts, collapse = ", ")))
  }
  
  # Step 7: combine records and clean species names
  message("\n[7/7] Combining all region records and cleaning species names...")
  region_info_complete <- dplyr::bind_rows(
    cds_records,
    igs_result$igs,
    intron_records
  ) %>%
    # Strip trailing "(modified)" suffixes from species names
    dplyr::mutate(
      species = stringr::str_replace(species, "\\s*\\(modified\\)\\s*$", "")
    ) %>%
    dplyr::arrange(species, region_start)
  
  message(sprintf("  Total records: %d", nrow(region_info_complete)))
  message(sprintf("    - CDS:        %d", nrow(cds_records)))
  message(sprintf("    - IGS:        %d", nrow(igs_result$igs)))
  message(sprintf("    - Intron:     %d", nrow(intron_records)))
  message("  Cleaned '(modified)' suffix from species names")
  
  # Save to file when an output path is provided
  if (!is.null(output_path)) {
    message(sprintf("\n[SAVE] Saving to: %s", output_path))
    readr::write_csv(region_info_complete, output_path)
    message("  Region info complete saved successfully!")
  }
  
  # Build a summary statistics report
  message("\n=======================================================")
  message("  PREPROCESSING STATISTICS                            ")
  message("=======================================================")
  
  stats <- region_info_complete %>%
    dplyr::group_by(species, region_type) %>%
    dplyr::summarise(
      count = dplyr::n(),
      total_length = sum(region_length, na.rm = TRUE),
      .groups = "drop"
    )
  print(stats)
  
  if (nrow(overlap_result$overlaps) > 0) {
    message("\nWARNING: Gene Overlaps Detected:")
    print(overlap_result$overlaps)
  }
  
  message(sprintf("\nIGS Filtering: %d overlaps + %d too short = %d filtered out",
                 igs_result$stats$overlaps,
                 igs_result$stats$too_short,
                 igs_result$stats$overlaps + igs_result$stats$too_short))
  
  message("\n=======================================================")
  message("  PREPROCESSING COMPLETE                              ")
  message("=======================================================")
  
  return(invisible(list(
    region_info = region_info_complete,
    gene_boundaries = gene_boundaries,
    overlaps = overlap_result$overlaps,
    igs_stats = igs_result$stats,
    statistics = stats
  )))
}


#' Load and Standardize Annotation Data
#'
#' @description
#' Load annotation data and standardise column names for Geneious Prime
#' and legacy exports. Shared with `pre_validate_region_info_v2()`.
#'
#' @param input_path Path to annotation CSV file
#' @return Standardized data frame
#' @keywords internal
load_and_standardize_annotation_data <- function(input_path) {
  
  # Detect file format
  if (grepl("\\.txt$", input_path)) {
    raw_data <- readr::read_tsv(input_path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  } else {
    raw_data <- readr::read_csv(input_path, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  }
  
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
                   grepl("tRNA", raw_data$Type, ignore.case = TRUE), na.rm = TRUE)
  
  raw_data <- raw_data %>%
    dplyr::filter(!grepl("^trn", Gene, ignore.case = TRUE)) %>%
    dplyr::filter(!grepl("tRNA", Type, ignore.case = TRUE))
  
  # Convert numeric columns
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

