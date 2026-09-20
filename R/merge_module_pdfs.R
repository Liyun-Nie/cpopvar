############################################################
#### merge_module_pdfs.R - Vectorized PDF Merger Tool ####
############################################################
#
# Merges individual vectorized PDF plots into module-based reports.
# Replaces the old PNG-based PDF generation system with a modern
# vectorized approach that provides publication-ready outputs.
#
# Key Features:
# - Module-based organization (M01, M02, M03)
# - Natural ordering of plots within modules
# - Error handling for missing modules
# - Publication-quality vector output
#
############################################################

# Load required libraries

#' Merge individual PDF plots into module-based reports
#' 
#' This function scans the session's plots directory for individual PDF files
#' and merges them by module (M01, M02, M03) into comprehensive reports.
#' 
#' @param session_path Character. Path to the session directory
#' @return List. Summary of created reports and any errors
#' 
#' @examples
#' \dontrun{
#' merge_pdfs_by_module("app_data/sessions/my_session_id")
#' }
#' @export
merge_pdfs_by_module <- function(session_path) {
  
  tryCatch({
    cat("[INFO] Starting vectorized PDF merger for session:", basename(session_path), "\n")
    
    # Validate session path
    if (!dir.exists(session_path)) {
      stop("Session path does not exist: ", session_path)
    }
    
    plots_dir <- file.path(session_path, "results", "plots")
    if (!dir.exists(plots_dir)) {
      stop("Plots directory does not exist: ", plots_dir)
    }
    
    results_dir <- file.path(session_path, "results")
    
    # Initialize results tracking
    merge_results <- list(
      created_reports = character(),
      errors = character(),
      total_pdfs_processed = 0
    )
    
    # Define module mapping with their display names
    modules <- list(
      "M01" = list(
        pattern = "M01_",
        subdir = "M01_distribution",
        report_name = "M01_vector_report.pdf",
        display_name = "Data Distribution Analysis"
      ),
      "M02" = list(
        pattern = "M02_",
        subdir = "M02_comparative", 
        report_name = "M02_vector_report.pdf",
        display_name = "Comparative Analysis"
      ),
      "M03" = list(
        pattern = "M03_",
        subdir = "M03_hotspot",
        report_name = "M03_vector_report.pdf",
        display_name = "Hotspot Analysis"
      )
    )
    
    # Process each module
    for (module_id in names(modules)) {
      module_info <- modules[[module_id]]
      
      tryCatch({
        cat(sprintf("[PROCESS] Processing %s (%s)...\n", module_id, module_info$display_name))
        
        # Find module directory
        module_dir <- file.path(plots_dir, module_info$subdir)
        if (!dir.exists(module_dir)) {
          warning(sprintf("Module directory not found: %s", module_dir))
          merge_results$errors <- c(merge_results$errors, 
                                  sprintf("%s: Directory not found", module_id))
          next
        }
        
        # Find all PDF files in module subdirectories
        pdf_files <- list.files(
          module_dir, 
          pattern = "\\.pdf$", 
          recursive = TRUE,
          full.names = TRUE
        )
        
        if (length(pdf_files) == 0) {
          warning(sprintf("No PDF files found in %s", module_dir))
          merge_results$errors <- c(merge_results$errors, 
                                  sprintf("%s: No PDF files found", module_id))
          next
        }
        
        # Sort files naturally to ensure logical order
        pdf_files <- pdf_files[order(basename(pdf_files))]
        
        cat(sprintf("   Found %d PDF files to merge\n", length(pdf_files)))
        
        # Validate PDF files before merging
        valid_pdfs <- character()
        for (pdf_file in pdf_files) {
          if (file.exists(pdf_file) && file.size(pdf_file) > 0) {
            # Quick validation - try to get page count
            tryCatch({
              page_count <- pdftools::pdf_length(pdf_file)
              if (page_count > 0) {
                valid_pdfs <- c(valid_pdfs, pdf_file)
                cat(sprintf("     [SUCCESS] %s (%d page%s)\n", 
                           basename(pdf_file), page_count, 
                           ifelse(page_count == 1, "", "s")))
              } else {
                warning(sprintf("PDF has no pages: %s", basename(pdf_file)))
              }
            }, error = function(e) {
              warning(sprintf("Invalid PDF file: %s - %s", basename(pdf_file), e$message))
            })
          } else {
            warning(sprintf("PDF file is empty or missing: %s", basename(pdf_file)))
          }
        }
        
        if (length(valid_pdfs) == 0) {
          warning(sprintf("No valid PDF files found for %s", module_id))
          merge_results$errors <- c(merge_results$errors, 
                                  sprintf("%s: No valid PDF files", module_id))
          next
        }
        
        # Create merged report
        output_path <- file.path(results_dir, module_info$report_name)
        
        cat(sprintf("   [MERGE] Merging %d PDFs into %s...\n", 
                   length(valid_pdfs), basename(output_path)))
        
        # Perform PDF merge
        pdftools::pdf_combine(valid_pdfs, output = output_path)
        
        # Verify the output
        if (file.exists(output_path) && file.size(output_path) > 0) {
          merged_pages <- pdftools::pdf_length(output_path)
          file_size <- round(file.size(output_path) / 1024 / 1024, 2)
          
          cat(sprintf("   [CREATED] Created: %s (%d pages, %.2f MB)\n", 
                     basename(output_path), merged_pages, file_size))
          
          merge_results$created_reports <- c(merge_results$created_reports, output_path)
          merge_results$total_pdfs_processed <- merge_results$total_pdfs_processed + length(valid_pdfs)
        } else {
          warning(sprintf("Failed to create merged PDF: %s", output_path))
          merge_results$errors <- c(merge_results$errors, 
                                  sprintf("%s: Merge failed", module_id))
        }
        
      }, error = function(e) {
        error_msg <- sprintf("%s processing failed: %s", module_id, e$message)
        warning(error_msg)
        merge_results$errors <- c(merge_results$errors, error_msg)
      })
    }
    
    # Generate summary
    cat("\n[SUMMARY] PDF Merger Summary:\n")
    cat(sprintf("   Reports created: %d\n", length(merge_results$created_reports)))
    cat(sprintf("   Individual PDFs processed: %d\n", merge_results$total_pdfs_processed))
    
    if (length(merge_results$errors) > 0) {
      cat(sprintf("   Errors encountered: %d\n", length(merge_results$errors)))
      for (error in merge_results$errors) {
        cat(sprintf("     [WARNING] %s\n", error))
      }
    }
    
    if (length(merge_results$created_reports) > 0) {
      cat("\n[REPORTS] Created Reports:\n")
      for (report in merge_results$created_reports) {
        cat(sprintf("   [SUCCESS] %s\n", basename(report)))
      }
    }
    
    cat("[SUCCESS] PDF merger completed successfully\n\n")
    
    return(merge_results)
    
  }, error = function(e) {
    error_msg <- sprintf("PDF merger failed: %s", e$message)
    cat(sprintf("[ERROR] %s\n", error_msg))
    return(list(
      created_reports = character(),
      errors = error_msg,
      total_pdfs_processed = 0
    ))
  })
}

#' Get list of available PDF reports for a session
#' 
#' @param session_path Character. Path to the session directory
#' @return Character vector. Paths to available PDF reports
get_available_pdf_reports <- function(session_path) {
  results_dir <- file.path(session_path, "results")
  
  if (!dir.exists(results_dir)) {
    return(character())
  }
  
  # Look for vector reports
  report_files <- list.files(
    results_dir,
    pattern = "_vector_report\\.pdf$",
    full.names = TRUE
  )
  
  # Filter to existing, non-empty files
  valid_reports <- report_files[file.exists(report_files) & file.size(report_files) > 0]
  
  return(valid_reports)
}

#' Create a manifest of all PDF files in the session
#' 
#' @param session_path Character. Path to the session directory
#' @return Data frame. Manifest of all PDF files with metadata
create_pdf_manifest <- function(session_path) {
  plots_dir <- file.path(session_path, "results", "plots")
  
  if (!dir.exists(plots_dir)) {
    return(data.frame())
  }
  
  # Find all PDF files
  all_pdfs <- list.files(
    plots_dir,
    pattern = "\\.pdf$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(all_pdfs) == 0) {
    return(data.frame())
  }
  
  # Create manifest
  manifest <- data.frame(
    file_path = all_pdfs,
    file_name = basename(all_pdfs),
    module = NA_character_,
    file_size_mb = NA_real_,
    pages = NA_integer_,
    created = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Extract module information and metadata
  for (i in seq_len(nrow(manifest))) {
    file_path <- manifest$file_path[i]
    file_name <- manifest$file_name[i]
    
    # Determine module
    if (grepl("^M01_", file_name)) {
      manifest$module[i] <- "M01"
    } else if (grepl("^M02_", file_name)) {
      manifest$module[i] <- "M02"
    } else if (grepl("^M03_", file_name)) {
      manifest$module[i] <- "M03"
    } else {
      manifest$module[i] <- "Other"
    }
    
    # Get file metadata
    if (file.exists(file_path)) {
      manifest$file_size_mb[i] <- round(file.size(file_path) / 1024 / 1024, 3)
      manifest$created[i] <- as.character(file.mtime(file_path))
      
      # Get page count
      tryCatch({
        manifest$pages[i] <- pdftools::pdf_length(file_path)
      }, error = function(e) {
        manifest$pages[i] <- NA_integer_
      })
    }
  }
  
  return(manifest)
}

# Export notification
# cat("[SUCCESS] Vectorized PDF merger system loaded successfully\n")
# cat("   Main function: merge_pdfs_by_module(session_path)\n")
# cat("   Features: Module-based merging, error handling, metadata tracking\n")