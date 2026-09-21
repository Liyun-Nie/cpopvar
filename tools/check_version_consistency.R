#!/usr/bin/env Rscript

check_version_consistency <- function(pkg_root = ".") {
  pkg_root <- normalizePath(pkg_root, winslash = "/", mustWork = TRUE)
  problems <- character(0)

  add_problem <- function(...) {
    problems <<- c(problems, paste0(...))
  }

  desc_path <- file.path(pkg_root, "DESCRIPTION")
  if (!file.exists(desc_path)) {
    stop("DESCRIPTION not found in ", pkg_root, call. = FALSE)
  }
  desc <- read.dcf(desc_path)
  version <- as.character(desc[1, "Version"])
  if (!nzchar(version)) {
    add_problem("DESCRIPTION Version is missing.")
  }

  version_parts <- strsplit(version, ".", fixed = TRUE)[[1]]
  is_dev_version <- length(version_parts) == 4L
  is_release_version <- length(version_parts) == 3L
  if (!is_dev_version && !is_release_version) {
    add_problem("DESCRIPTION Version must be three-part (release) or four-part (development): ", version)
  }

  news_path <- file.path(pkg_root, "NEWS.md")
  if (!file.exists(news_path)) {
    add_problem("NEWS.md is missing.")
  } else {
    news_lines <- readLines(news_path, warn = FALSE, encoding = "UTF-8")
    heading <- news_lines[grepl("^#\\s+", news_lines)][1]
    if (is.na(heading) || !nzchar(heading)) {
      add_problem("NEWS.md has no top-level heading.")
    } else {
      news_version <- sub("^#\\s+cpopvar\\s+([^[:space:]]+).*$", "\\1", heading)
      if (identical(news_version, heading) || !identical(news_version, version)) {
        add_problem(
          "NEWS.md first heading version (", heading, ") must match DESCRIPTION Version (", version, ")."
        )
      }
    }
  }

  cff_path <- file.path(pkg_root, "CITATION.cff")
  if (!file.exists(cff_path)) {
    add_problem("CITATION.cff is missing from the source tree.")
  } else {
    cff_lines <- readLines(cff_path, warn = FALSE, encoding = "UTF-8")
    version_line <- cff_lines[grepl("^version:\\s*", cff_lines)][1]
    if (is.na(version_line)) {
      add_problem("CITATION.cff has no version field.")
    } else {
      cff_version <- sub("^version:\\s*", "", version_line)
      cff_version <- gsub("^[\"']|[\"']$", "", cff_version)
      if (!identical(cff_version, version)) {
        add_problem(
          "CITATION.cff version (", cff_version, ") must match DESCRIPTION Version (", version, ")."
        )
      }
    }
    date_line <- cff_lines[grepl("^date-released:\\s*", cff_lines)][1]
    if (is_dev_version && !is.na(date_line)) {
      add_problem("Development versions must not set CITATION.cff:date-released.")
    }
  }

  citation_path <- file.path(pkg_root, "inst", "CITATION")
  if (!file.exists(citation_path)) {
    add_problem("inst/CITATION is missing.")
  } else {
    citation_text <- paste(readLines(citation_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    if (!grepl("meta$Version", citation_text, fixed = TRUE)) {
      add_problem("inst/CITATION must derive the package version from meta$Version.")
    }
    hardcoded <- regmatches(
      citation_text,
      gregexpr("\\b[0-9]+\\.[0-9]+\\.[0-9]+(?:\\.[0-9]+)?\\b", citation_text, perl = TRUE)
    )[[1]]
    hardcoded <- setdiff(hardcoded, "1.2.0")
    if (length(hardcoded) > 0) {
      add_problem(
        "inst/CITATION must not hard-code package versions (", paste(unique(hardcoded), collapse = ", "), ")."
      )
    }
  }

  expected_tarball <- paste0("cpopvar_", version, ".tar.gz")
  local_tarballs <- list.files(pkg_root, pattern = "^cpopvar_.*\\.tar\\.gz$", full.names = FALSE)
  unexpected <- setdiff(local_tarballs, expected_tarball)
  if (length(unexpected) > 0) {
    add_problem(
      "Package-root tarball names must be ", expected_tarball, "; found ", paste(unexpected, collapse = ", "), "."
    )
  }

  tarball_hint <- Sys.getenv("CPOPVAR_TARBALL", unset = "")
  if (nzchar(tarball_hint)) {
    if (!file.exists(tarball_hint)) {
      add_problem("CPOPVAR_TARBALL does not exist: ", tarball_hint)
    } else if (!identical(basename(tarball_hint), expected_tarball)) {
      add_problem(
        "Tarball file name (", basename(tarball_hint), ") must be ", expected_tarball, "."
      )
    }
  }

  git_dir <- file.path(pkg_root, ".git")
  if (dir.exists(git_dir) && nzchar(Sys.which("git"))) {
    tags <- tryCatch(
      system2("git", c("-C", pkg_root, "tag", "-l"), stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)
    )
    tags <- tags[nzchar(tags)]
    version_tag <- paste0("v", version)
    if (is_dev_version) {
      forbidden <- intersect(tags, c(version, version_tag))
      if (length(forbidden) > 0) {
        add_problem("Four-part development versions must not have git tags: ", paste(forbidden, collapse = ", "))
      }
    }
    if (is_release_version) {
      three_part_tags <- tags[grepl("^v[0-9]+\\.[0-9]+\\.[0-9]+$", tags)]
      mismatch <- setdiff(three_part_tags, version_tag)
      if (version_tag %in% tags && length(mismatch) > 0) {
        add_problem(
          "Release tags must equal v", version, "; extra three-part tags: ", paste(mismatch, collapse = ", ")
        )
      }
    }
  }

  namespace_path <- file.path(pkg_root, "NAMESPACE")
  if (!file.exists(namespace_path)) {
    add_problem("NAMESPACE is missing.")
  } else {
    ns_lines <- readLines(namespace_path, warn = FALSE)
    imported <- unique(c(
      sub("^import\\(([^),]+).*\\)$", "\\1", ns_lines[grepl("^import\\(", ns_lines)]),
      sub("^importFrom\\(([^),]+).*\\)$", "\\1", ns_lines[grepl("^importFrom\\(", ns_lines)])
    ))
    imported <- imported[!imported %in% ns_lines]
    declared <- character(0)
    for (field in c("Depends", "Imports")) {
      if (field %in% colnames(desc)) {
        raw <- gsub("\\s+", " ", desc[1, field])
        pkgs <- strsplit(raw, ",", fixed = TRUE)[[1]]
        pkgs <- trimws(sub("\\s*\\(.*$", "", pkgs))
        pkgs <- pkgs[nzchar(pkgs) & pkgs != "R"]
        declared <- c(declared, pkgs)
      }
    }
    base_pkgs <- c(
      "base", "compiler", "datasets", "graphics", "grDevices", "grid",
      "methods", "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils"
    )
    missing_decl <- setdiff(imported, c(declared, base_pkgs))
    if (length(missing_decl) > 0) {
      add_problem(
        "NAMESPACE imports not declared in DESCRIPTION: ", paste(missing_decl, collapse = ", ")
      )
    }
  }

  workflow_dir <- file.path(pkg_root, ".github", "workflows")
  if (dir.exists(workflow_dir)) {
    for (wf in list.files(workflow_dir, pattern = "\\.(yml|yaml)$", full.names = TRUE)) {
      wf_text <- paste(readLines(wf, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
      doc_starts <- gregexpr("(^|\\n)---[[:space:]]*(\\n|$)", wf_text, perl = TRUE)[[1]]
      if (length(doc_starts) > 1 && !identical(doc_starts, -1L)) {
        add_problem("Workflow has duplicate YAML documents: ", basename(wf))
      }
    }
  }
  citation_files <- list.files(file.path(pkg_root, "inst"), pattern = "^CITATION", full.names = FALSE)
  extra_citation <- setdiff(citation_files, "CITATION")
  if (length(extra_citation) > 0) {
    add_problem("inst/ has extra CITATION files: ", paste(extra_citation, collapse = ", "))
  }

  user_facing <- c(
    "DESCRIPTION",
    "README.md",
    "NEWS.md",
    "CITATION.cff",
    file.path("inst", "CITATION"),
    "LICENSE.md"
  )
  forbidden <- c(
    "CpPop-Glycine",
    "launch_workbench",
    "you@example.com",
    "your@email",
    "noreply@example"
  )
  whitelist_re <- "testthat|schema_version|pipeline_version|workflow_version"
  for (rel in user_facing) {
    path <- file.path(pkg_root, rel)
    if (!file.exists(path)) next
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    for (i in seq_along(lines)) {
      line <- lines[[i]]
      if (grepl(whitelist_re, line, ignore.case = TRUE)) next
      for (token in forbidden) {
        if (grepl(token, line, fixed = TRUE)) {
          add_problem(rel, ":", i, " contains forbidden token '", token, "'.")
        }
      }
      if (grepl("(?i)\\bv3\\.0\\.0\\b", line, perl = TRUE) && !grepl(whitelist_re, line)) {
        add_problem(rel, ":", i, " contains retired product version v3.0.0.")
      }
      if (grepl("github.com/Liyun-Nie/cpop(?!var)", line, perl = TRUE)) {
        add_problem(rel, ":", i, " contains old repository URL Liyun-Nie/cpop.")
      }
    }
  }

  if (length(problems) > 0) {
    cat("Version consistency check FAILED\n")
    cat(paste0(" - ", problems, collapse = "\n"), "\n", sep = "")
    stop("Version consistency check failed with ", length(problems), " problem(s).", call. = FALSE)
  }

  cat("Version consistency check PASSED\n")
  cat("DESCRIPTION Version:", version, "\n")
  cat("Expected tarball:", expected_tarball, "\n")
  invisible(list(version = version, tarball = expected_tarball, problems = problems))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  root <- if (length(args) >= 1) args[[1]] else "."
  check_version_consistency(root)
}
