test_that("installed NEWS and CITATION follow DESCRIPTION Version", {
  version <- as.character(utils::packageVersion("cpopvar"))
  news_path <- system.file("NEWS.md", package = "cpopvar")
  expect_true(nzchar(news_path) && file.exists(news_path))

  news_lines <- readLines(news_path, warn = FALSE, encoding = "UTF-8")
  heading <- news_lines[grepl("^#\\s+", news_lines)][1]
  expect_match(heading, paste0("^#\\s+cpopvar\\s+", gsub("\\.", "\\\\.", version)))

  citation_path <- system.file("CITATION", package = "cpopvar")
  expect_true(nzchar(citation_path) && file.exists(citation_path))
  citation_text <- paste(readLines(citation_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(citation_text, "meta\\$Version", perl = TRUE)
})

test_that("source-tree version consistency script passes when available", {
  script <- testthat::test_path("..", "..", "tools", "check_version_consistency.R")
  skip_if_not(file.exists(script), "source-tree version script is not available")

  source(script, local = TRUE)
  pkg_root <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
  expect_no_error(result <- check_version_consistency(pkg_root))
  expect_identical(result$version, as.character(utils::packageVersion("cpopvar")))
  expect_identical(result$tarball, paste0("cpopvar_", result$version, ".tar.gz"))
})
