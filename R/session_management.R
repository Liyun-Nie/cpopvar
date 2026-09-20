#' @title Session Management Utilities
#' @description Functions for creating, managing, and cleaning up user sessions.
#' @name session_management
NULL

#' Generate unique session ID for user isolation
#' @return A string representing the unique session ID.
#' @export
generate_session_id <- function() {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  random_string <- paste(sample(c(letters, LETTERS, 0:9), 8, replace = TRUE), collapse = "")
  session_id <- paste0("session_", timestamp, "_", random_string)
  return(session_id)
}


#' Initialize session management
#' @description Creates base directories and cleans up old sessions.
#' @export
init_session_management <- function() {
  sessions_dir <- file.path("app_data", "sessions")
  if (!dir.exists(sessions_dir)) {
    dir.create(sessions_dir, recursive = TRUE)
  }
  cat("[OK] Session management initialized\n")
}
