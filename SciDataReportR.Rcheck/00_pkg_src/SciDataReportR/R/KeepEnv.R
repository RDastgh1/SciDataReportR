#' Keep selected objects in an environment and remove everything else
#'
#' Keeps only the objects specified in `Keep` within `Env` and removes all other
#' objects in that environment. Optionally returns a summary of what was removed.
#'
#' @param Keep Character vector of object names to keep.
#' @param Env Environment to clean. Defaults to the calling environment.
#' @param DryRun If TRUE, does not remove anything and only reports what would be removed.
#' @param Invert If TRUE, removes only `Keep` and keeps everything else.
#' @param Quiet If TRUE, suppresses messages.
#'
#' @return Invisible list with `kept` and `removed` vectors.
#' @export
KeepEnv <- function(Keep,
                    #' @param Env Environment to clean. Defaults to the calling environment.
                    #' @param DryRun If TRUE, does not remove anything and only reports what would be removed.
                    #' @param Invert If TRUE, removes only `Keep` and keeps everything else.
                    #' @param Quiet If TRUE, suppresses messages.
                    Env = parent.frame(),
                    DryRun = FALSE,
                    Invert = FALSE,
                    Quiet = FALSE
) {
  if (!is.environment(Env)) stop("`Env` must be an environment.")
  if (missing(Keep) || is.null(Keep)) Keep <- character(0)
  if (!is.character(Keep)) stop("`Keep` must be a character vector of object names.")

  existing <- ls(envir = Env, all.names = TRUE)

  # Nothing to do
  if (length(existing) == 0) {
    if (!Quiet) message("Environment is already empty.")
    return(invisible(list(kept = character(0), removed = character(0))))
  }

  keep_set <- unique(Keep)

  # Decide what to remove
  to_remove <- if (!Invert) {
    setdiff(existing, keep_set)
  } else {
    intersect(existing, keep_set)
  }

  kept <- setdiff(existing, to_remove)

  if (!Quiet) {
    if (DryRun) {
      message("Dry run: would remove ", length(to_remove), " object(s).")
    } else {
      message("Removing ", length(to_remove), " object(s).")
    }
  }

  if (!DryRun && length(to_remove) > 0) {
    rm(list = to_remove, envir = Env)
  }

  invisible(list(kept = kept, removed = to_remove))
}
