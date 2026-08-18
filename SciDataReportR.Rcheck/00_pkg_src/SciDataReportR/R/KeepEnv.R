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
#' @details
#' A long analysis session accumulates intermediates - raw imports, temporary
#' merges, loop variables, half-finished plots - and only a few of those
#' objects are worth carrying forward. Saving the whole workspace preserves all
#' of it, producing `.RData` files that are large, slow to load, and full of
#' objects whose provenance nobody remembers.
#'
#' This is the inverse of writing out a long `rm()` call: name the handful of
#' objects that matter and everything else goes. Because the list is
#' allow-list rather than deny-list, objects created after the code was written
#' are cleaned up too, instead of surviving because nobody remembered to add
#' them to the `rm()`.
#'
#' Run it with `DryRun = TRUE` first. It reports exactly what would be removed
#' without touching anything, which is worth doing whenever the environment
#' holds something expensive to recompute.
#'
#' @return Invisible list with `kept` and `removed` vectors.
#'
#' @examples
#' # An analysis environment: a few results among many intermediates
#' env_Analysis <- new.env()
#' local({
#'   df_Raw <- data.frame(id = 1:5, value = rnorm(5))
#'   df_Clean <- df_Raw[!is.na(df_Raw$value), ]
#'   tmp_merge <- df_Clean
#'   i <- 3
#'   scratch_vector <- 1:100
#'   model_Final <- lm(value ~ id, data = df_Clean)
#'   df_Results <- data.frame(term = "id", estimate = coef(model_Final)[2])
#' }, envir = env_Analysis)
#'
#' ls(env_Analysis)
#'
#' # Check what would go, before anything is removed
#' preview <- KeepEnv(
#'   Keep = c("df_Results", "model_Final"),
#'   Env = env_Analysis,
#'   DryRun = TRUE
#' )
#' preview$removed
#'
#' # Then do it for real, and save a workspace holding only what matters
#' KeepEnv(c("df_Results", "model_Final"), Env = env_Analysis)
#' ls(env_Analysis)
#'
#' save(list = ls(env_Analysis), envir = env_Analysis,
#'      file = file.path(tempdir(), "analysis_results.RData"))
#'
#' # `Invert = TRUE`: drop the objects named, keep the rest
#' env_Other <- new.env()
#' assign("df_Huge", data.frame(x = 1:10), envir = env_Other)
#' assign("df_Small", data.frame(x = 1:2), envir = env_Other)
#' KeepEnv("df_Huge", Env = env_Other, Invert = TRUE)
#' ls(env_Other)
#'
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
