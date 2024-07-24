#' Insert Values into a String Array
#'
#' This function inserts one or more strings into an array of strings at specified locations.
#'
#' @param vec A character vector in which values will be inserted.
#' @param values A character vector or a single string to insert into the original vector.
#' @param after_value A character string indicating the element in `vec` before or after which the values will be inserted. Should have a length of one
#' @param location A character string specifying whether to insert `values` "before" or "after" the `after_value`. Defaults to "after".
#'
#' @return A character vector with the values inserted.
#' @export

InsertValues <- function(vec, values, after_value, location = "after") {
  # Check if after_value exists in the vector
  index <- which(vec == after_value)

  if (length(index) == 0) {
    stop(paste(after_value, "not found in the vector"))
  }

  # If location is "before", adjust the index to insert before the after_value
  if (location == "before") {
    index <- index - 1
  }

  # Insert the new values into the vector
  for (i in seq_along(values)) {
    vec <- append(vec, values[i], after = index + i - 1)
  }

  return(vec)
}
