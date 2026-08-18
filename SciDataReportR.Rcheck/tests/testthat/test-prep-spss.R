test_that("PrepSPSS truncates over-long factor levels for SPSS value labels", {
  skip_if_not_installed("janitor")

  long_a <- paste(rep("alpha", 40), collapse = " ")
  long_b <- paste(rep("bravo", 40), collapse = " ")

  df <- data.frame(group = factor(c(long_a, long_b, "short")))

  res <- PrepSPSS(df, quiet = TRUE)
  new_levels <- levels(res$data$group)

  expect_true(all(nchar(new_levels, type = "bytes") <= 120))
  expect_length(unique(new_levels), 3)
  expect_true("short" %in% new_levels)

  value_rows <- res$label_map[res$label_map$label_type == "value_label", ]
  expect_equal(nrow(value_rows), 2)
  expect_setequal(value_rows$original_label, c(long_a, long_b))
  expect_true(all(endsWith(value_rows$spss_label, "...")))
})

test_that("PrepSPSS de-duplicates factor levels that collide after truncation", {
  skip_if_not_installed("janitor")

  shared_prefix <- strrep("x", 150)
  lv1 <- paste0(shared_prefix, " ending one")
  lv2 <- paste0(shared_prefix, " ending two")

  df <- data.frame(group = factor(c(lv1, lv2)))

  res <- PrepSPSS(df, quiet = TRUE)
  new_levels <- levels(res$data$group)

  expect_true(all(nchar(new_levels, type = "bytes") <= 120))
  expect_length(unique(new_levels), 2)
  expect_false(anyNA(res$data$group))
})

test_that("PrepSPSS truncates over-long value labels on labelled vectors", {
  skip_if_not_installed("janitor")

  long_label <- strrep("a", 130)
  x <- structure(
    c(1, 2, 1),
    labels = stats::setNames(c(1, 2), c(long_label, "ok"))
  )

  df <- data.frame(score = 1:3)
  df$score <- x

  res <- PrepSPSS(df, quiet = TRUE)
  new_labels <- names(attr(res$data$score, "labels", exact = TRUE))

  expect_true(all(nchar(new_labels, type = "bytes") <= 120))
  expect_true("ok" %in% new_labels)
  expect_equal(sum(res$label_map$label_type == "value_label"), 1)
  expect_equal(res$label_map$original_label[[1]], long_label)
})

test_that("PrepSPSS caps variable labels at the SPSS 256-byte limit", {
  skip_if_not_installed("janitor")

  long_name <- paste(rep("question", 40), collapse = " ")
  df <- data.frame(x = 1:3)
  names(df) <- long_name

  res <- PrepSPSS(df, quiet = TRUE)
  new_label <- attr(res$data[[1]], "label", exact = TRUE)

  expect_true(nchar(new_label, type = "bytes") <= 256)
  expect_true(endsWith(new_label, "..."))

  variable_rows <- res$label_map[res$label_map$label_type == "variable_label", ]
  expect_equal(nrow(variable_rows), 1)
  expect_equal(variable_rows$original_label[[1]], long_name)
})

test_that("PrepSPSS returns an empty label map when nothing is truncated", {
  skip_if_not_installed("janitor")

  df <- data.frame(
    group = factor(c("A", "B")),
    value = c(1, 2)
  )

  res <- PrepSPSS(df, quiet = TRUE)

  expect_equal(nrow(res$label_map), 0)
  expect_named(
    res$label_map,
    c("spss_name", "label_type", "value", "original_label", "spss_label")
  )
})

test_that("PrepSPSS output with long labels writes successfully via write_sav", {
  skip_if_not_installed("janitor")
  skip_if_not_installed("haven")

  long_level <- paste(rep("chronic condition", 15), collapse = " ")
  df <- data.frame(
    diagnosis = factor(c(long_level, "none")),
    value = c(1, 2)
  )

  sav_path <- tempfile(fileext = ".sav")
  on.exit(unlink(sav_path), add = TRUE)

  expect_no_error(PrepSPSS(df, path = sav_path, quiet = TRUE))
  expect_true(file.exists(sav_path))

  reread <- haven::read_sav(sav_path)
  expect_equal(ncol(reread), 2)
})
