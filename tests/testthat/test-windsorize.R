test_that("windsorize applies IQR limits to the requested side", {
  x <- c(-100, 0:10, 100, NA_real_)

  both <- windsorize(x, method = "iqr", iqrlim = 1.5)
  right <- windsorize(x, method = "iqr", iqrlim = 1.5, side = "right")
  left <- windsorize(x, method = "iqr", iqrlim = 1.5, side = "left")

  expect_equal(both, windsorize(x, method = "iqr", iqrlim = 1.5, side = "both"))
  expect_equal(both, c(-7, 0:10, 17, NA_real_))
  expect_equal(right, c(-100, 0:10, 17, NA_real_))
  expect_equal(left, c(-7, 0:10, 100, NA_real_))
})

test_that("windsorize applies SD limits to the requested side", {
  x <- c(-100, 0:10, 100, NA_real_)
  valid_x <- x[!is.na(x)]
  lower <- mean(valid_x) - stats::sd(valid_x)
  upper <- mean(valid_x) + stats::sd(valid_x)

  both <- windsorize(x, method = "sd", sdlim = 1)
  right <- windsorize(x, method = "sd", sdlim = 1, side = "right")
  left <- windsorize(x, method = "sd", sdlim = 1, side = "left")

  expect_equal(both, pmin(pmax(x, lower), upper))
  expect_equal(right, ifelse(x > upper, upper, x))
  expect_equal(left, ifelse(x < lower, lower, x))
  expect_true(is.na(right[length(right)]))
  expect_true(is.na(left[length(left)]))
})

test_that("windsorize validates side", {
  expect_error(windsorize(1:10, side = "upper"), "side must be")
  expect_error(windsorize(1:10, side = NA_character_), "side must be")
  expect_error(windsorize(1:10, side = c("both", "right")), "side must be")
})
