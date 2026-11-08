test_that("LOO distances match a direct mean recalculation", {
  zmat <- as.matrix(EmbeddedTestData(variables = c("Var1", "Var2"))[1:3, ])
  nodes <- c(1L, 1L, 2L)
  codes <- rbind(colMeans(zmat[1:2, ]), zmat[3, ])
  expected <- vapply(1:3, function(i) {
    same <- which(nodes == nodes[i] & seq_along(nodes) != i)
    own <- if (length(same)) sqrt(sum((zmat[i, ] - colMeans(zmat[same, , drop = FALSE]))^2)) else Inf
    other <- sqrt(sum((zmat[i, ] - codes[3L - nodes[i], ])^2))
    min(own, other)
  }, numeric(1))
  expect_equal(.SOMLeaveOneOutDistances(zmat, codes, nodes), expected)
})
