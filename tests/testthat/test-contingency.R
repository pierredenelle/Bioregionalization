context("Create a contingency table from a data.frame")

test_that("Returned object is a matrix", {

  ex <- data.frame(sites = c(rep("A", 2), rep("B", 3), rep("C", 2)),
                   species = c("a", "b", "a", "c", "d", "b", "d"),
                   count = c(10, 100, 1, 20, 50, 10, 20))

  # Object is matrix
  expect_is(contingency(dat = ex, site = "sites", sp = "species",
                        ab = "count"), "matrix")
})
