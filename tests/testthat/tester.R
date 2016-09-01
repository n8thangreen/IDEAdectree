library(IDEAdectree)
context("IDEAdectree.simple")

test_that("basics", {
  expect_equal(nrow(IDEAdectree.simple(data=data, nsim = 10)$c), 10)
  expect_equal(nrow(IDEAdectree.simple(data=data, nsim = 10)$e), 10)

  expect_equal(ncol(IDEAdectree.simple(data=data, nsim = 10)$c), 2)
  expect_equal(ncol(IDEAdectree.simple(data=data, nsim = 10)$e), 2)
})

