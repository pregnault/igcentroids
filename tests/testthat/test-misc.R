# Binary search algorithm

test_that("binarysearch has three mandatory arguments", {
  tmp <- function(x) (x^2 - 2)
  expect_error(binarysearch(lb = NULL, ub = 2, fun = tmp))
  expect_error(binarysearch(lb = 0, ub = NULL, fun = tmp))
  expect_error(binarysearch(lb = 0, ub = 2, fun = NULL))
})

test_that("binarysearch expects an increasing function as input", {
  tmp <- function(x) (-x^2 + 2)
  expect_error(binarysearch(lb = 0, ub = 2, fun = tmp))
  tmp <- function(x) 1
  expect_error(binarysearch(lb = 0, ub = 2, fun = tmp))
})

test_that("binarysearch expects an 0 belongs to (fun(lb), fun(ub))", {
  tmp <- function(x) (x^2 + 2)
  expect_error(binarysearch(lb = 0, ub = 2, fun = tmp))
})

test_that("binarysearch provides an approximation of the antecedent of zero with tolerance 'tol'", {
  tmp <- function(x) (x^2 - 2)
  tol <- 1e-8
  bs <- binarysearch(lb = 0, ub = 2, fun = tmp, tol = tol)
  expect_equal(bs$sol, sqrt(2), tolerance = tol)
})
