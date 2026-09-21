# Binary search algorithm ----

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

test_that("binarysearch does not enter loop if solution is achieved for either the lower or upper bound", {
  tmp <- function(x) return(x-1)
  bs <- binarysearch(lb = 0, ub = 1, fun = tmp, tol = 1e-8)
  expect_equal(bs$iter, 0)
  tmp <- function(x) x
  bs <- binarysearch(lb = 0, ub = 1, fun = tmp, tol = 1e-8)
  expect_equal(bs$iter, 0)
})


# undo_igc ----

test_that("undo_igc expects input to be igc objects", {
  a <- igccatdist(1:3)
  b <- igccatdistn(rep(1/4, times = 4L))
  c <- igcmatcatdist(matrix(1:6, nrow = 3L))
  d <- igcmatcatdistn(matrix(c(dbinom(0:2, 2, prob = 0.2),
                               dbinom(0:2, 2, prob = 0.5)),
                             nrow =2, byrow = TRUE))
  expect_no_error(undo_igc(a))
  expect_no_error(undo_igc(b))
  expect_no_error(undo_igc(c))
  expect_no_error(undo_igc(d))
  e <- 1:4
  expect_error(undo_igc(e))
  f <- matrix(1:6, nrow = 3L)
  expect_error(undo_igc(f))
})
