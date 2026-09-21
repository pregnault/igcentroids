# Creating objetcs of classes igc... ----

## Objects igccatdist and igccatdistn ----

test_that("Function igccatdist returns an object of class igcccatdist inheriting from numeric vector", {
  p <- rep(1, times = 5L)
  pcd <- igccatdist(p)
  expect_equal(methods::S3Class(pcd), "igccatdist")
  expect_equal(mode(pcd), "numeric")
})

test_that("Function igccatdist returns an error message if input is not a positive numeric vector", {
  p <- letters[1:3]
  expect_error(igccatdist(p))
  p <- numeric()
  expect_error(igccatdist(p))
  p <- c(-1, 0, 2)
  expect_error(igccatdist(p))
})

test_that("Function igccatdist returns an error message if attribute mass is not a single positive value", {
  p <- igccatdist(rep(1, times= 3L))
  attributes(p) <- list(mass = "bug")
  expect_error(validate_igccatdist(p))
  attributes(p) <- list(mass = rep(1, times = 3L))
  expect_error(validate_igccatdist(p))
  attributes(p) <- list(mass = -3)
  expect_error(validate_igccatdist(p))
})

test_that("Function igccatdistn returns an object of class igcccatdist inheriting from numeric vector", {
  p <- rep(1/5, times = 5L)
  pcd <- igccatdist(p)
  expect_equal(methods::S3Class(pcd), "igccatdist")
  expect_equal(mode(pcd), "numeric")
})

test_that("Function igccatdistn returns an error message if input is not a positive numeric vector", {
  p <- matrix(rep(1/2, times = 4L), ncol = 2)
  expect_error(igccatdistn(p))
  p <- letters[1:3]
  expect_error(igccatdistn(p))
  p <- numeric()
  expect_error(igccatdistn(p))
  p <- c(-1, 0, 2)
  expect_error(igccatdistn(p))
})

test_that("Function igccatdistn returns an error message if attribute mass is not a single positive value", {
  p <- igccatdistn(rep(1/3, times= 3L))
  attributes(p) <- list(mass = "bug")
  expect_error(validate_igccatdist(p))
  attributes(p) <- list(mass = rep(1, times = 3L))
  expect_error(validate_igccatdist(p))
  attributes(p) <- list(mass = -3)
  expect_error(validate_igccatdist(p))
})

## Objects igcmatcatdist and igcmatcatdistn ----

test_that("Functions igcmatcatdist and igcmatcatdistn return objects of the eponyme class", {
  p <- matrix(rep(1/2, times = 4L), ncol = 2)
  pmcd <- igcmatcatdist(p)
  expect_equal(methods::S3Class(pmcd), "igcmatcatdist")
  pmcd <- igcmatcatdistn(p)
  expect_equal(methods::S3Class(pmcd), c("igcmatcatdistn", "igcmatcatdist"))
  expect_equal(mode(pmcd), "numeric")
})

test_that("Functions igcmatcatdist and igcmatcatdistn return an error message when input is not a positive matrix", {
  p <- matrix(letters[1:12], nrow = 3L, ncol = 4L)
  expect_error(igcmatcatdist(p))
  expect_error(igcmatcatdistn(p))
  p <- numeric()
  dim(p) <- c(0,0)
  expect_error(igcmatcatdist(p))
  expect_error(igcmatcatdistn(p))
  p <- matrix(c(-1, 1:11), nrow = 3L, ncol = 4L)
  expect_error(igcmatcatdist(p))
  expect_error(igcmatcatdistn(p))
})

test_that("Functions igcmatcatdist and igcmatcatdistn return an error message if attribute mass is not a positive vector", {
  p <- igcmatcatdistn(matrix(rep(1/4, times= 4L), nrow = 3L, ncol = 4L, byrow = TRUE))
  attributes(p) <- list(dim = c(3, 4),
                        class = c("igcmatcatdistn", "igcmatcatdist"),
                        mass = "bug")
  expect_error(validate_igcmatcatdist(p))
  expect_error(validate_igcmatcatdistn(p))
  attributes(p) <- list(dim = c(3, 4),
                        class = c("igcmatcatdistn", "igcmatcatdist"),
                        mass = rep(1, times = 4L))
  expect_error(validate_igcmatcatdist(p))
  expect_error(validate_igcmatcatdistn(p))
  attributes(p) <- list(dim = c(3, 4),
                        class = c("igcmatcatdistn", "igcmatcatdist"),
                        mass = c(1, 1, -1))
  expect_error(validate_igccatdist(p))
  expect_error(validate_igcmatcatdistn(p))
})

# Binding igc objects ----

## Converting catdist to matcatdist ----

test_that("Function cat2mat properly converts igccatdist objects to igcmatcatdist", {
  p <- igccatdist(rep(1, times = 5L))
  pcv <- cat2mat(p)
  expect_equal(dim(pcv), c(1, 5))
  expect_equal(class(pcv), "igcmatcatdist")
})

test_that("Function cat2mat properly converts igccatdistn objects to igcmatcatdistn", {
  p <- igccatdistn(rep(1/5, times = 5L))
  pcv <- cat2mat(p)
  expect_equal(dim(pcv), c(1, 5))
  expect_equal(class(pcv), c("igcmatcatdistn", "igcmatcatdist"))
})

## Binding igc objects ----

test_that("Function bind_igc returns a positive matrix when inputs are catdist or positive matrices", {
  A <- igcmatcatdist(matrix(1:6, nrow = 2))
  B <- igcmatcatdistn(matrix(c(dbinom(0:2, size = 2, prob = 0.2),
                               dbinom(0:2, size = 2, prob = 0.5)),
                             nrow = 2, byrow = TRUE))
  C <- igccatdist(7:9)
  resbindAB <- bind_igc(A,B)
  expect_equal(dim(resbindAB), c(4, 3))
  expect_equal(methods::S3Class(resbindAB), "igcmatcatdist")
  resbindAC <- bind_igc(A,C)
  expect_equal(dim(resbindAC), c(3, 3))
  expect_equal(methods::S3Class(resbindAC), "igcmatcatdist")
})

test_that("Function bind_igc returns an error message if inputs are not all of an igc class", {
  A <- igcmatcatdist(matrix(1:6, nrow = 2))
  B <-matrix(c(dbinom(0:2, size = 2, prob = 0.2),
               dbinom(0:2, size = 2, prob = 0.5)),
             nrow = 2, byrow = TRUE)
  expect_error(bind_igc(A, B))
  C <- 7:9
  expect_error(bind_igc(A, C))
})

test_that("Function bind_igc returns an error message if inputs have incompatible sizes", {
  B <-igcmatcatdistn(matrix(c(dbinom(0:2, size = 2, prob = 0.2),
                             dbinom(0:2, size = 2, prob = 0.5)),
                           nrow = 2, byrow = TRUE))
  E <- igcmatcatdistn(matrix(rep(1/4, times = 4L),
                             nrow = 1L, ncol = 4L))
  expect_error(bind_igc(B, E))
})


