# Methods for computing divergences ----

## KL divergence ----

test_that("kldiv returns an error message if inputs p and q are not consistent", {
  p <- igccatdist(1:5)
  q <- igccatdist(6:11)
  expect_error(kldiv(p, q))
})

test_that("kldiv returns a single positive or null value when inputs
          are positive vectors of same size", {
  p <- igccatdist(1:5)
  q <- igccatdist(6:10)
  div <- kldiv(p, q)
  expect_true(is.numeric(div) & length(div) == 1L)
  expect_gte(div, 0)
})

test_that("kldiv returns log(N) - S(p) when inputs are a probability distribution
          and the uniform distribution", {
  p <- igccatdistn(dbinom(0:4, size = 4, prob = 0.5))
  q <- igccatdistn(rep(1 / 5, times = 5))
  div <- kldiv(p, q)
  expect_equal(div, log(5) + sum(p * log(p)))
})

test_that("kldiv with base different from exp(1) returns consistent results", {
  p <- igccatdist(1:5)
  q <- igccatdist(6:10)
  div1 <- kldiv(p, q)
  div2 <- kldiv(p, q, base = 2)
  expect_equal(div2, div1 / log(2) + (sum(p) - sum(q)) * (1 / log(2) - 1))
})

## Jeffreys divergence ----

test_that("jdiv returns an error message if inputs p and q are not consistent", {
  p <- igccatdist(1:5)
  q <- igccatdist(6:11)
  expect_error(jdiv(p, q))
})

test_that("jdiv returns a single positive or null value when
            inputs are positive vectors of same size", {
  p <- igccatdist(1:5)
  q <- igccatdist(6:10)
  div <- jdiv(p, q)
  expect_true(is.numeric(div) & length(div) == 1L)
  expect_gte(div, 0)
})

test_that("jdiv is a symmetric functional", {
  p <- igccatdist(1:5)
  q <- igccatdist(6:10)
  divpq <- jdiv(p, q)
  divqp <- jdiv(q, p)
  expect_equal(divpq, divqp)
})

test_that("jdiv with base different from exp(1) returns consistent results", {
  p <- igccatdist(1:5)
  q <- igccatdist(6:10)
  div1 <- jdiv(p, q)
  div2 <- jdiv(p, q, base = 2)
  expect_equal(div2, div1 / log(2))
})

# Methods for computing centroids ----

## Arithmetic centroid ----

test_that("arithcent returns an error message if input is not an object
          of class igcmatcatdist", {
  p <- matrix(1:6, nrow = 2)
  expect_error(arithcent(p))
  p <- igccatdist(1:6)
  expect_error(arithcent(p))
})

test_that("arithcent returns an object of class igccatdist if input
          is of class icgmatcatdist", {
  p <- igcmatcatdist(matrix(1:6, nrow = 2))
  expect_true(class(arithcent(p)) == "igccatdist")
})

test_that("arithcent returns an object of class igccatdist with mass 1
          if input is of class icgmatcatdistn", {
  p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                               dbinom(0:4, size = 4, prob = 0.3)),
                             nrow = 2, byrow = TRUE))
  expect_equal(sum(arithcent(p)), 1)
})

test_that("Argument weight of arithcent is a positive vector", {
  p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                               dbinom(0:4, size = 4, prob = 0.3)),
                             nrow = 2, byrow = TRUE))
  weight <- c(-1, 1)
  expect_error(arithcent(p, weight = weight))
  weight <- 3
  expect_error(arithcent(p, weight = weight))
})

test_that("Argument weight of arithcent is invariant by scalar multiplication", {
  p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                               dbinom(0:4, size = 4, prob = 0.3)),
                             nrow = 2, byrow = TRUE))
  weight1 <- c(2, 2)
  weight2 <- c(10, 10)
  expect_equal(arithcent(p, weight = weight1), arithcent(p, weight = weight2))
})


## geometric centroid ----

test_that("geomcent returns an error message if input is not an object
          of class igcmatcatdist or igcmatcatdist", {
  p <- matrix(1:6, nrow = 2)
  expect_error(geomcent(p))
  p <- igccatdist(1:6)
  expect_error(geomcent(p))
})

test_that("geomcent returns an object of class igccatdist
          if input is of class icgmatcatdist", {
  p <- igcmatcatdist(matrix(1:6, nrow = 2))
  expect_true(class(geomcent(p)) == "igccatdist")
})

## Jeffreys centroid ----

test_that("jcent returns an error message if input is not an object
          of class igcmatcatdist or igcmatcatdist", {
  p <- matrix(1:6, nrow = 2)
  expect_error(jcent(p))
  p <- igccatdist(1:6)
  expect_error(jcent(p))
})

test_that("jcent returns an object of class igccatdist if input is
          of class icgmatcatdist, with adequate length", {
  p <- igcmatcatdist(matrix(1:6, nrow = 2))
  expect_true(class(jcent(p)) == "igccatdist")
  expect_true(length(jcent(p)) == ncol(p))
})

test_that("jcent returns an object of class igccatdistn if input is
          of class icgmatcatdistn, with adequate length", {
  p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                               dbinom(0:4, size = 4, prob = 0.3)),
                             nrow = 2, byrow = TRUE))
  expect_identical(class(jcent(p)), c("igccatdistn", "igccatdist"))
  expect_true(length(jcent(p)) == ncol(p))
})

test_that("jcent method for igcmatcatdist returns an igccatdist object whose
          sum is decreasing with argument lambda", {
  p <- igcmatcatdist(matrix(1:6, nrow = 2))
  lambda <- seq(-1, 0, by = 0.01)
  sumaslmbd <- apply(sapply(lambda, function(lmbd) jcent(p, lambda = lmbd)), 2, sum)
  expect_true(all(diff(sumaslmbd, lag = 1)< 0))
})

test_that("jcent method for igcmatcatdistn has a mandatory argument 'method'
          which is expected to be a character string among five possible values", {
  p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                               dbinom(0:4, size = 4, prob = 0.3)),
                             nrow = 2, byrow = TRUE))
  expect_error(jcent(p = p, method = NULL))
  expect_error(jcent(p = p, method = c("mag", "bsa")))
  expect_contains(class(jcent(p = p, method = "mag")), "igccatdistn")
  expect_contains(class(jcent(p = p, method = "njc")), "igccatdistn")
  expect_contains(class(jcent(p = p, method = "jfr")), "igccatdistn")
  expect_contains(class(jcent(p = p, method = "bsa")), "igccatdistn")
  expect_contains(class(jcent(p = p, method = "gbic")), "igccatdistn")
})

test_that("jcent method for igcmatcatdistn has a mandatory argument 'tol'
          which is expected to be a numeric value", {
  p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                              dbinom(0:4, size = 4, prob = 0.3)),
                            nrow = 2, byrow = TRUE))
  expect_error(jcent(p = p, tol = NULL))
  expect_contains(class(jcent(p = p, tol = 1e-16)), "igccatdistn")
})

test_that("Approximations returned by jcent are quite similar", {
  u <- rep(1/3, 3)
  prob <- c(0.9, 0.05, 0.05)
  p <- matrix(c(u, prob),
              nrow = 2, ncol =3, byrow = TRUE)
  p <- igcmatcatdistn(p)
  # Available methods for computing approximated values of jeffreys centroid
  meths <- c("mag", "njc", "jfr", "bsa", "gbic")
  # Approximations obtained from these methods
  jcents <- t(sapply(meths, function(m) jcent(p, method = m, tol = 1e-16)))
  dist_jcents <- dist(jcents, method = "manhattan")
  expect_all_true(dist_jcents < 1e-2)
})

