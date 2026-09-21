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

test_that("jvid applied to a couple of positive matrices returns the matrix of jdiv between rows of the first and second matrices", {
  p <- igcmatcatdist(matrix(1:6, nrow = 2L, ncol = 3L))
  q <- igcmatcatdist(matrix(1:12, nrow = 4L, ncol = 3L))
  r <- jdiv(p, q)
  expect_equal(dim(r), c(2, 4))
  expect_equal(r[1,2], jdiv(igccatdist(c(1, 3, 5)), igccatdist(c(2, 6, 10))))
})

test_that("jdiv applied to a single positive matrix returns the matrix of jdiv between rows of the input matrix", {
  p <- igcmatcatdist(matrix(1:6, nrow = 2L, ncol = 3L))
  r <- jdiv(p)
  expect_equal(dim(r), c(2, 2))
  expect_equal(r[1,2], jdiv(igccatdist(c(1, 3, 5)), igccatdist(c(2, 4, 6))))
  expect_equal(r[2,1], r[1,2])
})

## Jensen-Shannon divergence ----

test_that("jsdiv returns an error message if inputs p and q are not consistent", {
  p <- igccatdist(1:5)
  q <- igccatdist(6:11)
  expect_error(jsdiv(p, q))
})

test_that("jsdiv returns a single positive or null value when
            inputs are positive vectors of same size", {
              p <- igccatdist(1:5)
              q <- igccatdist(6:10)
              div <- jsdiv(p, q)
              expect_true(is.numeric(div) & length(div) == 1L)
              expect_gte(div, 0)
            })

test_that("jsdiv is a symmetric functional", {
  p <- igccatdist(1:5)
  q <- igccatdist(6:10)
  divpq <- jsdiv(p, q)
  divqp <- jsdiv(q, p)
  expect_equal(divpq, divqp)
})

test_that("jsdiv with base different from exp(1) returns consistent results", {
  p <- igccatdist(1:5)
  q <- igccatdist(6:10)
  div1 <- jsdiv(p, q)
  div2 <- jsdiv(p, q, base = 2)
  expect_equal(div2, div1 / log(2))
})


test_that("jsdiv returns values between 0 and log(2).", {
  set.seed(847)
  p <- matrix(sample(1:10, size = 500, replace = TRUE),
              nrow = 100, ncol = 5)
  p <- t(apply(p, 1, function(row) row / sum(row)))
  q <- matrix(sample(1:10, size = 500, replace = TRUE),
              nrow = 100, ncol = 5)
  q <- t(apply(q, 1, function(row) row / sum(row)))
  divs <- sapply(1:nrow(p), function(idrow) jsdiv(igccatdist(p[idrow,]), igccatdist(q[idrow,])))
  expect_all_true(divs >= 0 & divs <= log(2))
})

test_that("jsdiv equals the difference between the entropy of the mixture and the average of entropies.", {
  set.seed(847)
  p <- matrix(sample(1:10, size = 500, replace = TRUE),
              nrow = 100, ncol = 5)
  p <- t(apply(p, 1, function(row) row / sum(row)))
  q <- matrix(sample(1:10, size = 500, replace = TRUE),
              nrow = 100, ncol = 5)
  q <- t(apply(q, 1, function(row) row / sum(row)))
  divs <- sapply(1:nrow(p), function(idrow) jsdiv(igccatdist(p[idrow,]), igccatdist(q[idrow,])))
  entsp <- sapply(1:nrow(p), function(idrow) -sum(p[idrow,]*log(p[idrow,])))
  entsq <- sapply(1:nrow(q), function(idrow) -sum(q[idrow,]*log(q[idrow,])))
  entsmix <- sapply(1:nrow(p), function(idrow) -sum((p[idrow,]*0.5+q[idrow,]*0.5)*log((p[idrow,]*0.5+q[idrow,]*0.5))))
  divsalt <- entsmix - entsp*0.5 - entsq*0.5
  expect_true(all.equal(divs, divsalt))
})

test_that("jsdiv applied to a couple of positive matrices returns the matrix of jdiv between rows of the first and second matrices", {
  p <- igcmatcatdist(matrix(1:6, nrow = 2L, ncol = 3L))
  q <- igcmatcatdist(matrix(1:12, nrow = 4L, ncol = 3L))
  r <- jsdiv(p, q)
  expect_equal(dim(r), c(2, 4))
  expect_equal(r[1,2], jsdiv(igccatdist(c(1, 3, 5)), igccatdist(c(2, 6, 10))))
})

test_that("jsdiv applied to a single positive matrix returns the matrix of jdiv between rows of the input matrix", {
  p <- igcmatcatdist(matrix(1:6, nrow = 2L, ncol = 3L))
  r <- jsdiv(p)
  expect_equal(dim(r), c(2, 2))
  expect_equal(r[1,2], jsdiv(igccatdist(c(1, 3, 5)), igccatdist(c(2, 4, 6))))
  expect_equal(r[2,1], r[1,2])
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

test_that("Argument weight of geomcent is a positive vector", {
  p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                               dbinom(0:4, size = 4, prob = 0.3)),
                             nrow = 2, byrow = TRUE))
  weight <- c(-1, 1)
  expect_error(geomcent(p, weight = weight))
  weight <- 3
  expect_error(geomcent(p, weight = weight))
})

test_that("Argument weight of geomcent is invariant by scalar multiplication", {
  p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                               dbinom(0:4, size = 4, prob = 0.3)),
                             nrow = 2, byrow = TRUE))
  weight1 <- c(2, 2)
  weight2 <- c(10, 10)
  expect_equal(geomcent(p, weight = weight1), geomcent(p, weight = weight2))
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

## Jensen-Shannon centroid ----

test_that("jscent returns an error message if input is not an object
          of class igcmatcatdistn", {
            p <- matrix(1:6, nrow = 2)
            expect_error(jscent(p))
            p <- igccatdist(1:6)
            expect_error(jscent(p))
            p <- igcmatcatdist(matrix(1:6, nrow = 2))
            expect_error(jscent(p))
          })

test_that("jscent returns an object of class igccatdistn if input is
          of class icgmatcatdistn, with adequate length", {
            p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                                         dbinom(0:4, size = 4, prob = 0.3)),
                                       nrow = 2, byrow = TRUE))
            expect_identical(class(jscent(p)), c("igccatdistn", "igccatdist"))
            expect_true(length(jscent(p)) == ncol(p))
          })

test_that("jscent method for igcmatcatdistn has a mandatory argument 'tol'
          which is expected to be a numeric value", {
            p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                                         dbinom(0:4, size = 4, prob = 0.3)),
                                       nrow = 2, byrow = TRUE))
            expect_error(jscent(p = p, tol = NULL))
            expect_contains(class(jscent(p = p, tol = 1e-16)), "igccatdistn")
          })

test_that("Argument weight of jscent is a positive vector", {
  p <- igcmatcatdistn(matrix(c(dbinom(0:4, size = 4, prob = 0.5),
                               dbinom(0:4, size = 4, prob = 0.3)),
                             nrow = 2, byrow = TRUE))
  weight <- c(-1, 1)
  expect_error(jscent(p, weight = weight))
  weight <- 3
  expect_error(jscent(p, weight = weight))
})


# Plot methods ----

## Plot objects of class igcmatcatdistn ----

test_that("Method plot for class igcmatcatdistn unable to perform plot
          for matrices with only one column", {
  p <- igcmatcatdistn(matrix(rep(1, 10L), nrow = 10L))
  expect_error(plot(p))
})

test_that("Method plot for class igcmatcatdistn returns invisibly an object", {
  set.seed(987)
  p <- igcmatcatdistn(matrix(sample(1:10, size = 20L, replace = TRUE), nrow = 10L))
  expect_invisible(plot(p))
  p <- igcmatcatdistn(matrix(sample(1:10, size = 30L, replace = TRUE), nrow = 10L))
  expect_invisible(plot(p))
  p <- igcmatcatdistn(matrix(sample(1:10, size = 40L, replace = TRUE), nrow = 10L))
  expect_invisible(plot(p))
})

test_that("Two-state categorical distributions are drawn as points on a segment", {
  set.seed(987)
  p <- igcmatcatdistn(matrix(sample(1:10, size = 20L, replace = TRUE), nrow = 10L))
  res_plot <- plot(p)
  expect_match(tolower(res_plot$type), regexp = "segment")
})

test_that("Three-state categorical distributions are drawn as points in an equilateral triangle", {
  set.seed(987)
  p <- igcmatcatdistn(matrix(sample(1:10, size = 30L, replace = TRUE), nrow = 10L))
  res_plot <- plot(p)
  expect_match(tolower(res_plot$type), regexp = "triangle")
})

test_that("Three-state categorical distributions are properly transformed
          into coordinates in the equilateral triangle", {
  set.seed(987)
  matp <- matrix(c(1/3, 1/3, 1/3,
                   0.9, 0.05, 0.05,
                   0.1, 0.45, 0.45,
                   0.6, 0.25, 0.15),
                 nrow = 4L, byrow = TRUE)
  p <- igcmatcatdistn(matp)
  res_plot <- plot(p)
  correct_coords <- t(matrix(c(1, 0, 1/2, sqrt(3)/2), ncol = 2) %*% t(matp[, 1:2, drop = FALSE]))
  expect_true(all.equal(res_plot$coords, correct_coords))
})

test_that("Categorical distributions with more than three states are drawn as polylines.", {
  set.seed(987)
  p <- igcmatcatdistn(matrix(sample(1:10, size = 50L, replace = TRUE), nrow = 10L))
  res_plot <- plot(p)
  expect_match(tolower(res_plot$type), regexp = "polyline")
})


test_that("Methods points for igccatdistn and igcmatcatdistn add points to a pre-existing chart", {
            set.seed(987)
            matp <- matrix(c(1/3, 1/3, 1/3,
                             0.9, 0.05, 0.05,
                             0.1, 0.45, 0.45,
                             0.6, 0.25, 0.15),
                           nrow = 4L, byrow = TRUE)
            p <- igcmatcatdistn(matp)
            res_plot <- plot(p)
            matq <- bind_igc(arithcent(p),
                             geomcent(p))
            r <- jscent(p)
            expect_no_error(points(matq))
            expect_no_error(points(r))
          })

