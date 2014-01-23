source("helper-Revolve.R")

context("Mutation")

test_that("Invalid VCV matrices", {
  expect_that(make_mutation(matrix(nrow=0, ncol=0)), throws_error())
  expect_that(make_mutation(matrix(nrow=1, ncol=0)), throws_error())
  expect_that(make_mutation(matrix(nrow=2, ncol=3)), throws_error())
})

test_that("One dimensional mutation", {
  sd <- 3 * sqrt(2)
  mut <- make_mutation(sd * sd)

  # One resident type:
  set.seed(1)
  xm <- pi
  n <- 1000
  x <- mut(xm, n)
  expect_that(length(x), equals(n))
  expect_that(ks.test(x, pnorm, xm, sd)$p.value, is_greater_than(1/20))

  # Two resident types:
  xm <- c(pi, -pi)
  n  <- c(1000, 501)
  x <- mut(xm, n)

  expect_that(length(x), equals(sum(n)))

  type <- rep(seq_along(n), n)
  for (i in seq_along(n))
    expect_that(ks.test(x[type == i], pnorm, xm[i], sd)$p.value,
                is_greater_than(1/20))

  # Produce zero mutants:
  expect_that(mut(xm, rep(0, length(xm))),
              equals(numeric(0)))
})

test_that("Multidimensional mutation", {
  k <- 3

  # Matrix that is not positive definite:
  expect_that(make_mutation(matrix(2, k, k) - diag(k)), throws_error())

  # First, start in the covariance-free case
  set.seed(1)
  vcv <- diag(k) * runif(k)
  mut <- make_mutation(vcv)

  xm <- matrix(rnorm(k, seq_len(k) * 10), nrow=1)
  n  <- 1000

  # Orient the matrix correctly:
  expect_that(mut(t(xm), n), throws_error())
  expect_that(mut(t(xm), rep(n, k)), throws_error())

  x <- mut(xm, n)
  expect_that(x, is_a("matrix"))
  expect_that(dim(x), equals(c(n, k)))

  # Easy to check the conditional distributions conform to the
  # normals that we expect in this case...
  sd <- sqrt(diag(vcv))
  for (i in seq_len(k))
    expect_that(ks.test(x[,i], pnorm, xm[i], sd[i])$p.value,
                is_greater_than(1/20))

  # Add some covariance for extra good:
  set.seed(1)
  vcv[lower.tri(vcv, FALSE)] <- vcv[upper.tri(vcv, FALSE)] <-
    runif(k * (k - 1) / 2, -.5, .5)
  mut <- make_mutation(vcv)

  x <- mut(xm, n)
  expect_that(x, is_a("matrix"))
  expect_that(dim(x), equals(c(n, k)))

  sd <- sqrt(diag(vcv))
  for (i in seq_len(k))
    expect_that(ks.test(x[,i], pnorm, xm[i], sd[i])$p.value,
                is_greater_than(1/20))

  expect_that(max(abs(cov(x) - vcv)), is_less_than(1/sqrt(n)))

  # Two resident types:
  xm <- rbind(xm, -xm)
  n  <- c(1000, 501)
  x <- mut(xm, n)
  expect_that(x, is_a("matrix"))
  expect_that(dim(x), equals(c(sum(n), k)))

  type <- rep(seq_along(n), n)
  for (i in seq_len(k))
    for (j in seq_along(n))
      expect_that(ks.test(x[type == j,i], pnorm, xm[j,i], sd[i])$p.value,
                  is_greater_than(1/20))

  # Produce zero mutants:
  expect_that(mut(xm, rep(0, nrow(xm))),
              equals(matrix(NA_real_, nrow=0, ncol=k)))
})
