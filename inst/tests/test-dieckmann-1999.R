source("helper-Revolve.R")

context("Dieckmann & Dobelli 1999")

test_that("Model parameters", {
  m <- make_dieckmann_1999()

  defaults <- list(r=1, K_0=500, s2_C=0.16, s2_K=1)
  expect_that(names(m$parameters), equals(names(defaults)))
  expect_that(m$parameters$get(), equals(defaults))

  p_new <- list(K_0=1000)
  m$parameters$set(p_new)
  expect_that(m$parameters$get(), equals(modifyList(defaults, p_new)))

  expect_that(get("K_0", environment(m$fitness)),
              equals(p_new[["K_0"]]))
})

test_that("Model components", {
  m <- make_dieckmann_1999()
  expect_that(names(m), equals(c("fitness", "competition",
                                 "capacity", "parameters")))
  expect_that(m$fitness,     is_a("function"))
  expect_that(m$competition, is_a("function"))
  expect_that(m$capacity,    is_a("function"))
  expect_that(m$parameters,  is_a("parameters"))
})

test_that("Components work", {
  m <- make_dieckmann_1999()
  p <- m$parameters$get()

  ## Gaussian distribution with a peak at 1, in terms of variance not
  ## SD.
  gaussian <- function(x, mean, variance)
    dnorm(x, mean, sqrt(variance)) / dnorm(0, 0, sqrt(variance))

  ## Carrying capacity:
  xx <- seq(-5, 5, length.out=101)
  yy <- m$capacity(xx)
  expect_that(yy, equals(p$K * gaussian(xx, 0, p$s2_K)))

  expect_that(m$competition(xx, 0),
              equals(gaussian(xx, 0, p$s2_C)))
  expect_that(m$competition(0, xx),
              equals(gaussian(0, xx, p$s2_C)))
  set.seed(1)
  xp <- rnorm(1)
  expect_that(m$competition(xp, xx),
              equals(gaussian(xp, xx, p$s2_C)))

  expect_that(m$competition(xx, xx), throws_error())

  ## Now, try this on an established population:
  set.seed(1)
  x <- rnorm(5)
  y <- runif(length(x)) * m$capacity(x)

  w.cmp <- p$r * (1 - sum(y * gaussian(xp, x, p$s2_C)) /
                  (p$K_0 * gaussian(0, xp, p$s2_K)))
  expect_that(m$fitness(xp, x, y), equals(w.cmp))
})

## For a single case at equilibrium density, the growth rate should be
## zero:
test_that("Single species equilibrium", {
  m <- make_dieckmann_1999()
  p <- m$parameters$get()

  ## Per capita rate of growth in empty environment is zero:
  expect_that(m$fitness(0, 0, 0), equals(p$r))

  ## Growth rate at capacity is zero
  f <- function(x) m$fitness(x, x, m$capacity(x))
  n <- 101
  expect_that(sapply(seq(-5, 5, length.out=n), f), equals(rep(0, n)))
})
