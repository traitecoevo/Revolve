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
                                 "capacity", "equilibrium",
                                 "single_equilibrium",
                                 "parameters")))
  expect_that(m$fitness,     is_a("function"))
  expect_that(m$competition, is_a("function"))
  expect_that(m$capacity,    is_a("function"))
  expect_that(m$equilibrium, is_a("function"))
  expect_that(m$single_equilibrium, is_a("function"))
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
              equals(rbind(gaussian(xx, 0, p$s2_C))))
  expect_that(m$competition(0, xx),
              equals(cbind(gaussian(0, xx, p$s2_C))))
  set.seed(10)
  xp <- rnorm(1)
  expect_that(m$competition(xp, xx),
              equals(cbind(gaussian(xp, xx, p$s2_C))))

  ## Multiple response case:
  expect_that(m$competition(xx, xx),
              equals(outer(xx, xx, gaussian, p$s2_C)))
  zz <- sample(xx, 10)
  expect_that(m$competition(xx, zz),
              equals(outer(zz, xx, gaussian, p$s2_C)))
  expect_that(m$competition(zz, xx),
              equals(outer(xx, zz, gaussian, p$s2_C)))

  ## Now, try this on an established population:
  set.seed(100)
  x <- rnorm(5)
  y <- runif(length(x)) * m$capacity(x)

  w.cmp <- p$r * (1 - sum(y * m$competition(xp, x)) / m$capacity(xp))
  expect_that(m$fitness(xp, x, y), equals(w.cmp))

  ## Multiple mutants at once:
  set.seed(100)
  xp2 <- rnorm(2)

  w2.cmp1 <- sapply(xp2, m$fitness, x, y)
  w2.cmp2 <- p$r * (1 - colSums(y * m$competition(xp2, x)) / m$capacity(xp2))

  expect_that(w2.cmp1, equals(w2.cmp2))
  expect_that(m$fitness(xp2, x, y), equals(w2.cmp1))
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

  ## Single equlibrium
  eq <- m$single_equilibrium()
  expect_that(m$fitness(eq$x, eq$x, eq$y), equals(0))

  ## We can also get there using the equilibrium function:
  eq2 <- m$equilibrium(list(x=eq$x, y=1))
  expect_that(eq2, equals(eq))

  ## Multispecies equilibrium
  eq2 <- m$equilibrium(sys_split(eq, 0.1))
  expect_that(m$fitness(eq2$x, eq2$x, eq2$y), equals(c(0, 0)))

  eq3 <- m$equilibrium(sys_split(eq, c(-.1, .2)))
  expect_that(m$fitness(eq3$x, eq3$x, eq3$y), equals(c(0, 0)))
})
