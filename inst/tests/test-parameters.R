source("helper-Revolve.R")

context("Parameter helper")

test_that("Corner cases fail", {
  # Missing names
  expect_that(make_parameters(list(1)), throws_error())
  # Blank names
  tmp <- list(1)
  names(tmp) <- ""
  expect_that(make_parameters(tmp), throws_error())
  # Duplicated names
  expect_that(make_parameters(list(a=1, a=2)), throws_error())
})

test_that("Basic usage", {
  pars <- list(a=1, b=pi)
  p <- make_parameters(pars)
  expect_that(names(p), equals(names(pars)))
  expect_that(p$get(), is_identical_to(pars))
})

test_that("Setting parameters corner cases", {
  pars <- list(a=1, b=pi)
  p <- make_parameters(pars)

  # Empty set does not change anything:
  p$set()
  expect_that(p$get(), is_identical_to(pars))
  p$set(list())
  expect_that(p$get(), is_identical_to(pars))

  # Non-existant parameters will fail to set:
  expect_that(p$set(list(c=1)), throws_error())
  # Un-named parameters
  expect_that(p$set(list(1)), throws_error())
})

test_that("Setting parameters", {
  pars <- list(a=1, b=pi)
  p <- make_parameters(pars)

  pars1 <- list(a=exp(1))
  pars2 <- list(b=sqrt(2))
  pars3 <- list(b=runif(1), a=runif(1))

  expect_that(p$get(), is_identical_to(pars))

  p$set(pars1)
  expect_that(p$get(), is_identical_to(modifyList(pars, pars1)))

  p$set(pars2)
  expect_that(p$get(),
              is_identical_to(modifyList(modifyList(pars, pars1), pars2)))

  p$set(pars3)
  expect_that(p$get(),
              is_identical_to(modifyList(pars, pars3)))
})

test_that("New environment", {
  e <- new.env()
  expect_that(ls(e), is_identical_to(character(0)))

  pars <- list(a=1, b=pi)
  p <- make_parameters(pars, e)
  expect_that(ls(e), is_identical_to(names(pars)))

  expect_that(get("a", e), is_identical_to(pars$a))

  p$set(list(a=pi))
  expect_that(get("a", e), is_identical_to(pi))
})
