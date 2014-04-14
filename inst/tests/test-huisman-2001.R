source("helper-Revolve.R")

context("Huisman and Weissing 1999")

test_that("Model parameters", {
  m <- make_huisman_2001()

  defaults <- list(r=1, m=1/4, D=1/4, S=c(1, 1))
  expect_that(names(m$parameters), equals(names(defaults)))
  expect_that(m$parameters$get(),  equals(defaults))
})

test_that("Model components", {
  m <- make_huisman_2001()
  expect_that(names(m), equals(c("fitness",
                                 "equilibrium", "run", "parameters",
                                 "n", "k", "K", "C", "p", "Rstar",
                                 "single_equilibrium",
                                 "single_equilibrium_R",
                                 "run_fixed_density",
                                 "equilibrium_R")))

  expect_that(m$fitness,     is_a("function"))
  expect_that(m$equilibrium, is_a("function"))
  expect_that(m$run,         is_a("function"))
  expect_that(m$parameters,  is_a("parameters"))
  expect_that(m$n,           equals(NA))
  expect_that(m$k,           is_identical_to(2L))
  expect_that(m$K,           is_identical_to(huisman_mat_2))
  expect_that(m$C,           is_identical_to(huisman_mat_2))
  expect_that(m$p,           is_a("function"))
  expect_that(m$Rstar,       is_a("function"))
  expect_that(m$single_equilibrium,
              is_a("function"))
  expect_that(m$single_equilibrium_R,
              is_a("function"))
  expect_that(m$run_fixed_density,
              is_a("function"))
  expect_that(m$equilibrium_R,
              is_a("function"))
})

# Results to check:

# Several species, 1 resource.
#   - species with the lowest requirement for the limiting resource
#   (lowest R*) will displace all other species.
