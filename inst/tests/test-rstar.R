source("helper-Revolve.R")

context("R star")

mat2 <- rstar_matrices(rstar_mat_2, rstar_mat_2)

test_that("Model parameters", {
  m <- make_rstar(mat2)
  defaults <- list(r=1, m=1/4, D=1/4, S=c(1, 1))
  expect_that(names(m$parameters), equals(names(defaults)))
  expect_that(m$parameters$get(),  equals(defaults))
})

test_that("Model components", {
  m <- make_rstar(mat2)
  expect_that(names(m), equals(c("fitness",
                                 "equilibrium", "run", "parameters",
                                 "n", "k", "K", "C", "p", "Rstar",
                                 "single_equilibrium",
                                 "single_equilibrium_R",
                                 "run_fixed_density",
                                 "equilibrium_R",
                                 "derivs")))

  expect_that(m$fitness,     is_a("function"))
  expect_that(m$equilibrium, is_a("function"))
  expect_that(m$run,         is_a("function"))
  expect_that(m$parameters,  is_a("parameters"))
  expect_that(m$n,           equals(NA))
  expect_that(m$k,           is_identical_to(2L))
  expect_that(m$K,           is_a("function"))
  expect_that(m$C,           is_a("function"))
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
  expect_that(m$derivs,
              is_a("function"))
})

# Matrix handling:

test_that("Constant matrices", {
  K1 <- matrix(0.5, 1, 1)
  C1 <- matrix(0.1, 1, 1)
  x1 <- matrix(nrow=0, ncol=5)

  mat <- rstar_matrices(K1, C1)

  expect_that(mat$k, equals(1L))
  expect_that(mat$n, equals(NA))
  expect_that(mat$i.K, equals(integer(0)))
  expect_that(mat$i.C, equals(integer(0)))

  expect_that(mat$K(x1), equals(matrix(K1, 1, ncol(x1))))
  expect_that(mat$C(x1), equals(matrix(C1, 1, ncol(x1))))

  K2 <- matrix(runif(2), 2, 1)
  C2 <- matrix(runif(2), 2, 1)
  x2 <- matrix(nrow=0, ncol=5)

  mat <- rstar_matrices(K2, C2)

  expect_that(mat$k, equals(2L))
  expect_that(mat$n, equals(NA))
  expect_that(mat$i.K, equals(integer(0)))
  expect_that(mat$i.C, equals(integer(0)))

  expect_that(mat$K(x2), equals(matrix(K2, 2, ncol(x2))))
  expect_that(mat$C(x2), equals(matrix(C2, 2, ncol(x2))))
})

test_that("Fixed matries", {
  K <- matrix(runif(4), 2, 2)
  C <- matrix(runif(4), 2, 2)

  mat <- rstar_matrices_fixed(K, C)

  expect_that(mat$k, equals(2L))
  expect_that(mat$n, equals(2L))
  expect_that(mat$i.K, equals(integer(0)))
  expect_that(mat$i.C, equals(integer(0)))

  # With no arguments
  expect_that(mat$K(), equals(K))
  expect_that(mat$C(), equals(C))

  # With nonexistant arguments
  expect_that(mat$K(foo), equals(K))
  expect_that(mat$C(foo), equals(C))

  # With arguements that do exist but make no sense:
  expect_that(mat$K(runif(10)), equals(K))
  expect_that(mat$C(runif(10)), equals(C))
})

test_that("Identity matrices", {
  mat <- rstar_matrices(rstar_mat_1, rstar_mat_1)

  expect_that(mat$k, equals(1L))
  expect_that(mat$n, equals(NA))
  expect_that(mat$i.K, equals(1L))
  expect_that(mat$i.C, equals(2L))

  x1 <- matrix(runif(10), 2)
  expect_that(mat$K(x1), equals(x1[mat$i.K,,drop=FALSE]))
  expect_that(mat$C(x1), equals(x1[mat$i.C,,drop=FALSE]))

  mat <- rstar_matrices(rstar_mat_2, rstar_mat_2)

  expect_that(mat$k, equals(2L))
  expect_that(mat$n, equals(NA))
  expect_that(mat$i.K, equals(1:2))
  expect_that(mat$i.C, equals(3:4))

  x2 <- matrix(runif(20), 4)
  expect_that(mat$K(x2), equals(x2[mat$i.K,,drop=FALSE]))
  expect_that(mat$C(x2), equals(x2[mat$i.C,,drop=FALSE]))
})

test_that("Mixed identity and constant", {
  K1 <- matrix(0.5, 1, 1)
  C1 <- matrix(0.1, 1, 1)

  mat.K <- rstar_matrices(rstar_mat_1, C1)
  mat.C <- rstar_matrices(K1, rstar_mat_1)

  expect_that(mat.K$k, equals(1L))
  expect_that(mat.K$n, equals(NA))
  expect_that(mat.K$i.K, equals(1L))
  expect_that(mat.K$i.C, equals(numeric(0)))

  expect_that(mat.C$k, equals(1L))
  expect_that(mat.C$n, equals(NA))
  expect_that(mat.C$i.K, equals(numeric(0)))
  expect_that(mat.C$i.C, equals(1L))

  x <- matrix(runif(4), 1)

  expect_that(mat.K$K(x), equals(x))
  expect_that(mat.C$C(x), equals(x))

  expect_that(mat.K$C(x), equals(matrix(C1, 1, ncol(x))))
  expect_that(mat.C$K(x), equals(matrix(K1, 1, ncol(x))))
})

test_that("One resource", {
  m <- make_rstar(rstar_matrices(rstar_mat_1, rstar_mat_1), S=1)
  sys1 <- sys(matrix(0.5, nrow=2), y=1)

  ## Find the equilibrium with a single species (note that this is
  ## quite different to the D+D & Kisdi models because this is not the
  ## singular position for the model, but the demographic/resource
  ## equilibrium given a particular species trait).
  eq <- m$single_equilibrium(sys1$x)
  expect_that(m$derivs(0, c(eq$R, eq$y), sys1), equals(c(0, 0)))

  ## Computed numerically:
  eq2 <- m$equilibrium(sys1)
  expect_that(eq2, equals(eq[names(eq2)]))
})

# Results to check:

# Several species, 1 resource.
#   - species with the lowest requirement for the limiting resource
#   (lowest R*) will displace all other species.
