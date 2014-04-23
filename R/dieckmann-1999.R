##' Model implementing "On the origin of species by sympatric
##' speciation" by Dieckmann & Doebeli, Nature 1999.
##'
##' Still a test case only; much to come.
##' @title Dieckmann & Doebeli 1999
##' @param r Intrinsic birth rate
##' @param K_0 Scaling factor giving maximum population density
##' @param s2_C Width of the competition kernel
##' @param s2_K Width of the resource distribution kernel
##' @return TBD...
##' @author Rich FitzJohn
##' @export
make_dieckmann_1999 <- function(r=1, K_0=500, s2_C=0.16, s2_K=1) {
  defaults <- list(r=r, K_0=K_0, s2_C=s2_C, s2_K=s2_K)
  parameters <- make_parameters(defaults, environment())

  ## Fitness calculation:
  fitness <- function(x_new, x, y)
    r * (1 - colSums(y * competition(x_new, x)) / capacity(x_new))
  ## Model details:
  ##
  ## Returns a matrix with length(x) rows and length(n_new) columns;
  ## so reading down the column 'i' tells you about the competition
  ## experienced by mutant strategy x_new[i].
  competition <- function(x_new, x)
    exp(- outer(x, x_new, "-")^2 / (2 * s2_C))
  capacity <- function(x_new)
    K_0 * exp(- x_new^2 / (2 * s2_K))

  equilibrium <- function(sys, ...) {
    equilibrium_sys(sys, fitness, ...)
  }
  single_equilibrium <- function() {
    x <- 0 # fixed here
    sys(x, capacity(x))
  }

  ret <- list(fitness     = fitness,
              competition = competition,
              capacity    = capacity,
              equilibrium = equilibrium,
              single_equilibrium = single_equilibrium,
              parameters  = parameters)
  class(ret) <- "model" # TODO: less generic name!
  ret
}
