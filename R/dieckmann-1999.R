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
    r * (1 - sum(y * competition(x_new, x)) / capacity(x_new))
  ## Model details:
  competition <- function(x_new, x) {
    if (length(x_new) > 1 && length(x) > 1) # (TODO)
      stop("Only one of x_new or x can be vector (at present)")
    exp(- (x_new - x)^2 / (2 * s2_C))
  }
  capacity <- function(x_new)
    K_0 * exp(- x_new^2 / (2 * s2_K))

  ret <- list(fitness     = fitness,
              competition = competition,
              capacity    = capacity,
              parameters  = parameters)
  class(ret) <- "model" # TODO: less generic name!
  ret
}
