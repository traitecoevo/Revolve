##' Model implementing "Evolutionary branching under asymmetric
##' competition" by Kisdi, Journal of Theoretical Biology, 1999.
##'
##' Test case...
##' @title Kisdi 1999
##' @param intrinsic_growth_type String, indicating which type of
##' intrinsic growth to use; can be "linear", "gaussian" or
##' "convex" (or any contraction of these). Intrinsic growth is the
##' population growth rate in the absence of competition (including
##' self-competition). Units: inverse time
##' @param c Death rate, to be scaled by competition function.
##' Units: inverse time * inverse population
##' @param v Positioning parameter for competition function
##' @param k Trait scaling parameter for competition function.
##' Units: inverse trait.
##' @param beta Intercept of linear intrinsic growth function
##' @param b Slope of linear intrinsic growth function (also used in
##' convex intrinsic growth function)
##' @param a Scaling factor for gaussian intrinsic growth function
##' (also used in convex intrinsic growth function)
##' @param m Mean of gaussian intrinsic growth function
##' @param s2 Variance of gaussian intrinsic growth function
##' @param d Parameter of convex intrinsic growth function
##' @author Rich FitzJohn
##' @export
##' @importFrom numDeriv grad
make_kisdi_1999 <- function(intrinsic_growth_type="linear", c=2, v=1.2, k=4,
                            beta=1, b=1, a=1, m=0, s2=1, d=3.5) {
  defaults <- list(c=c, v=v, k=k, beta=beta, b=b, a=a, m=m, s2=s2, d=d)
  parameters <- make_parameters(defaults, environment())

  fitness <- function(x_new, x, y)
    intrinsic_growth(x_new) - colSums(y * competition(x_new, x))
  ## Equation 1:
  competition <- function(x_new, x)
    c * (1 - 1 / (1 + v * exp(- k * (-outer(x, x_new, "-")))))

  equilibrium <- function(sys, ...) {
    equilibrium_sys(sys, fitness, ...)
  }
  single_equilibrium <- function() {
    switch(intrinsic_growth_type,
           linear=single_equilibrium_linear,
           gaussian=single_equilibrium_gaussian,
           stop("Not implemented for type ", intrinsic_growth_type))()
  }
  single_equilibrium_linear <- function() {
    # Linear: monomorphic singular strategy at
    # beta / b + a(0) /  a'(0)
    alpha <- function(xp) {
      drop(competition(0, 0 - xp))
    }
    x <- beta / b + alpha(0) / grad(alpha, 0)
    y <- intrinsic_growth(x) / alpha(0)
    sys(x, y)
  }
  single_equilibrium_gaussian <- function() {
    # Gaussian: monomorphic singular strategy at
    #   m - a'(0) / a(0) * s2
    alpha <- function(xp) {
      drop(competition(0, 0 - xp))
    }
    x <- m - grad(alpha, 0) / alpha(0) * s2
    y <- intrinsic_growth(x) / alpha(0)
    sys(x, y)
  }

  ## Equation 12:
  intrinsic_growth_linear <- function(x)
    beta - b * x
  ## Equation 14:
  intrinsic_growth_gaussian <- function(x)
    a * exp(-(x - m)^2 / (2 * s2))
  ## Equation 17:
  intrinsic_growth_convex <- function(x)
    -a - b * (x - sqrt(x*x + d))
  intrinsic_growth_type <- match.arg(intrinsic_growth_type,
                          c("linear", "gaussian", "convex"))
  intrinsic_growth <- switch(intrinsic_growth_type,
                     linear=intrinsic_growth_linear,
                     gaussian=intrinsic_growth_gaussian,
                     convex=intrinsic_growth_convex)

  ret <- list(fitness            = fitness,
              competition        = competition,
              intrinsic_growth   = intrinsic_growth,
              parameters         = parameters,
              equilibrium        = equilibrium,
              single_equilibrium = single_equilibrium)
  class(ret) <- "model"
  ret
}
