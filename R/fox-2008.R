## Extending this model to multiple types (i.e., more than two) is
## going to be hard, because the 'y' matrix (here, 'Y') is as much of
## a property of the species as 'u' (the trait that can evolve).  But
## 'Y' needs expanding out to be a length(R) * length(u) matrix,
## rather than just being 2 * 2.  Georges is working on this.
##
## At present, parameters follow the paper, and defaults are as given
## in the paper.  The number of resources is fixed at two.
##
## The death rates 'd' is also a property of the species, at least in
## the general case.

##' Model implementing "Character Convergence under Competition for
##' Nutritionally Essential Resources.", Fox and Vasseur, The American
##' Naturalist, 2008.
##'
##' Test case...
##' @title Fox and Vasseur 2008
##' @param Y See paper
##' @param d See paper
##' @param D See paper
##' @param S See paper
##' @author Georges Kunstler and Rich FitzJohn
##' @export
make_fox_2008 <- function(Y=rbind(c(0.5, 1), c(1, 0.5)),
                          d=0.1, D=0.1, S=c(1, 1)) {
  defaults <- list(Y=Y, d=d, D=D, S=S)
  parameters <- make_parameters(defaults, environment())

  dydt <- function(t, ode.y, x, ...) {
    y <- ode.y[1:2]
    R <- ode.y[3:4]
    gj <- g(x, R)
    dRdt <- D * (S - R) - colSums(y * gj / Y)
    dydt <- y * (gj - d)
    c(dydt, dRdt)
  }

  ## Given traits 'x' and initial densities 'y', compute the
  ## equilibrium.  Return the equilibrium resource availability at the
  ## same time.
  equilibrium <- function(x, y, method="runsteady",
                          init_time=200,
                          max_time=1e5) {
    # Initial resources are set to 'S'
    ans <- equilibrium_(dydt, x, c(S, y), method, init_time, max_time)
    list(y=ans[1:2], R=ans[3:4])
  }

  run <- function(x, y, times) {
    y0 <- c(y, S)
    res <- lsoda_nolist(y0, times, dydt, x)
    list(t=res[,1],
         y=res[,2:3],
         R=res[,4:5])
  }

  g <- function(x, R) {
    pmin(Y[1,] * x       * R[1],
         Y[2,] * (1 - x) * R[2])
  }

  # Analytic solutions:
  # Equilibrium resource density when limiting
  analytic.resource.limiting <- function(x, sp=1) {
    cbind(d / (Y[1,sp] * x),       # B1a, B1c
          d / (Y[2,sp] * (1 - x))) # B1b, B1d
  }
  # Equilibrium resource density when not limiting
  analytic.resource.nonlimiting <- function(x, sp=1) {
    cbind(1 - Y[2,sp]/Y[1,sp] + d / (Y[1,sp] * (1 - x)),
          1 - Y[1,sp]/Y[2,sp] + d / (Y[2,sp] * x))
  }
  analytic <- list(resource=list(
                     limiting=analytic.resource.limiting,
                     nonlimiting=analytic.resource.nonlimiting))

  ret <- list(equilibrium  = equilibrium,
              run          = run,
              parameters   = parameters,
              analytic     = analytic)
  class(ret) <- "model"
  ret
}
