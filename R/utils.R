##' Utility for taking care of getting/setting/validating parameters.
##'
##' Eventually I'll add a hook here for checking that parameters are
##' allowable, for setting ranges, etc.
##' @title Utility for Parameter Handling
##' @param defaults List of default parameter values.  Allowable
##' parameter names are also taken from this list
##' @param where Environment into which parameters will be
##' set/retrieved.  If ommited a new environment is generated (so
##' access is only via get/set).
##' @return An object of class \code{parameters}
##' @author Rich FitzJohn
##' @export
make_parameters <- function(defaults, where=new.env()) {
  check.names <- function(x) {
    names <- names(x)
    if (length(x) > 0 &&
        (is.null(names) || any(is.na(names)) || any(names == "") ||
         any(duplicated(names))))
      stop("Names must be present, non-empty and non duplicated")
  }
  check <- function(parameters) {
    check.names(parameters)
    if (!all(names(parameters) %in% names))
      stop("Unknown parameters: ", setdiff(names(parameters), names))
    invisible(TRUE)
  }
  set <- function(parameters=list()) {
    check(parameters)
    for (v in names(parameters))
      assign(v, parameters[[v]], where)
  }
  get <- function()
    structure(lapply(names, base::get, where), names=names)

  check.names(defaults)
  names <- names(defaults)
  set(defaults)
  rm(defaults)

  ret <- list(set=set, get=get, names=function() names)
  class(ret) <- "parameters"
  ret
}

##' @S3method names parameters
names.parameters <- function(x, ...)
  x$names()
##' @S3method names<- parameters
`names<-.parameters` <- function(x, value)
  stop("Cannot set names of a parameters object")

##' Bin x into a grid, where x may have different densities.
##'
##' @title Bin Values
##' @param x Vector of locations
##' @param y Vector of densities/counts along x
##' @param grid Grid over which to bin the values.
##' @return A vector of length \code{length(grid) - 1} being the
##' number of things in each interval.
##' @author Rich FitzJohn
##' @export
to_grid <- function(x, y, grid) {
  i <- findInterval(x, grid, rightmost.closed=TRUE)
  unname(sapply(split(y, factor(i, seq_len(length(grid)-1))), sum))
}

##' Build a function that steps the system through one
##' generation. Order of events is: (1) resolve fitness (2) introduce
##' mutants.
##'
##' @title Step System Forward in Time
##' @param fitness Function for computing fitness.  Must take
##' arguments \code{x_new} (phenotypes to compute fitness for),
##' \code{x} (resident phenotypes) and \code{y} (resident abundances)
##' and return the per-capita growth rate (g).  The population will step
##' forward with an Euler step, so that the new \code{y} will be
##' \code{y + y * g * dt}.
##' @param mutation Function for generating mutants.  Must take
##' arguments \code{traits} (resident phenotypes) and \code{mutants}
##' (the number of mutants to generate, for each phenotype).  See
##' \code{\link{make_mutation}}, which generates a useful function.
##' @param dt Size of the time step
##' @param mu Mutation rate.  On average there will be \code{mu*dt*y}
##' mutations, with the actual number drawn from Poisson
##' distribution.
##' @param y_initial Abundance at which to introduce new mutants.
##' @return A Function that moves the system forward in time; it takes
##' a list with elements \code{x} (traits), \code{y} (abundances) and
##' \code{t} (time) and returns a list with these updated to the next
##' steps.
##' @author Rich FitzJohn
##' @export
make_step <- function(fitness, mutation, dt, mu, y_initial) {
  function(sys) {
    # births/time/individual * individuals * time -> births
    dy <- fitness(sys$x, sys$x, sys$y) * dt * sys$y
    # mutations/time/individual * individuals * time -> mutation
    mutants <- mutation(sys$x, rpois(length(sys$y), mu * dt * sys$y))

    sys$y <- sys$y + dy
    if (length(mutants) > 0) {
      sys$x <- c(sys$x, mutants)
      sys$y <- c(sys$y, rep_len(y_initial, length(mutants)))
    }
    sys$t <- sys$t + dt

    sys
  }
}

##' @export
##' @rdname make_step
make_step_continuous <- function(fitness, mutation, dt, mu, y_initial) {
  step_deterministic <- make_step_deterministic(fitness)
  function(sys) {
    y0 <- sys$y
    sys <- step_deterministic(sys, dt)

    # Work out mutational input from the integrated number of
    # individuals (quick and dirty) via trapezium rule times the
    # mutation rate (per time).
    # mutations/time/individual * individuals * time -> mutation
    mutants <- mutation(sys$x,
                        rpois(length(y0), mu * dt * (y0 + sys$y)/2))
    if (length(mutants) > 0) {
      sys$x <- c(sys$x, mutants)
      sys$y <- c(sys$y, rep_len(y_initial, length(mutants)))
    }
    sys
  }
}

##' @export
##' @rdname make_step
##' @importFrom deSolve lsoda
make_step_deterministic <- function(fitness) {
  ## NOTE: This will actually work with non-scalar dt, which is
  ## probably not desirable given that we write to sys directly.
  function(sys, dt) {
    derivs <- function(t, y, ...)
      list(fitness(sys$x, sys$x, y) * y)
    t0 <- sys$t
    t1 <- t0 + dt
    tmp <- lsoda(sys$y, c(t0, t1), derivs)[-1,,drop=FALSE]
    sys$t <- unname(tmp[,1])
    sys$y <- unname(tmp[,-1])
    sys
  }
}

##' Cleanup phenotypes that have dropped below a threshold abundance
##'
##' @title Cleanup Rare Phenotypes
##' @param eps Abundance below which species are considered extinct.
##' @author Rich FitzJohn
##' @export
make_cleanup <- function(eps) {
  function(sys) {
    keep <- sys$y >= eps
    if (!all(keep)) {
      sys$x <- sys$x[keep]
      sys$y <- sys$y[keep]
    }
    sys
  }
}

##' Run a system forward in time.
##'
##' @title Run System
##' @param sys a list with elements \code{x} (traits), \code{y} (abundances) and
##' \code{t} (time)
##' @param n_steps The number of steps to run the system for.
##' @param step Step function (see \code{\link{make_step}})
##' @param cleanup Cleanup function (optional - by default do no cleanup)
##' @param print_every Interval at which to print some progress information
##' @author Rich FitzJohn
##' @export
run <- function(sys, n_steps, step, cleanup=identity, print_every=0) {
  res <- vector("list", length(n_steps+1))
  res[[1]] <- sys
  for (i in seq_len(n_steps)) {
    if (print_every > 0 && i %% print_every == 0)
      message(i)
    sys <- step(sys)
    sys <- cleanup(sys)
    res[[i+1]] <- sys
  }
  res
}
