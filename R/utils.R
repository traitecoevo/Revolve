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

##' Turn model output into data to be plotted by \code{image}.
##'
##' @title Generate Heat Map of Community
##' @param res List of model states over time
##' @param xr Range of x values (if ommited will be fit to the data)
##' @param nx Number of x points to discretise over
##' @param sample Thin the time axis by this amount (10 means take one
##' point every 10 time steps).  This may not be ideal where time does
##' not proceed in even sized steps.
##' @author Rich FitzJohn
##' @export
discretise <- function(res, xr=NULL, nx=101, sample=10) {
  t <- sapply(res, "[[", "t")
  x <- lapply(res, "[[", "x")
  y <- lapply(res, "[[", "y")

  if (is.null(xr))
    xr <- unlist(range(x))
  xx <- seq(xr[1], xr[2], length.out=nx)

  j <- seq(1, length(res), by=sample)
  mat <- sapply(j, function(i)
                 to_grid(x[[i]], y[[i]], xx))
  mat[mat == 0] <- NA

  list(x=xx, t=t[j], y=mat)
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
make_step_continuous <- function(fitness, mutation, dt, mu, y_initial,
                                 ...) {
  step_deterministic <- make_step_deterministic(fitness, ...)
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
make_step_deterministic <- function(fitness, ...) {
  ## NOTE: This will actually work with non-scalar dt, which is
  ## probably not desirable given that we write to sys directly.
  ##
  ## NOTE: The time scaling argument here is useful in times where the
  ## rates of change are very very low so that the system wants to
  ## equilibrate over millions of time units rather than tens.
  function(sys, dt, hini=0, time_scale=1) {
    derivs <- function(t, y, ...)
      list(fitness(sys$x, sys$x, y) * y * time_scale)

    t0 <- sys$t
    t1 <- t0 + dt
    tmp <- lsoda(sys$y, c(t0, t1), derivs,
                 hmax=0, maxsteps=10000, hini=hini, ...)
    if (tmp[-1,1] < t1)
      stop("Did not complete integration")

    sys$t <- unname(tmp[-1,1]) * time_scale
    sys$y <- unname(tmp[-1,-1])
    sys$y[sys$y < 0] <- 0.0
    attr(sys, "step_size") <- attr(tmp, "rstate")[[1]] * time_scale
    sys
  }
}

##' \code{make_step_equilibrium} does a very simple minded attempt to
##' converge on the equilibrium of a system by running it for a long
##' time.  It advances the system by \code{dt} (should be a largeish
##' number) a number of times until the population size stabilises.
##' This could be replaced by a multidimensional root finding
##' approach.
##'
##' @export
##' @param method method to find equilibrium.  Must be one of
##' "runsteady" or "nleqslv" (or a contraction).
##' @param ... Additional arguments passed through to
##' either \code{rootSolve::runsteady} or \code{nleqslv::nleqslv}.
##' @rdname make_step
make_step_equilibrium <- function(fitness, ..., method) {
  switch(match.arg(method, c("runsteady", "nleqslv")),
         runsteady= make_step_equilibrium_runsteady,
         nleqslv  = make_step_equilibrium_nleqslv)(fitness, ...)
}

make_step_equilibrium_naive <- function(fitness, eps, dt=100, max_iter=100,
                                  ...) {
  step <- make_step_deterministic(fitness, ...)
  function(sys) {
    .sys <- sys
    time_scale <- 1
    for (i in seq_len(max_iter)) {
      sys$t <- 0
      y0 <- sys$y
      sys <- step(sys, dt, time_scale=time_scale)
      dy <- sys$y - y0

      if (all(dy < eps | sys$y < eps))
        return(sys)

      time_scale <- max(1, attr(sys, "step_size"))
    }
    stop("Failed to converge on equilibrium (maximum iterations reached)")
  }
}

##' @importFrom rootSolve runsteady
make_step_equilibrium_runsteady <- function(fitness, max_time=1e5, ...) {
  function(sys) {
    derivs <- function(t, y, pars)
      list(fitness(sys$x, sys$x, y) * y)
    sys$y <- runsteady(sys$y, derivs, times=c(0, max_time), ...)$y
    sys$y[sys$y < 0] <- 0
    sys
  }
}

##' @importFrom nleqslv nleqslv
make_step_equilibrium_nleqslv <- function(fitness, ...) {
  function(sys) {
    target <- function(y)
      fitness(sys$x, sys$x, y) * y
    sys$y <- nleqslv(sys$y, target, ...)$x
    sys$y[sys$y < 0] <- 0
    sys
  }
}

## Calculate the rate w_k = mu * \hat n * b(s_k) for the emergence of
## a mutant.  We are going to drop b(s_k), partly because it is
## constant for all genotypes and partly because it seems daft to
## invoke separation of timescales while mixing timescales.
##
## In theory, here we should compute the fitness of the individual and
## accept only if runif(1) < f(s') / b(s') but we're skipping this
## step for now.

##' @param step_equilibrium Function generated by
##' \code{make_step_equilibrium}.
##' @export
##' @rdname make_step
make_step_mutation_limited <- function(step_equilibrium,
                                       mutation, mu, y_initial) {
  force(step_equilibrium)
  function(sys) {
    t0 <- sys$t
    sys <- step_equilibrium(sys)

    ## NOTE: Notation from Ito (2007) -- w is not a good choice here.
    w <- sys$y * mu

    ## Update evolutionary time; ignore the ecological time.
    sys$t <- t0 + rexp(1, sum(w))

    ## Create a single mutation:
    i <- sample(length(w), 1, prob=w)
    mutant <- mutation(sys$x, as.integer(seq_along(w) == i))
    sys$x <- c(sys$x, mutant)
    sys$y <- c(sys$y, y_initial)

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
      message(sprintf("%d: t = %2.2f, n = %d", i, sys$t, length(sys$y)))
    sys <- step(sys)
    sys <- cleanup(sys)
    res[[i+1]] <- sys
  }
  res
}

sys <- function(x, y, t=0) {
  list(x=x, y=y, t=t)
}

# Excuse the slightly odd name: this is designed to be used from
# within functions called 'equilibrium', so we need to be able to
# distinguish it somehow.  May tidy that up later...
equilibrium_ <- function(dydt, x, y, method="runsteady",
                        init_time=200, max_time=1e5) {
  method <- match.arg(method, c("runsteady", "nleqslv", "simple"))

  if (init_time > 0) {
    y <- unname(lsoda_nolist(y, c(0, init_time), dydt, x)[2,-1])
  }

  if (method == "runsteady") {
    dydt.deSolve <- function(...) list(dydt(...))
    ans <- runsteady(y, dydt.deSolve, x, times=c(0, Inf), hmax=1)$y
  } else if (method == "nleqslv") {
    ans <- nleqslv(y, function(y) dydt(0, y, x), global="none")$x
  } else if (method == "simple") {
    ans <- unname(lsoda_nolist(y, c(init_time, max_time), dydt, x)[2,-1])
  }
  ans
}

lsoda_nolist <- function(y, times, func, ...) {
  lsoda(y, times, function(...) list(func(...)), ...)
}

colMins <- function(x) {
  apply(x, 2, min)
}
rowMins <- function(x) {
  apply(x, 1, min)
}

quadratic_roots <- function(a, b, c) {
  (-b + c(-1, 1) * sqrt(b*b - 4*a*c))/(2 * a)
}
