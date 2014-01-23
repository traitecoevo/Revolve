library(Revolve)

x.initial <- -0.3  # initial location
y.initial <- 1e-3  # initial mutant frequency
y.eps     <- 1e-10 # kill threshold
dt        <- 1     # Size of the time step (discrete gen == Euler?)
mutation.sd <- 1/20

make.deaths <- function(eps) {
  function(sys) {
    keep <- sys$y >= eps
    if (!all(keep)) {
      sys$x <- sys$x[keep]
      sys$y <- sys$y[keep]
    }
    sys
  }
}

make.births <- function(mutation, rate, y.initial) {
  function(sys) {
    n <- rpois(length(sys$x), sys$y * rate)
    x_new <- mutation(sys$x, n)
    sys$x <- c(sys$x, x_new)
    sys$y <- c(sys$y, rep(y.initial, length(x_new)))
    sys
  }
}

## NOTE: I think this will fail when we have more than one
## resident...
make.step <- function(fitness, dt) {
  function(sys) {
    dydt <- fitness(sys$x, sys$x, sys$y)
    sys$y <- sys$y * (1 + dydt * dt)
    sys$t <- sys$t + dt
    sys
  }
}

m <- make_dieckmann_1999()
deaths <- make.deaths(y.eps)
births <- make.births(make_mutation(mutation.sd),
                      m$parameters$get()$r * dt,
                      y.initial)
step <- make.step(m$fitness, dt)

sys <- list(x=x.initial, y=y.initial, t=0)
res <- list(sys)
for (i in 1:30) {
  sys <- step(sys)
  sys <- deaths(sys)
  # sys <- births(sys)
  res <- c(res, list(sys))
}

res <- as.data.frame(do.call(rbind, res))
res[] <- sapply(res, as.numeric)

plot(y ~ t, res, type="o")
