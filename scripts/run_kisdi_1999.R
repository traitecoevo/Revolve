library(Revolve)

## Some extra shared parameters:
y_initial <- 1e-3     # initial mutant frequency
dt        <- 0.1      # Size of the time step
mu        <- 1.0      # Mutation rate
s2_mu     <- (1/10)^2 # Mutational variance

## Here, we need to decide to run in always-equilibrium mode
## (separation of timescales) or do a stochastic assembly.  Both
## should give pretty similar results.

m <- make_kisdi_1999()

## Initial system state; a single individual of a single strategy with
## trait -0.3:
sys0 <- list(x=-1, y=y_initial, t=0)

## Remove species with < 1 individual.
cleanup <- make_cleanup(1e-3)
## Mutations normally distributed with variance s2_mu.
mutation <- make_mutation(s2_mu)

step.d <- make_step(m$fitness, mutation, dt, mu, y_initial)
step.c <- make_step_continuous(m$fitness, mutation, 1, mu, y_initial)

step.e <- make_step_equilibrium(m$fitness, 1e-5)
step.e(sys0)

set.seed(1)
res.d <- run(sys0, 10000, step.d, cleanup, print_every=1000)
res.c <- run(sys0, 1000,  step.c, cleanup, print_every=100)

discretise <- function(res, xr, nx=101, sample=10) {
  xx <- seq(xr[1], xr[2], length.out=nx)

  j <- seq(1, length(res), by=sample)
  mat <- sapply(j, function(i)
                 to_grid(res[[i]]$x, res[[i]]$y, xx))
  mat[mat == 0] <- NA

  list(x=xx, t=sapply(res[j], "[[", "t"), y=mat)
}

## Plot the community over time; binned into
xr <- c(-8, 1)
cols <- grey(((32:0)/32)^2)
col <- "#00000055"

sys <- res.d[[length(res.d)]]
img <- discretise(res.d, xr)

op <- par(mfrow=c(3, 1), mar=c(2.6, 4.6, .5, .5), oma=c(2.1, 0, 0, 0))
## Traits and population density through time:
image(img$x, img$t, img$y, col=cols, ylab="Time", las=1); box()
## Resident types and their densities:
plot(sys, xlim=xr, pch=19, col=col, ylab="Number", las=1)
## Fitness of the resident types:
plot(img$x, m$fitness(img$x, sys$x, sys$y), type="l", ylab="Fitness", las=1)
points(sys$x, m$fitness(sys$x, sys$x, sys$y), col=col)
abline(h=0)
mtext("Trait", 1, 3, xpd=NA)
par(op)

sys <- res.c[[length(res.c)]]
img <- discretise(res.c, xr)

op <- par(mfrow=c(4, 1), mar=c(2.6, 4.6, .5, .5), oma=c(2.1, 0, 0, 0))
## Traits and population density through time:
image(img$x, img$t, img$y, col=cols, ylab="Time", las=1); box()
## Resident types and their densities:
plot(sys, xlim=xr, pch=19, col=col, ylab="Number", las=1)
## Fitness of the resident types:
plot(img$x, m$fitness(img$x, sys$x, sys$y), type="l", ylab="Fitness", las=1)
points(sys$x, m$fitness(sys$x, sys$x, sys$y), col=col)
abline(h=0)
plot(img$x, m$fitness(img$x, sys$x, sys$y), type="l", ylab="Fitness",
     las=1, ylim=c(-.01, .01))
points(sys$x, m$fitness(sys$x, sys$x, sys$y), col=col)
abline(h=0)
mtext("Trait", 1, 3, xpd=NA)
par(op)
