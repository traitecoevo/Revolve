## As described in the paper, modelling waiting times for mutation:
library(Revolve)

eps_ext   <- 1e-6
y_initial <- eps_ext
s2_mu     <- (1/10)^2
mu        <- 1e-3

m <- make_kisdi_1999()
mutation <- make_mutation(s2_mu)
step_eq <- make_step_equilibrium(m$fitness, method="runsteady")
step <- make_step_mutation_limited(step_eq, mutation, mu, y_initial)
cleanup <- make_cleanup(eps_ext)

sys0 <- list(x=-1, y=y_initial, t=0)

set.seed(1)
res <- run(sys0, 300, step, cleanup, 50)

## Plot the community over time; binned into
cols <- grey(((32:0)/32)^2)
col <- "#ff000055"

sys <- res[[length(res)]]
img <- discretise(res, sample=1)
xr <- range(img$x)

op <- par(mfrow=c(3, 1), mar=c(2.6, 4.6, .5, .5), oma=c(2.1, 0, 0, 0))
## Traits and population density through time:
image(img$x, img$t, img$y, col=cols, ylab="Time", las=1, xaxs="r"); box()
## Resident types and their densities:
plot(sys, xlim=xr, pch=19, col=col, ylab="Number", las=1)
## Fitness of the resident types:
plot(img$x, m$fitness(img$x, sys$x, sys$y), type="l", ylab="Fitness",
     xlim=xr, las=1)
points(sys$x, m$fitness(sys$x, sys$x, sys$y), col=col, pch=19)
abline(h=0)
mtext("Trait", 1, 3, xpd=NA)
par(op)

sys <- res[[length(res)]]
set.seed(1)

## NOTE: Notation from Ito (2007) -- w is not a good choice here.
w <- sys$y * mu

repeat {
  ## Update evolutionary time; ignore the ecological time.
  sys$t <- sys$t + rexp(1, sum(w))

  ## Create a single mutation:
  i <- sample(length(w), 1, prob=w)
  mutant <- mutation(sys$x, as.integer(seq_along(w) == i))

  ## Check that the mutant has positive fitness at the current point:
  if (m$fitness(mutant, sys$x, sys$y) >= 0)
    break
}

## Add the mutant to the population and run to equilibrium:
sys$x <- c(sys$x, mutant)
sys$y <- c(sys$y, y_initial)
tmp <- step_eq(sys)

## Check that the parent of the mutant could still invade at
## equilibrium:
sys$x[i]
