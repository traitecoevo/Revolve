library(Revolve)

## Some extra shared parameters:
y_initial <- 1        # initial mutant frequency
dt        <- 1        # Size of the time step
mu        <- 0.001    # Mutation rate
s2_mu     <- (1/10)^2 # Mutational variance

## Two instances of the model; one has the competition kernel wider
## than the resource kernel (m1), which should prevent coexistance and
## the other has the competition kernel narrower than the resource
## kernel (m2), which should allow coexistance.
m1 <- make_dieckmann_1999(s2_C = 1.5, s2_K = 1)
m2 <- make_dieckmann_1999(s2_C = 0.4, s2_K = 1)

## Initial system state; a single individual of a single strategy with
## trait -0.3:
sys0 <- list(x=-0.3, y=y_initial, t=0)

## Remove species with < 1 individual.
cleanup <- make_cleanup(1)
## Mutations normally distributed with variance s2_mu.
mutation <- make_mutation(s2_mu)

## Functions to step the model forward in time
step1 <- make_step(m1$fitness, mutation, dt, mu, y_initial)
step2 <- make_step(m2$fitness, mutation, dt, mu, y_initial)

## Run the system for 3000 steps:
# set.seed(1)
res1 <- run(sys0, 3000, step1, cleanup)
res2 <- run(sys0, 10000, step2, cleanup)

## Plot the community over time; binned into
cols <- grey((32:0)/32)
col <- "#00000055"

img1 <- discretise(res1)
img2 <- discretise(res2)
xr <- range(img1$x, img2$x)

sys1 <- res1[[length(res1)]]
sys2 <- res2[[length(res2)]]

op <- par(mfrow=c(3, 2), mar=c(2.6, 3.1, .5, .5), oma=c(2.1, 2.1, 0, 0))
## Traits and population density through time:
image(img1$x, img1$t, img1$y, col=cols, ylab="", xlim=xr,
      xaxs="r", las=1); box()
mtext("Time", 2, 3, xpd=NA)
image(img2$x, img2$t, img2$y, col=cols, ylab="",
      xaxs="r", las=1); box()

## Resident types and their densities:
plot(sys1, xlim=xr, pch=19, col=col, ylab="", las=1)
mtext("Number", 2, 3, xpd=NA)
plot(sys2, xlim=xr, pch=19, col=col, ylab="", las=1)

## Fitness of the resident types:
plot(img1$x, m1$fitness(img1$x, sys1$x, sys1$y), type="l", ylab="",
     xlim=xr, las=1)
points(sys1$x, m1$fitness(sys1$x, sys1$x, sys1$y), col=col)
abline(h=0)
mtext("Fitness", 2, 3, xpd=NA)
mtext("Trait", 1, 3, xpd=NA)

plot(img2$x, m2$fitness(img2$x, sys2$x, sys2$y), type="l", ylab="",
     xlim=xr, las=1)
points(sys2$x, m2$fitness(sys2$x, sys2$x, sys2$y), col=col)
abline(h=0)
mtext("Trait", 1, 3, xpd=NA)
par(op)

# ## Here is the fitness as a single phenotype approaches the branching point.
# op <- par(mfrow=c(1,5), oma=c(2.5, 2.5, 2, 0), mar=c(2.1, 2.1, .5, .5))
# for (x_res in seq(-2, 0, length.out=5)) {
#   curve(m2$fitness(x, x_res, m2$capacity(x_res)),
#         type="l", xlim=c(-2,1), ann=FALSE, ylim=c(-1,1), las=1)
#   points(x_res,0, pch=16, col="red")
#   abline(h=0, col="grey", lty="dashed")
# }
# title("Change in fitness landscape approaching branching point", outer=TRUE)
# mtext("Trait value",1, outer=TRUE, line=.5)
# mtext("Fitness",2, outer=TRUE, line=.5)
# par(op)
