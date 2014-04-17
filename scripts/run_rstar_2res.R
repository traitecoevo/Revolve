# Fundamental result that we need: Species can coexist *at
# equilibrium* only if each species becomes limited by the resource
# for which it has, compared to its competitors, the highest
# requirement.
library(Revolve)
library(plyr)

## First, trajectories in the two resource model for one and two
## species.

mat <- rstar_matrices(rstar_mat_2_tradeoff, rstar_mat_2_tradeoff)
m <- make_rstar(mat)
col_died <- make_transparent("grey", .2)

## One species:
x1 <- matrix(0.2, nrow=2)
col1 <- "blue"
col1_tr <- make_transparent(col1, .2)

rstar_plot(m, x1, col=col1)
for (i in 1:100) {
  rstar_trajectory(m, x1, col=col1_tr, col_died=col_died, S=runif(2))
}

## Two species:
x2 <- cbind(x1, 0.7, deparse.level=0)
col2 <- c(col1, "red")
col2_tr <- make_transparent(col2, .2)

rstar_plot(m, x2, col=col2)
for (i in 1:100) {
  rstar_trajectory(m, x2, col=col2_tr, col_died=col_died, S=runif(2))
}

## Now, displace the single species system from equilibrium and look
## at the new level of resources:
m$parameters$set(list(S=c(1,1)))
eq <- m$single_equilibrium(x1)
dy <- eq$y * 0.1
t <- seq(0, 30, length=201)
res1 <- m$run_fixed_density(x1, eq$y - dy, t, eq$R)
res2 <- m$run_fixed_density(x1, eq$y + dy, t, eq$R)
eq1  <- m$equilibrium_R(x1, eq$y - dy, eq$R)
eq2  <- m$equilibrium_R(x1, eq$y + dy, eq$R)

matplot(res1$t, cbind(res1$R, res2$R), type="l",
        col=c(1,1,2,2), lty=c(1,2,1,2))
abline(h=eq1$R, lty=3, col=1)
abline(h=eq2$R, lty=3, col=2)

# Solving this analytically is going to be harder than for the one
# resource case.  However, it's possible that we can actually skip the
# hard work and solve this for each resource separately and then see
# which ones are plausible solutions (i.e., *assume* that the first or
# the second resource is limiting and then see if that makes sense).
# Realistically I'm going to need an pen and paper for this though.
#
# The alternative approach is to solve numerically.  That should be
# pretty easy to do, either with nleqslv or with runsteady, or
# possibly a uniroot approach on the equation itself.
#
# Even if a semianalytic solution is found, we're going to need the
# same trick before to work out where it lands using the supply curve.

# Look at the fitness landscape: how does the instantaneous growth
# rate look with respect to K (& C):
x.mutant <- rbind(seq(0, 1, length=301), x1[2])
w.mutant <- m$fitness(x.mutant, x1, eq$y, eq$R)

plot(x.mutant[1,], w.mutant, type="l")
abline(h=0, col="grey", lty=3)
abline(v=x1, lty=2)
