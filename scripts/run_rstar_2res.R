# Fundamental result that we need: Species can coexist *at
# equilibrium* only if each species becomes limited by the resource
# for which it has, compared to its competitors, the highest
# requirement.
library(Revolve)
library(testthat)

## First, trajectories in the two resource model for one and two
## species.

mat <- rstar_matrices(rstar_mat_2_tradeoff, rstar_mat_2_tradeoff)
m <- rstar(mat, 1)
col_died <- make_transparent("grey", .2)

## One species:
col1 <- "blue"
col1_tr <- make_transparent(col1, .2)

sys1 <- sys(matrix(0.2, nrow=2), y=1)

rstar_plot(m, sys1, col=col1)
for (i in 1:100) {
  rstar_trajectory(m, sys1, col=col1_tr, col_died=col_died, S=runif(2))
}

## Two species:
sys2 <- sys(cbind(sys1$x, 0.7, deparse.level=0), y=c(1, 1))
col2 <- c(col1, "red")
col2_tr <- make_transparent(col2, .2)

rstar_plot(m, sys2, col=col2)
for (i in 1:100) {
  rstar_trajectory(m, sys2, col=col2_tr, col_died=col_died, S=runif(2))
}

## Now, displace the single species system from equilibrium and look
## at the new level of resources:
m$S <- rep(1, 1)

rstar_plot(m, sys1, col=col2)

eq <- m$single_equilibrium(sys1$x)
rstar_trajectory(m, sys1, S=c(1, 1))

x1 <- sys1$x
x2 <- cbind(x1, .3, deparse.level=0)

R1 <- m$single_equilibrium(x2[,1,drop=FALSE])$R # .8166, .2666
R2 <- m$single_equilibrium(x2[,2,drop=FALSE])$R # .6714, .2333
expect_that(m$single_equilibrium(x2)$R,
            equals(cbind(R1, R2, deparse.level=0)))

xx <- matrix(rep(seq(0, 1, by=.1), each=2), 2)
m$single_equilibrium(xx)

dy <- eq$y * 0.1
t <- seq(0, 30, length=201)

sys.y1 <- modifyList(eq, list(y=eq$y - dy))
sys.y2 <- modifyList(eq, list(y=eq$y + dy))

res.y1 <- m$run_fixed_density(sys.y1, t)
res.y2 <- m$run_fixed_density(sys.y2, t)

eq.y1 <- m$equilibrium_R(sys.y1)
eq.y2 <- m$equilibrium_R(sys.y2)

matplot(res.y1$t, cbind(res.y1$R, res.y2$R), type="l",
        col=c(1,1,2,2), lty=c(1,2,1,2))
abline(h=eq.y1$R, lty=3, col=1)
abline(h=eq.y2$R, lty=3, col=2)

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
x.mutant <- rbind(seq(0, 1, length=301), eq$x[2])
w.mutant <- m$fitness(x.mutant, eq$x, eq$y, eq$R)

plot(x.mutant[1,], w.mutant, type="l")
abline(h=0, col="grey", lty=3)
abline(v=eq$x, lty=2)
