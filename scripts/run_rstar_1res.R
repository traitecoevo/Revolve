# Fundamental result that we need: Species can coexist *at
# equilibrium* only if each species becomes limited by the resource
# for which it has, compared to its competitors, the highest
# requirement.

# Several species competing for one resource (p. 2684)
#
#   The species that has the lowest resource requirement for the
#   limiting resouce (i.e., the species with the lowest R^*_{1i}) will
#   displace all other species.
#
# So let's look at the behaviour of the system with a single resource
# and move towards understanding how the competition function might
# look like.
library(Revolve)

mat <- rstar_matrices(rstar_mat_1, rstar_mat_1)

m <- rstar(mat, S=1)
sys0 <- sys(matrix(0.5, nrow=2), 1)
t <- seq(0, 30, length=201)

res <- m$run(sys0, t)
matplot(res$t, cbind(res$R, res$y), type="l", lty=1:2, xlab="Time",
        ylab="Abundance (red), Resource (black)")

eq <- m$single_equilibrium(sys0$x)
abline(h=unlist(eq[c("R", "y")]), col=1:2, lty=3)

# Displace the solution from equilibrium and look at the new level of
# resources:
dy <- eq$y * 0.1

sys.y1 <- modifyList(eq, list(y=eq$y - dy))
sys.y2 <- modifyList(eq, list(y=eq$y + dy))

res.y1 <- m$run_fixed_density(sys.y1, t)
res.y2 <- m$run_fixed_density(sys.y2, t)

eq.y1 <- m$equilibrium_R(sys.y1)
eq.y2 <- m$equilibrium_R(sys.y2)

matplot(res.y1$t, cbind(res.y1$R, res.y2$R), type="l", col=2:3)
abline(h=c(eq$R, eq.y1$R, eq.y2$R), col=1:3, lty=3)

# Next, we start working towards the instantaneous growth rate of a
# new type at this equilibrium

# Look at the fitness landscape: how does the instantaneous growth
# rate look with respect to K:
x.mutant <- seq(0, 1, length=101)
x.K <- rbind(x.mutant, eq$x[2], deparse.level=0) # K varying
x.C <- rbind(eq$x[2], x.mutant, deparse.level=0) # C varying

plot(x.mutant, m$fitness(x.K, eq$x, eq$y, eq$R), type="l",
     xlab="Trait (K)", ylab="Fitness")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)

# Fitness does not vary with respect to c when rare:
plot(x.mutant, m$fitness(x.C, eq$x, eq$y, eq$R), type="l",
     xlab="Trait (C)", ylab="Fitness")
abline(h=0, col="grey", lty=3)

# Fitness in an empty environment:
plot(x.mutant,
     m$fitness(x.K, eq$x, 0, m$S), type="l",
     xlab="Trait (K)", ylab="Fitness in an empty environment")
abline(h=0, col="grey", lty=3)
abline(v=eq$x[1], lty=2)
