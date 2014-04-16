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
library(plyr)
library(numDeriv)

mat <- huisman_matrices(huisman_mat_1, huisman_mat_1)

m <- make_huisman_2001(mat, S=1)
x <- matrix(0.5, nrow=2)
y0 <- 1
t <- seq(0, 30, length=201)

res <- m$run(x, y0, t)
matplot(res$t, cbind(res$R, res$y), type="l", lty=1:2, xlab="Time",
        ylab="Abundance (red), Resource (black)")

eq <- m$single_equilibrium(x)
abline(h=unlist(eq), col=1:2, lty=3)

# Displace the solution from equilibrium and look at the new level of
# resources:
dy <- eq$y * 0.1
res1 <- m$run_fixed_density(x, eq$y - dy, t, eq$R)
res2 <- m$run_fixed_density(x, eq$y + dy, t, eq$R)
eq1 <- m$equilibrium_R(x, eq$y - dy)
eq2 <- m$equilibrium_R(x, eq$y + dy)
matplot(res1$t, cbind(res1$R, res2$R), type="l", col=2:3)
abline(h=c(eq$R, eq1$R, eq2$R), col=1:3, lty=3)

# Next, we start working towards the instantaneous growth rate of a
# new type at this equilibrium

# Look at the fitness landscape: how does the instantaneous growth
# rate look with respect to K:
xx <- seq(0, 1, length=101)
x2 <- rbind(xx, x[2], deparse.level=0) # K varying
x3 <- rbind(x[2], xx, deparse.level=0) # C varying

plot(xx, m$fitness(x2, x, eq$y, eq$R), type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

# Fitness does not vary with respect to c when rare:
plot(xx, m$fitness(x3, x, eq$y, eq$R), type="l")
abline(h=0, col="grey", lty=3)

# Fitness in an empty environment:
plot(xx, m$fitness(x2, x, 0, m$parameters$get()[["S"]]), type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)
