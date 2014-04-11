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

m <- make_huisman_2001(S=1, matrices=huisman_matrices(huisman_mat_1))
x <- matrix(0.5, nrow=2)
y0 <- 1
t <- seq(0, 30, length=201)

res <- m$run(x, y0, t)
matplot(res$t, cbind(res$R, res$y), type="l", lty=1:2, xlab="Time",
        ylab="Abundance (red), Resource (black)")

eq <- m$single_equilibrium(x)
abline(h=unlist(eq), col=1:2, lty=3)

# Next, we start working towards the instantaneous growth rate of a
# new type at this equilibrium

# Look at the fitness landscape: how does the instantaneous growth
# rate look with respect to K:
xx <- seq(0, 1, length=51)
x2 <- rbind(xx, x[2], deparse.level=0) # K varying
x3 <- rbind(x[2], xx, deparse.level=0) # C varying

plot(xx, m$fitness(x2, x, eq$y, eq$R), type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

# Fitness does not vary with respect to c when rare, because species
# only consume resources when not rare.  Which doesn't mean that there
# will be no competition between species with different c values but
# the same K values of course: a different c value can drive a
# different species extinct!
plot(xx, m$fitness(x3, x, eq$y, eq$R), type="l")
abline(h=0, col="grey", lty=3)

# But then there's a broader view of what competition is here: The
# species with K above the resident species are all strongly competed
# against by the resident species, regardless of what the derivatives
# show -- because they can't grow there.  That's the insight that we'd
# get from the subtraction case:
plot(xx, m$fitness(x2, x, 0, m$parameters$get()[["S"]]), type="l")
abline(h=0, col="grey", lty=3)
abline(v=x[1], lty=2)

# So we can probably force a density of individuals and have the
# resources equilibrate (or go to zero).  That's not really the same
# as running the models (and violates the lack of separation of
# timescales) but might be a reasonable thing to do?

# We'll probably want to do this in a more formal way in a new model,
# I suspect.

# For now, because we're looking at derivatives only, we'll assume
# that we are sufficiently close to equilibrium that the change in
# resource is second order in the change in density so goes to zero
# fast enough.
jacobian(function(dy) m$fitness(x2, x, eq$y + dy, eq$R), 0)

# Except we *can't* do that; there is no dependence on R.  So we need
# to replumb the model a bit.
