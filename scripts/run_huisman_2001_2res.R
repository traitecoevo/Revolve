# Fundamental result that we need: Species can coexist *at
# equilibrium* only if each species becomes limited by the resource
# for which it has, compared to its competitors, the highest
# requirement.

library(Revolve)
library(plyr)

source("util.R")

mat <- huisman_matrices(huisman_mat_2_tradeoff, huisman_mat_2_tradeoff)
m <- make_huisman_2001(mat)

xx <- matrix(0.3, nrow=2)

rs <- m$Rstar(xx)
col <- "blue"
plot(NA, xlim=c(0, 1), ylim=c(0, 1), xlab="R1", ylab="R2")
abline(v=rs[1,], h=rs[2,], lty=3, col=col)
R.eq <- c(max(rs[1,]), max(rs[2,]))
R.eq <- c(max(rs[1,]), max(rs[2,]))
points(R.eq[1], R.eq[2], pch=19)
segments(rs[1,], rs[2,], rs[1,], par("usr")[4], col=col)
segments(rs[1,], rs[2,], par("usr")[2], rs[2,], col=col)
# Note: slope of these lines is (1 - xx) / xx, and in the single
# species case they pass through the origin.
abcline(R.eq[1], R.eq[2], m$C(xx)[2,1] / m$C(xx)[1,1], col=col[1], lty=2)
points(0, 0)

cols <- make.transparent(c(col, "grey"), 0.2)

set.seed(1)
for (i in seq_len(50)) {
  yy <- 1 # starting density
  tt <- seq(0, 300, length=301)
  S <- runif(2)
  m$parameters$set(list(S=S))
  obj.eq <- m$single_equilibrium(xx)
  obj.tr <- m$run(xx, yy, tt)

  survived <- obj.eq$y > 1e-6
  if (sum(survived) == 1) {
    col <- cols[1]
  } else {
    col <- cols[2]
  }

  lines(obj.tr$R[,1], obj.tr$R[,2], lty=2, col=col)
  points(c(S[1], obj.eq$R[1]), c(S[2], obj.eq$R[2]), pch=19, col=col, cex=.5)
}

# Displace the solution from equilibrium and look at the new level of
# resources:
m$parameters$set(list(S=c(1,1)))
eq <- m$single_equilibrium(xx)
dy <- eq$y * 0.1
t <- seq(0, 30, length=201)
res1 <- m$run_fixed_density(xx, eq$y - dy, t, eq$R)
res2 <- m$run_fixed_density(xx, eq$y + dy, t, eq$R)
eq1  <- m$equilibrium_R(xx, eq$y - dy, eq$R)
eq2  <- m$equilibrium_R(xx, eq$y + dy, eq$R)

matplot(res1$t, cbind(res1$R, res2$R), type="l", col=
        c(1,1,2,2), lty=c(1,2,1,2))
abline(h=eq1$R, lty=3, col=1)
abline(h=eq2$R, lty=3, col=2)

# Solving this analytically is going to be harder than for the two
# species case.  However, it's possible that we can actually skip the
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
x2 <- rbind(seq(0, 1, length=301), xx[2])

plot(x2[1,], m$fitness(x2, xx, eq$y, eq$R), type="l")
abline(h=0, col="grey", lty=3)
abline(v=xx, lty=2)
