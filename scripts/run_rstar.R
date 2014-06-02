library(Revolve)

mat <- rstar_matrices(rstar_mat_2_tradeoff, rstar_mat_2_tradeoff)
m <- rstar(mat, 1)

tt <- seq(0, 300, length=301)
sys0 <- sys(x=rbind(c(0.2, 0.7),
              c(0.2, 0.7)),
            y=c(0.3, 0.7))

## 1. At an arbitrarily chosen set of initial densities and states,
## here is the approach to equilibrium, with the values from runsteady
## added.
obj.eq <- m$equilibrium(sys0)
obj.tr <- m$run(sys0, tt)

op <- par(mfrow=c(2, 1), mar=c(4.1, 4.1, .5, .5))
matplot(obj.tr$t, obj.tr$y, type="l", xlab="", ylab="Abundance")
points(rep(max(tt), 2), obj.eq$y, col=1:2)
matplot(obj.tr$t, obj.tr$R, type="l", xlab="Time", ylab="Resource")
points(rep(max(tt), 2), obj.eq$R, col=1:2)
par(op)

## Now, look at the equilibrium of the single species setup.  For this
## we'll just start the second species at a density of 0 and with an
## arbitrary 0.5 for their state.
sys0$y[2] <- 0
sys0$x[,2] <- 0.5
sys0$x[,1] <- 0.4

obj.eq <- m$equilibrium(sys0)
obj.tr <- m$run(sys0, tt)

op <- par(mfrow=c(2, 1), mar=c(4.1, 4.1, .5, .5))
matplot(obj.tr$t, obj.tr$y, type="l", xlab="", ylab="Abundance",
        ylim=c(0, 2))
points(rep(max(tt), 2), obj.eq$y, col=1:2)
matplot(obj.tr$t, obj.tr$R, type="l", xlab="Time", ylab="Resource",
        ylim=c(0, 2))
points(rep(max(tt), 2), obj.eq$R, col=1:2)
par(op)
