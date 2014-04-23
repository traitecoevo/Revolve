library(rootSolve)
library(Revolve)
library(deSolve)
library(plyr)

xx <- c(0.2, 0.7) # same as Georges
yy <- c(0.5, 0.5)
tt <- seq(0, 300, length=301)

m <- make_fox_2008()

## 1. At an arbitrarily chosen set of initial densities and states,
## here is the approach to equilibrium, with the values from runsteady
## added.
obj.eq <- m$equilibrium(xx, yy)
obj.tr <- m$run(xx, yy, tt)

op <- par(mfrow=c(2, 1), mar=c(4.1, 4.1, .5, .5))
matplot(obj.tr$t, obj.tr$y, type="l", xlab="", ylab="Abundance")
points(rep(max(tt), 2), obj.eq$y, col=1:2)
matplot(obj.tr$t, obj.tr$R, type="l", xlab="Time", ylab="Resource")
points(rep(max(tt), 2), obj.eq$R, col=1:2)
par(op)

## Now, look at the equilibrium of the single species setup.  For this
## we'll just start the second species at a density of 0 and with an
## arbitrary 0.5 for their state.
yy[2] <- 0
xx[2] <- 0.5
xx[1] <- 0.4

obj.eq <- m$equilibrium(xx, yy)
obj.tr <- m$run(xx, yy, tt)

op <- par(mfrow=c(2, 1), mar=c(4.1, 4.1, .5, .5))
matplot(obj.tr$t, obj.tr$y, type="l", xlab="", ylab="Abundance",
        ylim=c(0, 1))
points(rep(max(tt), 2), obj.eq$y, col=1:2)
matplot(obj.tr$t, obj.tr$R, type="l", xlab="Time", ylab="Resource",
        ylim=c(0, 1))
points(rep(max(tt), 2), obj.eq$R, col=1:2)
par(op)

## Single species equilibrium helper function.
eq1 <- function(x1, m, y1=0.5, ...) {
  x <- c(x1, 0.5) # arbitrarily 0.5 for absent species
  y <- c(y1, 0)
  m$equilibrium(x, y, ...)
}
## Over a range of trait values
xx <- seq(0, 1, length=101)
obj <- lapply(xx, eq1, m)
eq.R <- laply(obj, "[[", "R")
eq.y <- laply(obj, "[[", "y")

## Colours for resources 1 and 2:
col <- c("red", "blue")
col2 <- c("#ff000055", "#0000ff55")

matplot(xx, m$analytic$resource$limiting(xx), type="l",
        ylim=c(0, 1), lty=1, col=col)
matlines(xx, m$analytic$resource$nonlimiting(xx), lty=2, col=col)
alive <- aaply(eq.y, 1, max) > 1e-5
points(xx, eq.R[,1], pch=ifelse(alive, 19, 1), col=col2[1])
points(xx, eq.R[,2], pch=ifelse(alive, 19, 1), col=col2[2])

## Two species equilibrium helper function
eq2 <- function(x, m, y=c(.5, .5), ...) {
  m$equilibrium(x, y, ...)
}

# This is terribly slow at the moment.  No effort taken to speed it
# up.
if (interactive()) {
  x2 <- as.matrix(expand.grid(x1=xx, x2=xx))
  res <- alply(x2, 1, eq2, m, .progress="text")

  y2 <- laply(res, "[[", "y")

  # Different types at equilibrium.
  z2 <- y2 > 1e-6
  type <- z2[,1] * 2 + z2[,2]
  image(xx, xx, matrix(type, nrow=length(xx)))
  abline(1, -1)

  # Lines from the paper -- they don't match up exactly because of
  # things like numerical issues.  Still have to get to the bottom of
  # that.
  extra <- function(x, m) {
    p <- m$parameters$get()
    d <- p$d
    Y <- p$Y
    S <- p$S

    data.frame(B3a=d * Y[1,2] / (d * Y[1,1]) * x,
               B3b=1 - d * Y[2,1] / (d * Y[2,2]) * (1 - x),
               B5a=d * Y[1,2] / (Y[1,1] * (Y[1,2] * S[1] -
                 Y[2,2] * S[2] + d / (1 - x))),
               B5b=1 - d * Y[2,1] / (Y[2,2] * (Y[2,1] * S[2] -
                 Y[1,1] * S[1] + d / x)))
  }

  tmp <- extra(xx, m)

  lines(tmp$B3a, xx)
  lines(xx, tmp$B3b)
  lines(tmp$B5a, xx)
  lines(xx, tmp$B5b)
}
