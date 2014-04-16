# Fundamental result that we need: Species can coexist *at
# equilibrium* only if each species becomes limited by the resource
# for which it has, compared to its competitors, the highest
# requirement.

library(Revolve)
library(plyr)

source("util.R")

## 1. At an arbitrarily chosen set of initial densities and states,
## here is the approach to equilibrium, with the values from runsteady
## added.

mat <- huisman_matrices(huisman_mat_2_tradeoff, huisman_mat_2_tradeoff)
m <- make_huisman_2001(mat)
xx <- rbind(c(0.2, 0.7),
            c(0.2, 0.7))

rs <- m$Rstar(xx)
col <- c("blue", "red")
plot(NA, xlim=c(0, 1), ylim=c(0, 1), xlab="R1", ylab="R2")
abline(v=rs[1,], h=rs[2,], lty=3, col=col)
R.eq <- c(max(rs[1,]), max(rs[2,]))
points(R.eq[1], R.eq[2], pch=19)
segments(rs[1,], rs[2,], rs[1,], par("usr")[4], col=col)
segments(rs[1,], rs[2,], par("usr")[2], rs[2,], col=col)
# Note: slope of these lines is (1 - xx) / xx
abcline(R.eq[1], R.eq[2], m$C(xx)[2,1] / m$C(xx)[1,1], col=col[1], lty=2)
abcline(R.eq[1], R.eq[2], m$C(xx)[2,2] / m$C(xx)[1,2], col=col[2], lty=2)

cols <- make.transparent(c(col, "purple", "grey"), 0.2)

for (i in seq_len(200)) {
  yy <- c(1, 1) # starting density
  tt <- seq(0, 1000, length=301)
  S <- runif(2)
  m$parameters$set(list(S=S))
  obj.eq <- m$equilibrium(xx, yy)
  obj.tr <- m$run(xx, yy, tt)

  survived <- obj.eq$y > 1e-6
  if (sum(survived) == 1) {
    col <- cols[which(survived)]
  } else if (sum(survived) == 2) {
    col <- cols[3]
  } else {
    col <- cols[4]
  }
  
  lines(obj.tr$R[,1], obj.tr$R[,2], lty=2, col=col)
  points(c(S[1], obj.eq$R[1]), c(S[2], obj.eq$R[2]), pch=19, col=col, cex=.5)
}
