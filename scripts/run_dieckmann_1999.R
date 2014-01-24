library(Revolve)

make_step <- function(fitness, mutation, dt, mu, y_initial) {
  function(sys) {
    # births/time/individual * individuals * time -> births
    dy <- fitness(sys$x, sys$x, sys$y) * dt * sys$y
    # mutations/time/individual * individuals * time -> mutation
    mutants <- mutation(sys$x, rpois(length(sys$y), mu * dt * sys$y))

    sys$y <- sys$y + dy
    if (length(mutants) > 0) {
      sys$x <- c(sys$x, mutants)
      sys$y <- c(sys$y, rep_len(y_initial, length(mutants)))
    }
    sys$t <- sys$t + dt

    sys
  }
}

make_cleanup <- function(eps) {
  function(sys) {
    keep <- sys$y >= eps
    if (!all(keep)) {
      sys$x <- sys$x[keep]
      sys$y <- sys$y[keep]
    }
    sys
  }
}

run <- function(sys, n.steps, step, cleanup, print.every=0) {
  res <- vector("list", length(n.steps+1))
  res[[1]] <- sys
  for (i in seq_len(n.steps)) {
    if (print.every > 0 && i %% print.every == 0)
      message(i)
    sys <- step(sys)
    sys <- cleanup(sys)
    res[[i+1]] <- sys
  }
  res
}

x_initial <- -0.3  # initial location
y_initial <- 1     # initial mutant frequency (ibm: 1)
y_eps     <- 1     # kill threshold (ibm: 1)
dt        <- 1     # Size of the time step (discrete gen == Euler?)
u         <- 0.001 # Fraction of births that are mutants
mutation_var <- (1/20)^2

m1 <- make_dieckmann_1999(s2_C = 1.5)
m2 <- make_dieckmann_1999()

mu <- m1$parameters$get()$r * u
step1 <- make_step(m1$fitness, make_mutation(mutation_var),
                   dt, u, y_initial)
step2 <- make_step(m2$fitness, make_mutation(mutation_var),
                   dt, u, y_initial)
cleanup <- make_cleanup(y_eps)

set.seed(1)
res1 <- run(list(x=x_initial, y=y_initial, t=0), 1000,
            step1, cleanup)
res2 <- run(list(x=x_initial, y=y_initial, t=0), 1000,
            step2, cleanup)

xr <- c(-1, 1)
ng <- 101
xx <- seq(xr[1], xr[2], length.out=ng)

j <- seq(1, length(res1), by=10)
mat1 <- sapply(j, function(i)
               to_grid(res1[[i]]$x, res1[[i]]$y, xx))
mat1[mat1 == 0] <- NA

mat2 <- sapply(j, function(i)
               to_grid(res2[[i]]$x, res2[[i]]$y, xx))
mat2[mat2 == 0] <- NA

t <- sapply(res1, function(x) x$t)

cols <- grey((32:0)/32)

last <- function(x) x[[length(x)]]
sys1 <- last(res1)
sys2 <- last(res2)

op <- par(mfrow=c(3, 2), mar=c(4.1, 4.5, .5, .5))
image(xx, t[j], mat1, col=cols, xlab="Trait", ylab="Time"); box()
image(xx, t[j], mat2, col=cols, xlab="Trait", ylab="Time"); box()

plot(sys1, xlab="Trait", ylab="Number", xlim=xr)
plot(sys2, xlab="Trait", ylab="Number", xlim=xr)

xx1 <- sort(c(seq(xr[1], xr[2], length=301), sort(sys1$x)))
xx2 <- sort(c(seq(xr[1], xr[2], length=301), sort(sys2$x)))

plot(xx1, m1$fitness(xx1, sys1$x, sys1$y), type="l", xlab="Trait",
     ylab="Fitness")
points(sys1$x, m1$fitness(sys1$x, sys1$x, sys1$y))
abline(h=0)

plot(xx2, m2$fitness(xx2, sys2$x, sys2$y), type="l", xlab="Trait",
     ylab="Fitness")
points(sys2$x, m2$fitness(sys2$x, sys2$x, sys2$y))
abline(h=0)
par(op)

## Look at the fitness landscape:
x <- seq(-2,2, length.out=50)
m <- make_dieckmann_1999()
plot(x, m$capacity(x), type="o", xlab = "trait")

par(mfrow=c(1,5), oma=rep(4,4))
for (x_res in seq(-2,0, length.out=5)) {
  plot(x, m$fitness(x, x_res,  m$capacity(x_res)), type="l",
       ann=FALSE, ylim=c(-1,1))
  points(x_res,0, pch=16, col="red")
  abline(h=0, col="grey", lty="dashed")
}
title("Change in fitness landscape approaching branching point", outer=TRUE)
mtext("Trait value",1, outer=TRUE, line=0)
mtext("Fitness",2, outer=TRUE, line=0)
