library(Revolve)

plot_huisman <- function(obj, ylim=NULL, lty=1:n, col=1:n) {
  np <- length(obj)
  n <- ncol(obj[[1]]$y)
  if (is.null(ylim)) {
    ylim <- range(sapply(obj, function(x) range(x$y)))
  }
  op <- par(mfrow=c(np, 1), mar=rep(.5, 4), oma=c(3, 4, 0, 0), cex=1)
  on.exit(par(op))

  for (i in seq_len(np)) {
    matplot(obj[[i]]$t, obj[[i]]$y, type="l", lty=lty, col=col,
            ylim=ylim, las=1, xaxt="n", xlab="", ylab="")
    axis(1, labels=i == np)
  }
  mtext("Time (d)", 1, line=2.1, xpd=NA)
  mtext("Population abundance", 2, line=2, outer=TRUE)
}

# First, results from the paper.  These are the matrices given in
# Appendix C.

# Figure2:
K.fig2 <- rbind(c(1.0, 0.6, 0.3),
                c(0.3, 1.0, 0.6),
                c(0.6, 0.3, 1.0))

C.fig2a <- rbind(c(0.07, 0.04, 0.04),
                 c(0.08, 0.10, 0.08),
                 c(0.10, 0.10, 0.14))

C.fig2b <- rbind(c(0.04, 0.07, 0.04),
                 c(0.08, 0.08, 0.10),
                 c(0.14, 0.10, 0.10))

C.fig2c <- rbind(c(0.04, 0.04, 0.07),
                 c(0.10, 0.08, 0.08),
                 c(0.10, 0.14, 0.10))

mat.fig2a <- huisman_matrices_fixed(K.fig2, C.fig2a)
mat.fig2b <- huisman_matrices_fixed(K.fig2, C.fig2b)
mat.fig2c <- huisman_matrices_fixed(K.fig2, C.fig2c)

S.fig2 <- c(6, 10, 14)
y0.fig2 <- rep(1, 3)
t.fig2 <- seq(0, 100, length=201)

m.fig2 <- list(make_huisman_2001(mat.fig2a, S=S.fig2),
               make_huisman_2001(mat.fig2b, S=S.fig2),
               make_huisman_2001(mat.fig2c, S=S.fig2))
res.fig2 <- lapply(m.fig2, function(m) m$run(NULL, y0.fig2, t.fig2))

plot_huisman(res.fig2, ylim=c(0, 100), lty=c(1, 2, 4))

# Figure 3
K.fig3a <- rbind(c(1.0, 0.9, 0.3),
                 c(0.3, 1.0, 0.9),
                 c(0.9, 0.3, 1.0))

K.fig3b <- K.fig3a
K.fig3b[rbind(c(1,2), c(2,3), c(3,1))] <- 0.5

K.fig3c <- K.fig3a
K.fig3c[rbind(c(1,2), c(2,3), c(3,1))] <- 0.4

C.fig3 <- rbind(c(0.04, 0.07, 0.04),
                c(0.08, 0.08, 0.10),
                c(0.14, 0.10, 0.10))

S.fig3 <- c(6, 10, 14)
y0.fig3 <- rep(1, 3)
t.fig3 <- seq(0, 400, length=201)

mat.fig3a <- huisman_matrices_fixed(K.fig3a, C.fig3)
mat.fig3b <- huisman_matrices_fixed(K.fig3b, C.fig3)
mat.fig3c <- huisman_matrices_fixed(K.fig3c, C.fig3)

m.fig3 <- list(make_huisman_2001(mat.fig3a, S=S.fig3),
               make_huisman_2001(mat.fig3b, S=S.fig3),
               make_huisman_2001(mat.fig3c, S=S.fig3))
res.fig3 <- lapply(m.fig3, function(m) m$run(NULL, y0.fig3, t.fig3))

# Note that this is different to the paper, because they run the last
# case out to to 2000, not 400.  But close enough for now.
plot_huisman(res.fig3, ylim=c(0, 150), lty=c(1, 2, 4))

# Figure 4:
K.fig4a <- rbind(c(1.0, 0.8, 0.4, 0.3, 0.2),
                 c(0.2, 1.0, 0.8, 0.4, 0.3),
                 c(0.3, 0.2, 1.0, 0.8, 0.4),
                 c(0.4, 0.3, 0.2, 1.0, 0.8),
                 c(0.8, 0.4, 0.3, 0.2, 1.0))

C.fig4 <- rbind(c(0.04, 0.07, 0.04, 0.04, 0.04),
                c(0.08, 0.08, 0.10, 0.08, 0.08),
                c(0.10, 0.10, 0.10, 0.14, 0.10),
                c(0.03, 0.03, 0.03, 0.03, 0.05),
                c(0.09, 0.07, 0.07, 0.07, 0.07))

K.fig4b <- K.fig4a
K.fig4b[rbind(c(1,2), c(2,3), c(3,4), c(4,5), c(5,1))] <- 0.6

S.fig4 <- c(6, 10, 14, 4, 9)
y0.fig4 <- rep(1, 5)
t.fig4a <- seq(0, 500, length=201)
t.fig4b <- seq(0, 3000, length=201)

mat.fig4a <- huisman_matrices_fixed(K.fig4a, C.fig4)
mat.fig4b <- huisman_matrices_fixed(K.fig4b, C.fig4)

m.fig4 <- list(make_huisman_2001(mat.fig4a, S=S.fig4),
               make_huisman_2001(mat.fig4b, S=S.fig4))

res.fig4 <- list(m.fig4[[1]]$run(NULL, y0.fig4, t.fig4a),
                 m.fig4[[2]]$run(NULL, y0.fig4, t.fig4b))
# Looks OK except that we start overstepping at some point and drive
# numbers negative, but apart from that we're all good really.
plot_huisman(res.fig4, ylim=c(0, 100))
