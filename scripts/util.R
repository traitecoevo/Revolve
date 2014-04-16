# Because we use abline, we can't recycle here, which is annoying.
abcline <- function(x, y, m, ...) {
  abline(y - x * m, m, ...)
}

make.transparent <- function(col, opacity=.5) {
  alpha <- opacity
  if ( length(alpha) > 1 && any(is.na(alpha)) ) {
    n <- max(length(col), length(alpha))
    alpha <- rep(alpha, length.out=n)
    col <- rep(col, length.out=n)
    ok <- !is.na(alpha)
    ret <- rep(NA, length(col))
    ret[ok] <- add.alpha(col[ok], alpha[ok])
    ret
  } else {
    tmp <- col2rgb(col)/255
    rgb(tmp[1,], tmp[2,], tmp[3,], alpha=alpha)
  }
}
