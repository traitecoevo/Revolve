##' Model implementing "Evolutionary branching under asymmetric
##' competition" by Kisdi, Journal of Theoretical Biology, 1999.
##'
##' Test case...
##' @title Kisdi 1999
##' @param capacity_type String, indicating which type of carrying
##' capacity to use; can be "linear", "gaussian" or "convex" (or any
##' contraction of these).
##' @param c Scaling factor for competition function
##' @param v Positioning parameter for competition function
##' @param k Rate parameter for competition function
##' @param beta Intercept of linear carrying capacity function
##' @param b Slope of linear carrying capacity function (also used in
##' convex carrying capacity function)
##' @param a Scaling factor for gaussian carrying capacity function
##' (also used in convex carrying capacity function)
##' @param m Mean of gaussian carrying capacity function
##' @param s2 Variance of gaussian carrying capacity function
##' @param d Parameter of convex carrying capacity function
##' @author Rich FitzJohn
##' @export
make_kisdi_1999 <- function(capacity_type="linear", c=2, v=1.2, k=4,
                            beta=1, b=1, a=1, m=0, s2=1, d=3.5) {
  defaults <- list(c=c, v=v, k=k, beta=beta, b=b, a=a, m=m, s2=s2, d=d)
  parameters <- make_parameters(defaults, environment())

  fitness <- function(x_new, x, y)
    capacity(x_new) - colSums(y * competition(x_new, x))
  ## Equation 1:
  competition <- function(x_new, x)
    c * (1 - 1 / (1 + v * exp(- k * (-outer(x, x_new, "-")))))

  ## Equation 12:
  capacity_linear <- function(x)
    beta - b * x
  ## Equation 14:
  capacity_gaussian <- function(x)
    a * exp(-(x - m)^2 / (2 * s2))
  ## Equation 17:
  capacity_convex <- function(x)
    -a - b * (x - sqrt(x*x + d))
  capacity_type <- match.arg(capacity_type, c("linear", "gaussian", "convex"))
  capacity <- switch(capacity_type,
                     linear=capacity_linear,
                     gaussian=capacity_gaussian,
                     convex=capacity_convex)

  ret <- list(fitness      = fitness,
              competition  = competition,
              capacity     = capacity,
              parameters   = parameters)
  class(ret) <- "model"
  ret
}
