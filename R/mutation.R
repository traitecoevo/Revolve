##' Really basic mutation of traits, assuming normally distributed
##' mutants for real valued traits.
##'
##' For now, this just does mulivariate-normal mutation (continuous
##' traits), but may expand to deal with diallelic loci etc.
##'
##' For some uses (e.g., what we do in \code{tree}) we may want to
##' enable mutation to be on the log scale.
##'
##' @title Mutation
##' @param vcv A variance-covariance matrix or a variance value.
##' @return A generating function, taking \code{traits} and
##' \code{mutants} as arguments.  For the 1d case, \code{traits} is a
##' vector over types/species, while for the 2d case it is a matrix
##' with rows representing types/species and columns traits.  The
##' vector \code{mutants} is always over types/species and is the
##' number of mutants to generate for each type.
##' @author Rich FitzJohn
##' @export
##' @importFrom mvtnorm rmvnorm
make_mutation <- function(vcv) {
  if (length(vcv) == 1)
    make_mutation_normal_uni(vcv)
  else
    make_mutation_normal_multi(vcv)
}

## Very similar, but different enough that dealing with both types in
## one function is tedious.
make_mutation_normal_uni <- function(s2) {
  sd <- sqrt(s2)
  function(traits, mutants) {
    if (length(traits) != length(mutants))
      stop("Traits and mutants must be the same length")
    i <- rep.int(seq_along(mutants), mutants)
    traits[i] + rnorm(length(i), 0, sd)
  }
}

make_mutation_normal_multi <- function(vcv) {
  n_traits <- nrow(vcv)
  if (ncol(vcv) != n_traits)
    stop("Mutational variance covariance matrix must be square")
  if (ncol(vcv) < 2)
    stop("Matrix must have at least two rows/columns")
  ## NOTE: This is stricter than that used by rmvnorm()
  if (any(eigen(vcv, symmetric=TRUE)$values < 0))
    stop("Variance matrix not positive-definite")
  zero <- rep(0, n_traits)
  function(traits, mutants) {
    if (nrow(traits) != length(mutants))
      stop("Traits and mutants must be the same length")
    if (ncol(traits) != n_traits)
      stop("Traits must have ", n_traits, " columns")
    i <- rep.int(seq_along(mutants), mutants)
    parent <- traits[i,,drop=FALSE]
    if (nrow(parent) > 0)
      parent + rmvnorm(nrow(parent), zero, vcv)
    else
      parent
  }
}
