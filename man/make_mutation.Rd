\name{make_mutation}
\alias{make_mutation}
\title{Mutation}
\usage{
make_mutation(vcv)
}
\arguments{
  \item{vcv}{A variance-covariance matrix or a variance
  value.}
}
\value{
A generating function, taking \code{traits} and
\code{mutants} as arguments.  For the 1d case,
\code{traits} is a vector over types/species, while for the
2d case it is a matrix with rows representing types/species
and columns traits.  The vector \code{mutants} is always
over types/species and is the number of mutants to generate
for each type.
}
\description{
Really basic mutation of traits, assuming normally
distributed mutants for real valued traits.
}
\details{
For now, this just does mulivariate-normal mutation
(continuous traits), but may expand to deal with diallelic
loci etc.

For some uses (e.g., what we do in \code{tree}) we may want
to enable mutation to be on the log scale.
}
\author{
Rich FitzJohn
}

