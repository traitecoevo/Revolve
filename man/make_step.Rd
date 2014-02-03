\name{make_step}
\alias{make_step}
\alias{make_step_continuous}
\alias{make_step_deterministic}
\title{Step System Forward in Time}
\usage{
make_step(fitness, mutation, dt, mu, y_initial)

make_step_continuous(fitness, mutation, dt, mu, y_initial)

make_step_deterministic(fitness)
}
\arguments{
  \item{fitness}{Function for computing fitness.  Must take
  arguments \code{x_new} (phenotypes to compute fitness
  for), \code{x} (resident phenotypes) and \code{y}
  (resident abundances) and return the per-capita growth
  rate (g).  The population will step forward with an Euler
  step, so that the new \code{y} will be \code{y + y * g *
  dt}.}

  \item{mutation}{Function for generating mutants.  Must
  take arguments \code{traits} (resident phenotypes) and
  \code{mutants} (the number of mutants to generate, for
  each phenotype).  See \code{\link{make_mutation}}, which
  generates a useful function.}

  \item{dt}{Size of the time step}

  \item{mu}{Mutation rate.  On average there will be
  \code{mu*dt*y} mutations, with the actual number drawn
  from Poisson distribution.}

  \item{y_initial}{Abundance at which to introduce new
  mutants.}
}
\value{
A Function that moves the system forward in time; it takes
a list with elements \code{x} (traits), \code{y}
(abundances) and \code{t} (time) and returns a list with
these updated to the next steps.
}
\description{
Build a function that steps the system through one
generation. Order of events is: (1) resolve fitness (2)
introduce mutants.
}
\author{
Rich FitzJohn
}

