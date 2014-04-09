\name{make_kisdi_1999}
\alias{make_kisdi_1999}
\title{Kisdi 1999}
\usage{
make_kisdi_1999(intrinsic_growth_type = "linear", c = 2, v = 1.2, k = 4,
  beta = 1, b = 1, a = 1, m = 0, s2 = 1, d = 3.5)
}
\arguments{
  \item{intrinsic_growth_type}{String, indicating which
  type of intrinsic growth to use; can be "linear",
  "gaussian" or "convex" (or any contraction of these).
  Intrinsic growth is the population growth rate in the
  absence of competition (including self-competition).
  Units: inverse time}

  \item{c}{Death rate, to be scaled by competition
  function. Units: inverse time * inverse population}

  \item{v}{Positioning parameter for competition function}

  \item{k}{Trait scaling parameter for competition
  function. Units: inverse trait.}

  \item{beta}{Intercept of linear intrinsic growth
  function}

  \item{b}{Slope of linear intrinsic growth function (also
  used in convex intrinsic growth function)}

  \item{a}{Scaling factor for gaussian intrinsic growth
  function (also used in convex intrinsic growth function)}

  \item{m}{Mean of gaussian intrinsic growth function}

  \item{s2}{Variance of gaussian intrinsic growth function}

  \item{d}{Parameter of convex intrinsic growth function}
}
\description{
Model implementing "Evolutionary branching under asymmetric
competition" by Kisdi, Journal of Theoretical Biology,
1999.
}
\details{
Test case...
}
\author{
Rich FitzJohn
}

