\name{make_kisdi_1999}
\alias{make_kisdi_1999}
\title{Kisdi 1999}
\usage{
make_kisdi_1999(capacity_type = "linear", c = 2, v = 1.2, k = 4,
  beta = 1, b = 1, a = 1, m = 0, s2 = 1, d = 3.5)
}
\arguments{
  \item{capacity_type}{String, indicating which type of
  carrying capacity to use; can be "linear", "gaussian" or
  "convex" (or any contraction of these).}

  \item{c}{Scaling factor for competition function}

  \item{v}{Positioning parameter for competition function}

  \item{k}{Rate parameter for competition function}

  \item{beta}{Intercept of linear carrying capacity
  function}

  \item{b}{Slope of linear carrying capacity function (also
  used in convex carrying capacity function)}

  \item{a}{Scaling factor for gaussian carrying capacity
  function (also used in convex carrying capacity
  function)}

  \item{m}{Mean of gaussian carrying capacity function}

  \item{s2}{Variance of gaussian carrying capacity
  function}

  \item{d}{Parameter of convex carrying capacity function}
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

