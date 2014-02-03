\name{to_grid}
\alias{to_grid}
\title{Bin Values}
\usage{
to_grid(x, y, grid)
}
\arguments{
  \item{x}{Vector of locations}

  \item{y}{Vector of densities/counts along x}

  \item{grid}{Grid over which to bin the values.}
}
\value{
A vector of length \code{length(grid) - 1} being the number
of things in each interval.
}
\description{
Bin x into a grid, where x may have different densities.
}
\author{
Rich FitzJohn
}

