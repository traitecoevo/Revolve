\name{rstar_plot}
\alias{rstar_plot}
\alias{rstar_trajectory}
\title{Plot Rstar ZNGIs}
\usage{
rstar_plot(m, sys, xlim = c(0, 1), ylim = c(0, 1), col = seq_len(ncol(x)))

rstar_trajectory(m, sys, col = seq_len(ncol(sys$x)), eps = 1e-06,
  t_max = 100, t_len = 101, S = NULL, col_died = "grey")
}
\arguments{
  \item{m}{Rstar model, made by \code{make_rstar}}

  \item{sys}{Model state (x is used)}

  \item{xlim}{X axis limits (resource 1)}

  \item{ylim}{Y axis limits (resource 2)}

  \item{col}{Vector of colours along species}

  \item{eps}{Density at which species are considered
  extinct}

  \item{t_max}{Maximum time to run simulation until}

  \item{t_len}{Number of time intervals to run simulation
  over}

  \item{S}{Resource supply rates (and initial resource
  levels)}

  \item{col_died}{Colour of communities that go extinct.}
}
\description{
Plot Rstar ZNGIs
}
\details{
Will change
}
\author{
Rich FitzJohn
}

